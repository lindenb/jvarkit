/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.samtools.util;

import java.math.BigInteger;
import java.util.Optional;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.function.UnaryOperator;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;


/**
 * parser of intervel
 * @author lindenb
 *
 */
public class IntervalParser implements Function<String,Optional<SimpleInterval>> {
static final Logger LOG = Logger.of(IntervalParser.class);

public static final String OPT_DESC= 
			"An interval as the following syntax : \"chrom:start-end\". Some jvarkit programs also allow the following syntax : \"chrom:middle+extend\"  or \"chrom:start-end+extend\" or \"chrom:start-end+extend-percent%\"."
			+ "A program might use a Reference sequence to fix the chromosome name (e.g: 1->chr1)" ;
private final int MAX_CONTIG_LENGTH = Integer.MAX_VALUE - 1000;
private final SAMSequenceDictionary _dict;
private boolean enableWholeContig=false;
private boolean trimToContigLength = false;
private boolean enable_single_point = false;
private boolean enable_extend = false;
private boolean enable_size_suffix = false;
private final DistanceParser distanceParser=new DistanceParser();
private boolean thow_on_error=false;
private UnaryOperator<String> ctgNameConverter;
public IntervalParser() {
	this(null);
	}

public IntervalParser(final SAMSequenceDictionary dict) {
	this._dict=dict;
	this.ctgNameConverter=(dict==null?null:ContigNameConverter.fromOneDictionary(dict));
	}


private BigInteger parseBigInteger(final String s) {
	if(enable_size_suffix) return this.distanceParser.parseAsBigInteger(s);
	return new BigInteger(s);
	}

public IntervalParser enableWholeContig(boolean b) {
	this.enableWholeContig=b;
	return this;
	}

public IntervalParser enableSinglePoint(boolean b) {
	this.enable_single_point = b;
	return this;
	}

public IntervalParser raiseExeceptionOnError(boolean b) {
	this.thow_on_error = b;
	return this;
	}

public IntervalParser raiseExeceptionOnError() {
	return raiseExeceptionOnError(true);
	}


public IntervalParser enableSinglePoint() {
	return enableSinglePoint(true);
	}

public IntervalParser enableWholeContig() {
	return enableWholeContig(true);
	}

public IntervalParser enableExtend(boolean b) {
	this.enable_extend = b;
	return this;
	}

public IntervalParser enableExtend() {
	return enableExtend(true);
	}


public IntervalParser enableSizeSuffix(boolean b) {
	this.enable_size_suffix = b;
	return this;
	}

public IntervalParser enableSizeSuffix() {
	return enableSizeSuffix(true);
	}

public boolean hasDictionary() {
	return getDictionary()!=null;
	}
private SAMSequenceDictionary getDictionary() {
	return this._dict;
	}

private <T> Optional<T> returnErrorOrNullInterval(final String msg) {
	if(!this.thow_on_error) return Optional.empty();
	throw new IllegalArgumentException(String.valueOf(msg));
	}

public static Supplier<IllegalArgumentException> exception(final String s) {
	return ()->new IllegalArgumentException(String.valueOf(s));
	}

private static int BigIntegerToInt(final BigInteger vbi)
	{
	final BigInteger INT_MAX=new BigInteger(String.valueOf(Integer.MAX_VALUE));
	final BigInteger INT_MIN=new BigInteger(String.valueOf(Integer.MIN_VALUE));

	if(INT_MAX.compareTo(vbi)<=0)
		{
		throw new IllegalArgumentException("value "+vbi+" cannot be greater or equal to "+INT_MAX);
		}
	if(vbi.compareTo(INT_MIN)<=0)
		{
		throw new IllegalArgumentException("value "+vbi+" cannot be lower or equal to "+INT_MIN);
		}
	return vbi.intValue();
	}
		


@Override
public Optional<SimpleInterval> apply(final String s)
	{
	SAMSequenceRecord ssr = null;
	final int colon=s.lastIndexOf(':');
	if(colon<1 || colon+1==s.length()) {
		/* chromosome alone */
		if(colon==-1 && this.enableWholeContig)
			{
			if(this._dict!=null) {
				ssr= this._dict.getSequence(this.ctgNameConverter.apply(s.trim()));
				if(ssr==null) return returnErrorOrNullInterval(JvarkitException.ContigNotFoundInDictionary.getMessage(s, this._dict));
				return Optional.of(new SimpleInterval(ssr));
				}
			else
				{
				return Optional.of(new SimpleInterval(s.trim(), 1, MAX_CONTIG_LENGTH));
				}
			}
		else
			{
			return returnErrorOrNullInterval("Expect CHR:START-END. Cannot find colon in "+s);
			}
		}
	final int hyphen = s.indexOf('-',colon+1);
	final int plus = this.enable_extend ? s.indexOf('+',colon+1) : -1;
	if(hyphen==-1 && plus==-1) {
		if(!this.enable_single_point) return returnErrorOrNullInterval("Cannot find hyphen or plus in "+s+". Single point position like 'chr1:234' are not allowed");
		/** single point mutation */
		String contig = s.substring(0,colon).trim();
		if(this._dict!=null) {
			ssr= this._dict.getSequence(this.ctgNameConverter.apply(contig));
			if(ssr==null) return returnErrorOrNullInterval(JvarkitException.ContigNotFoundInDictionary.getMessage(contig, this._dict));
			contig = ssr.getContig();
			}
		final BigInteger bPos = parseBigInteger(s.substring(colon+1).trim());
		final int pos= BigIntegerToInt(bPos);
		return Optional.of(new SimpleInterval(contig, pos, pos)); 
		}
	if(hyphen!=-1 && plus!=-1)
		{
		// chr1:123-+
		if(plus <= hyphen+1) return returnErrorOrNullInterval("both hyphen and plus in "+s);
		}
		
	BigInteger start=null,end=null;
	try {
		
		String contig = s.substring(0,colon).trim();
		if(this._dict!=null) {
			ssr = this._dict.getSequence(this.ctgNameConverter.apply(contig));
			if(ssr==null) return returnErrorOrNullInterval(JvarkitException.ContigNotFoundInDictionary.getMessage(contig, this._dict));
			contig = ssr.getContig();
			}
		
		if(hyphen!=-1)
			{
			final BigInteger extend;
			start = parseBigInteger(s.substring(colon+1,hyphen));
			
			// chr12:345-678
			if(plus==-1)
				{
				end = parseBigInteger(s.substring(hyphen+1));
				extend = BigInteger.ZERO;
				}
			// chr12:345-678+10
			else
				{
				end = parseBigInteger(s.substring(hyphen+1,plus));
				
				final String extendString = s.substring(plus+1).trim();
				// chr12:345-678+100%
				if(extendString.endsWith("%"))
					{
					final long fragLen =  end.subtract(start).longValueExact();
					final double percent = Double.parseDouble(extendString.substring(0,extendString.length()-1).trim())/100.0;
					if(percent<0) throw new IllegalArgumentException("Negative percent in "+s);
					extend = BigInteger.valueOf(Math.max(0L,(long)(fragLen*percent)));
					}
				// chr12:345-678+100
				else
					{
					extend = parseBigInteger(extendString);
					}
				}
			start = start.subtract(extend);
			if(start.compareTo(BigInteger.ZERO)<0) start=BigInteger.ZERO;
			end = end.add(extend);
			}
		else if(plus!=-1) // chr1:2+1
			{
			final BigInteger mid = parseBigInteger(s.substring(colon+1,plus));
			final BigInteger extend = parseBigInteger(s.substring(plus+1));
			start = mid.subtract(extend);
			if(start.compareTo(BigInteger.ZERO)<0) start=BigInteger.ZERO;
			end = mid.add(extend);
			}
		else
			{
			throw new IllegalArgumentException("boum");
			}
		
		if(this.trimToContigLength && ssr!=null )
			{
			final BigInteger ssrL = new BigInteger(String.valueOf(ssr.getSequenceLength()));
			if(ssrL.compareTo(start)<0) start=ssrL;
			if(ssrL.compareTo(end)<0) end=ssrL;
			}

		
		int istart = BigIntegerToInt(start);
		int iend = BigIntegerToInt(end);

		if(iend<istart) return returnErrorOrNullInterval("end < start in "+s);
		return Optional.of(new SimpleInterval(contig, istart, iend));
	} catch (final Exception err)
		{
		LOG.error(err);
		return returnErrorOrNullInterval("Cannot parse "+s +" : "+err.getMessage());
		}
	}
}
