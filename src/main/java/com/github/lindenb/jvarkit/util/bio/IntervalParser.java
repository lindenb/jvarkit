/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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


History:
* 2017 creation

*/
package com.github.lindenb.jvarkit.util.bio;

import java.math.BigInteger;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;

public class IntervalParser {
	private static final Logger LOG = Logger.build(IntervalParser.class).make();
	private static final int MAX_CONTIG_LENGTH = Integer.MAX_VALUE - 1000;
	public static final String OPT_DESC= 
			"An interval as the following syntax : \"chrom:start-end\" or \"chrom:middle+extend\"  or \"chrom:start-end+extend\" or \"chrom:start-end+extend-percent%\"."
			+ "A program might use a Reference sequence to fix the chromosome name (e.g: 1->chr1)" ;
	private SAMSequenceDictionary dict=null;
	private boolean raiseExceptionOnError=true;
	private boolean tryToFixContigName=true;
	private boolean trimToContigLength=false;
	private boolean contigAloneIsWholeContig=false;
	
	public IntervalParser()
		{
		this(null);
		}
	public IntervalParser(final SAMSequenceDictionary dict)
		{
		this.dict = dict;
		}

	/** there is a dictionary and we cannot find the SAMSequenceRecord
	 * try to find an equivalent chromosome e.g: 1->chr1 , chr2 -> 2, MT->chrM
	 * @param tryToFixContigName
	 * @return
	 */
	public IntervalParser setFixContigName(boolean tryToFixContigName) {
		this.tryToFixContigName = tryToFixContigName;
		return this;
		}
	
	public boolean isFixContigName() {
		return tryToFixContigName;
		}
	
	/** if dict is set and interval is larger to contig, adjust to contig's length */
	public void setTrimToContigLength(boolean trimToContig) {
		this.trimToContigLength = trimToContig;
		}
	
	public boolean isTrimToContigLength() {
		return trimToContigLength;
	}
	
	/** specifying 'chr1' would return this whole contig */
	public IntervalParser setContigNameIsWholeContig(final boolean contigAloneIsWholeContig) {
		this.contigAloneIsWholeContig = contigAloneIsWholeContig;
		return this;
		}

	public boolean isContigAloneIsWholeContig() {
		return contigAloneIsWholeContig;
		}
	
	
	public IntervalParser setDictionary(final SAMSequenceDictionary dict) {
		this.dict = dict;
		return this;
		}
	public boolean hasDictionary() {
		return getDictionary()!=null;
	}
	public SAMSequenceDictionary getDictionary() {
		return this.dict;
		}
	public IntervalParser setRaiseExceptionOnError(boolean raiseExceptionOnError) {
		this.raiseExceptionOnError = raiseExceptionOnError;
		return this;
		}
	
	private Interval returnErrorOrNullInterval(final String message) {
		if(raiseExceptionOnError)
			{
			throw new IllegalArgumentException(message);
			}
		return null;
		}
	
	private Set<String> findAlternativeContigNames(String contig)
		{
		if(!this.tryToFixContigName) return Collections.emptySet();
		final Set<String> set = new LinkedHashSet<>();
		if(contig.startsWith("chr"))
			{
			set.add(contig.substring(3));
			}
		else if(!contig.startsWith("chr"))
			{
			set.add("chr"+contig);
			}
			
		if(contig.equals("chrM")) set.add("MT");
		if(contig.equals("MT")) set.add("chrM");
		
		return set;
		}
	
	private SAMSequenceRecord getSAMSequenceRecord(final String contig)
		{
		if(!hasDictionary()) return null;
		SAMSequenceRecord rec= getDictionary().getSequence(contig);
		if(rec!=null) return rec;
		//case problem ?
		if(this.isFixContigName() )
			{
			for(final SAMSequenceRecord ssr: getDictionary().getSequences())
				{
				if(ssr.getSequenceName().equalsIgnoreCase(contig)) return ssr;
				}
			}
		for(final String c: findAlternativeContigNames(contig)) {
			rec= getDictionary().getSequence(c);
			if(rec!=null) return rec;
			}
		return null;
		}
	private BigInteger parseBigInteger( String s)
		{
		BigInteger factor = BigInteger.ONE;
		s=s.replace(",", "");
		
		
		 if(s.endsWith("bp"))
			{
			s=s.substring(0, s.length()-2).trim();
			}
		else if(s.toLowerCase().endsWith("kb"))
			{
			s=s.substring(0, s.length()-2).trim();
			factor = new BigInteger("1000");
			}
		else if(s.toLowerCase().endsWith("mb"))
			{
			s=s.substring(0, s.length()-2).trim();
			factor = new BigInteger("1000000");
			}
		else if(s.endsWith("b"))
			{
			s=s.substring(0, s.length()-1).trim();
			}
		else if(s.endsWith("k"))
			{
			s=s.substring(0, s.length()-1).trim();
			factor = new BigInteger("1000");
			}
		else if(s.endsWith("m"))
			{
			s=s.substring(0, s.length()-1).trim();
			factor = new BigInteger("1000000");
			}
		final BigInteger vbi = new BigInteger(String.valueOf(s)).multiply(factor);
		
		return vbi;
		}
	
	private int BigIntegerToInt(final BigInteger vbi)
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
	public Interval parse(final String s)
		{
		final int colon=s.indexOf(':');
		if(colon<1 || colon+1==s.length()) {
			/** chromosome alone */
			if(colon==-1 && isContigAloneIsWholeContig())
				{
				String chrom=s;
				if(hasDictionary())
					{
					final SAMSequenceRecord ssr= getSAMSequenceRecord(chrom);
					if(ssr==null)
						{
						return returnErrorOrNullInterval(
							"Cannot find chromosome \""+chrom+"\" in dictionary. Available chromosomes are : "+
							this.dict.getSequences().stream().
							map(S->"\""+S.getSequenceName()+"\"").
							collect(Collectors.joining(", ")));
						}
					return new Interval(ssr.getSequenceName(),1,ssr.getSequenceLength());
					}
				else
					{
					return new Interval(chrom, 1, MAX_CONTIG_LENGTH);
					}
				}
			
			return returnErrorOrNullInterval("Cannot find colon in "+s);
			}
		final int hyphen = s.indexOf('-',colon+1);
		final int plus = s.indexOf('+',colon+1);
		if(hyphen==-1 && plus==-1) return returnErrorOrNullInterval("Cannot find hyphen or plus in "+s);
		if(hyphen!=-1 && plus!=-1)
			{
			// chr1:123-+
			if(plus <= hyphen+1) return returnErrorOrNullInterval("both hyphen and plus in "+s);
			}
			
		BigInteger start=null,end=null;
		try {
			SAMSequenceRecord ssr=null;
			String chrom=s.substring(0,colon);
			if(hasDictionary())
				{
				ssr= getSAMSequenceRecord(chrom);
				if(ssr==null)
					{
					return returnErrorOrNullInterval(
						"Cannot find chromosome \""+chrom+"\" in dictionary. Available chromosomes are : "+
						this.dict.getSequences().stream().map(S->"\""+S.getSequenceName()+"\"").collect(Collectors.joining(", ")));
					}
				chrom = ssr.getSequenceName();
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
			else if(plus!=-1)
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
			
			if(isTrimToContigLength() && ssr!=null )
			if(isTrimToContigLength() && ssr!=null )
				{
				final BigInteger ssrL = new BigInteger(String.valueOf(ssr.getSequenceLength()));
				if(ssrL.compareTo(start)<0) start=ssrL;
				if(ssrL.compareTo(end)<0) end=ssrL;
				}

			
			int istart = BigIntegerToInt(start);
			int iend = BigIntegerToInt(end);

			if(iend<istart) return returnErrorOrNullInterval("end < start in "+s);
			return new Interval(chrom, istart, iend,false,s);
		} catch (final Exception err)
			{
			LOG.error(err);
			return returnErrorOrNullInterval("Cannot parse "+s +" : "+err.getMessage());
			}
		}
	
	}
