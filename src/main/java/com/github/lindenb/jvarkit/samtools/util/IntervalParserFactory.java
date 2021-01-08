/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.Maker;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;


public interface IntervalParserFactory extends Maker<Function<String,Optional<SimpleInterval>>> {
static final Logger LOG = Logger.build(IntervalParserFactory.class).make();

public static final String OPT_DESC= 
			"An interval as the following syntax : \"chrom:start-end\" or \"chrom:middle+extend\"  or \"chrom:start-end+extend\" or \"chrom:start-end+extend-percent%\"."
			+ "A program might use a Reference sequence to fix the chromosome name (e.g: 1->chr1)" ;
	
public static IntervalParserFactory newInstance() {
	return new IntervalParserFactoryImpl();
	}

public static IntervalParserFactory newInstance(final SAMSequenceDictionary dict) {
	return newInstance().dictionary(dict);
	}

public static Supplier<IllegalArgumentException> exception() {
	return exception(null);
}

public IntervalParserFactory dictionary(final SAMSequenceDictionary dic);
public IntervalParserFactory enableWholeContig();
public IntervalParserFactory throwOnError();
/** enable single point position like chr1:1334 */
public default IntervalParserFactory enableSinglePoint() {
	return enableSinglePoint(true);
	}
/** enable single point position like chr1:1334 */
public IntervalParserFactory enableSinglePoint(boolean b);


public static Supplier<IllegalArgumentException> exception(final String str) {
	return ()->new IllegalArgumentException("Cannot parse interval"+
			(StringUtils.isBlank(str)?"":" \""+str+"\"")
			);
}


class IntervalParserFactoryImpl implements IntervalParserFactory {
	
	private class IntervalParser implements Function<String, Optional<SimpleInterval>> {
		SAMSequenceDictionary dict;
		boolean enableWholeContig=false;
		boolean raiseExceptionOnError = false;
		final boolean tryToFixContigName = true;
		boolean trimToContigLength = false;
		private final int MAX_CONTIG_LENGTH = Integer.MAX_VALUE - 1000;
		boolean enable_single_point = false;

		ContigNameConverter ctgNameConverter = null;
		@Override
		public Optional<SimpleInterval> apply(String s) {
			return this.parse(s);
			}

		
		private boolean hasDictionary() {
			return getDictionary()!=null;
			}
		private SAMSequenceDictionary getDictionary() {
			return this.dict;
			}
		
		private <T> Optional<T> returnErrorOrNullInterval(final String message) {
			if(raiseExceptionOnError)
				{
				throw new IllegalArgumentException(message);
				}
			return Optional.empty();
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
			
			//use ctgName converter
			if(this.ctgNameConverter!=null) {
				final String ctg = this.ctgNameConverter.apply(contig);
				if(!StringUtils.isBlank(ctg))
					{
					rec= getDictionary().getSequence(contig);
					if(rec!=null) return rec;
					}
				}
			
			//case problem ?
			
			for(final SAMSequenceRecord ssr: getDictionary().getSequences())
				{
				if(ssr.getSequenceName().equalsIgnoreCase(contig)) return ssr;
				}
				
			for(final String c: findAlternativeContigNames(contig)) {
				rec= getDictionary().getSequence(c);
				if(rec!=null) return rec;
				}
			return null;
			}
		private BigInteger parseBigInteger( String s)
			{
			final DistanceParser distanceParser=new DistanceParser();
			return distanceParser.parseAsBigInteger(s);
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
		
		private Optional<String> parseChrom(final String chrom,final String context) {
			if(hasDictionary())
				{
				final SAMSequenceRecord ssr= getSAMSequenceRecord(chrom);
				if(ssr==null)
					{
					return returnErrorOrNullInterval(
						"Cannot find chromosome \""+chrom+"\" from \""+ context +"\" in dictionary. Available chromosomes are : "+
						this.dict.getSequences().stream().
						map(S->"\""+S.getSequenceName()+"\"").
						collect(Collectors.joining(", ")));
					}
				return Optional.of(ssr.getSequenceName());
				}
			else if(StringUtils.isBlank(chrom)) {
				return returnErrorOrNullInterval("empty chromosome in \""+context+"\".");
				}
			else
				{
				return Optional.of(chrom);
				}
			}
		
		private Optional<SimpleInterval> parse(final String s)
			{
			final int colon=s.lastIndexOf(':');
			if(colon<1 || colon+1==s.length()) {
				/* chromosome alone */
				if(colon==-1 && this.enableWholeContig)
					{
					final Optional<String> optChrom = parseChrom(s,s);
					if(!optChrom.isPresent()) return Optional.empty();
					if(hasDictionary())
						{
						final SAMSequenceRecord ssr= getSAMSequenceRecord(optChrom.get());
						if(ssr==null) throw new IllegalStateException("cannot get chrom in "+s);
						return Optional.of(new SimpleInterval(ssr));
						}
					else
						{
						return Optional.of(new SimpleInterval(optChrom.get(), 1, MAX_CONTIG_LENGTH));
						}
					}
				
				return returnErrorOrNullInterval("Cannot find colon in "+s);
				}
			final int hyphen = s.indexOf('-',colon+1);
			final int plus = s.indexOf('+',colon+1);
			if(hyphen==-1 && plus==-1) {
				if(!this.enable_single_point) return returnErrorOrNullInterval("Cannot find hyphen or plus in "+s+". Single point position like 'chr1:234' are not allowed");
				/** single point mutation */
				final Optional<String> optChrom = parseChrom(s.substring(0,colon).trim(),s);
				if(!optChrom.isPresent()) return Optional.empty();
				final BigInteger bPos = parseBigInteger(s.substring(colon+1).trim());
				final int pos= BigIntegerToInt(bPos);
				return Optional.of(new SimpleInterval(optChrom.get(), pos, pos)); 
				}
			if(hyphen!=-1 && plus!=-1)
				{
				// chr1:123-+
				if(plus <= hyphen+1) return returnErrorOrNullInterval("both hyphen and plus in "+s);
				}
				
			BigInteger start=null,end=null;
			try {
				
				final Optional<String> optChrom = parseChrom(s.substring(0,colon).trim(),s);
				if(!optChrom.isPresent()) return Optional.empty();
				final String chrom = optChrom.get();
				final SAMSequenceRecord ssr=getSAMSequenceRecord(chrom);
				
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
				
				if(this.trimToContigLength && ssr!=null )
					{
					final BigInteger ssrL = new BigInteger(String.valueOf(ssr.getSequenceLength()));
					if(ssrL.compareTo(start)<0) start=ssrL;
					if(ssrL.compareTo(end)<0) end=ssrL;
					}

				
				int istart = BigIntegerToInt(start);
				int iend = BigIntegerToInt(end);

				if(iend<istart) return returnErrorOrNullInterval("end < start in "+s);
				return Optional.of(new SimpleInterval(chrom, istart, iend));
			} catch (final Exception err)
				{
				LOG.error(err);
				return returnErrorOrNullInterval("Cannot parse "+s +" : "+err.getMessage());
				}
			}
		
		}
	private SAMSequenceDictionary dict;
	private boolean enableWholeContig=false;
	private boolean raiseExceptionOnError = false;
	private boolean trimToContigLength = false;
	private boolean enable_single_point = false;
	@Override
	public IntervalParserFactory dictionary(final SAMSequenceDictionary dict) {
		this.dict = dict;
		return this;
		}
	@Override
	public Function<String, Optional<SimpleInterval>> make() {
		final IntervalParser parser = new IntervalParser();
		parser.dict = this.dict;
		if(parser.dict!=null) {
			parser.ctgNameConverter=ContigNameConverter.fromOneDictionary(dict);
			}
		parser.enableWholeContig = this.enableWholeContig;
		parser.raiseExceptionOnError = this.raiseExceptionOnError;
		parser.trimToContigLength = this.trimToContigLength;
		parser.enable_single_point = this.enable_single_point;
		return parser;
		}
	@Override
	public IntervalParserFactory enableWholeContig() {
		this.enableWholeContig=true;
		return this;
		}
	@Override
	public IntervalParserFactory throwOnError() {
		this.raiseExceptionOnError = true;
		return this;
		}
	@Override
	public IntervalParserFactory enableSinglePoint(boolean b) {
		this.enable_single_point = b;
		return this;
		}
	}
}
