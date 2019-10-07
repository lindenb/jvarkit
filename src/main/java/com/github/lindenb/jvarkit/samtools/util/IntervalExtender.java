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

*/
package com.github.lindenb.jvarkit.samtools.util;

import java.util.function.Function;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;

/**
 * Utility to extend Intervals
 */
public interface IntervalExtender extends Function<Locatable,SimpleInterval> {
public static final String OPT_DESC="Extending interval.";

public static class StringConverter implements IStringConverter<IntervalExtender> {
	@Override
	public IntervalExtender convert(String arg0) {
		return IntervalExtender.of(arg0);
		}
}
	
/** creates an IntervalExtender of xtend bases */
public static IntervalExtender of(final int xtend) {
	return of(null,xtend);
	}

/** creates an IntervalExtender of xtend bases, validation using a samsequence dictionary */
public static IntervalExtender of(final SAMSequenceDictionary dict,final int xtend) {
	return new ExtendByDistance(dict,xtend);
}

/** creates an IntervalExtender of 'fraction' , validation using a samsequence dictionary */
public static IntervalExtender of(final double fraction) {
	return of(null,fraction);
}

/** creates an IntervalExtender of 'fraction', validation using a samsequence dictionary */
public static IntervalExtender of(final SAMSequenceDictionary dict,final double fraction) {
	return new ExtendByFraction(dict,fraction);
}

/** using a samsequence dictionary */
public static IntervalExtender of(final String decode) {
	return of(null,decode);
}

/** test whether interval will be reduced : fraction is < 1.0  of number of base is <0*/ 
public boolean isShriking();

/** using a samsequence dictionary */
public static IntervalExtender of(final SAMSequenceDictionary dict,String decode) {
	decode= decode.trim();
	if(decode.endsWith("%")) {
		final double fraction = Double.parseDouble(decode.substring(0,decode.length()-1).trim());
		return of(dict,fraction/100.0);
	} else if(decode.contains(".")) {
		final double fraction = Double.parseDouble(decode);
		return of(dict,fraction);
	} else
	{
		final int xtend = new DistanceParser().applyAsInt(decode);
		return of(dict,xtend);
	}
	
}


static abstract class AbstractExtender implements IntervalExtender {
	private final SAMSequenceDictionary dict;
	protected AbstractExtender(final SAMSequenceDictionary dict) {
		this.dict = dict;
		}
	protected SimpleInterval extend(final Locatable t,final int n) {
		final int mid = (int)(((long)t.getStart()+(long)t.getEnd())/2L);
		int beg = mid - n; /* not /2 : headeach if n==1 */
		int end = mid + n;
		// shrinking
		if(beg>end) {
			beg = mid;
			end = mid;
			}
		if(this.dict!=null) {
			final SAMSequenceRecord ssr = this.dict.getSequence(t.getContig());
			if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(t.getContig(), this.dict);
			if(t.getStart()> ssr.getSequenceLength()) {
				throw new IllegalArgumentException("in "+t+". start is greater than contig length "+ssr.getSequenceLength());
				}
			beg = Math.min(ssr.getSequenceLength(), Math.max(1, beg));
			end = Math.min(ssr.getSequenceLength(), Math.max(1, end));
			}
		return new SimpleInterval(t.getContig(), beg, end);
		}
	}


static class ExtendByFraction extends AbstractExtender {
	final double fract;
	ExtendByFraction(final SAMSequenceDictionary dict,double fract) {
		super(dict);
		this.fract = fract;
		}
	@Override
	public boolean isShriking() {
		return fract < 1.0;
		}
	@Override
	public SimpleInterval apply(final Locatable t) {
		final int oldLen = t.getLengthOnReference();
		final int newLen = (int)Math.ceil(oldLen* this.fract);
		return extend(t,newLen - oldLen);
		}
	@Override
	public String toString() {
		return getClass().getName()+" fraction : "+this.fract;
		}
	}

static class ExtendByDistance extends AbstractExtender {
	final int xtend;
	ExtendByDistance(final SAMSequenceDictionary dict,int xtend) {
		super(dict);
		this.xtend = xtend;
		}
	@Override
	public boolean isShriking() {
		return xtend < 0;
		}
	
	@Override
	public SimpleInterval apply(final Locatable t) {
		return extend(t,this.xtend);
		}
	@Override
	public String toString() {
		return getClass().getName()+" bases : "+this.xtend;
		}

	}
}
