
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
package com.github.lindenb.jvarkit.date;

import java.time.Duration;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;

/** simple duration parser */
public class DurationParser {
	public static final String OPT_DESC="format: <integer>(years|week|days|hours|minutes|seconds)";
	
	private static class SuffixFun {
		final String suffix;
		final int factor ;
		final Function<Long,Duration> fun;
		SuffixFun(final String suffix,int factor,Function<Long,Duration> fun) {
			this.suffix = suffix;
			this.factor = factor;
			this.fun = fun;
		}
	}
	
	public Duration convert(String s) {
		if(StringUtils.isBlank(s)) throw new IllegalArgumentException("cannot convert blank string to duration");
		final List<SuffixFun> suffixes = new ArrayList<>();
		suffixes.add(new SuffixFun("years",365,Duration::ofDays));
		suffixes.add(new SuffixFun("year",365,Duration::ofDays));
		suffixes.add(new SuffixFun("y",365,Duration::ofDays));
		suffixes.add(new SuffixFun("weeks",7,Duration::ofDays));
		suffixes.add(new SuffixFun("week",7,Duration::ofDays));
		suffixes.add(new SuffixFun("w",7,Duration::ofDays));
		suffixes.add(new SuffixFun("days",1,Duration::ofDays));
		suffixes.add(new SuffixFun("day",1,Duration::ofDays));
		suffixes.add(new SuffixFun("d",1,Duration::ofDays));
		suffixes.add(new SuffixFun("hours",1,Duration::ofHours));
		suffixes.add(new SuffixFun("hour",1,Duration::ofHours));
		suffixes.add(new SuffixFun("h",1,Duration::ofHours));
		suffixes.add(new SuffixFun("minutes",1,Duration::ofMinutes));
		suffixes.add(new SuffixFun("minute",1,Duration::ofMinutes));
		suffixes.add(new SuffixFun("min",1,Duration::ofMinutes));
		suffixes.add(new SuffixFun("m",1,Duration::ofMinutes));
		suffixes.add(new SuffixFun("secondes",1,Duration::ofSeconds));
		suffixes.add(new SuffixFun("seconde",1,Duration::ofSeconds));
		suffixes.add(new SuffixFun("secs",1,Duration::ofSeconds));
		suffixes.add(new SuffixFun("sec",1,Duration::ofSeconds));
		suffixes.add(new SuffixFun("s",1,Duration::ofSeconds));
		suffixes.add(new SuffixFun("milliseconds",1,Duration::ofMillis));
		suffixes.add(new SuffixFun("millisecond",1,Duration::ofMillis));
		suffixes.add(new SuffixFun("millisec",1,Duration::ofMillis));
		suffixes.add(new SuffixFun("ms",1,Duration::ofMillis));

		
		s=s.toLowerCase().trim();
		for(SuffixFun sf:suffixes) {
			if(!s.endsWith(sf.suffix)) continue;
			s = s.substring(0, s.length()-sf.suffix.length()).trim();
			final double d =  Double.parseDouble(s);
			if(d< 0) throw new IllegalArgumentException("cannot parse negative duration");
			return sf.fun.apply((long)(sf.factor *d));
			}
		
		throw new IllegalArgumentException("illegal suffix in "+s+" should be one of "+suffixes.stream().
			map(S->S.suffix).
			collect(Collectors.joining(" ")));
		}

}
