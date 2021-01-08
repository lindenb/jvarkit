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
package com.github.lindenb.jvarkit.net;

import java.util.function.Function;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.igv.IgvConstants;

import htsjdk.samtools.util.Locatable;

/** get a genomic URL hyperlink for a given locatable */
public interface Hyperlink extends Function<Locatable,String>{
	public static final String OPT_DESC="creates a hyperlink when 'click' in an area. "
			+ "The URL must contains __CHROM__, __START__ and __END__ that will be replaced by their values. Predefined values are 'hg19','hg38','igv'. "
			+ "IGV : \"http://localhost:60151/goto?locus=__CHROM__%3A__START__-__END__\" , "
			+ "UCSC: \"http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__\" "
			;

	public static class StringConverter implements IStringConverter<Hyperlink> {
		@Override
		public Hyperlink convert(final String s) {
			return Hyperlink.compile(s);
			}
		}
	
	static class HyperkinkImpl implements Hyperlink {
		private final String pattern ;
		private HyperkinkImpl(final String pattern) {
			this.pattern = pattern;
		}
		
		@Override
		public String apply(final Locatable loc) {
			if(loc==null) return null;
			if(isEmpty()) return null;
			return 	getPattern().
					replaceAll("__CHROM__", StringUtils.escapeHttp(loc.getContig())).
					replaceAll("__START__", String.valueOf(loc.getStart())).
					replaceAll("__END__", String.valueOf(loc.getEnd()))
					;
			}
		@Override
		public int hashCode() {
			return pattern.hashCode();
			}
		@Override
		public boolean equals(final Object obj) {
			return obj!=null && obj instanceof Hyperlink && Hyperlink.class.cast(obj).getPattern().equals(this.getPattern());
		}
		@Override
		public String getPattern() {
			return pattern;
		}
		@Override
		public String toString() {
			return getPattern();
			}
	}
	
	/** return the url pattern . never null, may be empty */
	public String getPattern();

	
	/** return true if the pattern doesn't contain any of __CHROM__ __START__ __END__ */
	public default boolean isEmpty() {
		final String s = this.getPattern();
		if(StringUtils.isBlank(s)) return true;
		
		if(!s.contains("__CHROM__"))  return true;
		if(!s.contains("__START__"))  return true;
		if(!s.contains("__END__"))  return true;
		return false;
		}

	
	/** return defaultHyperlink, always return a null url */
	public static Hyperlink empty() {
			return compile("");
		}

	
	public static Hyperlink compile(String s) {
		if(StringUtils.isBlank(s))  new HyperkinkImpl(s);
		if(s.equals("hg19") || s.equals("hg38")) {
			s=  "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db="+s+"&position=__CHROM__%3A__START__-__END__";
			}
		else if(s.equals("igv")) {
			s="http://"+IgvConstants.DEFAULT_HOST+":"+IgvConstants.DEFAULT_PORT+"/goto?locus=__CHROM__%3A__START__-__END__";
			}
		return new HyperkinkImpl(s);
		}
}
