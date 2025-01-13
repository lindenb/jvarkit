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
package com.github.lindenb.jvarkit.net;

import java.util.Optional;
import java.util.function.Function;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.igv.IgvConstants;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;

/** get a genomic URL hyperlink for a given locatable */
public interface Hyperlink extends Function<Locatable,Optional<String>> {
	static final String UCSC_PATTERN_HG19 = "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__";
	static final String UCSC_PATTERN = UCSC_PATTERN_HG19;
	public static final String OPT_DESC="creates a hyperlink when 'click' in an area. "
			+ "The URL must contains __CHROM__, __START__ and __END__ that will be replaced by their values. Predefined values are 'hg19','hg38','igv'. "
			+ "IGV : \"http://localhost:60151/goto?locus=__CHROM__%3A__START__-__END__\" , "
			+ "UCSC: \""+ UCSC_PATTERN+"\" "
			;

	public static class StringConverter implements IStringConverter<Hyperlink> {
		@Override
		public Hyperlink convert(final String s) {
			return Hyperlink.compile(s);
			}
		}
	
	
		static class HyperkinkImpl implements Hyperlink {
			private final String pattern ;
			Function<String, String> convertContig = S->S;;
			private HyperkinkImpl(final String pattern) {
				this.pattern = pattern;
			}
		
		@Override
		public Optional<String> apply(final Locatable loc) {
			if(loc==null) return Optional.empty();
			if(isEmpty()) return Optional.empty();
			return 	Optional.of(getPattern().
					replaceAll("__CHROM__", StringUtils.escapeHttp(convertContig.apply(loc.getContig()))).
					replaceAll("__START__", String.valueOf(loc.getStart())).
					replaceAll("__END__", String.valueOf(loc.getEnd()))
					);
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
		Function<String,String> convert = S->S;
		if(s.equals("hg19")) {
			s=  UCSC_PATTERN_HG19;
			convert = S->(S.startsWith("chr")?S:"chr"+S);
			}
		else if(s.equals("hg38")) {
			s=  UCSC_PATTERN_HG19.replace("hg19", "hg38");;
			convert = S->(S.startsWith("chr")?S:"chr"+S);
			}
		else if(s.equals("igv")) {
			s="http://"+IgvConstants.DEFAULT_HOST+":"+IgvConstants.DEFAULT_PORT+"/goto?locus=__CHROM__%3A__START__-__END__";
			}
		/* for testing purpose only */
		else if(s.equals("rotavirus_rf")) {
			s="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10933&locus=__CHROM__%3A__START__-__END__";
			}
		final HyperkinkImpl h= new HyperkinkImpl(s);
		h.convertContig = convert;
		return h;
		}
	
	/** provide a generic hyperlink for the given dictionary */
	public static Hyperlink compile(final SAMSequenceDictionary dict) {
		if(dict==null || dict.isEmpty()) return empty();
		if(SequenceDictionaryUtils.isGRCh37(dict)) return compile("hg19");
		if(SequenceDictionaryUtils.isGRCh38(dict)) return compile("hg38");
		if(SequenceDictionaryUtils.isRotavirusRF(dict)) return compile("rotavirus_rf");
		return empty();
		}

	
}
