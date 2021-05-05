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
package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.CharSplitter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/** see https://pcingola.github.io/SnpEff/adds/snpEff_lof_nmd.pdf */
public class SnpEffLofNmdParser {
public static final String LOF_TAG="LOF"; 
public static final String NMD_TAG="NMD"; 
private final boolean valid;
private final String tag;

public static interface Prediction {
	public String getGeneName();
	public String getGeneId();
	public int getNumberOfTranscripts();
	public float getPercentOfTranscriptsAffected();
	/** return NMD or LOF */
	public String getTag();
	}

private class PredictionImpl implements Prediction {
	String geneName;
	String geneId;
	int numberOfTranscripts;
	float percentOfTranscriptsAffected;
	@Override
	public String getGeneName() {
		return geneName;
		}
	@Override
	public String getGeneId() {
		return geneId;
		}
	@Override
	public int getNumberOfTranscripts() {
		return numberOfTranscripts;
		}
	@Override
	public float getPercentOfTranscriptsAffected() {
		return percentOfTranscriptsAffected;
		};
	@Override
	public int hashCode() {
		int i= geneName.hashCode();
		i = i*31 + geneId.hashCode();
		i= i*31 + Integer.hashCode(numberOfTranscripts);
		i= i*31 + Float.hashCode(percentOfTranscriptsAffected);
		return i;
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this)return true;
		if(obj==null || !(obj instanceof PredictionImpl)) return false;
		PredictionImpl o = PredictionImpl.class.cast(obj);
		return this.geneName.equals(o.geneName) && 
				this.geneId.equals(o.geneId) && 
				this.numberOfTranscripts == o.numberOfTranscripts &&
				this.percentOfTranscriptsAffected == o.percentOfTranscriptsAffected &&
				this.getTag().equals(o.getTag())
				;
		}
	@Override
	public String getTag() {
		return SnpEffLofNmdParser.this.getTag();
	}
	@Override
	public String toString() {
		return getTag()+"="+String.join("|",geneName,geneId,String.valueOf(numberOfTranscripts),String.valueOf(percentOfTranscriptsAffected));
		}
	}


public SnpEffLofNmdParser(final String tag,final VCFHeader header) {
	if(!tag.equals(NMD_TAG) && !tag.equals(LOF_TAG)) {
		throw new IllegalArgumentException("tag should be "+NMD_TAG+" or "+ LOF_TAG+" got "+tag);
		}
	this.tag = tag;
	final VCFInfoHeaderLine h = header.getInfoHeaderLine(getTag());
	String[] tokens = null;
	if(h!=null) {
		String s= h.getDescription();
		final String t = "Format:";
		int i = s.indexOf(t);
		if(i!=-1)
			{
			s = s.substring(i+t.length()).trim(); 
			if(s.startsWith("'") && s.endsWith("'")) {
				s= s.substring(1,s.length()-1).trim();
				}
			tokens = CharSplitter.PIPE.split(s);
			for(i=0;i< tokens.length;i++) {
				tokens[i]=tokens[i].trim();
				}
			}
		}
	//Gene_Name | Gene_ID | Number_of_transcripts_in_gene | 
	this.valid= tokens!=null && tokens.length ==4 &&
			tokens[0].equals("Gene_Name") &&
			tokens[1].equals("Gene_ID") &&
			tokens[2].equals("Number_of_transcripts_in_gene") &&
			tokens[3].equals("Percent_of_transcripts_affected")
			;
	}

public boolean isValid() {
	return valid;
	}

public List<Prediction> parse(final VariantContext ctx) {
	if(!isValid() || ctx==null || !ctx.hasAttribute(getTag())) return Collections.emptyList();
	return ctx.getAttributeAsStringList(getTag(), null).
			stream().
			map(S->parseOne(S).orElse(null)).
			filter(S->S!=null).
			collect(Collectors.toList())
			;
	}

public Optional<Prediction> parseOne(String s) {
	if(!isValid()) return Optional.empty();
	if(s.startsWith("(") && s.endsWith(")")) {
		s= s.substring(1,s.length()-1);
		}
	final String[] tokens = CharSplitter.PIPE.split(s);
	if(tokens.length!=4) return Optional.empty();
	final PredictionImpl lof = new PredictionImpl();
	lof.geneName = tokens[0];
	lof.geneId = tokens[1];
	lof.numberOfTranscripts = Integer.parseInt(tokens[2]);
	lof.percentOfTranscriptsAffected = Float.parseFloat(tokens[3]);
	return Optional.of(lof);
	}

public String getTag() {
	return this.tag;
	}
public static SnpEffLofNmdParser createLofParser(final VCFHeader header) {
	return new SnpEffLofNmdParser(LOF_TAG,header);
	}
public static SnpEffLofNmdParser createNmdParser(final VCFHeader header) {
	return new SnpEffLofNmdParser(NMD_TAG,header);
	}

}
