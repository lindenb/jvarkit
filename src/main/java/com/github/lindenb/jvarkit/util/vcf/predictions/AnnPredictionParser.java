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
package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;


/**
 * Parser for
 * http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf


     1	 Allele 
     2	 Annotation 
     3	 Annotation_Impact 
     4	 Gene_Name 
     5	 Gene_ID 
     6	 Feature_Type 
     7	 Feature_ID 
     8	 Transcript_BioType 
     9	 Rank 
    10	 HGVS.c 
    11	 HGVS.p 
    12	 cDNA.pos / cDNA.length 
    13	 CDS.pos / CDS.length 
    14	 AA.pos / AA.length 
    15	 Distance 
    16	 ERRORS / WARNINGS / INFO


 */
public class AnnPredictionParser
	implements PredictionParser
	{
	public static enum Impact 
		{
		// ordered from worst to lighter, keep this order
		 HIGH, MODERATE, MODIFIER,LOW, UNDEFINED
		}
	private static final Logger LOG=Logger.build(AnnPredictionParser.class).make();

	private final CharSplitter pipeRegex= CharSplitter.PIPE;
	private final CharSplitter ampRegex = CharSplitter.of('&');

	private final String tag;
	private final boolean valid;
	private SequenceOntologyTree soTree = SequenceOntologyTree.getInstance();
	
	AnnPredictionParser(final VCFHeader header)
		{		
		this(header,getDefaultTag());
		}
	
	public AnnPredictionParser sequenceOntologyTree( final SequenceOntologyTree soTree) {
		this.soTree = soTree;
		return this;
		}
	
	public static final String getDefaultTag()
		{
		return "ANN";
		}
	
	
	AnnPredictionParser(final VCFHeader header,final String tag)
		{	
		this.tag=(tag==null?getDefaultTag():tag);
		if(header!=null) {
			final VCFInfoHeaderLine info=(header==null?null:header.getInfoHeaderLine(this.tag));
			if(info==null || info.getDescription()==null)
				{
				LOG.warning("no INFO["+tag+"] or no description This VCF was probably NOT annotated with SnpEff(ANN version) . But it's not a problem if this tool doesn't need to access SnpEff Annotations.");
				//this.valid=false;
				//return;
				}
			}
		this.valid=true;
		}

	public String getTag()
		{
		return this.tag;
		}
	public boolean isValid() {
		return valid;
		}

	/** returns all the SO terms String found in this variant */
	public Set<String> getSOTermsStrings(final VariantContext ctx) {
		if(!isValid()) return Collections.emptySet();
		return getPredictions(ctx).stream().
				flatMap(P->P.getSOTermsStrings().stream()).
				collect(Collectors.toSet());
		}

	/** returns all the SO terms found in this variant */
	public Set<SequenceOntologyTree.Term> getSOTerms(final VariantContext ctx) {
		if(!isValid()) return Collections.emptySet();
		return getPredictions(ctx).stream().
				flatMap(P->P.getSOTerms().stream()).
				collect(Collectors.toSet());
		}

	
	public List<AnnPrediction> getPredictions(final VariantContext ctx)
		{
		if(!isValid())
			{
			return Collections.emptyList();
			}
		final List<? extends Object> L= ctx.getAttributeAsList(getTag());
		final ArrayList<AnnPrediction> preds= new ArrayList<AnnPrediction>(L.size());

		for(final Object o2:L)
			{
			final AnnPrediction pred= parseOnePrediction(o2);
			if(pred!=null) preds.add(pred);			
			}
		return preds;
		}
	
	public AnnPrediction  parseOnePrediction(final Object o)
		{
		if(o==null || !this.valid) return null;
		if(!(o instanceof String))
			{
			return parseOnePrediction( o.toString());
			}
		final String s=String.class.cast(o).trim();
		final List<CharSequence> tokens=this.pipeRegex.splitAsCharSequenceList(s);
		return new AnnPrediction(s,tokens);
		}
	
	/*
	private static class AAChange
		{
		String ref;
		int pos;
		String alt;
		}*/
			
	
	public class AnnPrediction
		implements Prediction
		{
		private final String originalStr;
		private final List<CharSequence> _tokens;
		private AnnPrediction(final String originalStr,final List<CharSequence> tokens)
			{
			this.originalStr = originalStr;
			this._tokens=tokens;
			}
		
		private String at(int i)
			{
			return(i<0 || i>=_tokens.size() ? null:_tokens.get(i).toString());
			}
		
		public String getAllele()
			{	
			/* dont user Allele.create because cancer alt don't allow this */
			return at(0);
			}
		
		/*
		@Override
		public String getGeneName()
			{
			return at(3);
			}
		
		public final String getHGNC()
			{
			return getGeneName();
			}
		
		@Override
		public String getEnsemblTranscript()
			{
			return null;
			}
		@Override
		public String getEnsemblProtein()
			{
			return null;
			}
		
		@Override
		public String getEnsemblGene() {
			return null;
			}
		
		private AAChange getAAChange()
			{
			String aa=null;//TODO
			if(aa==null || aa.isEmpty()) return null;
			if(!Character.isLetter(aa.charAt(0))) return null;
			if(!Character.isDigit(aa.charAt(1))) return null;
			AAChange change=new AAChange();
			change.ref=""+aa.charAt(0);
			change.pos=0;
			int i=1;
			while(i< aa.length() && Character.isDigit(aa.charAt(i)))
				{
				change.pos+=(aa.charAt(i)-'0')+change.pos*10;
				++i;
				}
			if(change.pos==0) return null;
			if(i<aa.length() && Character.isLetter(aa.charAt(i)))//occurs when synonymous-lony mutation e.g "M1"
				{
				change.alt=aa.substring(i);
				}
			return change;
			}
		
		@Override
		public String getAltAminoAcid()
			{
			AAChange aa=getAAChange();
			return aa==null?null:aa.alt;
			}
		
		@Override
		public Integer getAminoAcidPosition()
			{
			AAChange aa=getAAChange();
			return aa==null?null:aa.pos;
			}
		@Override
		public String getReferenceAminoAcid() {
			AAChange aa=getAAChange();
			return aa==null?null:aa.ref;
			}
	*/
		/** return true if SO-term-string is and is only equals to 'intergenic_region' */
		public boolean isIntergenicRegion()
			{
			return "intergenic_region".equals(this.getSOTermsString());
			}
		
		public String getSOTermsString() {
			final String so = at(1);
			return so==null?"":so;
			}
		
		public List<String> getSOTermsStrings() {
			final String soterms = getSOTermsString();
			if(StringUtil.isBlank(soterms)) return Collections.emptyList();
			return AnnPredictionParser.this.ampRegex.splitAsStringList(soterms);
			}
		
		//@Override
		public Set<SequenceOntologyTree.Term> getSOTerms()
			{
			final List<String> effects = getSOTermsStrings();
			if(effects.isEmpty()) return Collections.emptySet();
			final Set<SequenceOntologyTree.Term> set=new HashSet<>(effects.size());
			for(final String label:effects) {
				if(label.isEmpty()) continue;
				final SequenceOntologyTree.Term t = AnnPredictionParser.this.soTree.getTermByLabel(label);
				if(t==null) {
					LOG.warning("Current Sequence Ontology Tree doesn't contain \""+ label+"\". May be it's a deprecated term or the current version of this package is obsolete.");
					} 
				else
					{
					set.add(t);
					}
				}
			return set;
			}
		
		public Impact getPutativeImpact()
			{
			final String s=this.at(2);
			if(s==null || s.isEmpty() )return Impact.UNDEFINED;
			return Impact.valueOf(s.toUpperCase().trim());
			}
		
		public String getGeneName()
			{
			return at(3);
			}
		
		public String getGeneId()
			{
			return at(4);
			}
		public String getFeatureType()
			{
			return at(5);
			}
		
		public String getFeatureId()
			{
			return at(6);
			}
		
		public String getTranscriptBioType()
			{
			return this.at(7);
			}
		
		public String getRank()
			{
			return this.at(8);
			}

		public String getHGVSc()
			{
			return this.at(9);
			}
		public String getHGVSp()
			{
			return this.at(10);
			}
		public String getCDNAPos()
			{
			return this.at(11);
			}
		public String getCDSPos()
			{
			return this.at(12);
			}
		public String getAAPos()
			{
			return this.at(13);
			}
		public String getDistance()
			{
			return this.at(14);
			}
		public String getMessages()
			{
			return this.at(15);
			}

		
		private boolean _isEmpty(String s) {
			return s==null || s.trim().isEmpty();
		}
		
		/** returns some the posssible ways to name a gene */
		public Set<String> getGeneKeys() {
			final Set<String> keys= new HashSet<>();
			if(!_isEmpty(this.getFeatureType()) &&
				!_isEmpty(this.getFeatureId()) &&
				this.getFeatureType().equals("transcript")
				)
				{
				keys.add("ANN_FEATURE_TRANSCRIPT_"+this.getFeatureId());
				}
			if(!_isEmpty(this.getGeneName()))
				{
				keys.add("ANN_GENE_"+this.getGeneName());
				}
			if(!_isEmpty(this.getGeneId()))
				{
				keys.add("ANN_GENEID_"+this.getGeneId());
				}
			return keys;
			}
		
		/** return the prediction encoded in the original VariantContext */
		public String getOriginalAttributeAsString()
			{
			return this.originalStr;
			}
		
		@Override
		public String toString() {
			return this.originalStr;
			}
		}
	
	}
