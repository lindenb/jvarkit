/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import java.util.logging.Logger;
import java.util.regex.Pattern;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;


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
	{
	public static enum PutativeImpact 
		{
		UNDEFINED,HIGH,MODERATE,LOW,MODIFIER
		}
	private static final Logger LOG=Logger.getLogger("jvarkit");

	private final Pattern pipe=Pattern.compile("[\\|\\(\\)]");
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
	
	
	public AnnPredictionParser(final VCFHeader header,final String tag)
		{	
		this.tag=(tag==null?getDefaultTag():tag);
		final VCFInfoHeaderLine info=header.getInfoHeaderLine(tag);
		if(info==null || info.getDescription()==null)
			{
			LOG.warning("no INFO["+tag+"] or no description ");
			this.valid=false;
			return;
			}
		this.valid=true;
		}

	public String getTag()
		{
		return this.tag;
		}
	


	public List<AnnPrediction> getPredictions(final VariantContext ctx)
		{
		if(!this.valid)
			{
			return Collections.emptyList();
			}
		final Object o=ctx.getAttribute(getTag());
		final List<? extends Object> L= VCFUtils.attributeAsList(o);
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
		final String tokens[]=pipe.split(s);
		return new AnnPrediction(tokens);
		}
	
	/*
	private static class AAChange
		{
		String ref;
		int pos;
		String alt;
		}*/
			
	
	public class AnnPrediction
		//implements Prediction
		{
		private String tokens[];
		private AnnPrediction(final String tokens[])
			{
			this.tokens=tokens;
			}
		
		private String at(int i)
			{
			return(i<0 || i>=tokens.length ? null:tokens[i]);
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
		
		
		//@Override
		public Set<SequenceOntologyTree.Term> getSOTerms()
			{
			Set<SequenceOntologyTree.Term> set=new HashSet<SequenceOntologyTree.Term>();
			if(this.tokens.length<2) return set;
			String EFF=this.tokens[1];
			if(EFF==null) return set;
			for(String eff:EFF.split("[&]"))
				{
				if(eff.isEmpty()) continue;
				for(final SequenceOntologyTree.Term t: AnnPredictionParser.this.soTree.getTerms())
					{
					if(t.getLabel().equals(eff))
						{
						set.add(t);
						//break ?
						}
					}
				}
			return set;
			}
		
		public PutativeImpact getPutativeImpact()
			{
			if(this.tokens.length<3) return PutativeImpact.UNDEFINED;
			String s=this.tokens[2];
			return PutativeImpact.valueOf(s.toUpperCase().trim());
			}
		public String getGeneId()
			{
			if(this.tokens.length<5) return null;
			return this.tokens[4];
			}
		public String getFeatureType()
			{
			if(this.tokens.length<6) return null;
			return this.tokens[5];
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

		
		
		@Override
		public String toString() {
			return Arrays.toString(this.tokens);
			}
		}
	
	}
