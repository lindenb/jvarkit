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
* 2014 creation

*/
package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;


/**
 * ##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon  | GenotypeNum [ | ERRORS | WARNINGS ] )'">

 * @author lindenb
 *
 */
public class SnpEffPredictionParser implements PredictionParser
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");

	private enum COLS{ Effect , Effect_Impact , Functional_Class,Codon_Change, 
		Amino_Acid_change,Amino_Acid_length,Gene_Name , Gene_BioType , Coding , Transcript,
		Exon  , GenotypeNum , ERRORS , WARNINGS,Transcript_BioType,Gene_Coding,Transcript_ID,Exon_Rank,Genotype_Number};
	private Map<COLS, Integer> col2col=new HashMap<COLS, Integer>();
	private Pattern pipe=Pattern.compile("[\\|\\(\\)]");
	private String tag;
	public SnpEffPredictionParser(VCFHeader header)
		{		
		this(header,getDefaultTag());
		}
	
	public static final String getDefaultTag()
		{
		return "EFF";
		}
	
	
	public SnpEffPredictionParser(VCFHeader header,String tag)
		{	
		this.tag=(tag==null?getDefaultTag():tag);
		VCFInfoHeaderLine info=header.getInfoHeaderLine(tag);
		if(info==null || info.getDescription()==null)
			{
			LOG.warning("no INFO["+tag+"] or no description ");
			return;
			}
		String description=info.getDescription();
		String chunck="Format:";
		int i=description.indexOf(chunck);
		if(i==-1)
			{
			LOG.warning("Cannot find "+chunck+ " in "+description);

			return;
			}
		description=description.substring(i+chunck.length()).replace('(','|').replaceAll("[ \'\\.)\\[\\]]+","").trim();
		String tokens[]=pipe.split(description);
		for(i=0;i< tokens.length;++i)
			{
			if(tokens[i].isEmpty()) continue;
			COLS col=null;
			for(COLS c:COLS.values())
				{
				if(c.name().equalsIgnoreCase(tokens[i]))
					{
					col=c;
					}
				}
			if(col==null)
				{
				LOG.warning("Undefined SnpEff tag \""+tokens[i]+"\"");
				continue;
				}
			col2col.put(col, i);
			}
		}
	
	@Override
	public String getTag()
		{
		return this.tag;
		}
	

	@Override
	public List<SnpEffPrediction> getPredictions(VariantContext ctx)
		{
		 ArrayList<SnpEffPrediction> preds= new ArrayList<SnpEffPrediction>();
		if(col2col.isEmpty())
			{
			return preds;
			}
		Object o=ctx.getAttribute(getTag());
		List<? extends Object> L= VCFUtils.attributeAsList(o);
		for(Object o2:L)
			{
			 _predictions(preds,o2);
			}
		return preds;
		}
	
	private void _predictions( List<SnpEffPrediction> preds,Object o)
		{
		SnpEffPrediction pred= parseOnePrediction(o);
		if(pred==null) return;
		preds.add(pred);
		}
	
	public SnpEffPrediction  parseOnePrediction(Object o)
		{
		if(o==null) return null;
		if(!(o instanceof String))
			{
			return parseOnePrediction( o.toString());
			}
		String s=String.class.cast(o).trim();
		String tokens[]=pipe.split(s);
		return new SnpEffPrediction(tokens);
		}
	
	
	private static class AAChange
		{
		String ref;
		int pos;
		String alt;
		}
			
	
	public class SnpEffPrediction
		implements Prediction
		{
		private String tokens[];
		SnpEffPrediction(String tokens[])
			{
			this.tokens=tokens;
			}
		private String getByCol(COLS col)
			{
			Integer idx=col2col.get(col);
			if(idx==null || idx>=tokens.length || tokens[idx].isEmpty()) return null;
			return tokens[idx];
			}
		@Override
		public String getGeneName()
			{
			return getByCol(COLS.Gene_Name);
			}
		
		@Override
		public String getEnsemblTranscript()
			{
			String s=getByCol(COLS.Transcript);
			if(s==null || !s.startsWith("ENST")) return null;
			return s;
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
			String aa=getByCol(COLS.Amino_Acid_change);
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

		
		private Map<COLS,String> getMap()
			{
			Map<COLS, String> hash=new HashMap<COLS,String>();
			for(COLS c: col2col.keySet())
				{
				int idx=col2col.get(c);
				if(idx>=this.tokens.length) continue;
				hash.put(c, tokens[idx]);
				}
			return hash;
			}
		
		@Override
		public Set<SequenceOntologyTree.Term> getSOTerms()
			{
			Set<SequenceOntologyTree.Term> set=new HashSet<SequenceOntologyTree.Term>();
			String EFF=getByCol(COLS.Effect);
			if(EFF==null) return set;
			for(SequenceOntologyTree.Term t:SequenceOntologyTree.getInstance().getTerms())
				{
				if(t.getLabel().equals(EFF))
					{
					set.add(t);
					//break ?
					}
				}
			return set;
			}


		
		@Override
		public String toString() {
			return getMap().toString()+ " "+Arrays.asList(tokens);
			}
		}
	
	}
