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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;


/**
 * ##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon  | GenotypeNum [ | ERRORS | WARNINGS ] )'">

 * @author lindenb
 *
 */
public class SnpEffPredictionParser implements PredictionParser
	{
	private static final Logger LOG=Logger.build(SnpEffPredictionParser.class).make();

	private final Map<String, Integer> col2col=new HashMap<String, Integer>();
	private final Pattern pipe=Pattern.compile("[\\|\\(\\)]");
	private String tag;
	private SequenceOntologyTree soTree = SequenceOntologyTree.getInstance();
	private final boolean valid;
	
	SnpEffPredictionParser(final VCFHeader header)
		{		
		this(header,getDefaultTag());
		}
	
	public SnpEffPredictionParser sequenceOntologyTree( final SequenceOntologyTree soTree) {
		this.soTree = soTree;
		return this;
		}

	
	public static final String getDefaultTag()
		{
		return "EFF";
		}
	
	
	public SnpEffPredictionParser(final VCFHeader header,final String tag)
		{	
		this.tag=(tag==null?getDefaultTag():tag);
		final VCFInfoHeaderLine info=(header==null?null:header.getInfoHeaderLine(tag));
		if(info==null || info.getDescription()==null)
			{
			LOG.warning("no INFO["+tag+"] or no description. This VCF was probably NOT annotated with SnpEff (old version). But it's not a problem if this tool doesn't need to access SnpEff Annotations. ");
			this.valid=false;
			return;
			}
		String description=info.getDescription();
		String chunck="Format:";
		int i=description.indexOf(chunck);
		if(i==-1)
			{
			LOG.warning("Cannot find "+chunck+ " in "+description);
			this.valid=false;
			return;
			}
		description=description.substring(i+chunck.length()).replace('(','|').replaceAll("[ \'\\.)\\[\\]]+","").trim();
		String tokens[]=pipe.split(description);
		for(i=0;i< tokens.length;++i)
			{
			final String col=tokens[i];
			if(col.isEmpty()) continue;
			if(this.col2col.containsKey(col))
				{
				LOG.warn("duplicate SnpEff attribute:"+col);
				continue;
				}
			this.col2col.put(col, i);
			}
		this.valid=true;
		}
	
	@Override
	public String getTag()
		{
		return this.tag;
		}
	
	public boolean isValid() {
		return valid;
	}
	
	@Override
	public List<SnpEffPrediction> getPredictions(final VariantContext ctx)
		{
		if(!this.isValid() || this.col2col.isEmpty() || !ctx.hasAttribute(getTag()))
			{
			return Collections.emptyList();
			}
		final List<? extends Object> L= ctx.getAttributeAsList(getTag());
		final ArrayList<SnpEffPrediction> preds= new ArrayList<SnpEffPrediction>(L.size());
		for(final Object o2:L)
			{
			 _predictions(preds,o2);
			}
		return preds;
		}
	
	private void _predictions(final List<SnpEffPrediction> preds,final Object o)
		{
		final SnpEffPrediction pred= parseOnePrediction(o);
		if(pred==null) return;
		preds.add(pred);
		}
	
	public SnpEffPrediction  parseOnePrediction(final Object o)
		{
		if(!this.isValid() || this.col2col.isEmpty())
			{
			return null;
			}
		if(o==null) return null;
		if(!(o instanceof String))
			{
			return parseOnePrediction( o.toString());
			}
		final String tokens[]=pipe.split(String.class.cast(o).trim());
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
		private final String tokens[];
		SnpEffPrediction(final String tokens[])
			{
			this.tokens=tokens;
			}
		/** get column by name, may return null. Returns null if column is empty */
		private String getByCol(String col)
			{
			final Integer idx=col2col.get(col);
			if(idx==null || idx>=this.tokens.length || tokens[idx].isEmpty()) return null;
			return tokens[idx];
			}
		public String getGeneName()
			{
			return getByCol("Gene_Name");
			}
		
		public String getEnsemblTranscript()
			{
			final String s=getByCol("Transcript");
			if(s==null || !s.startsWith("ENST")) return null;
			return s;
			}
		
		
		private AAChange getAAChange()
			{
			String aa=getByCol("Amino_Acid_change");
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
		
		public String getAltAminoAcid()
			{
			AAChange aa=getAAChange();
			return aa==null?null:aa.alt;
			}
		
		public Integer getAminoAcidPosition()
			{
			AAChange aa=getAAChange();
			return aa==null?null:aa.pos;
			}
		public String getReferenceAminoAcid() {
			AAChange aa=getAAChange();
			return aa==null?null:aa.ref;
			}

		
		private Map<String,String> getMap()
			{
			final Map<String, String> hash = new HashMap<String,String>();
			for(final String c: col2col.keySet())
				{
				int idx=col2col.get(c);
				if(idx>=this.tokens.length) continue;
				hash.put(c, tokens[idx]);
				}
			return hash;
			}
		
		public Set<SequenceOntologyTree.Term> getSOTerms()
			{
			final Set<SequenceOntologyTree.Term> set=new HashSet<SequenceOntologyTree.Term>();
			String EFF=getByCol("Effect");
			if(EFF==null) return set;
			for(final SequenceOntologyTree.Term t: SnpEffPredictionParser.this.soTree.getTerms())
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
