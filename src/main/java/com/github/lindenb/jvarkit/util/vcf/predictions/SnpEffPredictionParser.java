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


import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;


/**
 * ##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon  | GenotypeNum [ | ERRORS | WARNINGS ] )'">

 * @author lindenb
 *
 */
public class SnpEffPredictionParser implements PredictionParser
	{
	private enum COLS{ Effect , Effect_Impact , Functional_Class,Codon_Change, 
		Amino_Acid_change,Amino_Acid_length,Gene_Name , Gene_BioType , Coding , Transcript,
		Exon  , GenotypeNum , ERRORS , WARNINGS};
	private Map<COLS, Integer> col2col=new HashMap<COLS, Integer>();
	private Pattern pipe=Pattern.compile("[\\|\\(\\)]");
	private String tag;
	public SnpEffPredictionParser(VCFHeader header)
		{		
		this(header,"EFF");
		
		}
	
	public SnpEffPredictionParser(VCFHeader header,String tag)
		{	
		this.tag=tag;
		VCFInfoHeaderLine info=header.getInfoHeaderLine(tag);
		if(info==null || info.getDescription()==null)
			{
			System.err.println("no INFO["+tag+"] or no description ");
			return;
			}
		String description=info.getDescription();
		String chunck="Format:";
		int i=description.indexOf(chunck);
		if(i==-1)
			{
			System.err.println("Cannot find "+chunck+ " in "+description);

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
				System.err.println("Undefined SnpEff tag \""+tokens[i]+"\"");
				}
			col2col.put(col, i);
			}
		}
	@Override
	public List<? extends Prediction> getPredictions(VariantContext ctx)
		{
		 ArrayList<Prediction> preds= new ArrayList<Prediction>();
		if(col2col.isEmpty())
			{
			return preds;
			}
		Object o=ctx.getAttribute(this.tag);
		if(o==null) return preds;
		if(o.getClass().isArray())
			{
			for(Object o2:(Object[])o) _predictions(preds,o2);
			}
		else if(o instanceof List)
			{
			for(Object o2:(List<?>)o)  _predictions(preds,o2);
			}
		else
			{
			_predictions(preds,Collections.singleton(o));
			}
		return preds;
		}
	
	private void _predictions( List<Prediction> preds,Object o)
		{
		if(o==null) return;
		if(!(o instanceof String))
			{
			_predictions(preds, o.toString());
			return;
			}
		String s=String.class.cast(o).trim();
		String tokens[]=pipe.split(s);
		preds.add(new MyPred(tokens));
		}
	
	private static class AAChange
		{
		String ref;
		int pos;
		String alt;
		}
			
	
	private class MyPred
		implements Prediction
		{
		private String tokens[];
		MyPred(String tokens[])
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

		
		public Map<COLS,String> getMap()
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
