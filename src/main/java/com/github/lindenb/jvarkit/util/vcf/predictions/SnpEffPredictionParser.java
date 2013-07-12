package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;


import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;


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
		public String toString() {
			return getMap().toString()+ " "+Arrays.asList(tokens);
			}
		}
	
	}
