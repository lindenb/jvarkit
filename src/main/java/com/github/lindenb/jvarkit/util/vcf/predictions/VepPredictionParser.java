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
 * ###INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_po
sition|Protein_position|Amino_acids|Codons|Existing_variation|HGNC|DISTANCE|SIFT|PolyPhen|CELL_TYPE"

 * @author lindenb
 *
 */
public class VepPredictionParser implements PredictionParser
	{

	private enum COLS{Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,HGNC,DISTANCE,SIFT,PolyPhen,CELL_TYPE};
	private Map<COLS, Integer> col2col=new HashMap<COLS, Integer>();
	private Pattern pipe=Pattern.compile("[\\|]");
	private String tag;
	public VepPredictionParser(VCFHeader header)
		{		
		this(header,"CSQ");
		}
	
	public VepPredictionParser(VCFHeader header,String tag)
		{	
		this.tag=tag;
		VCFInfoHeaderLine info=header.getInfoHeaderLine(tag);
		if(info==null || info.getDescription()==null)
			{
			System.err.println("NO "+tag+" found in header");
			return;
			}
		String description=info.getDescription();
		String chunck=" Format:";
		int i=description.indexOf(chunck);
		if(i==-1)
			{
			System.err.println("Cannot find "+chunck+ " in "+description);
	
			return;
			}
		description=description.substring(i+chunck.length()).replaceAll("[ \'\\.\\(\\)]+","").trim();
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
				System.err.println("Undefined VEP tag "+tokens[i]);
				}
			col2col.put(col, i);
			}
		}
	@Override
	public List<? extends Prediction> getPredictions(VariantContext ctx)
		{
		 ArrayList<Prediction> preds= new ArrayList<Prediction>();
		if(col2col.isEmpty()) return preds;
		Object o=ctx.getAttribute(this.tag);
		if(o==null)
			{
			return preds;
			}
		
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
			if(idx==null || idx>=tokens.length || tokens[idx].isEmpty())
				{
				return null;
				}
			return tokens[idx];
			}
		@Override
		public String getGeneName()
			{
			return getByCol(COLS.HGNC);
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
		public String getAltAminoAcid() {
			return null;
			}
		
		@Override
		public Integer getAminoAcidPosition()
			{
			return null;
			}
		@Override
		public String getReferenceAminoAcid() {
			return null;
			}
		
		@Override
		public String getEnsemblGene() {
			String s=getByCol(COLS.Gene);
			if(s.startsWith("ENSG")) return s;
			return s;
			}
		
		@Override
		public String getEnsemblProtein() {
			return null;
			}
		
		@Override
		public String getEnsemblTranscript() {
			String s=getByCol(COLS.Feature);
			if(s.startsWith("ENST")) return s;
			return null;
			}
	
	@Override
	public String toString() {
		return getMap().toString()+ " "+Arrays.asList(tokens);
		}
	}
		
	
	}
