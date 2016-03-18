package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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

import com.github.lindenb.jvarkit.tools.vcfannot.VCFPredictions;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;

/**
 * @author lindenb
 *
 */
public class MyPredictionParser implements PredictionParser
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	private Map<VCFPredictions.FORMAT1, Integer> col2col=new HashMap<VCFPredictions.FORMAT1, Integer>();

	private Pattern pipe=Pattern.compile("[\\|]");

	public final String getTag()
		{
		return VCFPredictions.TAG;
		}
	
	
	
	public MyPredictionParser(VCFHeader header)
		{	
		VCFInfoHeaderLine info=header.getInfoHeaderLine(getTag());
		if(info==null || info.getDescription()==null)
			{
			LOG.warning("NO "+getTag()+" found in header");
			return;
			}
		String description=info.getDescription();
		String chunck=" Format:";
		int i=description.indexOf(chunck);
		if(i==-1)
			{
			LOG.warning("Cannot find "+chunck+ " in "+description);
			return;
			}
		description=description.substring(i+chunck.length()).replaceAll("[ \'\\.\\(\\)]+","").trim();
		String tokens[]=pipe.split(description);

		for(i=0;i< tokens.length;++i)
			{
			if(tokens[i].isEmpty()) continue;
			VCFPredictions.FORMAT1 col=null;
			for(VCFPredictions.FORMAT1 c:VCFPredictions.FORMAT1.values())
				{
				if(c.name().equalsIgnoreCase(tokens[i]))
					{
					col=c;
					}
				}
			if(col==null)
				{
				LOG.warning("Undefined "+" tag "+tokens[i]);
				continue;
				}
			col2col.put(col, i);
			}
		}
	
	
	
	@Override
	public List<MyPrediction> getPredictions(VariantContext ctx)
		{
		ArrayList<MyPrediction> preds= new ArrayList<MyPrediction>();
		if(col2col.isEmpty()) return preds;
		Object o=ctx.getAttribute(getTag());
		if(o==null)
			{
			return preds;
			}
		
		if(o.getClass().isArray())
			{
			for(Object o2:(Object[])o) _predictions(preds,o2);
			}
		else if(o instanceof Collection)
			{
			for(Object o2:(Collection<?>)o)  _predictions(preds,o2);
			}
		else
			{
			_predictions(preds, o);
			}
		return preds;
		}
	
	private void _predictions( List<MyPrediction> preds,Object o)
		{
		if(o==null) return;
		if(!(o instanceof String))
			{
			_predictions(preds, o.toString());
			return;
			}
		String s=String.class.cast(o).trim();
		
		String tokens[]=pipe.split(s);
		preds.add(new MyPrediction(tokens));
		}
			
	
	public class MyPrediction
		implements Prediction
		{
		private String tokens[];
		MyPrediction(String tokens[])
			{
			this.tokens=tokens;
			}			
		private String getByCol(VCFPredictions.FORMAT1 col)
			{
			Integer idx=col2col.get(col);
			if(idx==null || idx>=tokens.length || tokens[idx].isEmpty())
				{
				return null;
				}
			return tokens[idx];
			}
		
		public String getTranscript()
			{
			return getByCol(VCFPredictions.FORMAT1.TRANSCRIPT);
			}
		
		
		@Override
		public String getGeneName()
			{
			return getTranscript();
			}
		
		
		private Map<VCFPredictions.FORMAT1,String> getMap()
			{
			Map<VCFPredictions.FORMAT1, String> hash=new HashMap<VCFPredictions.FORMAT1,String>();
			for(VCFPredictions.FORMAT1 c: col2col.keySet())
				{
				int idx=col2col.get(c);
				if(idx>=this.tokens.length) continue;
				hash.put(c, tokens[idx]);
				}
			return hash;
			}
		
		
		public String getCodonChange()
			{
			return getByCol(VCFPredictions.FORMAT1.CODON);
			}
		
		public String getAltCodon()
			{
			String s=getCodonChange();
			if(s==null || s.isEmpty()) return null;
			int slash=s.indexOf('/');
			return slash==-1?null:s.substring(slash+1);
			}
		
		public String getRefCodon()
			{
			String s=getCodonChange();
			if(s==null || s.isEmpty()) return null;
			int slash=s.indexOf('/');
			return slash==-1?null:s.substring(0,slash);
			}


		
		
		public String getAminoAcidChange()
			{
			return getByCol(VCFPredictions.FORMAT1.AA);
			}
		
		@Override
		public String getAltAminoAcid()
			{
			String s=getAminoAcidChange();
			if(s==null || s.isEmpty()) return null;
			int slash=s.indexOf('/');
			return slash==-1?null:s.substring(slash+1);
			}
		
		
		@Override
		public Integer getAminoAcidPosition()
			{
			String s= getByCol(VCFPredictions.FORMAT1.PROTPOS);
			if(s==null || s.isEmpty()) return null;
			try {
				return Integer.parseInt(s);
			} catch (Exception e) {
				return null;
				}
			}
		@Override
		public String getReferenceAminoAcid() {

			String s=getAminoAcidChange();
			if(s==null || s.isEmpty()) return null;
			int slash=s.indexOf('/');
			return slash==-1?null:s.substring(0,slash);
			}
		
		@Override
		public String getEnsemblGene() {
			return null;
			}
		
		@Override
		public String getEnsemblProtein() {
			return null;
			}
		
		@Override
		public String getEnsemblTranscript() {
			
			return null;
			}
		
		@Override
		public Set<SequenceOntologyTree.Term> getSOTerms()
			{
			Set<SequenceOntologyTree.Term> set=new HashSet<SequenceOntologyTree.Term>();
			String EFF=getByCol(VCFPredictions.FORMAT1.SEQONTOLOGY);
			if(EFF==null) return set;
			for(SequenceOntologyTree.Term t:SequenceOntologyTree.getInstance().getTerms())
				{
				for(String eff:EFF.split("[&]"))
					{
					if(t.getLabel().replace(' ', '_').equals(eff))
						{
						set.add(t);
						}
					else if(t.getAcn().equals(eff))
						{
						set.add(t);
						}
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
