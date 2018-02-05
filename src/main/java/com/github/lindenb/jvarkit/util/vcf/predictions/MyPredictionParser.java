package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
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

import com.github.lindenb.jvarkit.tools.vcfannot.VCFPredictions;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;

/**
 * @author lindenb
 *
 */
public class MyPredictionParser implements PredictionParser
	{
	private static final Logger LOG=Logger.build(MyPredictionParser.class).make();
	private Map<VCFPredictions.FORMAT1, Integer> col2col=new HashMap<VCFPredictions.FORMAT1, Integer>();
	private SequenceOntologyTree soTree = SequenceOntologyTree.getInstance();
	private final Pattern pipe=Pattern.compile("[\\|]");
	private final Pattern ampRegex = Pattern.compile("[&]");

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
	
	public MyPredictionParser sequenceOntologyTree( final SequenceOntologyTree soTree) {
	this.soTree = soTree;
	return this;
	}
	
	@Override
	public List<MyPrediction> getPredictions(final VariantContext ctx)
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
		final MyPrediction p= parseOnePrediction(o);
		if(p==null) return;
		preds.add(p);
		}
			
	public MyPrediction  parseOnePrediction(final Object o)
		{
		if(o==null) return null;
		if(!(o instanceof String))
			{
			return parseOnePrediction( o.toString());
			}
		final String s=String.class.cast(o).trim();
		final String tokens[]=pipe.split(s);
		return new MyPrediction(s,tokens);
		}
	
	
	public class MyPrediction
		implements Prediction
		{
		private final String originalInfoStr;
		private final String tokens[];
		MyPrediction(final String originalInfoStr,final String tokens[])
			{
			this.originalInfoStr = originalInfoStr;
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
		
		public String getAltAminoAcid()
			{
			String s=getAminoAcidChange();
			if(s==null || s.isEmpty()) return null;
			int slash=s.indexOf('/');
			return slash==-1?null:s.substring(slash+1);
			}
		
		
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

		public String getReferenceAminoAcid() {

			String s=getAminoAcidChange();
			if(s==null || s.isEmpty()) return null;
			int slash=s.indexOf('/');
			return slash==-1?null:s.substring(0,slash);
			}
		
		public String getSOTermsString()
			{
			final String EFF=getByCol(VCFPredictions.FORMAT1.SEQONTOLOGY);
			return EFF==null?"":EFF;
			}
		
		public List<String> getSOTermsStrings()
			{
			final String soterms=getSOTermsString();
			if(soterms==null || soterms.isEmpty()) return Collections.emptyList();
			return Arrays.asList(MyPredictionParser.this.ampRegex.split(soterms));
			}
		
		public Set<SequenceOntologyTree.Term> getSOTerms()
			{
			final Set<SequenceOntologyTree.Term> set=new HashSet<>();
			for(final String eff:getSOTermsStrings())
				{
				SequenceOntologyTree.Term t  = MyPredictionParser.this.soTree.getTermByLabel(eff);
				if(t!=null) set.add(t);
				}
			return set;
			}

		public String getOriginalAttributeAsString() {
			return this.originalInfoStr;
		}
		
		
	@Override
	public String toString() {
		return getMap().toString()+ " "+Arrays.asList(tokens);
		}
	}
		
	
	}
