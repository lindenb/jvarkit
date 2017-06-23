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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;

/**
 * ###INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_po
sition|Protein_position|Amino_acids|Codons|Existing_variation|HGNC|DISTANCE|SIFT|PolyPhen|CELL_TYPE"

 * @author lindenb
 *
 */
public class VepPredictionParser implements PredictionParser
	{
	private static final Logger LOG=Logger.build(VepPredictionParser.class).make();
	public static final String INDEL_SYMBOL_STR="<indel>";
	
	/* public, used in VcfBurdenFilterGene 
	public enum COLS{
		Allele,
		Gene,
		Feature,
		Feature_type,
		Consequence,
		cDNA_position,
		CDS_position,
		Protein_position,
		Amino_acids,
		Codons,
		Existing_variation,
		HGNC,
		DISTANCE,
		SIFT,
		PolyPhen,
		CELL_TYPE,
		RefSeq,
		EXON,
		STRAND,
		SYMBOL,
		SYMBOL_SOURCE,
		HGNC_ID,
		IMPACT,BIOTYPE,INTRON,HGVSc,
		HGVSp,ALLELE_NUM,CANONICAL,
		CCDS,ENSP,DOMAINS
		*/
	private final Map<String, Integer> col2col=new HashMap<String, Integer>();
	private final Pattern pipe=Pattern.compile("[\\|]");
	private final Pattern ampRegex = Pattern.compile("[&]");
	private final String tag;
	private SequenceOntologyTree soTree = SequenceOntologyTree.getInstance();
	private final boolean valid;
	
	VepPredictionParser(final VCFHeader header)
		{		
		this(header,getDefaultTag());
		}
	
	public VepPredictionParser sequenceOntologyTree( final SequenceOntologyTree soTree) {
		this.soTree = soTree;
		return this;
		}
	
	@Override
	public String getTag()
		{
		return this.tag;
		}
	
	public static final String getDefaultTag()
		{
		return "CSQ";
		}
	
	VepPredictionParser(final VCFHeader header,final String tag)
		{	
		this.tag=(tag==null?getDefaultTag():tag);
		final VCFInfoHeaderLine info=(header==null?null:header.getInfoHeaderLine(tag));
		if(info==null || info.getDescription()==null)
			{
			LOG.warning("NO INFO["+tag+"] found in header. This VCF was probably NOT annotated with VEP. But it's not a problem if this tool doesn't need to access VEP Annotations.");
			this.valid = false;
			return;
			}
		String description=info.getDescription();
		final String chunck=" Format:";
		int i=description.indexOf(chunck);
		if(i==-1)
			{
			this.valid = false;
			LOG.warning("Cannot find "+chunck+ " in "+description);
			return;
			}
		description=description.substring(i+chunck.length()).replaceAll("[ \'\\.\\(\\)]+","").trim();
		final String tokens[]= this.pipe.split(description);

		for(i=0;i< tokens.length;++i)
			{
			if(tokens[i].isEmpty()) continue;
			if(this.col2col.containsKey(tokens[i]))
				{
				LOG.warning("Column  "+tokens[i]+" defined twice in "+description);;
				continue;
				}
			this.col2col.put(tokens[i], i);
			}
		this.valid=true;
		}
	
	public boolean isValid() {
		return valid;
		}
	
	public Set<String> getCategories() {
		return Collections.unmodifiableSet(this.col2col.keySet());
	}
	
	@Override
	public List<VepPrediction> getPredictions(final VariantContext ctx)
		{
		if(!isValid() || this.col2col.isEmpty()) return Collections.emptyList();
		final List<? extends Object> L =ctx.getAttributeAsList(this.tag);
		ArrayList<VepPrediction> preds= new ArrayList<VepPrediction>(L.size());
		for(final Object o2:L)  _predictions(preds,o2,ctx);
		return preds;
		}
	
	public VepPrediction parseOnePrediction(final VariantContext ctx,final Object o)
		{
		if(o==null || !isValid()) return null;
		if(!(o instanceof String))
			{
			return parseOnePrediction(ctx,o.toString());
			}
		final String s=String.class.cast(o).trim();
		final String tokens[]= this.pipe.split(s);
		return new VepPrediction(tokens,s,ctx);
		}
	
	private void _predictions(final List<VepPrediction> preds,final Object o,final VariantContext ctx)
		{
		final VepPrediction pred= parseOnePrediction(ctx,o);
		if(pred!=null) preds.add(pred);
		}
			
	
	public class VepPrediction
		implements Prediction
		{
		private final String source;
		private final String tokens[];
		VepPrediction(final String tokens[],final String source,final VariantContext ctx)
			{
			this.source=source;
			this.tokens=tokens;
			/** special case for ALT, can be '-' */
			Integer idx_allele = VepPredictionParser.this.col2col.get("Allele");
			if(	idx_allele!=null && 
				idx_allele<tokens.length &&
				tokens[idx_allele].equals("-"))
				{
				if(ctx.getAlternateAlleles().size()==1)
					{
					this.tokens[idx_allele]=ctx.getAlternateAlleles().get(0).getDisplayString();
					}
				else
					{
					this.tokens[idx_allele]=INDEL_SYMBOL_STR;
					}
				}
			}
		public String getByCol(final String col)
			{
			if(col==null || col.isEmpty()) return null;
			final Integer idx= VepPredictionParser.this.col2col.get(col);
			if(idx==null || idx>=tokens.length || tokens[idx].isEmpty())
				{
				return null;
				}
			return this.tokens[idx];
			}
		
		/** alias of getByColl */
		public String get(final String col)
			{
			return this.getByCol(col);
			}
		
		public boolean hasAttribute(final String col) {
			String s= this.get(col);
			return s!=null && !s.isEmpty();
			}
		
		/** return the 'Allele' String or null if not found */
		public String getAlleleStr()
			{
			return getByCol("Allele");
			}
	
		
		public Allele getAllele()
			{
			final String s= getAlleleStr();
			return s==null?null:Allele.create(s, false);
			}
		
		
		public String getExon()
			{
			return getByCol("EXON");
			}
		
		
		/** return -1 negative strand, 1 position strand , else 0 */
		public int getStrand()
			{
			final String s=getByCol("STRAND");
			if(s==null) return 0;
			if(s.equals("1")) return 1;
			if(s.equals("-1")) return -1;
			return 0;
			}
		
		
		/** alias of getHGNC */
		public String getGeneName()
			{
			return getByCol("HGNC");
			}
		
		public String getHGNC()
			{
			return getByCol("HGNC");
			}
		public String getHgncId()
			{
			return getByCol("HGNC_ID");
			}
		
		public String getSymbol()
			{
			return getByCol("SYMBOL");
			}
		
		public String getRefSeq()
			{
			return getByCol("RefSeq");
			}
		
		public String getFeature()
			{
			return getByCol("Feature");
			}
		public String getFeatureType()
			{
			return getByCol("Feature_type");
			}

		public String getGene()
			{
			return getByCol("Gene");
			}
		
		public String getENSP()
			{
			return getByCol("ENSP");
			}
		
		public String getSymbolSource()
			{
			return getByCol("SYMBOL_SOURCE");
			}
		
		public Map<String,String> getMap()
			{
			final Map<String, String> hash=new LinkedHashMap<>();
			for(final String c: col2col.keySet())
				{
				int idx=col2col.get(c);
				if(idx>=this.tokens.length) continue;
				hash.put(c, tokens[idx]);
				}
			return hash;
			}
		
		
		
		public String getEnsemblGene() {
			String s=getByCol("Gene");
			if(s==null) return null;
			return s;
			}
		
		
		public Double getSift() {
			String s=getByCol("SIFT");
			if(s==null) return null;
			try
				{
				return new Double(s);
				}
			catch(Exception err)
				{
				return null;
				}
			}
		
		public Double getPolyphen() {
			String s=getByCol("PolyPhen");
			if(s==null) return null;
			try
				{
				return new Double(s);
				}
			catch(Exception err)
				{
				return null;
				}
			}
		
		/** return the "Consequence" String or NULL if not found */
		public String getSOTermsString()
		{
			return getByCol("Consequence");
		}
		
		/** return the "Consequence" splitted, as an array of String, empty if consequence is not found */
		public List<String> getSOTermsStrings()
		{
			final String EFF=getSOTermsString();
			if(EFF==null || EFF.isEmpty()) return Collections.emptyList();
			return Arrays.asList(VepPredictionParser.this.ampRegex.split(EFF));
		}
	
		/** convert the list of getConsequences() to a list of SequenceOntology Terms */
		public Set<SequenceOntologyTree.Term> getSOTerms()
			{
			final List<String> effects = getSOTermsStrings();
			if(effects.isEmpty()) return Collections.emptySet();
			final Set<SequenceOntologyTree.Term> set=new HashSet<>(effects.size());
			for(final String label:effects) {
				if(label.isEmpty()) continue;
				final SequenceOntologyTree.Term t =VepPredictionParser.this.soTree.getTermByLabel(label);
				if(t==null) {
					LOG.warning("Current Sequence Ontology Tree doesn't contain "+ label);
					} 
				else
					{
					set.add(t);
					}
				}
			return set;
			}
		
		public Integer getPositionInCDna()
		{
			final String s= getByCol("cDNA_position");
			if(s==null || s.trim().isEmpty()) return null;
			try {
				return Integer.parseInt(s);
			} catch (Exception e) {
				return null;
			}
		}
		
		public Integer getPositionInCDS()
		{
			final String s= getByCol("CDS_position");
			if(s==null || s.trim().isEmpty()) return null;
			try {
				return Integer.parseInt(s);
			} catch (Exception e) {
				return null;
			}
		}
		
		
		private boolean _isEmpty(String s) {
			return s==null || s.trim().isEmpty();
		}
		
		/** returns some the posssible ways to name a gene */
		public Set<String> getGeneKeys() {
			final Set<String> keys= new HashSet<>();
			if(!_isEmpty(this.getFeature()))
				{
				keys.add("VEP_FEATURE_"+this.getFeature());
				}
			if(!_isEmpty(this.getGene()))
				{
				keys.add("VEP_GENE_"+this.getGene());
				}
			if(!_isEmpty(this.getSymbol()))
				{
				keys.add("VEP_SYMBOL_"+this.getSymbol());
				}
			if(!_isEmpty(this.getRefSeq()))
				{
				keys.add("VEP_REFSEQ_"+this.getRefSeq());
				}
			return keys;
			}
		
	/** return the prediction encoded in the original VariantContext */
	public String getOriginalAttributeAsString()
		{
		return this.source;
		}
	@Override
	public String toString() {
		return this.source;
		}
	}
		
	
	}
