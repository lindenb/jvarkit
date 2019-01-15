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

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
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
	private static final Map<String,String> EFF2SO= new HashMap<>() ;
	{{{
		// https://github.com/pcingola/SnpEff/blob/master/src/main/java/org/snpeff/snpEffect/EffectType.java
		EFF2SO.put("CHROMOSOME_LARGE_DELETION","chromosome_number_variation");
		EFF2SO.put("CHROMOSOME_LARGE_DUPLICATION","duplication");
		EFF2SO.put("CHROMOSOME_LARGE_INVERSION","inversion");
		EFF2SO.put("CHROMOSOME","chromosome");
		EFF2SO.put("CHROMOSOME_ELONGATION","feature_elongation");
		EFF2SO.put("CODON_CHANGE","coding_sequence_variant");
		EFF2SO.put("CODON_CHANGE_PLUS_CODON_INSERTION","disruptive_inframe_insertion");
		EFF2SO.put("CODON_CHANGE_PLUS_CODON_DELETION","disruptive_inframe_deletion");
		EFF2SO.put("CODON_DELETION","conservative_inframe_deletion");
		EFF2SO.put("CODON_INSERTION","conservative_inframe_insertion");
		EFF2SO.put("DOWNSTREAM","downstream_gene_variant");
		EFF2SO.put("EXON","exon_region");
		EFF2SO.put("EXON_DELETED","exon_loss_variant");
		EFF2SO.put("EXON_DELETED_PARTIAL","exon_loss_variant");
		EFF2SO.put("EXON_DUPLICATION","duplication");
		EFF2SO.put("EXON_DUPLICATION_PARTIAL","duplication");
		EFF2SO.put("EXON_INVERSION","inversion");
		EFF2SO.put("EXON_INVERSION_PARTIAL","inversion");
		EFF2SO.put("FEATURE_FUSION","feature_fusion");
		EFF2SO.put("FRAME_SHIFT","frameshift_variant");
		EFF2SO.put("GENE","gene_variant");
		EFF2SO.put("GENE_INVERSION","inversion");
		EFF2SO.put("GENE_DELETED","feature_ablation");
		EFF2SO.put("GENE_DUPLICATION","duplication");
		EFF2SO.put("GENE_FUSION","gene_fusion");
		EFF2SO.put("GENE_FUSION_HALF","transcript_ablation");
		EFF2SO.put("GENE_FUSION_REVERESE","bidirectional_gene_fusion");
		EFF2SO.put("GENE_REARRANGEMENT","rearranged_at_DNA_level");
		EFF2SO.put("INTERGENIC","intergenic_region");
		EFF2SO.put("INTERGENIC_CONSERVED","conserved_intergenic_variant");
		EFF2SO.put("INTRON","intron_variant");
		EFF2SO.put("INTRON_CONSERVED","conserved_intron_variant");
		EFF2SO.put("INTRAGENIC","intragenic_variant");
		EFF2SO.put("MICRO_RNA","miRNA");
		EFF2SO.put("MOTIF","TF_binding_site_variant");
		EFF2SO.put("MOTIF_DELETED","TFBS_ablation");
		EFF2SO.put("NEXT_PROT","sequence_feature");
		EFF2SO.put("NON_SYNONYMOUS_CODING","missense_variant");
		EFF2SO.put("NON_SYNONYMOUS_START","initiator_codon_variant");
		EFF2SO.put("NON_SYNONYMOUS_STOP","stop_retained_variant");
		EFF2SO.put("PROTEIN_PROTEIN_INTERACTION_LOCUS","protein_protein_contact");
		EFF2SO.put("PROTEIN_STRUCTURAL_INTERACTION_LOCUS","structural_interaction_variant");
		EFF2SO.put("RARE_AMINO_ACID","rare_amino_acid_variant");
		EFF2SO.put("REGULATION","regulatory_region_variant");
		EFF2SO.put("SPLICE_SITE_ACCEPTOR","splice_acceptor_variant");
		EFF2SO.put("SPLICE_SITE_DONOR","splice_donor_variant");
		EFF2SO.put("SPLICE_SITE_REGION","splice_region_variant");
		EFF2SO.put("SPLICE_SITE_BRANCH","splice_branch_variant");
		EFF2SO.put("SPLICE_SITE_BRANCH_U12","splice_branch_variant");
		EFF2SO.put("START_LOST","start_lost");
		EFF2SO.put("START_GAINED","5_prime_UTR_premature_start_codon_gain_variant");
		EFF2SO.put("STOP_GAINED","stop_gained");
		EFF2SO.put("STOP_LOST","stop_lost");
		EFF2SO.put("SYNONYMOUS_CODING","synonymous_variant");
		EFF2SO.put("SYNONYMOUS_STOP","stop_retained_variant");
		EFF2SO.put("SYNONYMOUS_START","initiator_codon_variant");
		EFF2SO.put("TRANSCRIPT","non_coding_transcript_variant");
		EFF2SO.put("TRANSCRIPT_DELETED","transcript_ablation");
		EFF2SO.put("TRANSCRIPT_DUPLICATION","duplication");
		EFF2SO.put("TRANSCRIPT_INVERSION","inversion");
		EFF2SO.put("UPSTREAM","upstream_gene_variant");
		EFF2SO.put("UTR_3_PRIME","3_prime_UTR_variant");
		EFF2SO.put("UTR_3_DELETED","3_prime_UTR_truncation");
		EFF2SO.put("UTR_5_PRIME","5_prime_UTR_variant");
		EFF2SO.put("UTR_5_DELETED","5_prime_UTR_truncation");		
	}}}
		
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
		final String tokens[]=pipe.split(description);
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
		final String ostr = String.class.cast(o).trim();
		final String tokens[]=pipe.split(ostr);
		return new SnpEffPrediction(ostr,tokens);
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
		private final String originalAttributeAsString;
		private final String tokens[];
		SnpEffPrediction(final String originalAttributeAsString,final String tokens[])
			{
			this.originalAttributeAsString = originalAttributeAsString;
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
		public String getOriginalAttributeAsString()
			{
			return originalAttributeAsString;
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
		public String getSOTermsString()
			{
			return	getByCol("Effect");
			}
		
		public Set<String> getSOTermsStrings()
			{
			final String EFF=getSOTermsString();
			if(StringUtil.isBlank(EFF)) return Collections.emptySet();
			
			final String EFF2;
			if(EFF2SO.containsKey(EFF)) {
				EFF2 = EFF2SO.get(EFF);
				}
			else
				{
				EFF2 =EFF.toLowerCase();
				}
			return Collections.singleton(EFF2);
			}
		
		public Set<SequenceOntologyTree.Term> getSOTerms()
			{
			final Set<String> EFFs=getSOTermsStrings();
			if(EFFs.isEmpty()) return Collections.emptySet();
			
			final Set<SequenceOntologyTree.Term> set=new HashSet<>(EFFs.size());

			for(final String EFF: EFFs) {
				final SequenceOntologyTree.Term t = SnpEffPredictionParser.this.soTree.getTermByLabel(EFF);
				if(t==null) {
					LOG.warn("Cannot get snpeff prediction \""+EFF+"\" in Sequence Ontology");
					}
				set.add(t);
				}
			return set;
			}
		
		/** form : Codon_Change e.g: Tga/Cga */
		public String getCodonChange() {
			return getByCol("Codon_Change");
		}
		
		/** return ALT allele if found in Codon_Change ( e.g: Tga/Cga), may be null */
		public Allele getAllele()  {
			String s = getCodonChange();
			if(StringUtil.isBlank(s)) return null;
			final int slash = s.indexOf("/");
			if(slash==-1) return null;
			s = s.substring(slash+1).replaceAll("[^A-Z]+","");// remove all lower case
			if(StringUtil.isBlank(s)) return null;
			return Allele.create(s,false);
			}

		
		@Override
		public String toString() {
			return getMap().toString()+ " "+Arrays.asList(tokens);
			}
		}
	
	}
