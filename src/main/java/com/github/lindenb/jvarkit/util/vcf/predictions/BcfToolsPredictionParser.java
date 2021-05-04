/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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


*/
package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.Set;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;

/**
 * ##INFO=<ID=BCSQ,Number=.,Type=String,Description="Haplotype-aware consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change">

 * @author lindenb
 *
 */
public class BcfToolsPredictionParser implements PredictionParser
	{
	private static final Logger LOG=Logger.build(BcfToolsPredictionParser.class).make();
	
	private final Map<String, Integer> col2colidx=new HashMap<String, Integer>();
	private final CharSplitter pipe= CharSplitter.PIPE;
	private final CharSplitter ampRegex = CharSplitter.of('&');
	private final String tag;
	private SequenceOntologyTree soTree = SequenceOntologyTree.getInstance();
	private final boolean valid;
	private final Map<String, String> bcftools2so = new HashMap<>();
	
	BcfToolsPredictionParser(final VCFHeader header)
		{		
		this(header,getDefaultTag());
		}
	
	public BcfToolsPredictionParser sequenceOntologyTree( final SequenceOntologyTree soTree) {
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
		return "BCSQ";
		}
	
	BcfToolsPredictionParser(final VCFHeader header,final String tag)
		{	
		this.tag=(tag==null?getDefaultTag():tag);
		final VCFInfoHeaderLine info=(header==null?null:header.getInfoHeaderLine(tag));
		if(info==null || info.getDescription()==null)
			{
			this.valid = false;
			return;
			}
		String description=info.getDescription();
		final String chunck="Format: ";
		int i=description.indexOf(chunck);
		if(i==-1)
			{
			this.valid = false;
			LOG.warning("Cannot find "+chunck+ " in "+description);
			return;
			}
		description=description.substring(i+chunck.length()).trim();
		final List<String> tokens= this.pipe.splitAsStringList(description);

		for(i=0;i< tokens.size();++i)
			{
			final String token= tokens.get(i);
			if(StringUtil.isBlank(token)) continue;
			if(this.col2colidx.containsKey(token))
				{
				LOG.warning("Column  "+token+" defined twice in "+description);;
				continue;
				}
			this.col2colidx.put(token, i);
			}
		this.valid=true;
		
		// https://github.com/samtools/bcftools/blob/master/csq.c#L229
		this.bcftools2so.put("3_prime_utr", "3_prime_UTR_variant");
		this.bcftools2so.put("5_prime_utr", "5_prime_UTR_variant");
		this.bcftools2so.put("non_coding","non_coding_transcript_variant");
		this.bcftools2so.put("missense","missense_variant");
		this.bcftools2so.put("splice_acceptor","splice_acceptor_variant");
		this.bcftools2so.put("splice_donor","splice_donor_variant");
		this.bcftools2so.put("synonymous", "synonymous_variant");
		this.bcftools2so.put("stop_retained", "stop_retained_variant");
		this.bcftools2so.put("intron", "intron_variant");
		this.bcftools2so.put("intergenic", "intergenic_variant");
		this.bcftools2so.put("inframe_altering", "inframe_variant");//?
		this.bcftools2so.put("start_retained", "start_retained_variant");
		this.bcftools2so.put("coding_sequence", "coding_sequence_variant");
		}
	
	public boolean isValid() {
		return valid;
		}
	
	public Set<String> getCategories() {
		return Collections.unmodifiableSet(this.col2colidx.keySet());
	}
	
	@Override
	public List<BcfToolsPrediction> getPredictions(final VariantContext ctx)
		{
		if(!isValid() || this.col2colidx.isEmpty()) return Collections.emptyList();
		final List<? extends Object> L =ctx.getAttributeAsList(this.tag);
		ArrayList<BcfToolsPrediction> preds= new ArrayList<>(L.size());
		for(final Object o2:L)  _predictions(preds,o2,ctx);
		return preds;
		}
	
	public BcfToolsPrediction parseOnePrediction(final VariantContext ctx,final Object o)
		{
		if(o==null || !isValid()) return null;
		if(!(o instanceof String))
			{
			return parseOnePrediction(ctx,o.toString());
			}
		final String s=String.class.cast(o).trim();
		final String tokens[]= this.pipe.split(s);
		return new BcfToolsPrediction(ctx,s,tokens);
		}
	
	private void _predictions(final List<BcfToolsPrediction> preds,final Object o,final VariantContext ctx)
		{
		final BcfToolsPrediction pred= parseOnePrediction(ctx,o);
		if(pred!=null) preds.add(pred);
		}
			
	
	public class BcfToolsPrediction
		implements Prediction
		{
		private final String originalAttributeAsString;
		private final String tokens[];
		private final VariantContext ctx;
		BcfToolsPrediction(final VariantContext ctx,final String originalAttributeAsString,final String tokens[])
			{
			this.originalAttributeAsString = originalAttributeAsString;
			this.tokens=tokens;
			this.ctx = ctx;
			}
		/** get column by name, may return null. Returns null if column is empty */
		public String getByCol(final String col)
			{
			final Integer idx=col2colidx.get(col);
			if(idx==null || idx>=this.tokens.length || tokens[idx].isEmpty()) return null;
			return tokens[idx];
			}
		
		public boolean isIntergenicRegion() {
			return StringUtil.isBlank(getGeneName()) && StringUtil.isBlank(getTranscript());
		}
		
		public String getGeneName()
			{
			return getByCol("gene");
			}
		public String getTranscript()
			{
			return getByCol("transcript");
			}
		public String getTranscriptBioType() {
			return getByCol("biotype");
		}
		
		public String getOriginalAttributeAsString()
			{
			return originalAttributeAsString;
			}
		
		public String getStrand() {
			return getByCol("strand");
		}
		public String getAminoAcidChange() {
			return getByCol("amino_acid_change");
		}
		public String getDnaChange() {
			return getByCol("dna_change");
		}
		
		
		private Map<String,String> getMap()
			{
			final Map<String, String> hash = new HashMap<String,String>(col2colidx.size());
			for(final String c: col2colidx.keySet())
				{
				int idx=col2colidx.get(c);
				if(idx>=this.tokens.length) continue;
				hash.put(c, tokens[idx]);
				}
			return hash;
			}
		
		public String getSOTermsString()
			{
			/* reference to another position */
			if(this.originalAttributeAsString.startsWith("@")) return null;
			final String s = getByCol("Consequence");
			// " The consequence can start with the asterisk '*' prefix indicating a consequence downstream from a stop"
			if(s==null ) return null;
			return s.startsWith("*")?s.substring(1):s;
			}
		
		/** The consequence can start with the asterisk '*' prefix indicating a consequence downstream from a stop" */
		public boolean isDownstreamAStop() {
			final String s = getByCol("Consequence");
			return s!=null && s.startsWith("*");
		}
		/** Consequences of compound variants which span multiple sites are printed in one record only, the remaining records link to it by '@position */
		public OptionalInt getReferencePosition() {
			if(!this.originalAttributeAsString.startsWith("@")) return OptionalInt.empty();
			return OptionalInt.of(Integer.parseInt(this.originalAttributeAsString.substring(1)));
		}
		
		/** BCFtools csq doesn't use SO !!! */
		private String mapSoTerm(final String s) {
			return bcftools2so.getOrDefault(s, s);
		}
		
		public Set<SequenceOntologyTree.Term> getSOTerms()
			{
			final String EFFs=getSOTermsString();
			if(StringUtil.isBlank(EFFs)) return Collections.emptySet();
			final String tokens[] = ampRegex.split(EFFs);
			final Set<SequenceOntologyTree.Term> set=new LinkedHashSet<>(tokens.length);

			for(final String EFF: tokens) {
				final String soTerm = mapSoTerm(EFF);
				final SequenceOntologyTree.Term t = BcfToolsPredictionParser.this.soTree.getTermByLabel(soTerm);
				if(t==null) {
					LOG.warn("Cannot get CSQ prediction \""+EFF+"\"/\""+soTerm+"\" in Sequence Ontology.");
					}
				set.add(t);
				}
			return set;
			}
		
		/** return ALT allele DNA change e.g: 48305542T>TGGGCCTGGGATC+48305543C>A */
		public String getAllele()  {
			final String s = getDnaChange();
			if(StringUtil.isBlank(s)) return null;
			final String prefix= String.valueOf(this.ctx.getStart())+ this.ctx.getReference().getDisplayString()+">";
			for(final String token: CharSplitter.of('+').split(s)) {
				if(!token.startsWith(prefix)) continue;
				String allele =  token.substring(prefix.length()).trim();
				if(StringUtil.isBlank(allele)) continue;
				// a priori , no need to reverse complement according to strand
				return allele;
				}
			return null;
			}

		
		@Override
		public String toString() {
			return getMap().toString()+ " "+Arrays.asList(tokens);
			}
		}
		
	
	}
