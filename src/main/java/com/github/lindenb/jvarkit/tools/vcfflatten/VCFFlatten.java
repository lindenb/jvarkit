/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfflatten;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

# Motivation

flatten all variants to one variant. At the end, a Sample will be set to 1/1 if one or more of his Genotypes
in the VCF is contains an NON-(REF|NO_CALL) allele
The idea is to use the VCF output of this tool
in order to pipe it into `bcftools contrast`: we'll get a p-value for the whole VCF


# Example

```
$ java -jar dist/jvarkit.jar vcfflatten --gene-extractor src/test/resources/rotavirus_rf.ann.vcf.gz 2> /dev/null 
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=MULTIPLE_CONTIG,Number=0,Type=Flag,Description="Record spans multiple chromosomes. Only first chromosome is reported">
##INFO=<ID=N_VARIANTS,Number=1,Type=Integer,Description="Number of variants">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF01	970	.	N	<FLATTEN_VARIANT>	.	.	MULTIPLE_CONTIG;N_VARIANTS=45	GT	1/1	1/1	1/1	1/1	1/1
RF01	970	.	N	<ANN/GeneName:Gene_18_3284>	.	.	N_VARIANTS=1	GT	0/0	0/0	0/0	0/0	1/1
RF01	970	.	N	<ANN/FeatureId:AAA47319.1>	.	.	N_VARIANTS=1	GT	0/0	0/0	0/0	0/0	1/1
RF01	970	.	N	<ANN/GeneId:Gene_18_3284>	.	.	N_VARIANTS=1	GT	0/0	0/0	0/0	0/0	1/1
RF02	1962	.	N	<ANN/GeneId:UniProtKB/Swiss-Prot:P12472>	.	.	END=251;N_VARIANTS=5	GT	1/1	1/1	1/1	1/1	0/0
RF02	1962	.	N	<ANN/FeatureId:CAA32215.1>	.	.	END=251;N_VARIANTS=5	GT	1/1	1/1	1/1	1/1	0/0
RF02	1962	.	N	<ANN/FeatureId:CAA32213.1>	.	.	END=251;N_VARIANTS=5	GT	1/1	1/1	1/1	1/1	0/0
RF02	1962	.	N	<ANN/GeneName:Gene_1621_1636>	.	.	END=251;N_VARIANTS=5	GT	1/1	1/1	1/1	1/1	0/0
```

END_DOC
 */
@Program(
	name="vcfflatten",
	description="Flatten variants to one variant",
	keywords={"vcf","burden","contrast"},
	creationDate = "20230222",
	modificationDate  = "20230222",
	menu="VCF Manipulation",
	jvarkit_amalgamion = true
	)
public class VCFFlatten extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VCFFlatten.class);
	@Parameter(names={"-i","--id"},description="Default Variant ID")
	protected String default_variant_id = "FLATTEN_VARIANT";
	@Parameter(names={"--gene-extractor"},description="Activate default gene extractors. Variant will be grouped by gene using snpeff/bcftools/vep annotations")
	protected boolean use_gene_extractors = false;

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	private Group default_group = null;
	private final List<GeneExtractorFactory.GeneExtractor> geneExtractors = new ArrayList<>();
	private final Map<GeneExtractorFactory.KeyAndGene,Group> keyAndGeneToGroup = new HashMap<>();
	
	private static class Group implements Locatable {
		final String name;
		final BitSet has_ALT;
		String contig=null;
		int start=0;
		int end=0;
		boolean multiple_contigs = false;
		int n_variants = 0;
		
		Group(final String name,int n_samples) {
			this.name = name;
			this.has_ALT = new BitSet(n_samples);
			}
		@Override
		public String getContig() {return this.contig;}
		@Override
		public int getStart() {return this.start;}
		@Override
		public int getEnd() {return this.end;}
		
		
		Allele getAlt() {
			return Allele.create("<"+name+">", false);
			}
		void visit(final VariantContext ctx) {
			this.n_variants++;
			if(contig==null) {
				contig = ctx.getContig();
				start = ctx.getStart();
				end = ctx.getEnd();
				}
			else if(contig.equals(ctx.getContig())) {
				start = Math.min(ctx.getStart(), start);
				end = Math.max(ctx.getEnd(), end);
				}
			else
				{
				multiple_contigs = true;
				}
			for(int i=0;i< ctx.getNSamples();i++) {
				if(this.has_ALT.get(i)) continue;
				final Genotype gt = ctx.getGenotype(i);
				if(gt.getAlleles().stream().anyMatch(A->!(A.isReference() || A.isNoCall()))) {
					this.has_ALT.set(i);
					}
				}
			}
		}
	
	private Collection<Group>  extractGroups(final VariantContext ctx) {
		if(default_group==null) {
			this.default_group = new Group(this.default_variant_id, ctx.getNSamples());
			}
		final List<Group> ret = new ArrayList<>();
		ret.add(this.default_group);

		for(GeneExtractorFactory.GeneExtractor extractor: this.geneExtractors) {
			Map<GeneExtractorFactory.KeyAndGene,Set<String>> gene2values = extractor.apply(ctx);
			if(gene2values.isEmpty()) continue;
			
			for(final GeneExtractorFactory.KeyAndGene keyAndGene :gene2values.keySet()) {
				Group g = this.keyAndGeneToGroup.get(keyAndGene);
				if(g==null) {
					g = new Group(keyAndGene.getMethod()+":"+keyAndGene.getKey(), ctx.getNSamples());
					this.keyAndGeneToGroup.put(keyAndGene,g);
					}
				ret.add(g);
				}
			}
		return ret;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		final VCFHeader headerin = iterin.getHeader();
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(headerin);
		final Set<VCFHeaderLine> metaData = new HashSet<>();
		VCFStandardHeaderLines.addStandardFormatLines(metaData, true,VCFConstants.GENOTYPE_KEY);
		VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.END_KEY);
		final VCFInfoHeaderLine info_multi_flag = new VCFInfoHeaderLine("MULTIPLE_CONTIG", 1, VCFHeaderLineType.Flag ,"Record spans multiple chromosomes. Only first chromosome is reported");
		metaData.add(info_multi_flag);
		final VCFInfoHeaderLine info_n_variant = new VCFInfoHeaderLine("N_VARIANTS", 1, VCFHeaderLineType.Integer ,"Number of variants");
		metaData.add(info_n_variant);

		final VCFHeader header = new VCFHeader(metaData,headerin.getGenotypeSamples());
		header.setSequenceDictionary(dict);
		final List<String> sampleNames = headerin.getGenotypeSamples();
		
		if(use_gene_extractors) {
			this.geneExtractors.addAll( new GeneExtractorFactory(headerin).getAllExtractors());
		}
		
		/* paranoid test */
		if( !headerin.getSampleNameToOffset().equals( header.getSampleNameToOffset()) ||
				!headerin.getGenotypeSamples().equals( header.getGenotypeSamples())) {
			LOG.error("internal test failed");
			return -1;
		}
		
		JVarkitVersion.getInstance().addMetaData(this, header);
		out.writeHeader(header);
		
		while(iterin.hasNext()) {
			final VariantContext ctx = iterin.next();
			for(final Group group: extractGroups(ctx)) {
				group.visit(ctx);
				}
			}
		
		final List<Group> all_groups = new ArrayList<>(1+this.keyAndGeneToGroup.size());
		all_groups.add(this.default_group);
		all_groups.addAll(this.keyAndGeneToGroup.values());
		
		Collections.sort(all_groups, new ContigDictComparator(dict).createLocatableComparator());
		
		for(Group g:all_groups) {
			final Allele ALT= g.getAlt();
			final Allele REF= Allele.REF_N;
			final List<Allele> HOM_VAR = Arrays.asList(ALT,ALT);
			final List<Allele> HOM_REF = Arrays.asList(REF,REF);
			final VariantContextBuilder vcb=new VariantContextBuilder(
				null,
				g.getContig(),
				g.getStart(),
				g.getEnd(),
				Arrays.asList(REF,ALT)
				);
			final List<Genotype> genotypes = new ArrayList<>(sampleNames.size());
			for(int i=0;i< sampleNames.size();i++) {
				genotypes.add(
					new GenotypeBuilder(
						sampleNames.get(i),
						g.has_ALT.get(i)?HOM_VAR:HOM_REF
						).make());
				}
			if(g.getStart()!=g.getEnd()) vcb.attribute(VCFConstants.END_KEY,g.getEnd());
			if(g.multiple_contigs) {
				vcb.attribute(info_multi_flag.getID(), true);
				}
			vcb.attribute(info_n_variant.getID(), g.n_variants);
			vcb.genotypes(genotypes);
			out.add(vcb.make());
			}
		this.keyAndGeneToGroup.clear();
		this.geneExtractors.clear();
		return 0;
		}
public static void main(final String[] args) {
	new VCFFlatten().instanceMainWithExit(args);
	}
}
