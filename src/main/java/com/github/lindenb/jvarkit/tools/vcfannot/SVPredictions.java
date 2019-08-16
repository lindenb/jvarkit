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


*/
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GftReader;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;


/**
BEGIN_DOC


END_DOC
*/

@Program(name="svpredictions",
	description="Basic Variant Effect prediction using gtf",
	keywords={"vcf","annotation","prediction","sv"},
	creationDate="20190815",
	modificationDate="20190815"
	)
public class SVPredictions extends Launcher
	{
	private static final Logger LOG = Logger.build(SVPredictions.class).make();
	private final IntervalTreeMap<Gene> all_gene = new IntervalTreeMap<>();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-g","--gtf"},description="GTF File",required=true)
	private Path gtfPath = null;
	@Parameter(names={"-u","--upstream"},description="Upstream size. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int upstream_size =  5_000;
	@Parameter(names={"--max-genes"},description="don't print the genes names if their count exceed 'x'")
	private int max_genes_count =  20;
	@Parameter(names={"-nti","--no-transcript-id"},description="don't print transcript id (reduce length of annotation)")
	private boolean ignore_transcript_id = false;

	@Parameter(names={"--tag"},description="VCF info attribute")
	private String info_tag = "SVCSQ";
		
	private String annotTranscript(final Transcript transcript, final VariantContext ctx) {
		final String suffix="|"+transcript.getGene().getId()+"|"+(ignore_transcript_id?"":transcript.getId()+"|")+transcript.getGene().getGeneName()+"|"+transcript.getGene().getGeneBiotype();
		
		if(!transcript.overlaps(ctx)) {
			return "upstream_transcript_variant"+suffix;
			}
		
		if(ctx.contains(transcript)) {
			return "transcript_ablation"+suffix;
			}
		
		if(transcript.getUTRs().stream().anyMatch(U->U.contains(ctx))) {
			return "utr"+suffix;
			}
		
		if(transcript.hasCDS() && transcript.getAllCds().stream().anyMatch(U->U.contains(ctx))) {
			return "cds"+suffix;
			}

		if(transcript.getExons().stream().anyMatch(U->ctx.contains(U))) {
			return "exon"+suffix;
			}
		
		if(transcript.getIntrons().stream().anyMatch(U->ctx.contains(U))) {
			return "intron"+suffix;
			}
		
		final Set<String> set=new HashSet<>();
		
		if(transcript.getUTRs().stream().anyMatch(U->U.overlaps(ctx))) {
			set.add("utr");
			}
		
		
		if(transcript.hasCDS() && transcript.getAllCds().stream().anyMatch(U->U.overlaps(ctx))) {
			set.add("cds");
			}
		

		if(transcript.getExons().stream().anyMatch(U->ctx.overlaps(U))) {
			set.add("exon");
			}
		
		if(transcript.getIntrons().stream().anyMatch(U->ctx.overlaps(U))) {
			set.add("intron");
			}
		return String.join("&", set) + suffix;
		}
	
	private Set<String> annotGene(final Gene gene, final VariantContext ctx) {
		final Set<String> annot = new HashSet<>();
		for(final Transcript transcript: gene.getTranscripts()) {
			annot.add(annotTranscript(transcript, ctx));
			}
		return annot;
		}
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator r, VariantContextWriter w)
		{
		try {
		final VCFHeader header= r.getHeader();
		final SAMSequenceDictionary dict=  header.getSequenceDictionary();
		try(final GftReader gtfReader=new GftReader(this.gtfPath)) {
			if(dict!=null) gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
			gtfReader.getAllGenes().stream().forEach(G->this.all_gene.put(new Interval(G), G));
			}
		
		final VCFHeader h2=new VCFHeader(header);
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				this.info_tag,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				"Structural variant consequence."
				));
		JVarkitVersion.getInstance().addMetaData(this, h2);
		w.writeHeader(h2);

		final ProgressFactory.Watcher<VariantContext> progress=ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
		while(r.hasNext())
			{
			final VariantContext ctx=progress.apply(r.next());
			
			final Collection<Gene> genes  = this.all_gene.getOverlapping(new Interval(
					ctx.getContig(),
					Math.max(1,ctx.getStart()-this.upstream_size),
					ctx.getEnd()+this.upstream_size
					));
			
			if(genes.isEmpty())
				{
				Gene leftGene = null;
				for(final Gene g:this.all_gene.getOverlapping(new Interval(
					ctx.getContig(),
					1,
					ctx.getStart()
					)))
					{
					if(leftGene==null || leftGene.getEnd() < g.getEnd()) leftGene=g;
					}

				final String leftId = (leftGene==null?"":leftGene.getId());
				final String leftName =  (leftGene==null?"":leftGene.getGeneName());
				
				
				Gene rightGene = null;
				
				for(final Gene g:this.all_gene.getOverlapping(new Interval(
						ctx.getContig(),
						ctx.getStart(),
						Integer.MAX_VALUE
						)))
						{
						if(rightGene==null || rightGene.getStart() > g.getStart()) rightGene=g;
						}
				final String rightId = (rightGene==null?"":rightGene.getId());
				final String rightName =  (rightGene==null?"":rightGene.getGeneName());
				
				//intergenic
				if(leftGene==null && rightGene==null) {
					w.add(ctx);
					} 
				else {
					w.add(new VariantContextBuilder(ctx).
						attribute(this.info_tag, "intergenic|"+leftId+"-"+rightId+"|"+leftName+"-"+rightName).
						make());
					}
				}
			else
				{
				if(genes.size()>this.max_genes_count) {
					w.add(new VariantContextBuilder(ctx).
							attribute(this.info_tag, "gene_abalation|multiple_genes|"+genes.size()).
							make());
					}
				else
					{
					final Set<String> annots = genes.stream().
							flatMap(G->annotGene(G, ctx).stream()).
							collect(Collectors.toSet());
					w.add(new VariantContextBuilder(ctx).
							attribute(this.info_tag,new ArrayList<>(annots)).
							make()
							);					
					}
				}				
			}
		progress.close();
		return 0;
		} catch(final Throwable err ) {
			LOG.error(err);
			return -1;
		} finally {
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		return doVcfToVcf(args,outputFile);
		}
	
	public static void main(final String[] args)
		{
		new SVPredictions().instanceMainWithExit(args);
		}
	}
