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
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;


/**
BEGIN_DOC

## Example

```
java -Xmx3g -Djava.io.tmpdir=. -jar dist/svpredictions.jar --max-genes 30  --gtf "human.gtf.gz" in.vcf >   out.vcf

more out.vcf
(...)
chr19	54672382	MantaBND:2392:1:3:1:0:0:0	G	[chr9:87618877[G	.	.	BND_PAIR_COUNT=7;CIPOS=-135,135;CLUSTER=CTX3378;IMPRECISE;MATEID=MantaBND:2392:1:3:1:0:0:1;PAIR_COUNT=7;SVCSQ=upstream_transcript_variant|ENSG00000167608|ENST00000416963|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000494594|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000468343|TMC4|protein_coding,exon|ENSG00000167608|ENST00000446291|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000453320|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000414665|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000437868|MBOAT7|protein_coding,intron|ENSG00000167608|ENST00000479750|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000494142|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000391754|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000465790|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000495398|TMC4|protein_coding,exon|ENSG00000167608|ENST00000476013|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000474910|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000449249|MBOAT7|protein_coding,cds|ENSG00000167608|ENST00000376591|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000338624|MBOAT7|protein_coding,cds|ENSG00000167608|ENST00000301187|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000495968|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000497518|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000491216|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000245615|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000449860|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000495279|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000464098|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000431666|MBOAT7|protein_coding;SVTYPE=BND
chr21	10475514	MantaINS:141610:0:0:0:1:0	AG	AAAAAAAAAAAAAAA	.	.	CIGAR=1M14I1D;CLUSTER=CTX3514;DOWNSTREAM_PAIR_COUNT=0;END=10475515;PAIR_COUNT=0;SVCSQ=exon|ENSG00000270533|ENST00000604687|bP-21201H5.1|pseudogene;SVLEN=14;SVTYPE=INS;UPSTREAM_PAIR_COUNT=0
chr22	23478420	MantaDEL:144501:0:1:0:0:0	T	<DEL>	.	.	CIEND=-160,160;CIPOS=-174,175;CLUSTER=CTX3616;DOWNSTREAM_PAIR_COUNT=16;END=23479619;IMPRECISE;PAIR_COUNT=16;SVCSQ=utr&cds&intron&exon|ENSG00000100218|ENST00000406876|RTDR1|protein_coding,intron|ENSG00000100218|ENST00000216036|RTDR1|protein_coding,upstream_transcript_variant|ENSG00000272019|ENST00000606537|Metazoa_SRP|misc_RNA,transcript_ablation|ENSG00000221069|ENST00000408142|AC000029.1|miRNA,intron|ENSG00000100218|ENST00000439064|RTDR1|protein_coding,upstream_transcript_variant|ENSG00000100218|ENST00000421213|RTDR1|protein_coding,utr&intron&exon|ENSG00000100218|ENST00000452757|RTDR1|protein_coding;SVLEN=-1199;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=16
(...)

```


END_DOC
*/

@Program(name="svpredictions",
	description="Basic Variant Effect prediction using gtf",
	keywords={"vcf","annotation","prediction","sv"},
	creationDate="20190815",
	modificationDate="20200623"
	)
public class SVPredictions extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(SVPredictions.class).make();
	private final IntervalTreeMap<Gene> all_gene = new IntervalTreeMap<>();

	private enum WhereInGene {
		intergenic,
		upstream,
		downstream,
		cds,
		utr,
		exon,
		intron,
		start,
		stop
	}
	
	private static class Annotation {
		final String label;
		final Set<WhereInGene> where = new TreeSet<>();
		
		Annotation(final String label) {
			this.label  = label;
		 	}
		}
	
	@Parameter(names={"-g","--gtf"},description=GtfReader.OPT_DESC,required=true)
	private Path gtfPath = null;
	@Parameter(names={"-u","--upstream"},description="Gene Upstream/Downstream length. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int upstream_size =  5_000;
	@Parameter(names={"--max-genes"},description="don't print the genes names if their count exceed 'x'. '-1' = ignore/unlimited")
	private int max_genes_count =  -1;
	@Parameter(names={"-nti","--no-transcript-id"},description="don't print transcript id (reduce length of annotation)")
	private boolean ignore_transcript_id = false;
	@Parameter(names={"--tag"},description="VCF info attribute")
	private String info_tag = "SVCSQ";
	@Parameter(names={"-w","--where"},description="where in gene should overlap the variant. Empty string: no limit/use all possible annotations. "
			+ "Should be a comma/space/semicolon string with the following items: 'intergenic|gene|transcript|intron|exon|utr|cds|downstream|upstream' ")
	private String whereStr="";
	@Parameter(names={"-F","--filter"},description="FILTER to set if variant failing prediction of option --where. Empty: no FILTER, discard variant.")
	private String filterStr="BAD_SV_PRED";
	@Parameter(names={"-r","--remove-attribute"},description="Do not print the annotations that don't contain the contraint for the argument  --where")
	private boolean remove_attributes =  false;
	@Parameter(names={"--bnd"},description="Ignore the INFO/END attribute for SVTYPE=BND, so it is just considered as a single point mutation.")
	private boolean ignore_bnd_end =  false;

	
	private final Set<WhereInGene> limitWhere = new HashSet<>();
	
	
	private Annotation annotTranscript(final Transcript transcript, final VariantContext ctx) {
		final  Annotation ret = new Annotation("|"+transcript.getGene().getId()+"|"+(ignore_transcript_id?"":transcript.getId()+"|")+transcript.getGene().getGeneName()+"|"+transcript.getGene().getGeneBiotype());
		if(!transcript.overlaps(ctx)) {
			
			if(
				(transcript.isNegativeStrand() && transcript.getEnd() < ctx.getStart()) ||
				(ctx.getEnd()  < transcript.getStart() && !transcript.isNegativeStrand() )
				)
				{
				ret.where.add(WhereInGene.upstream);
				return ret;
				}
			else
				{
				ret.where.add(WhereInGene.downstream);
				return ret;
				}
			}
		if(transcript.hasUTR() && transcript.getUTRs().stream().flatMap(U->U.getIntervals().stream()).anyMatch(U->U.overlaps(ctx))) {
			ret.where.add(WhereInGene.utr);
			ret.where.add(WhereInGene.exon);
			}

		if(transcript.hasCDS() && transcript.getAllCds().stream().anyMatch(U->U.overlaps(ctx))) {
			 ret.where.add(WhereInGene.cds);
			 ret.where.add(WhereInGene.exon);
			}
		
		if(transcript.hasExon() && transcript.getExons().stream().anyMatch(U->ctx.overlaps(U))) {
			ret.where.add(WhereInGene.exon);
			}
		if(transcript.hasIntron() && transcript.getIntrons().stream().anyMatch(U->ctx.overlaps(U))) {
			ret.where.add(WhereInGene.intron);
			}
		if(transcript.hasCodonStopDefined() && transcript.getCodonStop().get().getBlocks().stream().anyMatch(U->ctx.overlaps(U))) {
			ret.where.add(WhereInGene.stop);
			}
		if(transcript.hasCodonStartDefined() && transcript.getCodonStart().get().getBlocks().stream().anyMatch(U->ctx.overlaps(U))) {
			ret.where.add(WhereInGene.start);
			}
		return ret;
		}
	
	private List<Annotation> annotGene(final Gene gene, final VariantContext ctx) {
		return gene.getTranscripts().
				stream().
				map(T->annotTranscript(T,ctx)).
				filter(A->!this.remove_attributes || A.where.stream().anyMatch(X->this.limitWhere.contains(X))).
				collect(Collectors.toList())
				;
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator r, VariantContextWriter w) {
		try {
			final VCFHeader header= r.getHeader();
			final SAMSequenceDictionary dict=  header.getSequenceDictionary();
			try(final GtfReader gtfReader=new GtfReader(this.gtfPath)) {
				if(dict!=null) gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
				gtfReader.getAllGenes().stream().forEach(G->this.all_gene.put(new Interval(G), G));
				}
			
			final VCFHeader h2=new VCFHeader(header);
			
			
			/* split 'limit' string */
			if(StringUtils.isBlank(this.whereStr)) {
				//add all
				Arrays.stream(WhereInGene.values()).forEach(C->this.limitWhere.add(C));
			} else {
				for(final String ss: this.whereStr.split("[ \t,;]+")) {
					if(StringUtils.isBlank(ss)) continue;
					// gene and transcript expand to everything but intergenic
					if(ss.toLowerCase().equals("gene") || ss.toLowerCase().equals("transcript")) {
						Arrays.stream(WhereInGene.values()).
							filter(C->!C.equals(WhereInGene.intergenic)).
							forEach(C->this.limitWhere.add(C));
						} 
					else
						{
						final WhereInGene g = Arrays.stream(WhereInGene.values()).
								filter(A->A.name().equalsIgnoreCase(ss)).
								findFirst().
								orElseThrow(()->new IllegalArgumentException(
									"Bad identifier expected one of :" +
									Arrays.stream(WhereInGene.values()).
									map(X->X.name()).
									collect(Collectors.joining(", "))));
						this.limitWhere.add(g);
						}
				}
			if(this.limitWhere.isEmpty()) {
				LOG.error("--where option provided but no identifier was found.");
				return -1;
				}
			}
			
			final VCFFilterHeaderLine filterHeader;
			if(!StringUtils.isBlank(this.filterStr))
				{
				filterHeader = new VCFFilterHeaderLine(this.filterStr, "variant failing locations: " + 
						this.limitWhere.stream().map(V->V.name()).collect(Collectors.joining(","))
						);
				h2.addMetaDataLine(filterHeader);
				}
			else
				{
				filterHeader = null;
				}
	
			
			h2.addMetaDataLine(new VCFInfoHeaderLine(
					this.info_tag,
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Structural variant consequence."
					));
			JVarkitVersion.getInstance().addMetaData(this, h2);			
			w.writeHeader(h2);
	
			while(r.hasNext())
				{
				final VariantContext ctx= r.next();
				
				final int ctx_bnd_end;
				if(this.ignore_bnd_end &&
					ctx.hasAttribute(VCFConstants.SVTYPE) &&
					ctx.getAttributeAsString(VCFConstants.SVTYPE, "").equals("BND")) {
					ctx_bnd_end = ctx.getStart();
					}
				else
					{
					ctx_bnd_end = ctx.getEnd();
					}
				
				final Collection<Gene> genes  = this.all_gene.getOverlapping(new Interval(
						ctx.getContig(),
						Math.max(1,ctx.getStart()-this.upstream_size),
						ctx_bnd_end + this.upstream_size
						));
				
				
				if(genes.isEmpty()) // intergenic
					{
					// discard anyway
					if(!this.limitWhere.contains(WhereInGene.intergenic) && filterHeader==null) continue;
					
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
					
					final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
					// FILTER
					if(!this.limitWhere.contains(WhereInGene.intergenic) && filterHeader!=null) {
						vcb.filter(filterHeader.getID());
						}
					
					
					if(!(leftGene==null && rightGene==null)) {
						vcb.attribute(this.info_tag, "intergenic|"+leftId+"-"+rightId+"|"+leftName+"-"+rightName);
						}
					w.add(vcb.make());
					}
				else
					{
					final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
					final List<List<Annotation>> annotations = genes.stream().
							map(G->annotGene(G, ctx)).
							collect(Collectors.toList());
					
					final boolean match_user_filter = annotations.stream().
							flatMap(L->L.stream()).
							flatMap(A->A.where.stream()).
							anyMatch(A->this.limitWhere.contains(A))
							;
					
					if(!match_user_filter) {
						if(filterHeader==null) continue;
						vcb.filter(filterHeader.getID());
						}
					
					
					
					if(this.max_genes_count!=-1 && genes.size()>this.max_genes_count) {
						
						final String prefix= annotations.stream().
								flatMap(L->L.stream()).
								flatMap(A->A.where.stream()).
								map(A->A.name()).
								collect(Collectors.toCollection(TreeSet::new)).
								stream().
								collect(Collectors.joining("&"));
	
						vcb.attribute(this.info_tag, 
								(StringUtils.isBlank(prefix)?".":prefix)+
								"|multiple_genes|"+genes.size());
						}
					else
						{
						final List<String> annots = annotations.
								stream().
								flatMap(L->L.stream()).
								map(A->A.where.stream().
									map(X->X.name()).
									collect(Collectors.toCollection(TreeSet::new)).
									stream().
									collect(Collectors.joining("&")) + A.label ).
								collect(Collectors.toSet()).
								stream().
								collect(Collectors.toList());
						
						if(!annots.isEmpty()) vcb.attribute(this.info_tag,annots);
						}
					w.add(vcb.make());
					}
				}
			return 0;
		} catch(final Throwable err ) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(w);
			}
		
	}
	
	public static void main(final String[] args)
		{
		new SVPredictions().instanceMainWithExit(args);
		}
	}
