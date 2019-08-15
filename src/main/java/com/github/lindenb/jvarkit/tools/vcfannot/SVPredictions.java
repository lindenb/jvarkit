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
import htsjdk.variant.variantcontext.Allele;
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

First described in http://plindenbaum.blogspot.fr/2011/01/my-tool-to-annotate-vcf-files.html.

### Examples


#### Example 1

```
$java -jar  dist/vcfpredictions.jar  \
	-R human_g1k_v37.fasta \
	-k <(curl  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c | awk -F '   ' '{if($2 ~ ".*_.*") next; OFS="        "; gsub(/chr/,"",$2);print;}'  ) \
	I=~/WORK/variations.gatk.annotations.vcf.gz | bgzip > result.vcf.gz

```

#### Annotate 1000 genomes

```
$  curl "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" |\
gunzip -c |\
java -jar dist/vcfpredictions.jar \
  -R  human_g1k_v37.fasta \
  -k <(curl  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c | awk -F '      ' '{if($2 ~ ".*_.*") next; OFS=" "; gsub(/chr/,"",$2);print;}'  )  |\
cut -d '    ' -f 1-8 | grep -A 10 CHROM


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10583	rs58108140	G	A	100	PASS	AA=.;AC=314;AF=0.14;AFR_AF=0.04;AMR_AF=0.17;AN=2184;ASN_AF=0.13;AVGPOST=0.7707;ERATE=0.0161;EUR_AF=0.21;LDAF=0.2327;PRED=|||||intergenic_variant;RSQ=0.4319;SNPSOURCE=LOWCOV;THETA=0.0046;VT=SNP
1	10611	rs189107123	C	G	100	PASS	AA=.;AC=41;AF=0.02;AFR_AF=0.01;AMR_AF=0.03;AN=2184;ASN_AF=0.01;AVGPOST=0.9330;ERATE=0.0048;EUR_AF=0.02;LDAF=0.0479;PRED=|||||intergenic_variant;RSQ=0.3475;SNPSOURCE=LOWCOV;THETA=0.0077;VT=SNP
1	13302	rs180734498	C	T	100	PASS	AA=.;AC=249;AF=0.11;AFR_AF=0.21;AMR_AF=0.08;AN=2184;ASN_AF=0.02;AVGPOST=0.8895;ERATE=0.0058;EUR_AF=0.14;LDAF=0.1573;PRED=uc010nxq.1|||||intron_variant;RSQ=0.6281;SNPSOURCE=LOWCOV;THETA=0.0048;VT=SNP
1	13327	rs144762171	G	C	100	PASS	AA=.;AC=59;AF=0.03;AFR_AF=0.02;AMR_AF=0.03;AN=2184;ASN_AF=0.02;AVGPOST=0.9698;ERATE=0.0012;EUR_AF=0.04;LDAF=0.0359;PRED=uc010nxq.1|||||intron_variant;RSQ=0.6482;SNPSOURCE=LOWCOV;THETA=0.0204;VT=SNP
1	13957	rs201747181	TC	T	28	PASS	AA=TC;AC=35;AF=0.02;AFR_AF=0.02;AMR_AF=0.02;AN=2184;ASN_AF=0.01;AVGPOST=0.8711;ERATE=0.0065;EUR_AF=0.02;LDAF=0.0788;PRED=uc010nxq.1|||||;RSQ=0.2501;THETA=0.0100;VT=INDEL
1	13980	rs151276478	T	C	100	PASS	AA=.;AC=45;AF=0.02;AFR_AF=0.01;AMR_AF=0.02;AN=2184;ASN_AF=0.02;AVGPOST=0.9221;ERATE=0.0034;EUR_AF=0.02;LDAF=0.0525;PRED=uc010nxq.1|||||3_prime_UTR_variant;RSQ=0.3603;SNPSOURCE=LOWCOV;THETA=0.0139;VT=SNP
1	30923	rs140337953	G	T	100	PASS	AA=T;AC=1584;AF=0.73;AFR_AF=0.48;AMR_AF=0.80;AN=2184;ASN_AF=0.89;AVGPOST=0.7335;ERATE=0.0183;EUR_AF=0.73;LDAF=0.6576;PRED=|||||intergenic_variant;RSQ=0.5481;SNPSOURCE=LOWCOV;THETA=0.0162;VT=SNP
1	46402	rs199681827	C	CTGT	31	PASS	AA=.;AC=8;AF=0.0037;AFR_AF=0.01;AN=2184;ASN_AF=0.0017;AVGPOST=0.8325;ERATE=0.0072;LDAF=0.0903;PRED=|||||intergenic_variant;RSQ=0.0960;THETA=0.0121;VT=INDEL
1	47190	rs200430748	G	GA	192	PASS	AA=G;AC=29;AF=0.01;AFR_AF=0.06;AMR_AF=0.0028;AN=2184;AVGPOST=0.9041;ERATE=0.0041;LDAF=0.0628;PRED=|||||intergenic_variant;RSQ=0.2883;THETA=0.0153;VT=INDEL
1	51476	rs187298206	T	C	100	PASS	AA=C;AC=18;AF=0.01;AFR_AF=0.01;AMR_AF=0.01;AN=2184;ASN_AF=0.01;AVGPOST=0.9819;ERATE=0.0021;EUR_AF=0.01;LDAF=0.0157;PRED=|||||intergenic_variant;RSQ=0.5258;SNPSOURCE=LOWCOV;THETA=0.0103;VT=SNP
```



### See also

 *  http://plindenbaum.blogspot.fr/2013/10/yes-choice-of-transcripts-and-software.html
 *  VcfFilterSequenceOntology
 *  GroupByGene

## History

 *  2013-Dec : moved to a standard arg/argv command line







END_DOC
*/

@Program(name="vcfpredictions",
	description="Basic Variant Effect prediction using gtf",
	keywords={"vcf","annotation","prediction"},
	modificationDate="20190815"
	)
public class SVPredictions extends Launcher
	{
	private static final Logger LOG = Logger.build(SVPredictions.class).make();
	private IntervalTreeMap<Gene> all_gene=null;

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-g","--gtf"},description="GTF File",required=true)
	private Path gtfPath = null;
	@Parameter(names={"-u","--upstream"},description="Upstream size. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int upstream_size =  5_000;
	
	private class Annotation
		{
		Transcript transcript;
		Allele alt2;
		final Set<String> seqont=new HashSet<>();
		String intergenicName="";
		String intergenicID="";

		
		public String toString()
			{
			final StringBuilder b=new StringBuilder();
			
			b.append(alt2==null?"":alt2.getBaseString());
			b.append('|');
			b.append(this.seqont.stream().
					map(T->T.replaceAll("[ ]", "_")).
					collect(Collectors.joining("&")));
			b.append('|');
			b.append("MODIFIER");//impact 
			b.append(transcript==null?this.intergenicID:transcript.getGene().getId());// gene name
			b.append('|');
			b.append(transcript==null?this.intergenicName:transcript.getProperties().getOrDefault("gene_name","")); //gene id
			b.append('|');//
			b.append("transcript");//Feature_Type
			b.append(transcript==null?"":transcript.getId());//Feature_ID
			b.append('|');
			b.append(transcript==null?"":transcript.getProperties().getOrDefault("transcript_biotype",""));
			b.append('|');
			b.append("");//RAnk
			b.append('|');
			b.append("");//HGVS.c 
			b.append('|');
			b.append("");//HGVS.p 
			b.append('|');
			b.append("");//cDNA.pos / cDNA.length 
			b.append('|');
			b.append("");//CDS.pos / CDS.length 
			b.append('|');
			b.append("");//AA.pos / AA.length 
			b.append('|');
			b.append("");//Distance 
			b.append('|');
			b.append("");//ERRORS / WARNINGS / INFO
						
			return b.toString();
			}
		
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
		h2.addMetaDataLine(new VCFInfoHeaderLine("ANN",
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'"
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
			
			
			final List<Annotation> ctx_annotations=new ArrayList<Annotation>();
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
				
				Gene rightGene = null;
				
				for(final Gene g:this.all_gene.getOverlapping(new Interval(
						ctx.getContig(),
						ctx.getStart(),
						Integer.MAX_VALUE
						)))
						{
						if(rightGene==null || rightGene.getStart() > g.getStart()) rightGene=g;
						}
				
				//intergenic
				final Annotation a=new Annotation();
				if(leftGene==null && rightGene==null) {
					
				} else if(leftGene!=null && rightGene!=null) {
					a.intergenicID = leftGene.getId()+"-"+rightGene.getId();
				} else if(leftGene!=null && rightGene==null) {
					a.intergenicID = leftGene.getId()+"-";
					}
				else if(leftGene==null && rightGene!=null) {
					a.intergenicID = "-"+rightGene.getId();
					}
				
				a.seqont.add("intergenic_variant");
				ctx_annotations.add(a);
				}
			else
				{
				for(final Gene g: genes) {				
					for(final Transcript transcript: g.getTranscripts()) {
						final Annotation a=new Annotation();
						
						if(!transcript.overlaps(ctx)) {
							a.seqont.add("upstream_transcript_variant");
							continue;
							}
						
						if(ctx.contains(transcript)) {
							a.seqont.add("transcript_ablation");
							ctx_annotations.add(a);
							continue;
							}
						
						if(transcript.getUTRs().stream().anyMatch(U->U.overlaps(ctx))) {
							a.seqont.add("UTR");
							}
						if(transcript.getExons().stream().anyMatch(U->ctx.contains(U))) {
							a.seqont.add("exon_loss_variant");
							}
						if(transcript.getExons().stream().anyMatch(U->U.overlaps(ctx))) {
							a.seqont.add("UTR");
							}
						if(transcript.getIntrons().stream().anyMatch(U->U.overlaps(ctx))) {
							a.seqont.add("intron");
							}
						if(transcript.getAllCds().stream().anyMatch(U->U.overlaps(ctx))) {
							a.seqont.add("CDS_region");
							}
						}
					}
				
				
				
				}
				
			final Set<String> info= ctx_annotations.stream().map(S->S.toString()).collect(Collectors.toSet());
			final VariantContextBuilder vb=new VariantContextBuilder(ctx);
			vb.attribute("ANN", info.toArray());
			w.add(vb.make());
			}
		progress.close();
		return RETURN_OK;
		} catch(Exception err ) {
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
