/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dgv.DGVBedTabixVariantAnnotator;
import com.github.lindenb.jvarkit.gff3.GffEnsemblRegVariantAnnotator;
import com.github.lindenb.jvarkit.gnomad.GnomadSVBedTabixVariantAnnotator;
import com.github.lindenb.jvarkit.gtf.GtfTabixSVVariantAnnotator;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.regulomedb.RegulomeDBTabixAnnotator;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.vcf.AbstractOnePassVcfAnnotator;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;


/**
BEGIN_DOC

## Example

```
java -Xmx3g -Djava.io.tmpdir=. -jar dist/jvarkit.jar vcfsvannotator --gtf "human.gtf.gz" in.vcf >   out.vcf

more out.vcf
(...)
chr19	54672382	MantaBND:2392:1:3:1:0:0:0	G	[chr9:87618877[G	.	.	BND_PAIR_COUNT=7;CIPOS=-135,135;CLUSTER=CTX3378;IMPRECISE;MATEID=MantaBND:2392:1:3:1:0:0:1;PAIR_COUNT=7;SVCSQ=upstream_transcript_variant|ENSG00000167608|ENST00000416963|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000494594|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000468343|TMC4|protein_coding,exon|ENSG00000167608|ENST00000446291|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000453320|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000414665|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000437868|MBOAT7|protein_coding,intron|ENSG00000167608|ENST00000479750|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000494142|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000391754|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000465790|TMC4|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000495398|TMC4|protein_coding,exon|ENSG00000167608|ENST00000476013|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000474910|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000449249|MBOAT7|protein_coding,cds|ENSG00000167608|ENST00000376591|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000338624|MBOAT7|protein_coding,cds|ENSG00000167608|ENST00000301187|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000495968|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000497518|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000491216|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000245615|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000167608|ENST00000449860|TMC4|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000495279|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000464098|MBOAT7|protein_coding,upstream_transcript_variant|ENSG00000125505|ENST00000431666|MBOAT7|protein_coding;SVTYPE=BND
chr21	10475514	MantaINS:141610:0:0:0:1:0	AG	AAAAAAAAAAAAAAA	.	.	CIGAR=1M14I1D;CLUSTER=CTX3514;DOWNSTREAM_PAIR_COUNT=0;END=10475515;PAIR_COUNT=0;SVCSQ=exon|ENSG00000270533|ENST00000604687|bP-21201H5.1|pseudogene;SVLEN=14;SVTYPE=INS;UPSTREAM_PAIR_COUNT=0
chr22	23478420	MantaDEL:144501:0:1:0:0:0	T	<DEL>	.	.	CIEND=-160,160;CIPOS=-174,175;CLUSTER=CTX3616;DOWNSTREAM_PAIR_COUNT=16;END=23479619;IMPRECISE;PAIR_COUNT=16;SVCSQ=utr&cds&intron&exon|ENSG00000100218|ENST00000406876|RTDR1|protein_coding,intron|ENSG00000100218|ENST00000216036|RTDR1|protein_coding,upstream_transcript_variant|ENSG00000272019|ENST00000606537|Metazoa_SRP|misc_RNA,transcript_ablation|ENSG00000221069|ENST00000408142|AC000029.1|miRNA,intron|ENSG00000100218|ENST00000439064|RTDR1|protein_coding,upstream_transcript_variant|ENSG00000100218|ENST00000421213|RTDR1|protein_coding,utr&intron&exon|ENSG00000100218|ENST00000452757|RTDR1|protein_coding;SVLEN=-1199;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=16
(...)

```


END_DOC
*/

@Program(name="vcfsvannotator",
	description="SV Variant Effect prediction using gtf, gnomad, etc",
	keywords={"vcf","annotation","prediction","sv","gnomad","gtf","regulome"},
	creationDate="20190815",
	modificationDate="20230509",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
	)
public class VCFSVAnnotator extends AbstractOnePassVcfAnnotator
	{
	private static final Logger LOG = Logger.build(VCFSVAnnotator.class).make();

	@Parameter(names={"--gtf"},description=GtfTabixSVVariantAnnotator.OPT_DESC)
	private String gtfPath = null;
	@Parameter(names={"--gnomad"},description=GnomadSVBedTabixVariantAnnotator.OPT_DESC)
	private String gnomadPath = null;
	@Parameter(names={"--dgv"},description=DGVBedTabixVariantAnnotator.OPT_DESC)
	private String dgvPath = null;
	@Parameter(names={"--regulomedb"},description=RegulomeDBTabixAnnotator.OPT_DESC)
	private String regulomePath = null;
	@Parameter(names={"--ensemblreg"},description=GffEnsemblRegVariantAnnotator.OPT_DESC)
	private String ensemblRegGTFPath = null;

	
	

	@SuppressWarnings("serial")
	@DynamicParameter(names={"-D","--define"},description="Dynamic parameters -Dkey=value ."
			+ "'extends': GTF Gene Upstream/Downstream length. "
			+ "'fraction': min common Fraction between two SVs/CNVs."
			)
	private Map<String,String> __dynaParams = new HashMap<String,String>() {{{
		put("extend","1000");
		put("fraction","0.9");
		}}};
	
	
	
	private static class AddSvLen implements VariantAnnotator {
		private final String SVLEN="SVLEN";
		@Override
		public void fillHeader(VCFHeader header) {
			if(header.getInfoHeaderLine(SVLEN)==null) {
				header.addMetaDataLine(new VCFInfoHeaderLine(SVLEN, 1, VCFHeaderLineType.Integer,"SV length"));
				}
			}
		
		@Override
		public void close() {
			}
		
		@Override
		public List<VariantContext> annotate(VariantContext ctx) throws IOException {
			if(ctx.hasAttribute(VCFConstants.SVTYPE) && !ctx.hasAttribute(SVLEN)) {
				int svlen = ctx.getLengthOnReference();
				if(ctx.getAttribute(VCFConstants.SVTYPE,".").equals("DEL")) {
					svlen = -svlen;
					}
				ctx = new VariantContextBuilder(ctx).attribute(SVLEN, svlen).make();
				}
			return Collections.singletonList(ctx);
			}
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected List<VariantAnnotator> createVariantAnnotators() {
		final List<VariantAnnotator> annotators = new ArrayList<>();
		final AttributeMap dynaParams = AttributeMap.wrap(this.__dynaParams);
		final DistanceParser parser=  new DistanceParser();
		try {
			
			annotators.add(new AddSvLen());

			
			if(!StringUtils.isBlank(this.gtfPath)) {
				annotators.add(new GtfTabixSVVariantAnnotator(this.gtfPath,
						AttributeMap.fromPairs(
								GtfTabixSVVariantAnnotator.EXTEND_KEY,
								String.valueOf(parser.applyAsInt(dynaParams.getAttribute("extend", "0")))
						)));
				}
			if(!StringUtils.isBlank(this.gnomadPath)) {
				annotators.add(new GnomadSVBedTabixVariantAnnotator(this.gnomadPath,
						AttributeMap.fromPairs(
								GnomadSVBedTabixVariantAnnotator.FRACTION_KEY,
								String.valueOf(dynaParams.getDoubleAttribute("fraction").orElse(0.9))
						)));
				}
			if(!StringUtils.isBlank(this.dgvPath)) {
				annotators.add(new DGVBedTabixVariantAnnotator(this.dgvPath,
						AttributeMap.fromPairs(
								DGVBedTabixVariantAnnotator.FRACTION_KEY,
								String.valueOf(dynaParams.getDoubleAttribute("fraction").orElse(0.9))
						)));
				}
			if(!StringUtils.isBlank(this.regulomePath)) {
				annotators.add(new RegulomeDBTabixAnnotator(this.regulomePath,
						AttributeMap.fromPairs(
								RegulomeDBTabixAnnotator.EXTEND_KEY,
								String.valueOf("0")
						)));
				}
			if(!StringUtils.isBlank(this.ensemblRegGTFPath)) {
				annotators.add(new GffEnsemblRegVariantAnnotator(this.ensemblRegGTFPath,
						AttributeMap.empty()
						));
				}
			
			
			
			return annotators;
			} 
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	
	
	public static void main(final String[] args)
		{
		new VCFSVAnnotator().instanceMainWithExit(args);
		}
	}
