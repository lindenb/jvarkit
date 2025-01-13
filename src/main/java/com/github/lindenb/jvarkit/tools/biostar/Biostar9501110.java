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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.FindVariantInSamRecord;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFReader;

/** 

BEGIN_DOC

## Example

```bash
java -jar dist/biostar9501110.jar -V src/test/resources/rotavirus_rf.freebayes.vcf.gz src/test/resources/S*.bam


@HD	VN:1.6	GO:none	SO:coordinate
@SQ	SN:RF01	LN:3302
(...)
@CO	biostar9501110. compilation:20211210154516 githash:efc3d5165 htsjdk:2.24.1 date:20211210155226. cmd:-V src/test/resources/rotavirus_rf.freebayes.vcf.gz src/test/resources/S1.bam src/test/resources/S2.bam src/test/resources/S3.bam src/test/resources/S4.bam src/test/resources/S5.bam
RF01_188_587_2:0:0_2:0:0_af	99	RF01	188	60	70M	=	518	400	AGCTCTTAGTTGAATATAGCGATGTTATGGAGAATGCCACACTGTTGTCAATATTCTCGTACTCTTATGA	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S4	NM:i:2	AS:i:60	XS:i:0
RF01_198_667_4:0:0_1:0:0_49	99	RF01	198	60	70M	=	598	470	AGAATATAGCGATGTTATGGAGAATGCCACACTGTTGTCAATATTCTCGAACTCTTATGATAAATATAAG	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S2	NM:i:4	AS:i:58	XS:i:0
RF01_198_667_4:0:0_1:0:0_49	99	RF01	198	60	70M	=	598	470	AGAATATAGCGATGTTATGGAGAATGCCACACTGTTGTCAATATTCTCGAACTCTTATGATAAATATAAG	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S3	NM:i:4	AS:i:58	XS:i:0
RF01_201_685_1:0:0_0:0:0_a8	99	RF01	201	60	70M	=	616	485	ATATAGCGATGTTATGGAGAATGCCACACTGTTGTCAATATTCTCGTACTCTTATGATAAATATAACGCT	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S4	NM:i:1	AS:i:65	XS:i:0
RF01_244_811_1:0:0_2:0:0_4c	163	RF01	244	60	70M	=	742	568	TCGTACTCTTATGATAAATATAACGCTGTTGAAAGGCAATTAGTAAAATATGCAAAAGGTAAGCCGCTAG	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S5b	NM:i:1	AS:i:65	XS:i:0
RF01_257_807_2:0:0_0:0:0_2b	99	RF01	257	60	70M	=	738	551	ATAAATATAACGCTGTTGAAAGGCAATTAGTAAAATATGCAAAAGGTAAGCCGGTAGAAGCAGATTTGAC	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S5	NM:i:2	AS:i:60	XS:i:0
RF01_314_833_2:0:0_0:0:0_84	99	RF01	314	60	70M	=	764	520	AAGCAGATTTGAGAGTGAATGAGTTGGATTATGAAAATAACAAGATAACATCTGAACATTTCCCAACAGC	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S2	NM:i:2	AS:i:60	XS:i:0
RF01_314_833_2:0:0_0:0:0_84	99	RF01	314	60	70M	=	764	520	AAGCAGATTTGAGAGTGAATGAGTTGGATTATGAAAATAACAAGATAACATCTGAACATTTCCCAACAGC	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S3	NM:i:2	AS:i:60	XS:i:0
RF01_329_808_2:0:0_0:0:0_b0	99	RF01	329	60	70M	=	739	480	TGAATGAGTTGGATTATGAAAAAAACAAGATAACATCTGAACTTTTCCCAACAGCAGAGGAATAAACTGA	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S3	NM:i:2	AS:i:60	XS:i:0
RF01_329_808_2:0:0_0:0:0_b0	99	RF01	329	60	70M	=	739	480	TGAATGAGTTGGATTATGAAAAAAACAAGATAACATCTGAACTTTTCCCAACAGCAGAGGAATAAACTGA	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S2	NM:i:2	AS:i:60	XS:i:0
RF01_350_917_6:0:0_1:0:0_a9	99	RF01	350	60	70M	=	848	568	AAAACAAGAAAACATGTGAACTTTTCCGAACAGCAGAGGAATATACTGAATCATTTATGGATCCAGCAAT	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S2	NM:i:6	AS:i:43	XS:i:0
```


## See also

 * BamPhased01



END_DOC

*/

@Program(name="biostar9501110",
description="Keep reads including/excluding variants from VCF",
keywords= {"sam","bam","vcf"},
creationDate="20211210",
modificationDate="20211213",
biostars=9501110,
jvarkit_amalgamion =  true,
menu="Biostars"
)
public class Biostar9501110 extends OnePassBamLauncher
	{
	private static final Logger LOG = Logger.build(Biostar9501110.class).make();
	@Parameter(names={"-clip","--clip"},description="search variant in clipped section of reads")
	protected boolean use_clip =false;
	@Parameter(names={"-V","--variants","--vcf"},description="indexed vcf file",required=true)
	protected Path vcfFile = null;
	@Parameter(names={"--buffer-size"},description=BufferedVCFReader.OPT_BUFFER_DESC,splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int buffSizeInBp = 1_000;
	@Parameter(names={"--inverse"},description="inverse selection. Keep reads that does NOT contain any variant")
	private boolean inverse_selection = false;
	@Parameter(names={"-m","--min-variants"},description="Find a least 'x' variants in each read")
	private int min_num_variants = 1;
	@Parameter(names={"--tag"},description="add attribute with this tag containing the informations about the variants. Empty:ignore")
	private String attribute = "XV";


	private	VCFReader vcfReader = null;
	private BufferedVCFReader bufferedVCFReader = null;
	private final FindVariantInSamRecord findVariantInSamRecord = new FindVariantInSamRecord();
	
	private VariantContext simplify(final VariantContext vc) {
		return new VariantContextBuilder(vc).
				noID().
				noGenotypes().
				passFilters().
				log10PError(VariantContext.NO_LOG10_PERROR).
				attributes(Collections.emptyMap()).
				make();
		}

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeSam()
		{
		if(this.min_num_variants<1) {
			LOG.error("--min-variants < 1");
			return -1;
			}
		if(!StringUtils.isBlank(this.attribute)) {
			if(this.attribute.length()!=2 || !this.attribute.startsWith("X")) {
				LOG.error("attribute should have length==2 and start with X but got "+this.attribute+".");
				return -1;
			}
		}
		this.vcfReader = VCFReaderFactory.makeDefault().open(this.vcfFile, true);
		this.bufferedVCFReader = new BufferedVCFReader(this.vcfReader, this.buffSizeInBp);
		this.bufferedVCFReader.setSimplifier(V->simplify(V));
		this.findVariantInSamRecord.setUseClip(this.use_clip);
		return super.beforeSam();
		}
	
	@Override
	protected void afterSam() {
		CloserUtil.close(this.bufferedVCFReader);
		CloserUtil.close(this.vcfReader);
		super.afterSam();
		}
	
	private boolean findVariants(final SAMRecord record) {
		if (record.getReadUnmappedFlag()) {
			boolean keep =false;
			if(this.inverse_selection) keep = !keep;
			return keep;
			}
		final Locatable recloc = this.use_clip?
				new SimpleInterval(record.getContig(),record.getUnclippedStart(),record.getUnclippedEnd()):
				record
				;
		final Set<String> atts = new HashSet<>();
		int count_variant = 0;
		try(CloseableIterator<VariantContext> iter = this.bufferedVCFReader.query(recloc)) {
			while(iter.hasNext() && count_variant<this.min_num_variants) {
				final VariantContext ctx = iter.next();
				final FindVariantInSamRecord.Match match  = this.findVariantInSamRecord.apply(record, ctx);
				if(match.getAllele().isPresent() && !match.getAllele().get().isReference()) {
					count_variant++;
					if(!StringUtils.isBlank(this.attribute)) {
						char delim = '|';
						final StringBuilder sb=new StringBuilder();
						sb.append(ctx.getStart()).append(delim);
						if(ctx.hasID()) sb.append(ctx.getID()).append(delim);
						sb.append(ctx.getReference().getDisplayString()).append(delim);
						sb.append(match.getAllele().get().getDisplayString());
						atts.add(sb.toString());
						}
					}
				}
			}
		if(!atts.isEmpty()) {
			record.setAttribute(this.attribute, String.join(",", atts));
			}
		boolean keep = count_variant>=min_num_variants;
		if(this.inverse_selection) keep = !keep;
		return keep;
		}
	
	@Override
	protected Function<SAMRecord, List<SAMRecord>> createSAMRecordFunction()
		{
		return (record)->{
			return findVariants(record) ?Collections.singletonList(record):Collections.emptyList();
			};
		}
		
	public static void main(final String[] args) throws IOException
		{
		new Biostar9501110().instanceMainWithExit(args);
		}
		

	}
