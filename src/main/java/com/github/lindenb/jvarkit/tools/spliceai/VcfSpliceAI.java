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
package com.github.lindenb.jvarkit.tools.spliceai;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.spliceai.SpliceAI;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;


/**
BEGIN_DOC

# Example

```
java -jar dist/jvarkit.jar vcfspliceai --annot /path/to/spliceai_scores.masked.indel.hg38.vcf.gz  src/test/resources/test_vcf01.vcf 

(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
(...)
1	866893	.	T	C	431	PASS	AA=t;AC=7;AF=0.7;AN=10;SpliceAI=SAMD11|0.00|0.00|0.00|0.00|-13|24|-13|-48
1	870317	.	G	A	12	PASS	AC=11;AF=0.917;AN=12;SpliceAI=SAMD11|0.00|0.00|0.00|0.00|2|17|16|-12
1	875770	.	A	G	338	PASS	AA=a;AC=8;AF=0.8;AN=10;SpliceAI=SAMD11|0.00|0.00|0.01|0.00|-1|-45|-50|-46
1	903245	.	A	G	199	PASS	AA=a;AC=6;AF=0.6;AN=10;SpliceAI=PLEKHN1|0.00|0.00|0.00|0.00|48|-37|-22|1
1	905130	.	ATG	A	487	PASS	AC=3;AF=0.5;AN=6;CIGAR=1M2D;IDREP=1;REFREP=2;RU=TG;SpliceAI=PLEKHN1|0.00|0.00|0.00|0.00|-43|21|-33|-37
1	909238	.	G	C	229	PASS	AA=C;AC=8;AF=0.667;AN=12;SpliceAI=PLEKHN1|0.00|0.01|0.00|0.00|-43|-50|39|-7
1	912049	.	T	C	400	PASS	AA=T;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.01|0.01|0.00|-28|-14|-27|-23
1	913889	.	G	A	372	PASS	AA=G;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.01|0.00|0.00|-46|9|2|-45
1	914333	.	C	G	556	PASS	AA=G;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.00|0.00|0.00|-3|27|-3|-38
1	914852	.	G	C	525	PASS	AA=C;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.00|0.00|0.00|22|-22|48|49
1	914940	.	T	C	488	PASS	AA=C;AC=5;AF=0.625;AN=8;SpliceAI=PERM1|0.00|0.00|0.00|0.00|28|-30|-39|3
(...)
```

END_DOC
 */
@Program(name="vcfspliceai",
description="Annotate VCF with local spiceai vcf",
keywords={"vcf","splice","splicing","spliceai"},
creationDate="20201107",
modificationDate="20260602",
jvarkit_amalgamion = true
)
public class VcfSpliceAI  extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfSpliceAI.class);
	
	
	@Parameter(names={"--tag"},description="INFO tag")
	private String tag=SpliceAI.getTag();
	@Parameter(names={"--vcf","--annotation","--spliceai"},description="SpliceAI VCF.vcf.gz indexed with tabix",required = true)
	private Path spliceAiVCF = null;
	@Parameter(names={"--buffer"},description=BufferedVCFReader.OPT_BUFFER_DESC+" "+DistanceParser.OPT_DESCRIPTION,splitter = NoSplitter.class, converter = DistanceParser.StringConverter.class)
	private int distance=1000;

	private VCFFileReader localSpliceAIVcf=null;
	private BufferedVCFReader localSpliceAIBufferedVcf=null;
	
	@Override
	protected Logger getLogger()
		{
		return LOG;
		}

	@Override
	protected int beforeVcf()
		{
		try {
			this.localSpliceAIVcf = new VCFFileReader(this.spliceAiVCF,true);
			final VCFHeader annotHeader = this.localSpliceAIVcf.getHeader();
			final VCFInfoHeaderLine info = annotHeader.getInfoHeaderLine(this.tag);
			if(info==null) throw new IllegalArgumentException("INFO/"+this.tag+" missing in "+this.spliceAiVCF);

			this.localSpliceAIBufferedVcf = new BufferedVCFReader(this.localSpliceAIVcf, this.distance);
			}
		catch(final Throwable err) {
			LOG.error(err);
			}
		return super.beforeVcf();
		}
	
	@Override
	protected void afterVcf() {
		if(this.localSpliceAIBufferedVcf!=null) {
			try {
				this.localSpliceAIBufferedVcf.close();
				this.localSpliceAIBufferedVcf = null;
				}
			catch(Throwable err) {
				LOG.error(err);
				}
			}
		if(this.localSpliceAIVcf!=null) {
			try {
				this.localSpliceAIVcf.close();
				this.localSpliceAIVcf = null;
				}
			catch(Throwable err) {
				LOG.error(err);
				}
			}
		super.afterVcf();
		}
	
	
	
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iterin,
			final VariantContextWriter out)
		{
		final VCFHeader header = iterin.getHeader();
		final VCFHeader annotHeader = this.localSpliceAIBufferedVcf.getHeader();
		final VCFInfoHeaderLine info_hdr = annotHeader.getInfoHeaderLine(this.tag);
		if(info_hdr==null) throw new IllegalArgumentException("INFO/"+this.tag+" missing in "+this.spliceAiVCF);
		final SAMSequenceDictionary annotdict = SequenceDictionaryUtils.extractRequired(annotHeader);
		final ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(annotdict);
		
		
		
		header.addMetaDataLine(info_hdr);
		JVarkitVersion.getInstance().addMetaData(this, header);
		
		out.writeHeader(header);
		while(iterin.hasNext()) {
			final VariantContext ctx= iterin.next();
			final String ctg_annot = ctgConverter.apply(ctx.getContig());
			if(StringUtils.isBlank(ctg_annot)) {
				out.add(ctx);
				continue;
				}
			final Set<SpliceAI> predictions =new HashSet<>();
			// previous values
			if(ctx.hasAttribute(this.tag)) {
				for(SpliceAI prediction: SpliceAI.parse(ctx,tag)) {
					final Allele alt = Allele.create(prediction.getAllele(),false);
					if(!ctx.hasAlternateAllele(alt)) continue;
					predictions.add(prediction);
					}
				}
			try(CloseableIterator<VariantContext> iter2 =  this.localSpliceAIBufferedVcf.query(ctg_annot,ctx.getStart(),ctx.getEnd())) {
				while(iter2.hasNext()) {
					final VariantContext ctx2 = iter2.next();
					if(!ctx2.hasAttribute(tag)) continue;
					if(ctx.getStart()!=ctx2.getStart()) continue;
					if(!ctx.getReference().equals(ctx2.getReference())) continue;
					for(SpliceAI prediction: SpliceAI.parse(ctx2,tag)) {
						final Allele alt = Allele.create(prediction.getAllele(),false);
						if(!ctx.hasAlternateAllele(alt)) continue;
						predictions.add(prediction);
						}
					}
				}
			
			if(predictions.isEmpty()) {
				out.add(ctx);
				}
			else
				{
				out.add(
						new VariantContextBuilder(ctx)
							.attribute(this.tag,predictions.stream()
							.map(P->P.getAttribute())
							.collect(Collectors.toList()))
							.make()
					);
				}
			}
		return 0;
		}
	
	public static void main(final String[] args) {
		new VcfSpliceAI().instanceMainWithExit(args);
		}
	
	}
