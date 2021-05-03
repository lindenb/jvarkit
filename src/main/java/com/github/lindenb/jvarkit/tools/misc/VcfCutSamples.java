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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import htsjdk.variant.vcf.VCFIterator;


/**
 * 
BEGIN_DOC

##Example

```bash
$ curl  "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  java -jar dist/vcfcutsamples.jar  -S M10475 -S M128215 |\
   grep "CHROM" -A 2

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M128215	M10475
chr1	145273345	.	T	C	289.85	.	AC=3;AF=0.38;AN=8;BaseQRankSum=1.062;CSQ=missense_variant|Tct/Cct|S/P|ENSG00000213240|NOTCH2NL|ENST00000369340|4/6|benign(0.238)|tolerated(0.45),missense_variant|Tct/Cct|S/P|ENSG00000213240|NOTCH2NL|ENST00000362074|3/5|benign(0.238)|tolerated(0.45),missense_variant&NMD_transcript_variant|Tct/Cct|S/P|ENSG00000255168||ENST00000468030|3/23|benign(0.416)|tolerated(0.55),missense_variant|Tct/Cct|S/P|ENSG00000213240|NOTCH2NL|ENST00000344859|3/6|possibly_damaging(0.545)|tolerated(0.44);DP=1000;DS;Dels=0.00;EFF=EXON(MODIFIER|||||RP11-458D21.5|nonsense_mediated_decay|NON_CODING|ENST00000468030|3),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Tct/Cct|S67P|230|NOTCH2NL|protein_coding|CODING|ENST00000344859|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Tct/Cct|S67P|236|NOTCH2NL|protein_coding|CODING|ENST00000362074|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Tct/Cct|S67P|236|NOTCH2NL|protein_coding|CODING|ENST00000369340|);FS=3.974;HRun=1;HaplotypeScore=17.4275;MQ=29.25;MQ0=0;MQRankSum=-1.370;QD=0.39;ReadPosRankSum=-1.117	GT:AD:DP:GQ:PL	0/1:215,34:250:99:269,0,3796	0/0:226,22:250:99:0,158,4259
chr1	156011444	.	T	C	2523.46	.	AC=4;AF=0.50;AN=8;BaseQRankSum=-0.490;CSQ=missense_variant|atA/atG|I/M|ENSG00000160803|UBQLN4|ENST00000368309|10/11|benign(0.012)|tolerated(0.3),downstream_gene_variant|||ENSG00000160803|UBQLN4|ENST00000459954|||,missense_variant|Atc/Gtc|I/V|ENSG00000160803|UBQLN4|ENST00000368307|6/7|unknown(0)|tolerated(0.88);DP=204;Dels=0.00;EFF=DOWNSTREAM(MODIFIER|||||UBQLN4|processed_transcript|CODING|ENST00000459954|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Atc/Gtc|I148V|226|UBQLN4|protein_coding|CODING|ENST00000368307|),NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|atA/atG|I495M|601|UBQLN4|protein_coding|CODING|ENST00000368309|);FS=4.328;HRun=0;HaplotypeScore=4.3777;MQ=35.24;MQ0=0;MQRankSum=-0.101;QD=14.93;ReadPosRankSum=1.575	GT:AD:DP:GQ:PL	0/0:34,1:35:69:0,69,717	0/1:24,15:40:99:214,0,443
```

END_DOC
 *
 */
@Program(name="vcfcutsamples",
		description="Select/Exclude some samples from a VCF",
		keywords={"vcf","sample"},
		deprecatedMsg="use bcftools or gatk SelectVariants"
		)
public class VcfCutSamples extends OnePassVcfLauncher {
	private static final Logger LOG=Logger.build(VcfCutSamples.class).make();
	
	/** selected sample */
	@Parameter(names={"-S","--samples"},description="Sample name")
	private Set<String> user_samples=new HashSet<String>();

	@Parameter(names={"-f","--samplefile"},description="read file containing sample names")
	private Path sampleFile=null;
	
	@Parameter(names="--invert",description=" invert selection")
	private boolean invert=false;
	/** remove variant if no call at all */
	@Parameter(names="-r",description="remove variant if no call at all")
	private boolean removeCtxIfNoCall=false;
	
	@Parameter(names="-E",description=" a missing sample is an error")
	private boolean missing_sample_is_error=true;
	
	@ParametersDelegate
	private VariantAttributesRecalculator recalculator = new VariantAttributesRecalculator();

	
	public VcfCutSamples()
		{
		}
	
	@Override
	protected Logger getLogger()
		{
		return LOG;
		}

	@Override
	protected int beforeVcf()
		{
		if( this.sampleFile!=null)
			{
			try(BufferedReader r=IOUtils.openPathForBufferedReading(this.sampleFile)) {
				String line;
				while((line=r.readLine())!=null)
					{
					if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
					this.user_samples.add(line);
					}
				}
			catch(Exception err)
				{
				LOG.error(err);
				return -1;
				}
			}
		return super.beforeVcf();
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in,VariantContextWriter out) {
		VCFHeader header=in.getHeader();
		final Set<String> samples1=new HashSet<String>(header.getSampleNamesInOrder());
		
		for(String my:this.user_samples)
			{
			if(!samples1.contains(my))
				{
				String msg="user sample "+my+" is not present in VCF Header : "+samples1;
				if(this.missing_sample_is_error)
					{
					throw new RuntimeException(msg);
					}
				else
					{
					LOG.warning(msg);
					}
				}
			}
		
		final List<String> samples2=new ArrayList<String>();

		for(final String sample: header.getSampleNamesInOrder())
			{
			if(this.user_samples.contains(sample))
				{
				if(!invert)
					{
					samples2.add(sample);
					}
				}
			else
				{
				if(invert)
					{
					samples2.add(sample);
					}
				}
			}
		
		
		final VCFHeader header2=new VCFHeader(
				header.getMetaDataInInputOrder(),
				samples2
				);
		JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), header2);
		this.recalculator.setHeader(header2);
		out.writeHeader(header2);
		while(in.hasNext())
			{	
			final VariantContext ctx= in.next();
			
			final VariantContextBuilder vb=new VariantContextBuilder(ctx);
			final List<Genotype> genotypes=new ArrayList<Genotype>();
			final Set<Allele> alleles=new HashSet<Allele>();
			boolean only_no_call=true;
			for(final String sample:samples2)
				{
				final Genotype g=ctx.getGenotype(sample);
				if(g.isNoCall() || !g.isCalled()) continue;
				alleles.addAll(g.getAlleles());
				genotypes.add(g);
				if(g.isCalled()) only_no_call=false;
				}
			
			if(removeCtxIfNoCall && only_no_call) continue;
			
			alleles.add(ctx.getReference());
			vb.alleles(alleles);
			vb.genotypes(genotypes);
			out.add(this.recalculator.apply(vb.make()));
			}
		return 0;
		}
	
	
	public static void main(String[] args)
		{
		new VcfCutSamples().instanceMainWithExit(args);
		}
	}
