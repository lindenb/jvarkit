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


History:
* 2014 creation
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.misc;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.tests.TestIt;
import com.github.lindenb.jvarkit.tests.TestRunner;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
 BEGIN_DOC
 

```bash
$ curl -s "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
 java -jar dist/vcfhead.jar -n 2 | grep -v "##"

#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  1094PC0005  1094PC0009  1094PC0012  1094PC0013
chr1    30860   .   G   C   33.46   .   AC=2;AF=0.053;AN=38;BaseQRankSum=2.327;DP=49;Dels=0.00;EFF=DOWNSTREAM(MODIFIER||||85|FAM138A|protein_coding|COD
ING|ENST00000417324|),DOWNSTREAM(MODIFIER|||||FAM138A|processed_transcript|CODING|ENST00000461467|),DOWNSTREAM(MODIFIER|||||MIR1302-10|miRNA|NON_CODING
|ENST00000408384|),INTRON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST00000469289|),INTRON(MODIFIER|||||MIR1302-10|antisense|NON_CODING|ENST000004
73358|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000423562|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING
|ENST00000430492|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000438504|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene
|NON_CODING|ENST00000488147|),UPSTREAM(MODIFIER|||||WASH7P|unprocessed_pseudogene|NON_CODING|ENST00000538476|);FS=3.128;HRun=0;HaplotypeScore=0.6718;In
breedingCoeff=0.1005;MQ=36.55;MQ0=0;MQRankSum=0.217;QD=16.73;ReadPosRankSum=2.017 GT:AD:DP:GQ:PL  0/0:7,0:7:15.04:0,15,177    0/0:2,0:2:3.01:0,3,39   0
/0:6,0:6:12.02:0,12,143    0/0:4,0:4:9.03:0,9,119
chr1    69270   .   A   G   2694.18 .   AC=40;AF=1.000;AN=40;DP=83;Dels=0.00;EFF=SYNONYMOUS_CODING(LOW|SILENT|tcA/tcG|S60|305|OR4F5|protein_coding|CODI
NG|ENST00000335137|exon_1_69091_70008);FS=0.000;HRun=0;HaplotypeScore=0.0000;InbreedingCoeff=-0.0598;MQ=31.06;MQ0=0;QD=32.86 GT:AD:DP:GQ:PL  ./. ./. 1/
1:0,3:3:9.03:106,9,0  1/1:0,6:6:18.05:203,18,0
```
 
 END_DOC
 */
@Program(
		name="vcfhead",
		description="print the first variants of a vcf",
		keywords={"vcf"},
		modificationDate="20200518",
		creationDate="20131210"
		)
public class VcfHead extends OnePassVcfLauncher
	{
	private static final Logger LOG=Logger.build(VcfHead.class).make();
	
	@Parameter(names={"-n","-N","--count"},description="number of variants")
	private long count=10;
	@Parameter(names={"-c","--bycontig"},description="Print first variant for each contig; Implies VCF is sorted",order=1)
	private boolean by_contig=false;
	
	public VcfHead()
		{
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter w) {
		if(count<0) {
			LOG.error("count < 0");
			return -1;
			}
		try {
			final VCFHeader header = in.getHeader();
			
			JVarkitVersion.getInstance().addMetaData(this, header);
			String prev_contig=null;
			long n=0L;
			w.writeHeader(header);
			while(in.hasNext())
				{
				final VariantContext ctx = in.next();
				if(this.by_contig && 
						(prev_contig==null ||
						!prev_contig.equals(ctx.getContig())))
					{
					prev_contig = ctx.getContig();
					n=0L;
					};
					
				if(n< this.count) {
					w.add(ctx);
					++n;
					}
				else if(!this.by_contig )
					{
					break;
					}
				}
			return 0;
		} catch (final Throwable err) {
			LOG.error(err);
			return -1;
			} 
		}
	
	@TestIt
	private static void test01(final TestRunner tr) {
		
	}
	
	public static void main(final String[] args)
		{
		new VcfHead().instanceMainWithExit(args);
		}
	}
