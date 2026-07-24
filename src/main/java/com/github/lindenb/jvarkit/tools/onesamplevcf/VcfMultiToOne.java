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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.onesamplevcf;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.UnaryOperator;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/*
BEGIN_DOC

## Input

if there is only one input with the '.list' suffix, it is interpreted as a file containing the path to the vcf files



```bash
$ curl -s "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz" |\
gunzip -c |\
java -jar dist/jvarkit.jar vcfmulti2one  -c -r -a  |\
grep -v '##' |\
grep -E '(CHROM|SAMPLENAME)' | head | verticalize 


>>> 2
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00096;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 2

>>> 3
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00097;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 0|1
<<< 3

>>> 4
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00099;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 0|1
<<< 4

>>> 5
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00100;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 5

>>> 6
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00102;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 6

>>> 7
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00103;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 7

>>> 8
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00105;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 8

>>> 9
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00106;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 9

>>> 10
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00114;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 0|1
<<< 10
```



END_DOC
 */
@Program(name="vcfmulti2one",
	biostars=130456,
	description="Convert VCF with multiple samples to a VCF with one SAMPLE, duplicating variant and adding the sample name in the INFO column. ",
	keywords={"vcf","sample"},
	creationDate="20150312",
	modificationDate="20260724",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfMultiToOne extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.of(VcfMultiToOne.class);

	@Parameter(names={"-c","--nc","-nc","--discard_no_call"},description="discard if variant is no-call")
	private boolean discard_no_call = false;
	@Parameter(names={"-r","--hr","-hr","--discard_hom_ref"},description="discard if variant is hom-ref")
	private boolean discard_hom_ref = false;
	@Parameter(names={"-a","--discard_non_available"},description="discard if variant is not available (see htsjdk definition 'available if the type of this genotype is set')")
	private boolean discard_non_available = false;
	@Parameter(names={"--anonymize","-x"},description="anonymize samples")
	private boolean anonymize_flag = false;

	public static final String DEFAULT_VCF_SAMPLE_NAME="SAMPLE";
	public static final String DEFAULT_SAMPLE_TAGID="SAMPLENAME";
	public static final String SAMPLE_HEADER_DECLARATION="VcfMultiToOne.Sample";
	
	public VcfMultiToOne()
		{
		}
	@Override
	protected Logger getLogger() {
		return  LOG;
		}
		
	/** general utility for program using VCFMulti2One:
	 *  Extract SampleNames
	 */
	static Set<String> extractSampleNames(final VCFHeader header)
		{
		final List<String> sample_list =header.getSampleNamesInOrder();
		if(sample_list.size()!=1 || !sample_list.get(0).equals(DEFAULT_VCF_SAMPLE_NAME))
			{
			throw new IllegalArgumentException("Not a VCF produced by VcfMultiToOne");
			}
		final Set<String> samples = new TreeSet<String>();
		for(final VCFHeaderLine h:header.getMetaDataInInputOrder())
			{
			if(h.getKey().equals(SAMPLE_HEADER_DECLARATION))
				{
				sample_list.add(h.getValue());
				}
			}
		return samples;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter out) {
		final VCFHeader h0 = in.getHeader();	
		final Set<VCFHeaderLine> metaData = new LinkedHashSet<VCFHeaderLine>(h0.getMetaDataInInputOrder());
		final UnaryOperator<String> rename;
		
		if(this.anonymize_flag) {
			rename= S->StringUtils.md5(S);
			}
		else
			{
			rename= S->S;
			}
		//addMetaData(metaData);
		metaData.add(new VCFInfoHeaderLine(
				DEFAULT_SAMPLE_TAGID,1,VCFHeaderLineType.String,
				"Sample Name from multi-sample vcf"
				));
			
			for(final String sample:h0.getSampleNamesInOrder())
				{
				metaData.add(
					new VCFHeaderLine(
					SAMPLE_HEADER_DECLARATION,
					rename.apply(sample)));
				}
			
			final VCFHeader h2 =  new VCFHeader(
					metaData,
					Collections.singleton(DEFAULT_VCF_SAMPLE_NAME)
					);
			JVarkitVersion.getInstance().addMetaData(this, h2);
			out.writeHeader(h2);
			while(in.hasNext()) {
				final VariantContext ctx = in.next();
				// no genotype
				if(!ctx.hasGenotypes())
					{
					if(!this.discard_no_call)
						{
						final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
						vcb.genotypes(GenotypeBuilder.createMissing(DEFAULT_VCF_SAMPLE_NAME,2));
						out.add(vcb.make());
						}
					continue;
					}
				//loop over samples
				for(Genotype g: ctx.getGenotypes())
					{
					
					if(this.discard_no_call && g.getAlleles().stream().allMatch(A->A.isNoCall())) continue;
					if(!g.isAvailable() && this.discard_non_available) continue;
					if(this.discard_hom_ref && g.getAlleles().stream().allMatch(A->A.isReference())) continue;
					
					
					final GenotypeBuilder gb=new GenotypeBuilder(g);
					gb.name(DEFAULT_VCF_SAMPLE_NAME);
					
					
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.attribute(DEFAULT_SAMPLE_TAGID, rename.apply( g.getSampleName()));
					
					vcb.genotypes(Collections.singletonList( gb.make()));
					out.add(vcb.make());
					}
				} //end while vcfiterator
				
			
			return 0;
			}
	
	
	
	public static void main(final String[] args)
		{
		new VcfMultiToOne().instanceMainWithExit(args);
		}
	}
