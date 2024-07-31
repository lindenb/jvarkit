/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfmulti2oneallele;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.MultiToOneAllele;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
/*
BEGIN_DOC

## Deprecated

use bcftools norm

## Example

Exac contains multi-ALT  variants:

```bash
$ gunzip -c ExAC.r0.3.sites.vep.vcf.gz | grep rs3828049

1	889238	rs3828049	G	A,C	8422863.10	PASS	AC=6926,3;AC_AFR=220,0;AC_AMR=485,1;AC_Adj=6890,3;AC_EAS=746,0;AC_FIN=259,0;AC_Het=6442,3,0;AC_Hom=224,0;AC_NFE=3856,0;AC_OTH=41,0;AC_SAS=1283,2;AF=0.057,2.472e-05;AN=121358;AN_AFR=10148;AN_AMR=11522;AN_Adj=119272;AN_EAS=8582;AN_FIN=6358;AN_NFE=65282;AN_OTH=876;AN_SAS=16504;(...)

```

processed with this tools:
```
$ java -jar dist/jvarkit.jar vcfmulti2oneallele  ExAC.r0.3.sites.vep.vcf.gz   | grep rs3828049

1	889238	rs3828049	G	A	8422863.10	PASS	AC=6926;AC_AFR=220;AC_AMR=485;AC_Adj=6890;AC_EAS=746;AC_FIN=259;AC_Het=6442;AC_Hom=224;AC_NFE=3856;AC_OTH=41;AC_SAS=1283;AF=0.057;AN=121358;AN_AFR=10148;AN_AMR=11522;AN_Adj=119272;AN_EAS=8582;AN_FIN=6358;AN_NFE=65282;AN_OTH=876;AN_SAS=16504;BaseQRankSum=-2.170e-01;VCF_MULTIALLELIC_SRC=A|C;(...)
1	889238	rs3828049	G	C	8422863.10	PASS	AC=3;AC_AFR=0;AC_AMR=1;AC_Adj=3;AC_EAS=0;AC_FIN=0;AC_Het=3;AC_Hom=0;AC_NFE=0;AC_OTH=0;AC_SAS=2;AF=2.472e-05;AN=121358;AN_AFR=10148;AN_AMR=11522;AN_Adj=119272;AN_EAS=8582;AN_FIN=6358;AN_NFE=65282;AN_OTH=876;AN_SAS=16504;VCF_MULTIALLELIC_SRC=A|C;(....)
```

END_DOC
 */
@Program(
		name="vcfmulti2oneallele",
		description="'one variant with N ALT alleles' to 'N variants with one ALT'",
		deprecatedMsg = "Use bcftools norm",
		modificationDate = "20240731",
		keywords={"vcf"},
		menu="VCF Manipulation",
		jvarkit_amalgamion = true
		)
public class VcfMultiToOneAllele extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VcfMultiToOneAllele.class).make();
	@Parameter(names={"--keepSpanningDeletions"},description="Keep Alt Spanning deletion alleles "+Allele.SPAN_DEL_STRING)
	private boolean keepSpanningDeletions=false;
	@Parameter(names={"--most-frequent"},description="Keep only most frequent allele.")
	private boolean keepOnlyMostFrequentAltAllele = false;
	@Parameter(names={"-flag","--flag"},description="Info field name that will be added to recall the original alleles.")
	private String flagTag ="VCF_MULTIALLELIC_SRC";
	@Parameter(names={"--print-no-alt"},description="Print Variants without ALT allele")
	private boolean addNoAltVariant=false;


	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	public int doVcfToVcf(
			final String inputName,
			final VCFIterator r,
			final VariantContextWriter w
			)   
		{
		final VCFHeader header0 = r.getHeader();
		final VCFHeader header2 = new VCFHeader(header0);
		header2.addMetaDataLine(new VCFInfoHeaderLine(this.flagTag, 1, VCFHeaderLineType.Flag, "Variant was multiallelic"));
		
		JVarkitVersion.getInstance().addMetaData(this, header2);
		final MultiToOneAllele mapper = new MultiToOneAllele(header2).
			setKeepOnlyMostFrequentAltAllele(keepOnlyMostFrequentAltAllele);
		w.writeHeader(header2);
		while(r.hasNext())
			{
			final VariantContext ctx0 =r.next();
			for(VariantContext ctx:mapper.apply(ctx0)) {
				if(ctx.getNAlleles()==1 && !addNoAltVariant) continue;
				if(ctx.getNAlleles()>1) {
					if(!keepSpanningDeletions && ctx.getAlleles().get(1).equals(Allele.SPAN_DEL)) continue;
					}
				
				if(ctx.getNAlleles()!=ctx0.getNAlleles()  && !StringUtils.isDouble(flagTag)) {
					ctx = new VariantContextBuilder(ctx).attribute(this.flagTag, Boolean.TRUE).make();
					}
				
				w.add(ctx);
				}
			}
		return 0;
		}	
	
	public static void main(final String[] args)
		{
		new VcfMultiToOneAllele().instanceMainWithExit(args);
		}
	
	}
