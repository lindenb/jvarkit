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
package com.github.lindenb.jvarkit.tools.vcfmulti2oneinfo;

import htsjdk.samtools.util.StringUtil;

import java.util.Collections;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.vcf.AbstractOnePassVcfAnnotator;

/**
BEGIN_DOC

before:

```
RF11	74	.	CAAAAA	CAA	113	.	AC=1;AN=10;ANN=CAA|frameshift_variant&start_lost|HIGH|Gene_78_374|Gene_78_374|transcript|AAG15312.1|protein_coding|1/1|c.-2_1delAAA|p.Met1fs||1/297|1/98||,CAA|disruptive_inframe_deletion|MODERATE|Gene_20_616|Gene_20_616|transcript|AAG15311.1|protein_coding|1/1|c.57_59delAAA|p.Lys19del|57/597|57/597|19/198||,CAA|upstream_gene_variant|MODIFIER|Gene_78_374|Gene_78_374|transcript|AAG15312.1|protein_coding|1/1|c.-2_1delAAA|||||2|;DP=82;DP4=55,0,5,0;HOB=0.02;ICB=0.0439024;IDV=7;IMF=0.388889;INDEL;LOF=(Gene_78_374|Gene_78_374|1|1.00);MQ=60;MQ0F=0;SGB=5.5074;VDB=0.00911888	GT:PL	0/0:0,33,246	0/0:0,36,255	0/0:0,36,255	0/1:149,0,205	0/0:0,33,244
```

run:

```
java -jar dist/jvarkit.jar  vcfmulti2oneinfo -i ANN src/test/resources/rotavirus_rf.ann.vcf.gz
```

after:

```
RF11	74	.	CAAAAA	CAA	113	.	AC=1;AN=10;ANN=CAA|frameshift_variant&start_lost|HIGH|Gene_78_374|Gene_78_374|transcript|AAG15312.1|protein_coding|1/1|c.-2_1delAAA|p.Met1fs||1/297|1/98||;DP=82;DP4=55,0,5,0;HOB=0.02;ICB=0.0439024;IDV=7;IMF=0.388889;INDEL;LOF=(Gene_78_374|Gene_78_374|1|1.00);MQ=60;MQ0F=0;SGB=5.5074;VDB=0.00911888	GT:PL	0/0:0,33,246	0/0:0,36,255	0/0:0,36,255	0/1:149,0,205	0/0:0,33,244
RF11	74	.	CAAAAA	CAA	113	.	AC=1;AN=10;ANN=CAA|disruptive_inframe_deletion|MODERATE|Gene_20_616|Gene_20_616|transcript|AAG15311.1|protein_coding|1/1|c.57_59delAAA|p.Lys19del|57/597|57/597|19/198||;DP=82;DP4=55,0,5,0;HOB=0.02;ICB=0.0439024;IDV=7;IMF=0.388889;INDEL;LOF=(Gene_78_374|Gene_78_374|1|1.00);MQ=60;MQ0F=0;SGB=5.5074;VDB=0.00911888	GT:PL	0/0:0,33,246	0/0:0,36,255	0/0:0,36,255	0/1:149,0,205	0/0:0,33,244
RF11	74	.	CAAAAA	CAA	113	.	AC=1;AN=10;ANN=CAA|upstream_gene_variant|MODIFIER|Gene_78_374|Gene_78_374|transcript|AAG15312.1|protein_coding|1/1|c.-2_1delAAA|||||2|;DP=82;DP4=55,0,5,0;HOB=0.02;ICB=0.0439024;IDV=7;IMF=0.388889;INDEL;LOF=(Gene_78_374|Gene_78_374|1|1.00);MQ=60;MQ0F=0;SGB=5.5074;VDB=0.00911888	GT:PL	0/0:0,33,246	0/0:0,36,255	0/0:0,36,255	0/1:149,0,205	0/0:0,33,244
```

END_DOC
 */
@Program(
		name="vcfmulti2oneinfo",
		description="'one variant with INFO with N values' to 'N variants with one INFO'",
		keywords={"vcf"},
		creationDate = "20260106",
		modificationDate = "20230524",
		jvarkit_amalgamion = true,
		menu="VCF Manipulation"
		)
public class VcfMultiToOneInfo
	extends AbstractOnePassVcfAnnotator
	{
	private static final Logger LOG = Logger.of(VcfMultiToOneInfo.class);

	@Parameter(names={"-i","--info"},description="The INFO tag",required=true)
	private String infoTag = null;
	
	@Override
	protected List<VariantAnnotator> createVariantAnnotators() {
		return Collections.singletonList(new MultiToOneInfoVariantAnnotator(this.infoTag));
		}
	
	 @Override
	protected int beforeVcf() {
		 if(StringUtil.isBlank(this.infoTag)) {
			LOG.error("No info tag defined");
			return -1;
			}
		return super.beforeVcf();
		}
	
	 @Override
	protected Logger getLogger() {
		 return LOG;
		 }
	 
	public static void main(final String[] args)
		{
		new VcfMultiToOneInfo().instanceMainWithExit(args);
		}
	
	}
