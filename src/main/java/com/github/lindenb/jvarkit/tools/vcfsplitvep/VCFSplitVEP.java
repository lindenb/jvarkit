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
package com.github.lindenb.jvarkit.tools.vcfsplitvep;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.vcf.AbstractOnePassVcfAnnotator;

import htsjdk.samtools.util.RuntimeIOException;

/*
BEGIN_DOC

## Example

```
$ java  -jar dist/vcfpolyx.jar -R reference.fa input.vcf
(...)
2   1133956 .   A   G   2468.84 .   POLYX=23
2   1133956 .   A   AG  3604.25 .   POLYX=23
2   2981671 .   T   G   47.18   .   POLYX=24
(...)
```

## Cited in:

  * "Multiscale heterogeneity in gastric adenocarcinomaevolution is an obstacle to precision medicine" https://assets.researchsquare.com/files/rs-62554/v1/7883b5d6-a5e6-4d39-8554-e9fef719ac42.pdf
  * Maitena Tellaetxe-Abete, Borja Calvo, Charles Lawrie, Ideafix: a decision tree-based method for the refinement of variants in FFPE DNA sequencing data, NAR Genomics and Bioinformatics, Volume 3, Issue 4, December 2021, lqab092, https://doi.org/10.1093/nargab/lqab092
  * Pol32, an accessory subunit of DNA polymerase delta, plays an essential role in genome stability and pathogenesis of Candida albicans. Shraddheya Kumar Patel & al. Gut Microbes. https://doi.org/10.1080/19490976.2022.2163840 2023.
  * Heczko, L., Hlavac, V., Holy, P. et al. Prognostic potential of whole exome sequencing in the clinical management of metachronous colorectal cancer liver metastases. Cancer Cell Int 23, 295 (2023). https://doi.org/10.1186/s12935-023-03135-x
  * Heczko, L., Liska, V., Vycital, O. et al. Targeted panel sequencing of pharmacogenes and oncodrivers in colorectal cancer patients reveals genes with prognostic significance. Hum Genomics 18, 83 (2024). https://doi.org/10.1186/s40246-024-00644-2
  *  Critical roles of Dpb3-Dpb4 sub-complex of DNA polymerase epsilon in DNA replication, genome stability, and pathogenesis of Candida albicans. Bhabasha Gyanadeep Utkalaja, Shraddheya Kumar Patel, Satya Ranjan Sahu, Abinash Dutta, Narottam Acharya.  DOI: 10.1128/mbio.01227-24 . mBio. 2024 Aug 29:e0122724.



END_DOC
*/
@Program(name="vcfsplitvep",
	description="Split CSQ vep annotations",
	keywords={"vcf","vep","annotation"},
	creationDate="20250517",
	modificationDate="20250517",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
	)
public class VCFSplitVEP extends AbstractOnePassVcfAnnotator
	{
	private static final Logger LOG = Logger.of(VCFSplitVEP.class);

	@Parameter(names={"-t","--tag","--tags"},description="VEP tags name{:type{:aggregate}},name2{:type2{:aggregate2}},etc...  where type is a VCF info type Integer,String,Float an aggregate is one of none,min,max,uniq,first")
	private String tags = "";
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	@Override
	protected List<VariantAnnotator> createVariantAnnotators() {
		try {
			final VCFSplitVEPAnnotator ann = new VCFSplitVEPAnnotator();
			ann.setTags(this.tags);
			return Collections.singletonList(ann);
			}
		catch(Throwable err) {
			throw new RuntimeIOException(err);
			}
		}
	

	public static void main(final String[] args)
		{
		new VCFSplitVEP().instanceMainWithExit(args);
		}

}
