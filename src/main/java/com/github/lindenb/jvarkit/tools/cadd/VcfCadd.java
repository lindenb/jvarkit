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
package com.github.lindenb.jvarkit.tools.cadd;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.vcf.AbstractOnePassVcfAnnotator;

import htsjdk.samtools.util.RuntimeIOException;

/**
BEGIN_DOC
## Example

```bash
$ java -Dhttp.proxyHost=my.proxy.host.fr -Dhttp.proxyPort=1234 -jar dist/jvarkit.jar vcfcadd \
	-u "http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz"  \
	src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz 2> /dev/null | ~/package/bcftools/bcftools annotate -x '^INFO/CADD_SCORE,INFO/CADD_PHRED'

##fileformat=VCFv4.2
(...)
##INFO=<ID=CADD_PHRED,Number=A,Type=Float,Description="PHRED expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD scores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc.  URI was http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz">
##INFO=<ID=CADD_SCORE,Number=A,Type=Float,Description="Score suggests that that variant is likely to be  observed (negative values) vs simulated(positive values).However, raw values do have relative meaning, with higher values indicating that a variant is more likely to be simulated (or -not observed-) and therefore more likely to have deleterious effects. URI was http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz">
##VcfCaddCmdLine=-u http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	905606	rs540662886	G	C,A	41743.9	PASS	CADD_PHRED=3.426,.;CADD_SCORE=0.082875,.
(...)
1	905621	rs368876607	G	A	14291.5	PASS	CADD_PHRED=6.025;CADD_SCORE=0.334762
(...)
1	905669	rs111483874	C	G,T	86574.3	PASS	CADD_PHRED=12.77,.;CADD_SCORE=1.39614,.
(...)
1	905723	rs150703609	G	A	15622.1	PASS	CADD_PHRED=23.7;CADD_SCORE=4.05532
1	905726	rs751084833	C	T,A	8733.36	PASS	.
1	905727	rs761609807	G	A	12936.9	PASS	.
(..)
```
## Note to self

I got problem with the certificate. Fixed with `-Dcom.sun.security.enableAIAcaIssuers=true -Dcom.sun.net.ssl.checkRevocation=false `

END_DOC
*/
@Program(
	name="vcfcadd",
	description= "Annotate VCF with  Combined Annotation Dependent Depletion (CADD) (Kircher & al. "+
		"A general framework for estimating the relative pathogenicity of human genetic variants. "+
		"Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892." +
		"PubMed PMID: 24487276.",
	creationDate="20220119",
	modificationDate="20240524",
	keywords={"vcf","prediction","cadd","annotation"},
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)

public class VcfCadd extends AbstractOnePassVcfAnnotator
	{
	private static final Logger LOG = Logger.build(VcfCadd.class).make();
	@Parameter(names={"-u","--uri","--tabix","--cadd"},description="Combined Annotation Dependent Depletion (CADD) Tabix file URI ",required = true)
	private String ccaduri="";
	@Parameter(names={"-d","--buffer-size"},description="Buffer size / processing window size. " + DistanceParser.OPT_DESCRIPTION , converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int buffer_distance= 1_000;
	@Parameter(names={"-na"},description="value  used for 'allele-not-found'.")
	private float NA_value = -999f;

	
	public VcfCadd() {
		}
	
	@Override
	protected List<VariantAnnotator> createVariantAnnotators() {
		try {
			final CaddVariantAnnotator ann= new CaddVariantAnnotator(this.ccaduri);
			ann.setBufferDistance(this.buffer_distance);
			ann.setNAValue(this.NA_value);			
			return Collections.singletonList(ann);
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	@Override
	protected Logger getLogger() {
		return LOG;
		}

	public static void main(final String[] args)
		{
		new VcfCadd().instanceMainWithExit(args);
		}
	}
