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
package com.github.lindenb.jvarkit.tools.vcfregulomedb;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.regulomedb.RegulomeDBTabixAnnotator;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

Build the database, for grch38


```
# here, we use `head` to get a short example
$ wget -q  -O - "https://encode-public.s3.amazonaws.com/2023/03/03/d38f3202-3364-415d-86ed-8690330cb7a2/ENCFF250UJY.tsv" |\
	head -n 1000 | sed 's/^chrom/#chrom/' > regulome.bed
$ bgzip -f  regulome.bed 
$ tabix -f -p bed regulome.bed.gz 
```

for grch37, the database is available at http://legacy.regulomedb.org/downloads/RegulomeDB.dbSNP141.txt.gz .

## Example

```bash
$ wget -q -O - "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz" |\
	bcftools view --types snps |\
	java -jar dist/jvarkit.jar vcfregulomedb --regulomedb regulome.bed.gz |\
	bcftools query -i 'REGULOMEDB>0' -f '%CHROM\t%POS\t%REF\t%ALT\t%REGULOMEDB\n' 
chr1	10177	A	C	0.829
chr1	10177	A	G	0.829
chr1	10181	A	C	0.342
chr1	10181	A	G	0.342
chr1	10181	A	T	0.342
chr1	10248	A	C	0.609
chr1	10248	A	G	0.609
chr1	10248	A	T	0.609
chr1	10250	A	C	0.609
chr1	10250	A	T	0.609
chr1	10255	A	T	0.609
chr1	10257	A	C	0.609
chr1	10327	T	A	0.796
chr1	10327	T	C	0.796
chr1	10327	T	G	0.796


```


END_DOC

 */
@Program(name="vcfregulomedb",
description="Annotate a VCF with the Regulome2 data (https://regulomedb.org/)",
creationDate = "20140709",
modificationDate = "20230512",
keywords = {"vcf","regulomedb"},
menu="VCF Manipulation",
jvarkit_amalgamion = true
)
public class VcfRegulomeDB extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfRegulomeDB.class);

	@Parameter(names={"-b","--bed","--tabix","--regulomedb"},description="RegulomeDB bed sorted, bgzipped and indexed with tabix.",required=true)
	private String bedFile=null;
	@Parameter(names={"-x","--extends"},description="(int) base pairs. look.for data around the variation +/- 'x'. " + DistanceParser.OPT_DESCRIPTION, splitter=NoSplitter.class, converter = DistanceParser.StringConverter.class)
	private int extend = 0;
	@Parameter(names={"-r","--ranking-regex"},description="if defined, only accept the rank matching the regular expression. see https://regulomedb.org/regulome-help/ . For example: 1a	eQTL/caQTL + TF binding + matched TF motif + matched Footprint + chromatin accessibility peak")
	private String acceptRegexStr="";
	
	private RegulomeDBTabixAnnotator annotator=null;
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
 
	@Override
	protected int doVcfToVcf(
			final String inputName, 
			final VCFIterator in,
			final VariantContextWriter out)
		{
		try {
			final VCFHeader header=in.getHeader();
			this.annotator.fillHeader(header);
			JVarkitVersion.getInstance().addMetaData(this, header);
			out.writeHeader(header);
			
			while(in.hasNext())
				{
				for(VariantContext ctx: annotator.annotate(in.next())) {
					out.add(ctx);
					}
				}
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		return 0;
		}
	
		
	@Override
	protected int beforeVcf() {
		if(this.bedFile==null)
			{
			LOG.error("Bed file indexed with tabix is missing");
			return -1;
			}
		
		try
			{
			annotator=new RegulomeDBTabixAnnotator(this.bedFile, 
					AttributeMap.fromPairs(
							RegulomeDBTabixAnnotator.EXTEND_KEY,String.valueOf(this.extend),
							RegulomeDBTabixAnnotator.REGEX_KEY,String.valueOf(this.acceptRegexStr)
					));
			
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		return super.beforeVcf();
		}
	
	@Override
	protected void afterVcf() {
		if(this.annotator!=null) {
			this.annotator.close();
			}
		super.afterVcf();
		}
	
	
	public static void main(final String[] args)
		{
		new VcfRegulomeDB().instanceMainWithExit(args);
		}
	}
