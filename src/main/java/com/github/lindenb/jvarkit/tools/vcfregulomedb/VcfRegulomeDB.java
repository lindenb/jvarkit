/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;
import java.util.regex.Pattern;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
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
modificationDate = "20230505",
keywords = {"vcf","regulomedb"},
menu="VCF Manipulation",
jvarkit_amalgamion = true
)
public class VcfRegulomeDB extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(VcfRegulomeDB.class).make();

	@Parameter(names={"-b","--bed","--tabix","--regulomedb"},description="RegulomeDB bed sorted, bgzipped and indexed with tabix.",required=true)
	private String bedFile=null;
	@Parameter(names="-T",description="tag in vcf INFO.")
	private String infoTag="REGULOMEDB";
	@Parameter(names={"-x","--extends"},description="(int) base pairs. look.for data around the variation +/- 'x' ",hidden = true)
	private int extend = 0;
	@Parameter(names={"-r","--ranking-regex"},description="if defined, only accept the rank matching the regular expression. see https://regulomedb.org/regulome-help/ . For example: 1a	eQTL/caQTL + TF binding + matched TF motif + matched Footprint + chromatin accessibility peak")
	private String acceptRegexStr=null;
	
	private Pattern acceptRegex=null;
	private FileHeader fileHeader = null;
	private TabixReader regDataTabixFileReader=null;
	private int database_version=-1; /* 1 : grch37 , 2 : grch38 */
	
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
			final ContigNameConverter ctgConvert = ContigNameConverter.fromContigSet(this.regDataTabixFileReader.getChromosomes());
			final VCFHeader header=in.getHeader();
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict!=null && SequenceDictionaryUtils.isGRCh37(dict) && database_version==2) {
				throw new IllegalArgumentException("VCF Sequence dictionary looks like grch37 and database grch38");
				}
			else if(dict!=null && SequenceDictionaryUtils.isGRCh38(dict) && database_version==1) {
				throw new IllegalArgumentException("VCF Sequence dictionary looks like grch38 and database grch37");
				}
			
			final VCFInfoHeaderLine infoLine = new VCFInfoHeaderLine(
					this.infoTag,
					1,
					VCFHeaderLineType.Float,
					"Mean probability_score for regulomedb "+this.bedFile + " ranking-regex:"+(this.acceptRegexStr)+
					". The scoring scheme refers to the following supporting evidence for that particular location or variant id. "
					+ "In general, if more supporting data is available, the higher is its likelihood of being functional and hence receives a higher score (with 1 being higher and 7 being lower score)."
					);
			header.addMetaDataLine(infoLine);
			
			JVarkitVersion.getInstance().addMetaData(this, header);
			out.writeHeader(header);
			
			while(in.hasNext())
				{
				final VariantContext ctx=in.next();
				final String ctg = ctgConvert.apply(ctx.getContig());
				
				final int start=Math.max(0,ctx.getStart()-this.extend);
				final int end=ctx.getEnd()+this.extend;
				
				final TabixReader.Iterator iter = this.regDataTabixFileReader.query(StringUtils.isBlank(ctg)?ctx.getContig():ctg, Math.max(start-1,0), end);
				
				double sum_probability_score = 0.0;
				int count_probability_score = 0;
				for(;;)
					{
					final String line = iter.next();
					if(line==null) break;
					final Map<String,String> curr = this.fileHeader.toMap(CharSplitter.TAB.split(line));
					if(this.acceptRegex!=null && 
					   !this.acceptRegex.matcher(curr.get("ranking")).matches()
					   )
						{
						continue;
						}
					final int x0 = Integer.parseInt(curr.get("start")+1);
					final int x1 = Integer.parseInt(curr.get("end"));
					if(!CoordMath.overlaps(start, end, x0, x1)) continue;
					
					final double probability_score = Double.parseDouble(curr.get("probability_score"));
					sum_probability_score+=probability_score;
					count_probability_score++;
					}
				if(count_probability_score==0)
					{
					out.add(ctx);
					continue;
					}
				out.add(new VariantContextBuilder(ctx).
						attribute(this.infoTag,sum_probability_score/count_probability_score).
						make());
				}
			}
		catch(IOException err) {
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
			if(!StringUtils.isBlank(this.acceptRegexStr))
				{
				this.acceptRegex=Pattern.compile(this.acceptRegexStr);
				}
			
			this.regDataTabixFileReader=new TabixReader(bedFile);
			final String line = this.regDataTabixFileReader.readLine();
			if(line==null) throw new IOException("Cannot read header of "+ this.bedFile);
			if(line.startsWith("chr1")) {
				this.fileHeader = new FileHeader(Arrays.asList(
						"col1", "start", "end", "col4", "col5", "col6",
						"col7", "col8", "col9", "col10",
						"col11", "col12", "col13", "col14",
						"col15", "col16", "col17", "col18",
						"col19", "col20", "col21", "col22",
						"col23"));
				this.database_version = 1;
				/* ok, let's ignore this for now, I cannot find an equivalent description of regulomedb2.2 */
				LOG.error("cannot decode regulomedb header :"+line);
				return -1;
				}
			else
				{				
				this.fileHeader = new FileHeader(CharSplitter.TAB.split(line));
				this.fileHeader.getColumnIndex("ranking");
				this.fileHeader.getColumnIndex("probability_score");
				this.database_version = 2;
				}
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
		if(this.regDataTabixFileReader!=null) {
			this.regDataTabixFileReader.close();
			}
		super.afterVcf();
		}
	
	
	public static void main(final String[] args)
		{
		new VcfRegulomeDB().instanceMainWithExit(args);
		}
	}
