/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Building the Tabix File for regulomeDB:
```bash
(for S in 1 2 3 4 5 6 ; do curl -s "http://regulome.stanford.edu/downloads/RegulomeDB.dbSNP132.Category${S}.txt.gz"  | gunzip -c |  cut -f 1,2,5 | sed -e 's/^chrX/23/'  -e 's/^chr//'  | awk -F ' ' '{printf("%s\t%d\t%s\t%s\n",$1,int($2)-1,$2,$3);}' | uniq ; done)| LC_ALL=C sort -t ' ' -k1,1n -k2,2n -k3,3n |  sed 's/^23/X/' | bgzip -c > regulomeDB.bed.gz && tabix  -p bed -f regulomeDB.bed.gz 
```

## Example

```bash
$   curl -kLs "https://raw.githubusercontent.com/arq5x/gemini/master/test/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.snippet.snpEff.vcf" |\
   java -jar dist/vcfstripannot.jar -k '*' |\
   java -jar dist/vcfregulomedb.jar -b regulomeDB.bed.gz -x 10  |\
   grep -i REGU | head

##INFO=<ID=REGULOMEDB,Number=.,Type=String,Description="Format: Position|Distance|Rank">
##VcfRegulomeDBCmdLine=-b regulomeDB.bed.gz -x 10
##VcfRegulomeDBVersion=235a18b083ea15c4ad94060de512f1edc74cec42
1	10583	rs58108140	G	A	100	PASS	REGULOMEDB=10582|0|5
1	13327	rs144762171	G	C	100	PASS	REGULOMEDB=13326|0|6
1	13980	rs151276478	T	C	100	PASS	REGULOMEDB=13971|8|6,13979|0|6,13980|1|6
1	46402	.	C	CTGT	31	PASS	REGULOMEDB=46402|1|6
1	55164	rs3091274	C	A	100	PASS	REGULOMEDB=55163|0|6
1	55299	rs10399749	C	T	100	PASS	REGULOMEDB=55298|0|6
1	55313	rs182462964	A	T	100	PASS	REGULOMEDB=55321|9|6
1	55326	rs3107975	T	C	100	PASS	REGULOMEDB=55321|4|6,55325|0|6
1	55330	rs185215913	G	A	100	PASS	REGULOMEDB=55321|8|6,55325|4|6

```

Those SNPs can be seen at:

* http://regulome.stanford.edu/snp/chr1/10582
* http://regulome.stanford.edu/snp/chr1/13326
* http://regulome.stanford.edu/snp/chr1/13971
* http://regulome.stanford.edu/snp/chr1/46402
* etc...

END_DOC

 */
@Program(name="vcfregulomedb",description="Annotate a VCF with the Regulome data (http://regulome.stanford.edu/")
public class VcfRegulomeDB extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfRegulomeDB.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names="-b",description=" bed indexed with tabix. Format: chrom(tab)start(tab)end(tab)rank",required=true)
	private String bedFile=null;
	@Parameter(names="-T",description="tag in vcf INFO.")
	private String infoTag="REGULOMEDB";
	@Parameter(names="-x",description="(int) base pairs. look.for data around the variation +/- 'x' ")
	private int extend=5;
	@Parameter(names="-r",description="if defined, only accept the rank matching the regular expression ")
	private String acceptRegexStr=null;	
	private Pattern acceptRegex=null;
	
	private RegDataTabixFileReader regDataTabixFileReader=null;
	
	
	private static class RegData
		{
		@SuppressWarnings("unused")
		String chrom;
		int chromSart;
		@SuppressWarnings("unused")
		int chromEnd;
		String rank;
		RegData(String tokens[])
			{
			this.chrom=tokens[0];
			this.chromSart=Integer.parseInt(tokens[1]);
			this.chromEnd=Integer.parseInt(tokens[2]);
			this.rank=tokens[3];
			}
		}
	
	private static class RegDataTabixFileReader
		extends AbstractTabixObjectReader<RegData>
		{
		RegDataTabixFileReader(String uri) throws IOException
			{
			super(uri);
			}
		@Override
		protected Iterator<RegData> iterator(Iterator<String> delegate) {
			return new MyIterator(delegate);
			}
		private class MyIterator
	    	extends AbstractMyIterator
	    	{
	    	private Pattern tab=Pattern.compile("[\t]");
	    	MyIterator(Iterator<String> delegate)
	    		{
	    		super(delegate);
	    		}
	    	
	    	@Override
	    	public RegData next() {
	    		return new RegData(this.tab.split(delegate.next(),5));
	    		}
	    	}

		}
	
	private VcfRegulomeDB()
		{
		
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in,
			VariantContextWriter out)
		{
		
		VCFHeader header=in.getHeader();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		header.addMetaDataLine(new VCFInfoHeaderLine(
				this.infoTag,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				"Format: Position|Distance|Rank"
				));
		
		
		out.writeHeader(header);
		
		while(in.hasNext())
			{
			List<String>  regDataList=new ArrayList<String>();
			VariantContext ctx=in.next();
			
			progress.watch(ctx.getContig(),ctx.getStart());
			
			int start=Math.max(0,ctx.getStart()-this.extend);
			int end=ctx.getEnd()+this.extend;
			
			for(Iterator<RegData> iter=this.regDataTabixFileReader.iterator(ctx.getContig(), start, end);
					iter.hasNext();
					)
				{
				RegData curr=iter.next();
				if(this.acceptRegex!=null && 
				   !this.acceptRegex.matcher(curr.rank).matches()
				   )
					{
					continue;
					}
				String str=
						String.valueOf(curr.chromSart)+"|"+
						String.valueOf(Math.abs(curr.chromSart-(ctx.getStart()-1)))+"|"+
						curr.rank;
				regDataList.add(str);
				}
			if(regDataList.isEmpty())
				{
				out.add(ctx);
				continue;
				}
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.attribute(this.infoTag, regDataList.toArray());
			out.add(vcb.make());
			}
		progress.finish();
		return 0;
		}
	
		

	@Override
	public int doWork(List<String> args)
		{
		
		
		if(bedFile==null)
			{
			LOG.error("Bed file indexed with tabix is missing");
			return -1;
			}
		
		try
			{
			if(this.acceptRegexStr!=null)
				{
				this.acceptRegex=Pattern.compile(this.acceptRegexStr);
				}
			
			LOG.info("Opening "+bedFile);
			this.regDataTabixFileReader=new RegDataTabixFileReader(bedFile);
			return doVcfToVcf(args, outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.regDataTabixFileReader);
			}
		}
	public static void main(String[] args)
		{
		new VcfRegulomeDB().instanceMainWithExit(args);
		}
	}
