/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
BEGIN_DOC
## Example

```bash
$  curl -s  "https://raw2.github.com/arq5x/gemini/master/test/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.snippet.vcf" | \
 java -jar dist/vcfcadd.jar \
      -u "http://krishna.gs.washington.edu/download/CADD/v1.2/whole_genome_SNVs.tsv.gz" |\
 grep -E '(CADD|#)'

(...)
##INFO=<ID=CADD,Number=.,Type=String,Description="(Allele|Score|Phred) Score suggests that that variant is likely to be  observed (negative values) vs simulated(positive 
values).However, raw values do have relative meaning, with higher values indicating that a variant is more likely to be simulated (or -not observed-) and therefore more l
ikely to have deleterious effects.PHRED expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD sc
ores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc">
(..)
##VcfCaddCmdLine=-u http://krishna.gs.washington.edu/download/CADD/v1.0/1000G.tsv.gz
##VcfCaddVersion=6123910f68df940c1f3986d142f9b0414f76a43a
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	1308871	.	A	T	6.20	.	AC1=2;AF1=1;CADD=T|-0.549120|0.089;...
1	1308982	.	A	G	6.98	.	AC1=2;AF1=1;CADD=G|-0.000088|3.329;...
1	1657021	.	T	C	3.02	.	AC1=2;AF1=1;CADD=C|-0.271229|0.740;...
(..)
```
END_DOC
*/
@Program(name="vcfcadd",description= "Annotate VCF with  Combined Annotation Dependent Depletion (CADD) (Kircher & al. "+
		"A general framework for estimating the relative pathogenicity of human genetic variants. "+
		"Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892." +
		"PubMed PMID: 24487276.")

public class VcfCadd extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfCadd.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	public static final String DEFAULT_URI="http://krishna.gs.washington.edu/download/CADD/v1.2/whole_genome_SNVs.tsv.gz";
	private TabixFileReader tabix=null;
	@Parameter(names="-u",description="Combined Annotation Dependent Depletion (CADD) Tabix file URI ")
	private String ccaduri=DEFAULT_URI;
	private Pattern TAB=Pattern.compile("[\t]");
	private static final String CADD_FLAG="CADD";
	@Parameter(names="-d",description="processing window size")
	private int buffer_distance=1000;
	
	
	private static class Record
		{
		int pos;
		Allele ref;
		Allele alt;
		float score;
		float phred;
		}
	
	public VcfCadd()
		{
		}
	
	
	private void runTabix(List<VariantContext> buffer)
			throws IOException
		{
		if(buffer.isEmpty()) return;
		String chrom=buffer.get(0).getContig();
		int start=Integer.MAX_VALUE;
		int end=0;
		for(VariantContext ctx:buffer)
			{
			start= Math.min(start, ctx.getStart());
			end= Math.max(end, ctx.getEnd());
			}
		LOG.info(chrom+":"+start+"-"+end);
		List<Record> caddList=new ArrayList<Record>(buffer.size());
		
		Iterator<String> iter = tabix.iterator(
				chrom,
				(int)Math.max(1,start),
				(int)Math.min(Integer.MAX_VALUE,end)
				);
		while(iter.hasNext() )
			{
			String line=iter.next();
			String tokens[]=TAB.split(line);
			if(tokens.length!=6) throw new IOException("Bad CADD line . Expected 6 fields:"+line);
			Record rec=new Record();
			rec.pos= Integer.parseInt(tokens[1]);
			rec.ref=Allele.create(tokens[2],true);
			rec.alt=Allele.create(tokens[3],false);
			rec.score = Float.parseFloat(tokens[4]);
			rec.phred = Float.parseFloat(tokens[5]);
			caddList.add(rec);
			}
		CloserUtil.close(iter);
		
		
		for(int i=0;i< buffer.size();++i)
			{
			List<String> cadd_array=new ArrayList<>();
			VariantContext ctx=buffer.get(i);
			for(Record rec:caddList)
				{
				if(rec.pos!=ctx.getStart()) continue;
				if(!ctx.getReference().equals(rec.ref)) continue;

				
				for(Allele alt:ctx.getAlternateAlleles())
					{
					if(alt.isSymbolic() || !alt.equals(rec.alt)) continue;
					cadd_array.add(alt.getDisplayString()+"|"+rec.score+"|"+rec.phred);
					}
				}
			
			if(!cadd_array.isEmpty())
				{
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(CADD_FLAG, cadd_array);
				buffer.set(i, vcb.make());
				}
			}

		
		}
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		
			try {
			VCFHeader header=in.getHeader();
			if(header.getSequenceDictionary()!=null)
				{
				SAMSequenceDictionary dict=header.getSequenceDictionary();
				Set<String> vcfchr=new HashSet<String>();
				for(SAMSequenceRecord ssr:dict.getSequences()) vcfchr.add(ssr.getSequenceName());
				if(!vcfchr.retainAll(this.tabix.getChromosomes()))//nothing changed
					{
					LOG.warning("#### !!!! NO common chromosomes between tabix and vcf file. Check chromosome 'chr' prefix ? tabix chroms:"+this.tabix.getChromosomes());
					}
				}
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			
			header.addMetaDataLine(new VCFInfoHeaderLine(CADD_FLAG,VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					"(Allele|Score|Phred) Score suggests that that variant is likely to be  observed (negative values) vs simulated(positive values)."+
					"However, raw values do have relative meaning, with higher values indicating that a variant is more likely to be simulated (or -not observed-) and therefore more likely to have deleterious effects."
					+ "PHRED expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD scores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc"
							));
			 
			
			out.writeHeader(header);
			List<VariantContext> buffer= new ArrayList<>();
			for(;;)
				{	
				VariantContext ctx=null;
				if(in.hasNext())
					{
					ctx=progress.watch(in.next());
					}
				if( ctx==null ||
					(!buffer.isEmpty() &&
					(buffer.get(0).getContig().equals(ctx.getContig())) && (ctx.getEnd()-buffer.get(0).getStart())>buffer_distance))
					{
					if(!buffer.isEmpty())
						{
						runTabix(buffer);
						for(VariantContext c:buffer)
							{
							out.add(c);
							}
						}
				
					if(ctx==null) break;
					buffer.clear();
					}
			
				buffer.add(ctx);
				}
			progress.finish();
			return 0;
			} 
		catch(Exception err) {
			LOG.error(err);
			return -1;
			}

		}
	@Override
	public int doWork(List<String> args) {
		try
			{
			
			if(this.ccaduri==null || !this.ccaduri.endsWith(".gz"))
				{
				LOG.error("CCAD uri should end with gz. got "+this.ccaduri);
				return -1;
				}
			
			LOG.info("Loading index for "+this.ccaduri+". Please wait...");
			this.tabix=new TabixFileReader(this.ccaduri);
			LOG.info("End loading index");
			
			return doVcfToVcf(args,outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(tabix);
			tabix=null;
			}
		}

	public static void main(String[] args)
		{
		new VcfCadd().instanceMainWithExit(args);
		}
	}
