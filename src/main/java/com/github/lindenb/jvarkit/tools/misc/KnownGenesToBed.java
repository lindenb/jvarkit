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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.List;
import java.util.function.Predicate;

import org.apache.commons.jexl2.JexlContext;

import htsjdk.samtools.util.CloserUtil;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jexl.JexlPredicate;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**

BEGIN_DOC

### Example

```
$ curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" |\
  gunzip -c |\
  java -jar dist/kg2bed.jar
chr1	11873	14409	+	uc001aaa.3	TRANSCRIPT	uc001aaa.3
chr1	11873	12227	+	uc001aaa.3	EXON	Exon 1
chr1	12227	12612	+	uc001aaa.3	INTRON	Intron 1
chr1	11873	12227	+	uc001aaa.3	UTR	UTR3
chr1	12612	12721	+	uc001aaa.3	EXON	Exon 2
chr1	12721	13220	+	uc001aaa.3	INTRON	Intron 2
chr1	12612	12721	+	uc001aaa.3	UTR	UTR3
chr1	13220	14409	+	uc001aaa.3	EXON	Exon 3
chr1	13220	14409	+	uc001aaa.3	UTR	UTR3
chr1	11873	14409	+	uc010nxr.1	TRANSCRIPT	uc010nxr.1
chr1	11873	12227	+	uc010nxr.1	EXON	Exon 1
chr1	12227	12645	+	uc010nxr.1	INTRON	Intron 1
chr1	11873	12227	+	uc010nxr.1	UTR	UTR3
chr1	12645	12697	+	uc010nxr.1	EXON	Exon 2
chr1	12697	13220	+	uc010nxr.1	INTRON	Intron 2
```

END_DOC
*/


@Program(name="kg2bed",
	description="converts UCSC knownGenes file to BED.",
	keywords={"ucsc","bed","knownGenes"},
	biostars=151628,
	modificationDate="20190427"
	)
public class KnownGenesToBed extends Launcher
	{
	private static final Logger LOG = Logger.build(KnownGenesToBed.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-hide","--hide"},description="don't show the following items (comma separated, one of 'INTRON,UTR,CDS,EXON,TRANSCRIPT,NON_CODING,CODING'). Empty don't hide anything")
	private String hideStr="";
	@Parameter(names={"-s","--select"},description="JEXL select expression. Object 'kg' is an instance of KnownGene (https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/ucsc/KnownGene.java)." +JexlPredicate.OPT_WHAT_IS_JEXL)
	private String selectExpr="";

	private PrintStream out;
	
	private static class KgContext implements JexlContext
		{
		private final KnownGene kg;
		KgContext(final KnownGene kg) {
			this.kg=kg;
			}
		@Override
		public Object get(final String s) {
			if(s.equals("kg") || s.equals("gene") || s.equals("transcript")) return this.kg;
			return null;
			}
		@Override
		public boolean has(final String s) {
			if(s.equals("kg") || s.equals("gene") || s.equals("transcript")) return true;
			return false;
			}
		@Override
		public void set(String arg0, Object arg1) {
			throw new UnsupportedOperationException();
			}
		@Override
		public String toString() {
			return kg.toString();
			}
		}
	
	private void print(final KnownGene kg,final int start,final int end,final String type,final String name)
		{
		if(start>=end) return;
		out.print(kg.getChromosome());
		out.print('\t');
		out.print(start);
		out.print('\t');
		out.print(end);
		out.print('\t');
		out.print(kg.isPositiveStrand()?'+':'-');
		out.print('\t');
		out.print(kg.getName());
		out.print('\t');
		out.print(type);
		out.print('\t');
		out.print(name);
		out.println();
		}
	// 
	private void scan(final BufferedReader r) throws IOException
		{
		boolean hide_introns = false;
		boolean hide_utrs = false;
		boolean hide_cds = false;
		boolean hide_exons = false;
		boolean hide_transcripts = false;
		boolean hide_non_coding = false;
		boolean hide_coding = false;
		final Predicate<JexlContext> predicate = StringUtils.isBlank(selectExpr)?
				KG->true:
				new JexlPredicate(this.selectExpr)
				; 
		
		for(String str : CharSplitter.COMMA.split(this.hideStr)) {
			if(StringUtils.isBlank(str)) continue;
			str = str.trim().toUpperCase();
			if(str.equals("INSTRON") || str.equals("INSTRONS")) hide_introns=true;
			if(str.equals("UTR") || str.equals("UTRs")) hide_utrs=true;
			if(str.equals("CDS")) hide_cds=true;
			if(str.equals("EXON") || str.equals("EXONS")) hide_exons=true;
			if(str.equals("TRANSCRIPT") || str.equals("TRANSCRIPTS")) hide_transcripts=true;
			if(str.equals("NON_CODING")) hide_non_coding=true;
			if(str.equals("CODING")) hide_coding=true;
			}		
		
		String line;
		final CharSplitter tab=CharSplitter.TAB;
		while((line=r.readLine())!=null)
			{
			if(out.checkError()) break;
			final String tokens[]=tab.split(line);
			final KnownGene kg=new KnownGene(tokens);
			
			if(hide_coding && !kg.isNonCoding()) continue;
			if(hide_non_coding && kg.isNonCoding()) continue;
			
			if(!predicate.test(new KgContext(kg))) continue;
			
			if(!hide_transcripts) print(kg,kg.getTxStart(),kg.getTxEnd(),"TRANSCRIPT",kg.getName());
			for(int i=0;i< kg.getExonCount();++i)
				{
				final KnownGene.Exon exon=kg.getExon(i);
				if(!hide_exons) print(kg,exon.getStart(),exon.getEnd(),"EXON",exon.getName());
				
				if(!hide_utrs && kg.getCdsStart()>exon.getStart())
					{
					print(kg,exon.getStart(),
							Math.min(kg.getCdsStart(),exon.getEnd()),"UTR","UTR"+(kg.isPositiveStrand()?"5":"3"));
					}
				
				if(!hide_cds && !(kg.getCdsStart()>=exon.getEnd() || kg.getCdsEnd()<exon.getStart()))
					{
					print(kg,
							Math.max(kg.getCdsStart(),exon.getStart()),
							Math.min(kg.getCdsEnd(),exon.getEnd()),
							"CDS",exon.getName()
							);
					}
				
				final KnownGene.Intron intron=exon.getNextIntron();
				if(!hide_introns && intron!=null)
					{
					print(kg,intron.getStart(),intron.getEnd(),"INTRON",intron.getName());
					}
				
				if(!hide_utrs && kg.getCdsEnd()<exon.getEnd())
					{
					print(kg,Math.max(kg.getCdsEnd(),exon.getStart()),
							exon.getEnd(),
							"UTR","UTR"+(kg.isPositiveStrand()?"3":"5"));
					}
				
				}
			}
		}

	@Override
	public int doWork(final List<String> args) {
		BufferedReader r=null;
		try
			{
			this.out = super.openPathOrStdoutAsPrintStream(this.outputFile);
			if(args.isEmpty())
				{
				r=IOUtils.openStreamForBufferedReader(super.stdin());
				scan(r);
				r.close();r=null;
				}
			else
				{
				for(final String filename:args)
					{
					r=IOUtils.openURIForBufferedReading(filename);
					scan(r);
					r.close();r=null;
					}
				}
			this.out.flush();
			this.out.close();
			this.out=null;
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(this.out);
			this.out=null;
			}
		}

	
	public static void main(final String[] args) {
	new KnownGenesToBed().instanceMainWithExit(args);

	}

}
