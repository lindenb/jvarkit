/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloserUtil;


import com.github.lindenb.jvarkit.io.IOUtils;
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
	biostars=151628
	)
public class KnownGenesToBed extends Launcher
	{
	private static final Logger LOG = Logger.build(KnownGenesToBed.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-i","--intron"},description="Hide Introns")
	private boolean hide_introns = false;

	@Parameter(names={"-u","--utr"},description="Hide UTRs")
	private boolean hide_utr = false;

	@Parameter(names={"-c","--cds"},description="Hide CDSs")
	private boolean hide_cds = false;

	@Parameter(names={"-x","--exon"},description="Hide Exons")
	private boolean hide_exons = false;

	@Parameter(names={"-t","--transcript"},description="Hide Transcript")
	private boolean hide_transcripts = false;

	private PrintStream out;
	
	
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
	
	private void scan(final BufferedReader r) throws IOException
		{
		String line;
		final Pattern tab=Pattern.compile("[\t]");
		while((line=r.readLine())!=null)
			{
			if(out.checkError()) break;
			final String tokens[]=tab.split(line);
			final KnownGene kg=new KnownGene(tokens);
			if(!this.hide_transcripts) print(kg,kg.getTxStart(),kg.getTxEnd(),"TRANSCRIPT",kg.getName());
			for(int i=0;i< kg.getExonCount();++i)
				{
				final KnownGene.Exon exon=kg.getExon(i);
				if(!this.hide_exons) print(kg,exon.getStart(),exon.getEnd(),"EXON",exon.getName());
				
				if(!this.hide_utr && kg.getCdsStart()>exon.getStart())
					{
					print(kg,exon.getStart(),
							Math.min(kg.getCdsStart(),exon.getEnd()),"UTR","UTR"+(kg.isPositiveStrand()?"5":"3"));
					}
				
				if(!this.hide_cds && !(kg.getCdsStart()>=exon.getEnd() || kg.getCdsEnd()<exon.getStart()))
					{
					print(kg,
							Math.max(kg.getCdsStart(),exon.getStart()),
							Math.min(kg.getCdsEnd(),exon.getEnd()),
							"CDS",exon.getName()
							);
					}
				
				KnownGene.Intron intron=exon.getNextIntron();
				if(!this.hide_introns && intron!=null)
					{
					print(kg,intron.getStart(),intron.getEnd(),"INTRON",intron.getName());
					}
				
				if(!this.hide_utr && kg.getCdsEnd()<exon.getEnd())
					{
					print(kg,Math.max(kg.getCdsEnd(),exon.getStart()),
							exon.getEnd(),
							"UTR","UTR"+(kg.isPositiveStrand()?"3":"5"));
					}
				
				}
			}
		}

	@Override
	public int doWork(List<String> args) {
		BufferedReader r=null;
		try
			{
			this.out = super.openFileOrStdoutAsPrintStream(this.outputFile);
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				r=IOUtils.openStreamForBufferedReader(super.stdin());
				scan(r);
				CloserUtil.close(r);
				}
			else
				{
				for(final String filename:args)
					{
					LOG.info("Reading from "+filename);
					r=IOUtils.openURIForBufferedReading(filename);
					scan(r);
					CloserUtil.close(r);
					}
				}
			this.out.flush();
			return RETURN_OK;
			}
		catch(final Exception err)
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
