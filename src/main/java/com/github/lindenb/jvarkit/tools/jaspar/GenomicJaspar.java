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
package com.github.lindenb.jvarkit.tools.jaspar;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.SubSequence;
import com.github.lindenb.jvarkit.util.bio.RevCompCharSequence;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
## Example

```bash
$ java -jar dist/genomicjaspar.jar  -J pfm_vertebrates.txt human_g1k_v37.fasta

1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	10819	10825	MA0130.1 ZNF354C	978	-	6	ATCCAC	CTCCAC
1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	10895	10901	MA0130.1 ZNF354C	978	-	6	ATCCAC	CTCCAC
1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	10971	10977	MA0130.1 ZNF354C	978	-	6	ATCCAC	CTCCAC
1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	11088	11094	MA0006.1 Arnt::Ahr	957	-	6	TGCGTG	CGCGTG
1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	11104	11112	MA0067.1 Pax2	951	+	8	AGTCACGG	CGTCACGG
1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	11421	11427	MA0056.1 MZF1_1-4	1000	+	6	TGGGGA	TGGGGA
1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	11550	11558	MA0033.1 FOXL1	959	-	8	TATACATA	TAAACATA
1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	11554	11560	MA0151.1 ARID3A	1000	-	6	ATTAAA	ATTAAA
1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	11556	11561	MA0075.1 Prrx2	1000	-	5	AATTA	AATTA
1 dna:chromosome chromosome:GRCh37:1:1:249250621:1	11629	11635	MA0130.1 ZNF354C	978	+	6	ATCCAC	CTCCAC

```

## See also

* [[VcfJaspar]]
* http://www.biostars.org/p/90823/
* http://jaspar.genereg.net/

 */
@Program(name="genomicjaspar",
	description="Find jaspar patterns in FASTA sequences. Reports a BED file.",
	keywords={"jaspar","genomic","pattern"}
	)
public class GenomicJaspar extends Launcher
	{
	private static final Logger LOG = Logger.build(GenomicJaspar.class).make();
	@Parameter(names="-o",description=OPT_OUPUT_FILE_OR_STDOUT)
	private File OUT=null;
	@Parameter(names="-J",description=" jaspar PFM uri. required. example: http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt",required=true)
	private String jasparUri=null;
	@Parameter(names="-f",description="(0<ratio<1) fraction of best score")
	private double fraction_of_max=0.95;

	
	private final List<Matrix> jasparDb=new ArrayList<Matrix>();
	GenomicJaspar()
		{
		}
	
	
	private void digest(
			final PrintWriter out,
			final String seqName,
			final int position0,
			final StringBuilder sequence
			)
		{
		
		for(final Matrix matrix:this.jasparDb)
			{
			if(matrix.length()>sequence.length()) continue;
			
			final CharSequence forward=new SubSequence(sequence,0,matrix.length());
			final CharSequence revcomp=new RevCompCharSequence(forward);
			

			
			for(int strand=0;strand<2;++strand)
				{
				double score= matrix.score(strand==0?forward:revcomp);
				if(score<=0) continue;
				
				if(score>= matrix.max()*this.fraction_of_max)
					{
					out.print(seqName);
					out.print('\t');
					out.print(position0);
					out.print('\t');
					out.print(position0+matrix.length());
					out.print('\t');
					out.print(matrix.getName());
					out.print('\t');
					out.print((int)(1000.0*(score/matrix.max())));
					out.print('\t');
					out.print(strand==1?'-':'+');
					out.print('\t');
					out.print(matrix.length());
					out.print('\t');
					out.print(matrix.getArchetype());
					out.print('\t');
					out.print(strand==0?forward:revcomp);
					out.println();
					
					break;
					}
				}
			}
		}
	
	private void run(final PrintWriter out,final Reader in) throws IOException
		{
		int longest=0;
		for(final Matrix m:this.jasparDb)
			{
			longest=Math.max(m.length(), longest);
			}
		LOG.info("longest:"+longest);
		String seqName="";
		int position0=0;
		final StringBuilder sequences=new StringBuilder(longest);
		for(;;)
			{
			int c=in.read();
			if(c==-1 || c=='>')
				{
				while(sequences.length()!=0)
					{
					digest(out,seqName,position0,sequences);
					++position0;
					sequences.deleteCharAt(0);
					}
				if(c==-1) break;
				final StringBuilder b=new StringBuilder();
				while((c=in.read())!=-1 && c!='\n')
					{
					b.append((char)c);
					}
				seqName=b.toString();
				position0=0;
				}
			else if(!Character.isWhitespace(c))
				{
				sequences.append((char)Character.toUpperCase(c));
				if(sequences.length()==longest)
					{
					digest(out,seqName,position0,sequences);
					if(out.checkError())  return ;
					++position0;
					sequences.deleteCharAt(0);
					if(position0%1000000==0)
						{
						LOG.info(seqName+" "+position0);
						}
					}
				}
			}
		}

	@Override
	public int doWork(final List<String> args) {
		if(jasparUri==null)
			{
			LOG.error("Undefined jaspar-uri");
			return -1;
			}
		PrintWriter out=null;
		try
			{
			out = super.openFileOrStdoutAsPrintWriter(OUT);
			LOG.info("Reading "+jasparUri);
			LineReader lr= IOUtils.openURIForLineReader(jasparUri);
			LineIterator liter=new LineIteratorImpl(lr);
			Iterator<Matrix> miter=Matrix.iterator(liter);
			while(miter.hasNext())
				{
				Matrix matrix = miter.next();
				this.jasparDb.add(matrix.convertToPWM());
				}
			lr.close();
			LOG.info("JASPAR size: "+this.jasparDb.size());
			
			
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				run(out,new InputStreamReader(stdin()));
				}
			else
				{
				for(final String fname:args)
					{
					LOG.info("Opening "+fname);
					Reader in=IOUtils.openURIForBufferedReading(fname);
					run(out,in);
					in.close();
					}
					
				}
			out.flush();
			out.close();
			out=null;
			return 0;
			}
		catch(Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new GenomicJaspar().instanceMainWithExit(args);
		}

	}
