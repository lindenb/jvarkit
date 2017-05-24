/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.SynchronousLineReader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.SubSequence;
import com.github.lindenb.jvarkit.util.bio.RevCompCharSequence;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="genomicjaspar",description="Find jaspar patterns in FASTA sequences. Reports a BED file.")
public class GenomicJaspar extends Launcher
	{
	private static final Logger LOG = Logger.build(GenomicJaspar.class).make();
	@Parameter(names="-J",description=" jaspar PFM uri. required. example: http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt")
	private String jasparUri=null;
	@Parameter(names="-f",description="(0<ratio<1) fraction of best score")
	private double fraction_of_max=0.95;

	
	private List<Matrix> jasparDb=new ArrayList<Matrix>();
	private GenomicJaspar()
		{
		}
	
	
	private void digest(
			String seqName,
			int position0,
			final StringBuilder sequence
			)
		{
		
		for(final Matrix matrix:this.jasparDb)
			{
			if(matrix.length()>sequence.length()) continue;
			
			CharSequence forward=new SubSequence(sequence,0,matrix.length());
			CharSequence revcomp=new RevCompCharSequence(forward);
			

			
			for(int strand=0;strand<2;++strand)
				{
				double score= matrix.score(strand==0?forward:revcomp);
				if(score<=0) continue;
				
				if(score>= matrix.max()*this.fraction_of_max)
					{
					System.out.print(seqName);
					System.out.print('\t');
					System.out.print(position0);
					System.out.print('\t');
					System.out.print(position0+matrix.length());
					System.out.print('\t');
					System.out.print(matrix.getName());
					System.out.print('\t');
					System.out.print((int)(1000.0*(score/matrix.max())));
					System.out.print('\t');
					System.out.print(strand==1?'-':'+');
					System.out.print('\t');
					System.out.print(matrix.length());
					System.out.print('\t');
					System.out.print(matrix.getArchetype());
					System.out.print('\t');
					System.out.print(strand==0?forward:revcomp);
					System.out.println();
					
					break;
					}
				}
			}
		}
	
	private void run(Reader in) throws IOException
		{
		int longest=0;
		for(Matrix m:this.jasparDb)
			{
			longest=Math.max(m.length(), longest);
			}
		LOG.info("longest:"+longest);
		String seqName="";
		int position0=0;
		StringBuilder sequences=new StringBuilder(longest);
		for(;;)
			{
			int c=in.read();
			if(c==-1 || c=='>')
				{
				while(sequences.length()!=0)
					{
					digest(seqName,position0,sequences);
					++position0;
					sequences.deleteCharAt(0);
					}
				if(c==-1) break;
				StringBuilder b=new StringBuilder();
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
					digest(seqName,position0,sequences);
					if(System.out.checkError())  return ;
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
	public int doWork(List<String> args) {
		if(jasparUri==null)
			{
			LOG.error("Undefined jaspar-uri");
			return -1;
			}
		
		
		
		try
			{
			LOG.info("Reading "+jasparUri);
			LineReader lr=new SynchronousLineReader(IOUtils.openURIForReading(jasparUri));
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
				run(new InputStreamReader(stdin()));
				}
			else
				{
				for(final String fname:args)
					{
					LOG.info("Opening "+fname);
					Reader in=IOUtils.openURIForBufferedReading(fname);
					run(in);
					in.close();
					}
					
				}
			return 0;
			}
		catch(Throwable err)
			{
			LOG.error(err);
			return -1;
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
