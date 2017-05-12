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
import java.io.InputStream;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="bedrenamechr",description="Convert the names of the chromosomes in a Bed file",keywords={"bed","chromosome","contig","convert"})
public class ConvertBedChromosomes
	extends Launcher
	{
	private static final Logger LOG = Logger.build(ConvertBedChromosomes.class).make();
	
	@Parameter(names={"-convert","--convert"},description="What should I do when  a converstion is not found")
	private ContigNameConverter.OnNotFound onNotFound=ContigNameConverter.OnNotFound.RAISE_EXCEPTION;
	@Parameter(names={"-f","--mapping","-m"},description="load a custom name mapping. Format (chrom-source\\tchrom-dest\\n)+",required=true)
	private File mappingFile=null;
	@Parameter(names={"-c","--column"},description="1-based chromosome column")
	private int chromColumn1=1;
	@Parameter(names={"-o","--out"},description="output bed. Default stdout")
	private File outputFile= null;

	
	private ContigNameConverter customMapping=ContigNameConverter.getIdentity();
	private Set<String> unmappedChromosomes=new HashSet<String>();
	
	private ConvertBedChromosomes()
		{
		
		}
	
	private String convertName(final String chrom)throws IOException
		{
		if(chrom==null) throw new NullPointerException();
		String newname=customMapping.apply(chrom);
		if(newname==null)
			{
			if(!unmappedChromosomes.contains(chrom))
				{
				LOG.warning("unmapped chromosome "+chrom);
				unmappedChromosomes.add(chrom);
				}
			return null;
			}
		return newname;
		}
	
	@SuppressWarnings("resource")
	protected int doWork(InputStream in,PrintStream out)
			throws IOException
		{
		final int chromColumn0=chromColumn1-1;
	
		Pattern tab=Pattern.compile("[\t]");
		LineIterator lr=new LineIteratorImpl(new AsciiLineReader(in));
		
		
		while(lr.hasNext())
			{	
			String line=lr.next();
			if(BedLine.isBedHeader(line))
				{
				out.println(line);
				continue;
				}
			final String tokens[]=tab.split(line, (chromColumn0+2));
			if(chromColumn0 >=tokens.length) throw new IOException("Bad BED line : "+line+" extected at least "+(chromColumn0+2)+" columns");
			final String chrom=convertName(tokens[chromColumn0]);
			if(chrom==null) continue;
			for(int i=0;i< tokens.length;++i)
				{
				if(i>0) out.print("\t");
				out.print(i==chromColumn0?chrom:tokens[i]);
				}
			out.println();
			}
		out.flush();
		return 0;
		}
	
	@Override
	public int doWork(List<String> args) {
		if(this.chromColumn1<1)
			{
			LOG.error("bad chromosome index (<1): "+this.chromColumn1);
			return -1;
			}
		if(this.mappingFile!=null) {
			LOG.info("reading custom mapping "+mappingFile);
			this.customMapping=ContigNameConverter.fromFile(mappingFile);
			}
		
		PrintStream out=null;
		try
			{
			out = super.openFileOrStdoutAsPrintStream(this.outputFile);
			if(args.isEmpty())
				{
				LOG.info("reading stdin");
				doWork(stdin(), out);
				}
			else
				{
				for(final String filename:args)
					{
					InputStream in=IOUtils.openURIForReading(filename);
					doWork(in, out);
					CloserUtil.close(in);
					}
				}
			if(!unmappedChromosomes.isEmpty())
				{
				LOG.warning("Unmapped chromosomes:"+unmappedChromosomes);
				}
			out.flush();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	

	public static void main(String[] args)
		{
		new ConvertBedChromosomes().instanceMainWithExit(args);
		}
	}
