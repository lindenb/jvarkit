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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReaderUtil;

import htsjdk.samtools.util.CloserUtil;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;

public class ConvertBedChromosomes
	extends AbstractConvertBedChromosomes
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(ConvertBedChromosomes.class);

	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractConvertBedChromosomes.AbstractConvertBedChromosomesCommand
	 	{

	
		private Map<String,String> customMapping=new HashMap<String,String>();
		private Set<String> unmappedChromosomes=new HashSet<String>();
		
		
		private String convertName(String chrom)throws IOException
			{
			if(chrom==null) throw new NullPointerException();
			String newname=customMapping.get(chrom);
			if(newname==null)
				{
				if(!unmappedChromosomes.contains(chrom))
					{
					LOG.warn("unmapped chromosome "+chrom);
					unmappedChromosomes.add(chrom);
					}
				return null;
				}
			return newname;
			}
		
		@SuppressWarnings("resource")
		private int doWork(InputStream in,PrintStream out)
				throws IOException
			{
			Pattern tab=Pattern.compile("[\t]");
			LineIterator lr=new LineIteratorImpl(LineReaderUtil.fromBufferedStream(in));
			while(lr.hasNext())
				{	
				String line=lr.next();
				String tokens[]=tab.split(line, (chromColumn0+2));
				if(chromColumn0 >=tokens.length) throw new IOException("Bad BED line : "+line+" extected at least "+(chromColumn0+2)+" columns");
				String chrom=convertName(tokens[chromColumn0]);
				if(chrom==null) continue;
				for(int i=0;i< tokens.length;++i)
					{
					if(i>0) out.print("\t");
					out.print(i==chromColumn0?chrom:tokens[i]);
					}
				out.println();
				}
			
			return 0;
			}
		
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception
			{
			if(super.mappingFile==null)
				{
				return wrapException("undefined mapping file");
				}
			else
				{
				try
					{
					this.customMapping =  super.loadCustomChromosomeMapping(getMappingFile());					}
				catch(Exception err)
					{
					return wrapException(err);
					}
				
				}
				
				if(customMapping.isEmpty())
					{
					LOG.error("No custom mapping defined");
					}
				PrintStream out = null;
				try
					{
					out =  openFileOrStdoutAsPrintStream();
					if(inputName==null)
						{
						doWork(System.in, System.out);
						}
					else
						{
						InputStream in=IOUtils.openURIForReading(inputName);
						doWork(in, System.out);
						CloserUtil.close(in);
						}
					if(!unmappedChromosomes.isEmpty())
						{
						LOG.warn("Unmapped chromosomes:"+unmappedChromosomes);
						}
					return RETURN_OK;
					}
				catch(Exception err)
					{
					return wrapException(err);
					}
				finally{
					CloserUtil.close(out);
					}
				}
	 	}

	public static void main(String[] args)
		{
		new ConvertBedChromosomes().instanceMainWithExit(args);
		}
	}
