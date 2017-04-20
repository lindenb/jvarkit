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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;



@Program(name="howmanybamdict",description="finds if there's are some differences in the sequence dictionaries.")
public class HowManyBamDict extends Launcher {
	private static final Logger LOG = Logger.build(HowManyBamDict.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	public HowManyBamDict()
		{
		 
		}
	
	
	private class Dict
		{
		SAMSequenceDictionary ssd;
		String hash;
		File representative;
		Dict(SAMSequenceDictionary ssd,File representative)
			{
			this.ssd=ssd;
			this.representative=representative;
	    	this.hash =ssd.md5();
			}
		
		@Override
		public boolean equals(Object obj)
			{
			if(this==obj) return true;
			if(obj==null) return false;
			return this.ssd.equals(Dict.class.cast(obj).ssd);
			}
		
		@Override
		public int hashCode() {
			return ssd.hashCode();
			}
		void print()
			{
			System.out.print("DICT");
			System.out.print("\t");
			System.out.print(this.hash);
			System.out.print("\t");
			System.out.print(ssd.size());
			System.out.print("\t");
			System.out.print(ssd.getReferenceLength());
			System.out.print("\t");
			boolean first=true;
			for(SAMSequenceRecord ssr:ssd.getSequences())
				{
				if(!first) System.out.print(";");
				first=false;
				System.out.print(ssr.getSequenceName());
				System.out.print('=');
				System.out.print(ssr.getSequenceLength());
				}
			System.out.print("\t");
			System.out.print(this.representative);
			System.out.println();
			}
		}

	private Dict empty=null;
	private Set<Dict>  allditcs=new LinkedHashSet<Dict>();
	
	
 	
 	private void handle(PrintWriter out,File f) throws IOException
 		{
 		SamReader sfr=null;
 		try {
 			LOG.info(f.getPath());
			sfr=SamFileReaderFactory.mewInstance().open(f);
			SAMFileHeader header=sfr.getFileHeader();
			if(header==null || header.getSequenceDictionary()==null)
				{
				if(this.empty==null)
					{
					this.empty=new Dict(new SAMSequenceDictionary(),f);
					allditcs.add(this.empty);
					this.empty.print();
					}
				System.out.print("BAM\t");
				System.out.print(f.getPath());
				System.out.print("\t");
				System.out.print(this.empty.hash);
				System.out.println();
				}
			else
				{
				Dict d=new Dict(header.getSequenceDictionary(), f);
				if(this.allditcs.add(d))
					{
					d.print();
					}
				System.out.print("BAM\t");
				System.out.print(f.getPath());
				System.out.print("\t");
				System.out.print(d.hash);
				System.out.println();
				}
 			} 
 		catch (Exception e)
			{
			LOG.error(e.getMessage(),e);
			throw new IOException(e);
			}
 		finally
 			{
 			CloserUtil.close(sfr);
 			}
 		}
	
 	@Override
 	public int doWork(List<String> args) {
		PrintWriter out=null;
		try
			{
			out = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				String line;
				
					BufferedReader in=new BufferedReader(new InputStreamReader(System.in));
					while((line=in.readLine())!=null)
						{
						if(line.isEmpty() || line.endsWith(File.separator) || line.startsWith("#")) continue;
						handle(out,new File(line));
						}
					in.close();
					
				}
			else
				{
				for(String filename:args)
					{
					handle(out,new File(filename));
					}
				}
			out.flush();
			return RETURN_OK;
			}
		catch(IOException err)
			{
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(out);
			}
		}

	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new HowManyBamDict().instanceMainWithExit(args);
		}

}
