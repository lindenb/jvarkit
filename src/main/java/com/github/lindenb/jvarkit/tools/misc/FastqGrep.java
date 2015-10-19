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
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

public class FastqGrep
	extends AbstractFastqGrep
	{

	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(FastqGrep.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractFastqGrep.AbstractFastqGrepCommand
	 	{		
		private Map<String,Integer> readNames=new HashMap<String,Integer>(); 

		
		private String getReadName(FastqRecord r)
			{
			return getReadName(r.getReadHeader());
			}
		
		private String getReadName(String s)
			{
			int beg=(s.startsWith(FastqConstants.SEQUENCE_HEADER)?1:0);
			int end=s.indexOf(' ');
			if(end==-1) end=s.length();
			s= s.substring(beg, end);
			return s;
			}
		private void run(FastqReader r,FastqWriter out)
			{
			long nRec=0L;
			r.setValidationStringency(ValidationStringency.LENIENT);
			while(r.hasNext())
				{
				FastqRecord fastq=r.next();
				boolean keep=false;
				String readName=getReadName(fastq);
				Integer count=readNames.get(readName);
				if(count!=null)
					{
					keep=true;
					}
				if(super.inverse) keep=!keep;
				if(keep)
					{
					++nRec;
					out.write(fastq);
					}
				
				if(super.n_before_remove!=-1 && !super.inverse && keep)
					{
					count++;
					if(count>=n_before_remove)
						{
						this.readNames.remove(readName);
						if(this.readNames.isEmpty()) break;
						}
					else
						{
						this.readNames.put(readName,count);
						}
					}
				
				
				}
			LOG.info("Done. N-Reads:"+nRec);
			}
		
	@Override
	public Collection<Throwable> call() throws Exception
		{
		
		if(super.readNameFile!=null)
			{
			BufferedReader in=null;
			try {
				in=IOUtils.openFileForBufferedReading(super.readNameFile);
		    	String line;
		    	while((line=in.readLine())!=null)
		    		{
		    		line=line.trim();
		    		if(line.isEmpty()) continue;
		    		this.readNames.put(getReadName(line),0);
		    		}
			} catch (Exception e) {
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(in);
				}
			}
		for(String rn: super.readNameStrings)
			{
			this.readNames.put(getReadName(rn),0);
			}
		
		if(readNames.isEmpty())
    		{
    		LOG.warn("no read name found.");
    		}
		List<String> args = this.getInputFiles();
		FastqWriter out=null;
		try
			{
			out  =  openFastqWriter();
			
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				FastqReader fqR=new FourLinesFastqReader(stdin());
				run(fqR,out);
				fqR.close();
				}
			else for(String arg:args)
				{
				File f=new File(arg);
				LOG.info("Reading from "+f);
				FastqReader fqR=new FourLinesFastqReader(f);
				run(fqR,out);
				fqR.close();
				}
			
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	
	 	}
	
	public static void main(String[] args) {
		new FastqGrep().instanceMainWithExit(args);

	}

}
