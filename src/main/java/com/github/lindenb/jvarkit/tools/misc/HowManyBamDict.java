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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

public class HowManyBamDict extends AbstractHowManyBamDict
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(HowManyBamDict.class);

	
	
	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractHowManyBamDict.AbstractHowManyBamDictCommand
	 	{		

	private MessageDigest md5;
	
	private class Dict
		{
		SAMSequenceDictionary ssd;
		String hash;
		File representative;
		Dict(SAMSequenceDictionary ssd,File representative)
			{
			this.ssd=ssd;
			this.representative=representative;
		    
	    	
	    	 md5.reset();
	    	 md5.update(String.valueOf(ssd.size()).getBytes());
	    	 for(SAMSequenceRecord ssr:ssd.getSequences())
	    	 	{
	    		 md5.update(ssr.getSequenceName().getBytes());
	    		 md5.update(String.valueOf(ssr.getSequenceLength()).getBytes());
	    	 	}
	         
	         hash = new BigInteger(1, md5.digest()).toString(16);
	         if (hash.length() != 32) {
	             final String zeros = "00000000000000000000000000000000";
	             hash= zeros.substring(0, 32 - hash.length()) + hash;
	         }
	       

			
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
		void print(PrintStream out)
			{
			out.print("DICT");
			out.print("\t");
			out.print(this.hash);
			out.print("\t");
			out.print(ssd.size());
			out.print("\t");
			out.print(ssd.getReferenceLength());
			out.print("\t");
			boolean first=true;
			for(SAMSequenceRecord ssr:ssd.getSequences())
				{
				if(!first) out.print(";");
				first=false;
				out.print(ssr.getSequenceName());
				out.print('=');
				out.print(ssr.getSequenceLength());
				}
			out.print("\t");
			out.print(this.representative);
			out.println();
			}
		}

	private Dict empty=null;
	private Set<Dict>  allditcs=new LinkedHashSet<Dict>();
	
	
 	
 	
 	private void handle(PrintStream out,File f) throws IOException
 		{
 		SamReader sfr=null;
 		try {
 			LOG.info(f);
			sfr=SamFileReaderFactory.mewInstance().open(f);
			SAMFileHeader header=sfr.getFileHeader();
			if(header==null || header.getSequenceDictionary()==null)
				{
				if(this.empty==null)
					{
					this.empty=new Dict(new SAMSequenceDictionary(),f);
					allditcs.add(this.empty);
					this.empty.print(out);
					}
				out.print("BAM\t");
				out.print(f.getPath());
				out.print("\t");
				out.print(this.empty.hash);
				out.println();
				}
			else
				{
				Dict d=new Dict(header.getSequenceDictionary(), f);
				if(this.allditcs.add(d))
					{
					d.print(out);
					}
				out.print("BAM\t");
				out.print(f.getPath());
				out.print("\t");
				out.print(d.hash);
				out.println();
				}
 			} 
 		catch (Exception e)
			{
			LOG.error(e);
			throw new IOException(e);
			}
 		finally
 			{
 			CloserUtil.close(sfr);
 			}
 		}
	
 	@Override
 		public Collection<Throwable> call() throws Exception {
 		
		 try {
   		  this.md5 = MessageDigest.getInstance("MD5");
         } catch (NoSuchAlgorithmException e) {
        	 LOG.error("MD5 algorithm not found");
            return wrapException(e);
         }
		final List<String> args  = super.getInputFiles();
		PrintStream out = null;
		try
			{
			out = super.openFileOrStdoutAsPrintStream();
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				String line;
				
					BufferedReader in=new BufferedReader(new InputStreamReader(stdin()));
					while((line=in.readLine())!=null)
						{
						if(line.isEmpty() || line.endsWith(File.separator) || line.startsWith("#")) continue;
						handle(out,new File(line));
						}
					in.close();
					
				}
			else
				{
				for(String arg:args)
					{
					handle(out,new File(arg));
					}
				}
			return RETURN_OK;
			}
		catch(IOException err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(out);
			}
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
