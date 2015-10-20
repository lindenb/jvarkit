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

import java.util.Collection;
import java.util.List;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.command.Command;

public class PadEmptyFastq extends AbstractPadEmptyFastq
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(PadEmptyFastq.class);

	
	private static final int DEFAULT_LENGTH=50;
	
	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractPadEmptyFastq.AbstractPadEmptyFastqCommand
	 	{

	
	private void copyTo(FastqReader r,FastqWriter w)
		{
		int padLength=this.N;
		long nReads=0L;
		long nFill=0L;
		String fillN=null;
		String fillQ=null;
		while(r.hasNext())
			{
			FastqRecord rec=r.next();
			
			
			if(++nReads%1E6==0)
				{
				LOG.info("Read "+nReads +" reads. empty reads="+nFill);
				}
			if(rec.getReadString().isEmpty())
				{
				++nFill;
				if(padLength<1)
					{
					padLength=DEFAULT_LENGTH;
					}
				if(fillN==null)
					{
					StringBuilder b1=new StringBuilder();
					while(b1.length()< padLength) b1.append("N");
					fillN=b1.toString();
					fillQ=fillN.replace('N', '#');
					}
				
				rec=new FastqRecord(
						rec.getReadHeader(),
						fillN,
						rec.getBaseQualityHeader(),
						fillQ
						);
				}
			else if(padLength<1)
				{
				padLength=rec.getReadString().length();
				}
			w.write(rec);
			}
		LOG.info("Done. Read "+nReads +" reads. empty reads="+nFill);
		}
	

	
	@Override
		public Collection<Throwable> call() throws Exception {		
		FastqWriter fqw=null;		
		final List<String> args = getInputFiles();
		try
			{
			
			fqw = openFastqWriter();
			
			if(args.isEmpty())
				{
				FastqReader fqr= super.openFastqReader(null);
				copyTo(fqr,fqw);
				fqr.close();
				}
			else
				{
				for(String filename:args)
					{
					LOG.info("Reading from "+filename);
					FastqReader fqr= super.openFastqReader(filename);
					copyTo(fqr,fqw);
					fqr.close();
					}
				}
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(fqw);
			}
		}
	 	}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new PadEmptyFastq().instanceMainWithExit(args);
		}

}
