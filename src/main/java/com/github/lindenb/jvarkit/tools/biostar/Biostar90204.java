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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.Collections;

import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;

public class Biostar90204 extends AbstractBiostar90204
	{
	
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar90204.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBiostar90204.AbstractBiostar90204Command
		{    
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			if(suffix_length<=0)
				{
				return wrapException("Bad value of suffix_length:"+suffix_length);
				}
			if(record_per_file<1L)
				{
				return wrapException("Bad value of record_per_file:"+record_per_file);
				}
			
			final SAMFileWriterFactory swfactory=new SAMFileWriterFactory();
			
			SAMFileWriter sfw=null;
			SAMRecordIterator iter=null;
			SamReader samFileReader=null;
			PrintWriter manifest=new PrintWriter(new NullOuputStream());
			try
				{
				samFileReader  = openSamReader(inputName);
				SAMFileHeader header=samFileReader.getFileHeader();
				
				int split_file_number=0;
				long nReads=0L;
				iter=samFileReader.iterator();
				
				if(manifestFile!=null)
					{
					manifest.close();
					manifest=new PrintWriter(manifestFile);
					}
				SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
				while(iter.hasNext())
					{
					SAMRecord rec= progress.watch(iter.next());
					++nReads;
					if(sfw==null)
						{
						split_file_number++;
						String pathname=this.prefix+"."+String.format("%0"+suffix_length+"d", split_file_number)+".bam";
						File out=new File(pathname);
						LOG.info("Opening "+out);
						manifest.write(pathname);
						manifest.write("\t"+(nReads)+"\t");
						
						SAMFileHeader header2=header.clone();
						header2.addComment("SPLIT:"+split_file_number);
						header2.addComment("SPLIT:Starting from Read"+nReads);
						
						sfw=swfactory.makeSAMOrBAMWriter(header2, true, out);
						}
					sfw.addAlignment(rec);
					
					if(nReads%super.record_per_file==0)
						{
						LOG.info("Closing "+sfw);
						sfw.close();
						manifest.write((nReads)+"\n");
						sfw=null;
						}
					
					}
				if(sfw!=null)
					{
					sfw.close();
					manifest.write((nReads)+"\n");
					}
				manifest.flush();
				progress.finish();
				return Collections.emptyList();
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(manifest);
				CloserUtil.close(sfw);
				CloserUtil.close(iter);
				CloserUtil.close(samFileReader);
				}
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar90204().instanceMainWithExit(args);

	}

}
