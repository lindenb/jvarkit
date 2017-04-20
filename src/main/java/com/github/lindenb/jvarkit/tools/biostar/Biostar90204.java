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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;

@Program(name="biostar90204",description="Bam version of linux split. See also http://www.biostars.org/p/90204/",biostars=90204)
public class Biostar90204 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar90204.class).make();

	@Parameter(names="-p",description="(prefix) output file prefix.")
	private String prefix="_splitbam";
	@Parameter(names="-a",description="suffix length")
	private int suffix_length=2;
	@Parameter(names="-M",description=" manifest file. Optional")
	private File manifestFile=null;
	
	@Parameter(names="-n",description="Records per file")
	private long record_per_file=-1L;

	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	
	private  Biostar90204() {
		}
	
	@Override
	public int doWork(List<String> args) {
		
		
		if(suffix_length<0)
			{
			LOG.error("Bad value of suffix_length:"+suffix_length);
			return -1;
			}
		if(record_per_file<1L)
			{
			LOG.error("Bad value of record_per_file:"+record_per_file);
			return -1;
			}
		
		
		SAMFileWriter sfw=null;
		SAMRecordIterator iter=null;
		SamReader samFileReader=null;
		PrintWriter manifest=new PrintWriter(new NullOuputStream());
		try
			{
			samFileReader = super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader header=samFileReader.getFileHeader();
			
			
			int split_file_number=0;
			long nReads=0L;
			iter=samFileReader.iterator();
			
			if(manifestFile!=null)
				{
				manifest.close();
				manifest=new PrintWriter(manifestFile);
				}
			
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
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
					
					sfw=this.writingBamArgs.openSAMFileWriter(out,header2, true);
					}
				sfw.addAlignment(rec);
				
				if(nReads%record_per_file==0)
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
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(manifest);
			CloserUtil.close(sfw);
			CloserUtil.close(iter);
			CloserUtil.close(samFileReader);
			}
		return 0;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar90204().instanceMainWithExit(args);

	}

}
