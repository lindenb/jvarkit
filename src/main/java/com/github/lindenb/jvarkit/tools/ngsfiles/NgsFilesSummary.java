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
package com.github.lindenb.jvarkit.tools.ngsfiles;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Date;
import java.util.List;

import htsjdk.variant.vcf.VCFHeader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfIteratorImpl;

@Program(name="ngsfilessummary",description="Scan folders and generate a summary of the files (SAMPLE/BAM SAMPLE/VCF etc..)")
public class NgsFilesSummary extends AbstractScanNgsFilesProgram
	{
	private static final Logger LOG = Logger.build(NgsFilesSummary.class).make();
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	
    private NgsFilesSummary()
    	{
    	
    	}		
    
  
    
    private void print(String sample,InfoType it,File f)
    	{
		System.out.print(sample);
		System.out.print('\t');
		System.out.print(it);
		System.out.print('\t');
		System.out.print(f);
		if(f.isFile())
			{
			System.out.print('\t');
			System.out.print(f.length());
			System.out.print('\t');
			System.out.print(new Date(f.lastModified()));
			}
		System.out.println();
		}
    	
    
	@Override
    protected void readBam(File f)
    	{
    	if(!f.canRead()) return;
    	SamReader r=null;
    	try {
			r=SamFileReaderFactory.mewInstance().open(f);
			SAMFileHeader h=r.getFileHeader();
			if(h!=null && h.getReadGroups()!=null)
				{
				for(SAMReadGroupRecord rg: h.getReadGroups())
					{
					String sample=rg.getSample();
					if(sample==null || sample.isEmpty()) continue;
					print(sample,InfoType.BAM, f);
					}
				}
			} 
    	catch (Exception e)
    		{
    		LOG.warning(e);
			}
    	finally
    		{
    		CloserUtil.close(r);
    		}
    	}
   
   @Override
   protected void readFastq(File f)
		{
    	//File parent=f.getParentFile();
    	//if(parent==null || super.VERBOSITY==Log.LogLevel.) return;
    	FastQName fq=FastQName.parse(f);
    	
		
		if(!fq.isValid())
			{
			//bad name
			return;
			}
		
    	print(fq.getSample(),InfoType.FASTQ, f);
		}
    
   @Override
   protected void readVCF(File f)
		{
    	if(!f.canRead()) return;
    	LOG.debug("readVCF  "+f);
    	    	

    	VcfIterator r=null;
    	InputStream in=null;
    	try
    		{
    		in=IOUtils.openFileForReading(f);
    		
    		r=new VcfIteratorImpl(in);
        	VCFHeader header=r.getHeader();
        	for(String sample:header.getSampleNamesInOrder())
	        	{
	        	print(sample,InfoType.VCF, f);
	    		}
    		}
    	catch(Exception err)
    		{
    		LOG.error(err);
    		}
    	finally
    		{
    		CloserUtil.close(r);
    		CloserUtil.close(in);
    		}
    	
		}

    private void scan(BufferedReader in) throws IOException
    	{
    	String line;
    	while((line=in.readLine())!=null)
    			{
    			if(line.isEmpty() || line.startsWith("#")) continue;
    			analyze(new File(line));
    			}
    	}
    
	
    @Override
	public int doWork(List<String> args) {
		try
			{
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				scan(new BufferedReader(new InputStreamReader(stdin())));
				}
			else
				{
				for(final String filename:args)
					{
					LOG.info("Reading from "+filename);
					BufferedReader r=IOUtils.openURIForBufferedReading(filename);
					scan(r);
					r.close();
					}
				}
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new NgsFilesSummary().instanceMainWithExit(args);

	}

}
