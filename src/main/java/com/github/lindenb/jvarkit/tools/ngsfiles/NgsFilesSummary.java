/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.ngsfiles;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Date;
import java.util.List;

import htsjdk.variant.vcf.VCFHeader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

/**
BEGIN_DOC

## Example

```bash
find /projects/align01/ -type f |\
  java -jar dist/ngsfilessummary.jar 

SAMPLE1	BAM	/projects/align01/Samples/SAMPLE1/BAM/SAMPLE1_final.bam	321262321	Wed Jun 26 10:30:07 CEST 2013
SAMPLE1	FASTQ	/project/align01/fastq/SAMPLE1/SAMPLE1_CGATGT_L008_R1_002.fastq.gz	35828879	Fri Oct 18 16:15:58 CEST 2013
SAMPLE1	VCF	/projects/align01/Samples/SAMPLE1/VCF/SAMPLE1_variations.freebayes.vcf.gz	184191	Mon Jun 17 14:47:22 CEST 2013
SAMPLE1	VCF	/projects/align01/Samples/SAMPLE1/VCF/SAMPLE1_variations.gatk.vcf.gz	113341	Mon Jun 17 11:57:19 CEST 2013
SAMPLE1	VCF	/projects/align01/Samples/SAMPLE1/VCF/SAMPLE1_variations.samtools.vcf.gz	57518	Mon Jun 17 11:58:49 CEST 2013
SAMPLE2	BAM	/projects/align01/Samples/SAMPLE2/BAM/SAMPLE2_final.bam	286100773	Wed Jun 26 10:47:09 CEST 2013
SAMPLE2	FASTQ	/project/align01/fastq/SAMPLE2/SAMPLE2_CGATGT_L008_R1_002.fastq.gz	356828879	Fri Oct 18 16:15:58 CEST 2013
SAMPLE2	VCF	/projects/align01/Samples/SAMPLE2/VCF/SAMPLE2_variations.freebayes.vcf.gz	172970	Mon Jun 17 14:45:51 CEST 2013
SAMPLE2	VCF	/projects/align01/Samples/SAMPLE2/VCF/SAMPLE2_variations.gatk.vcf.gz	106390	Mon Jun 17 11:57:19 CEST 2013
SAMPLE2	VCF	/projects/align01/Samples/SAMPLE2/VCF/SAMPLE2_variations.samtools.vcf.gz	52709	Mon Jun 17 11:58:04 CEST 2013
```
END_DOC

 */
@Program(name="ngsfilessummary",
	description="Scan folders and generate a summary of the files (SAMPLE/BAM SAMPLE/VCF etc..). Useful to get a summary of your samples.",
	keywords= {"sam","bam","vcf","util"}
		)
public class NgsFilesSummary extends AbstractScanNgsFilesProgram
	{
	private static final Logger LOG = Logger.build(NgsFilesSummary.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-header","--header"},description="[20180725]print header")
	private boolean show_header= false;

	private PrintWriter printWriter=null;
	
    public NgsFilesSummary()
    	{
    	
    	}		
    
  
    
    private void print(final String sample,final InfoType it,final File f)
    	{
		this.printWriter.print(sample);
		this.printWriter.print('\t');
		this.printWriter.print(it);
		this.printWriter.print('\t');
		this.printWriter.print(f);
		if(f.isFile())
			{
			this.printWriter.print('\t');
			this.printWriter.print(f.length());
			this.printWriter.print('\t');
			this.printWriter.print(new Date(f.lastModified()));
			}
		this.printWriter.println();
		}
    	
    
	@Override
    protected void readBam(final File f)
    	{
    	if(!f.canRead()) return;
    	SamReader r=null;
    	try {
			r = super.openSamReader(f.getPath());
			SAMFileHeader h=r.getFileHeader();
			if(h!=null && 
					h.getReadGroups()!=null && 
					h.getReadGroups().isEmpty())
				{
				for(final SAMReadGroupRecord rg: h.getReadGroups())
					{
					String sample=rg.getSample();
					if(StringUtil.isBlank(sample)) {
						sample = "_NO_SAMPLE_RG_";
						}
					print(sample,InfoType.BAM, f);
					}
				}
			else
				{
				print("_NO_READ_GROUP_",InfoType.BAM, f);
				}
			} 
    	catch (final Exception e)
    		{
    		LOG.warning(e);
			}
    	finally
    		{
    		CloserUtil.close(r);
    		}
    	}
   
   @Override
   protected void readFastq(final File f)
		{
    	//File parent=f.getParentFile();
    	//if(parent==null || super.VERBOSITY==Log.LogLevel.) return;
    	final FastQName fq=FastQName.parse(f);
    	
		
		if(!fq.isValid())
			{
			//bad name
			return;
			}
		
    	print(fq.getSample(),InfoType.FASTQ, f);
		}
    
   @Override
   protected void readVCF(final File f)
		{
    	if(!f.canRead()) return;
    	LOG.debug("readVCF  "+f);
    	    	

    	VCFIterator r=null;
    	InputStream in=null;
    	try
    		{
    		in=IOUtils.openFileForReading(f);
    		
    		r= new VCFIteratorBuilder().open(in);
        	final VCFHeader header=r.getHeader();
        	for(final String sample:header.getSampleNamesInOrder())
	        	{
	        	print(sample,InfoType.VCF, f);
	    		}
    		}
    	catch(final Exception err)
    		{
    		LOG.error(err);
    		}
    	finally
    		{
    		CloserUtil.close(r);
    		CloserUtil.close(in);
    		}
    	
		}

    private void scan(final BufferedReader in) throws IOException
    	{
    	String line;
    	while((line=in.readLine())!=null)
    			{
    			if(StringUtil.isBlank(line) || line.startsWith("#")) continue;
    			analyze(new File(line));
    			}
    	}
    
	
    @Override
	public int doWork(final List<String> args) {
		try
			{
			this.printWriter = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			if(show_header)
				{
				this.printWriter.print("#SAMPLE");
				this.printWriter.print('\t');
				this.printWriter.print("TYPE");
				this.printWriter.print('\t');
				this.printWriter.print("FILE");
				this.printWriter.print('\t');
				this.printWriter.print("FILE_SIZE");
				this.printWriter.print('\t');
				this.printWriter.print("DATE");
				this.printWriter.println();
				}
			
			if(args.isEmpty())
				{
				scan(new BufferedReader(new InputStreamReader(stdin())));
				}
			else
				{
				for(final String filename:args)
					{
					final BufferedReader r=IOUtils.openURIForBufferedReading(filename);
					scan(r);
					r.close();
					}
				}
			this.printWriter.flush();
			this.printWriter.close();
			this.printWriter=null;
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.printWriter);
			}
		}
	
	public static void main(final String[] args) {
		new NgsFilesSummary().instanceMainWithExit(args);

	}

}
