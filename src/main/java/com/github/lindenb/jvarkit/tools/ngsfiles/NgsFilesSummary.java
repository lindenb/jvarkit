/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Date;
import java.util.List;
import java.util.Set;

import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

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
	keywords= {"sam","bam","vcf","util"},
	modificationDate="20190906"
	)
public class NgsFilesSummary extends Launcher
	{
	private static final Logger LOG = Logger.build(NgsFilesSummary.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-header","--header"},description="[20180725]print header")
	private boolean show_header= false;
	@Parameter(names={"-R","--reference"},description="[20190905]restrict to that reference. Also is used to read CRAM files")
	private Path faidxPath = null;
	@Parameter(names={"-i","--indexed"},description="[20190905]VCF or BAM must be indexed")
	private boolean must_be_indexed=false;
	@Parameter(names={"-p","--partition"},description="For BAM files: "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	@Parameter(names={"--no-read-group"},description="Flag form SAM/VCF without read group/ samples")
	private String noReadGroupFlag = "__NO_READ_GROUP__";

	
	private PrintWriter printWriter=null;
	private SAMSequenceDictionary dict = null;
	
    public NgsFilesSummary()
    	{
    	}		
    
  
    
    private void print(final String sample,final String type,final Path f,final String index)
    	{
		this.printWriter.print(StringUtils.isBlank(sample)?this.noReadGroupFlag:sample);
		this.printWriter.print('\t');
		this.printWriter.print(type);
		this.printWriter.print('\t');
		this.printWriter.print(f.toAbsolutePath());
		this.printWriter.print('\t');
		this.printWriter.print(index);

		if(Files.isRegularFile(f))
			{
			long size;
			Date modif;
			try {
				size=Files.size(f);
				}
			catch(IOException err) {
				size=-1L;
				}
			try {
				modif=new Date(Files.getLastModifiedTime(f).toMillis());
				}
			catch(final IOException err) {
				modif = null;
				}
			
			this.printWriter.print('\t');
			this.printWriter.print(size);
			this.printWriter.print('\t');
			this.printWriter.print(modif);
			}
		this.printWriter.println();
		}
    	
    
    private void readBam(final Path f)
    	{
		final SamReaderFactory srf = super.createSamReaderFactory();
		srf.validationStringency(ValidationStringency.SILENT);
		if(this.faidxPath!=null) srf.referenceSequence(this.faidxPath);
    	SamReader r=null;
    	try {
			r = srf.open(f);
			boolean indexed = r.hasIndex();
    		if(this.must_be_indexed && !indexed) return;

			
			final SAMFileHeader h=r.getFileHeader();
    		if(this.dict!=null) {
    			SAMSequenceDictionary dict2 = h.getSequenceDictionary();
    			if(!SequenceUtil.areSequenceDictionariesEqual(this.dict, dict2)) return ;
    			}

			
    		final Set<String> names = this.partition.getPartitions(h);
    		
			if(!names.isEmpty())
				{
				for(final String sample:names)
					{
					print(sample,"BAM", f,String.valueOf(r.hasIndex()));
					}
				}
			else
				{
				print(null,"BAM", f,String.valueOf(r.hasIndex()));
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
   
    private void readFastq(final Path f)
		{
    	//File parent=f.getParentFile();
    	//if(parent==null || super.VERBOSITY==Log.LogLevel.) return;
    	final FastQName fq=FastQName.parse(f.toFile());
    	
		
		if(!fq.isValid())
			{
			//bad name
			return;
			}
		
    	print(fq.getSample(),"FASTQ", f,".");
		}
    
    private void readVCF(final Path f)
		{
    	VCFReader r=null;
    	InputStream in=null;
    	try
    		{
    		boolean indexed= VCFUtils.isTabixVcfPath(f) || VCFUtils.isTribbleVcfPath(f);
    		if(this.must_be_indexed && !indexed) return;
    		
    		in=IOUtils.openPathForReading(f);
    		
    		r= VCFReaderFactory.makeDefault().open(f,false);
        	final VCFHeader header=r.getHeader();
        	
    		if(this.dict!=null) {
    			SAMSequenceDictionary dict2 = header.getSequenceDictionary();
    			if(dict2==null) return;
    			if(!SequenceUtil.areSequenceDictionariesEqual(this.dict, dict2)) return ;
    			}
    		final List<String> sns = header.getSampleNamesInOrder();
        	
    		if(sns.isEmpty())
    			{
    			print(null,"VCF", f,String.valueOf(indexed));
    			}
    		else
	    		{
	        	for(final String sample:sns)
		        	{
		        	print(sample,"VCF", f,String.valueOf(indexed));
		    		}
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
    	in.lines().
    		filter(L->!(StringUtil.isBlank(L) || L.startsWith("#"))).
    		map(L->Paths.get(L)).
    		filter(L->Files.isRegularFile(L)).
    		filter(L->Files.isReadable(L)).
    		forEach(L->scan(L));
    	}
    
    private void scan(final Path f)
		{
		if(f==null) return;
		if(!Files.isRegularFile(f)) return;		
		if(!Files.isReadable(f)) return;		
		
		final String name=f.getFileName().toString();
		if(name.endsWith(FileExtensions.SAM)  && !this.must_be_indexed)
			{
			readBam(f);
			}
		else if(name.endsWith(FileExtensions.BAM))
			{
			readBam(f);
			}
		else if(name.endsWith(FileExtensions.CRAM) && this.dict!=null) {
			final SAMSequenceDictionary dict2=SAMSequenceDictionaryExtractor.extractDictionary(f);
			if(dict2==null) return;
			if(!SequenceUtil.areSequenceDictionariesEqual(this.dict, dict2)) return;
			readBam(f);
			}
		else if(FileExtensions.VCF_LIST.stream().anyMatch(E->name.endsWith(E)))
			{
			readVCF(f);
			}
		else if( (name.endsWith(".fastq") || name.endsWith(".fastq.gz") ||
				name.endsWith(".fq") || name.endsWith(".fq.gz")))
			{
			readFastq(f);
			}
		}

    @Override
	public int doWork(final List<String> args) {
		try
			{
			if(this.faidxPath!=null) {
				this.dict = SequenceDictionaryUtils.extractRequired(this.faidxPath);
			}
			
			this.printWriter = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			if(show_header)
				{
				this.printWriter.print("#"+this.partition.name());
				this.printWriter.print('\t');
				this.printWriter.print("TYPE");
				this.printWriter.print('\t');
				this.printWriter.print("FILE");
				this.printWriter.print('\t');
				this.printWriter.print("INDEXED");
				this.printWriter.print('\t');
				this.printWriter.print("FILE_SIZE");
				this.printWriter.print('\t');
				this.printWriter.print("DATE");
				this.printWriter.println();
				}
			
			if(args.isEmpty())
				{
				try(BufferedReader r = new BufferedReader(new InputStreamReader(stdin()))) {
					scan(r);
					}
				}
			else
				{
				for(final String filename:args)
					{
					try(final BufferedReader r=IOUtils.openURIForBufferedReading(filename)) {
						scan(r);
						}
					}
				}
			this.printWriter.flush();
			this.printWriter.close();
			this.printWriter=null;
			return 0;
			}
		catch(final Throwable err)
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
