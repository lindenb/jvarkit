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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.zip.GZIPInputStream;

import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfIteratorImpl;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

/*

BEGIN_DOC

## Example

```
$ find  DIR1 DIR2 -type f |\
java -jar dist/findcorruptedfiles.jar \
	-V SILENT 2> /dev/null > redo.txt
```

END_DOC

 */

@Program(name="findcorruptedfiles",keywords={"vcf","bam","fastq","bed"},description="Reads filename from stdin and prints corrupted NGS files (VCF/BAM/FASTQ/BED/TBI/BAI)")
public class FindCorruptedFiles extends Launcher
	{
	private static final Logger LOG = Logger.build(FindCorruptedFiles.class).make();
	@Parameter(names={"-E","--noempty"},description="empty file is an error.")
	private boolean emptyIsError=false;

	@Parameter(names={"-N"},description="number of features (samrecord, variant) to read. -1= read everything.")
	private long NUM=100L;
    
	@Parameter(names={"-V","--stringency"},description="BAM ValidationStringency")
    private ValidationStringency validationStringency=ValidationStringency.LENIENT;
    
    
    private void emptyFile(final File f)
    	{
    	if(this.emptyIsError)
    		{
    		stdout().println(f);
    		}
    	else
    		{
    		LOG.warning("Empty content:"+f);
    		}
    	}
	
    private void testBam(File f)
    	{
    	LOG.fine("Test BAM for "+f);
    	
    	//Test BGZ-EOF
    	try
    		{
    		BlockCompressedInputStream.FileTermination type=BlockCompressedInputStream.checkTermination(f);
    		if(type!=BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK)
    			{
    			LOG.warning("bgz:"+type+" for "+f);
        		stdout().println(f);
        		return;
    			}
    		}
    	catch(IOException err)
    		{
    		LOG.warning("Error in "+f);
    		stdout().println(f);
    		return;
    		}
    	
    	
    	long n=0L;
    	SamReader r=null;
    	SAMRecordIterator iter=null;
    	try {
			r= super.createSamReaderFactory().validationStringency(this.validationStringency).open(f);
			r.getFileHeader();
			iter=r.iterator();
			while(iter.hasNext() && (NUM<0 || n<NUM))
				{
				iter.next();
				++n;
				}
			if(n==0)
				{
				emptyFile(f);
				}
			} 
    	catch (final Exception e)
    		{
    		LOG.warning( "Error in "+f);
    		stdout().println(f);
			}
    	finally
    		{
    		CloserUtil.close(iter);
    		CloserUtil.close(r);
    		}
    	}
    
    private void testBed(File f)
		{
    	LineIterator iter=null;
    	try
	    	{
    		java.util.regex.Pattern tab=java.util.regex.Pattern.compile("[\t]");
			long n=0;
		
			iter=IOUtils.openFileForLineIterator(f);
			while(iter.hasNext() && (NUM<0 || n<NUM))
				{
				String line=iter.next();
				String tokens[]=tab.split(line);
				if(tokens.length<3)
					{
					LOG.warning( "BED error. Line "+(n+1)+" not enough columns in "+f);
					stdout().println(f);
					break;
					}
				if(tokens[0].trim().isEmpty())
					{
					LOG.warning( "BED error. Line "+(n+1)+" Bad chrom in "+f);
					stdout().println(f);
					break;
					}
				int start=0;
				int end=0;
				
				try
					{
					start=Integer.parseInt(tokens[1]);
					}
				catch(NumberFormatException err)
					{
					LOG.warning( "BED error. Line "+(n+1)+" Bad start in "+f);
					stdout().println(f);
					break;
					}
				
				if(start<0)
					{
					LOG.warning( "BED error. Line "+(n+1)+" Bad start in "+f);
					stdout().println(f);
					break;
					}
				
				try
					{
					end=Integer.parseInt(tokens[2]);
					}
				catch(NumberFormatException err)
					{
					LOG.warning( "BED error. Line "+(n+1)+" Bad end in "+f);
					stdout().println(f);
					break;
					}
				if(end<start)
					{
					LOG.warning( "BED error. Line "+(n+1)+" end<start in "+f);
					stdout().println(f);
					break;
					}


				++n;
				}
			if(n==0)
				{
				emptyFile(f);
				}
	    	}
    	catch (Exception e)
			{
    		LOG.warning( "Error in "+f);
			stdout().println(f);
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}


    
    private void testVcf(final File f,InputStream in) throws IOException,TribbleException
    	{
    	long n=0;
    	final VCFIterator iter= VCFUtils.createVCFIteratorFromStream(in);
    	iter.getHeader();
    	while(iter.hasNext() &&  (NUM<0 || n<NUM))
    		{
    		iter.next();
    		++n;
    		}
    	if(n==0)
			{
    		emptyFile(f);
			}
    	iter.close();
    	}
    
    private void testFastq(File f)
		{
    	long n=0;
    	LOG.fine("Test Fastq for "+f);
    	FastqReader r=null;
    	try
    		{
    		r=new FourLinesFastqReader(f);
    		r.setValidationStringency(this.validationStringency);
    		while(r.hasNext() && (NUM<0 || n<NUM))
    			{
    			r.next();
    			++n;
    			}
    		if(n==0)
    			{
    			emptyFile(f);
    			}
    		}
    	catch(Exception err)
    		{
    		LOG.fine("Cannot read "+f);
    		stdout().println(f);
    		}
    	finally
    		{
    		CloserUtil.close(r);
    		}
		}
    
    private void testVcf(File f)
		{
    	LOG.fine("Test VCF for "+f);
    	Exception error1=null;
    	BlockCompressedInputStream in1=null;
    	try
    		{
    		in1=new BlockCompressedInputStream(f);
    		testVcf(f,in1);
    		}
    	catch(Exception err)
    		{
    		error1=err;
    		}
    	finally
    		{
    		if(in1!=null) try{ in1.close();in1=null;} catch(IOException e) {}
    		}
    	if(error1==null)
    		{
    		try
	    		{
	    		BlockCompressedInputStream.FileTermination type=BlockCompressedInputStream.checkTermination(f);
	    		if(type!=BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK)
	    			{
	    			LOG.warning("bgz:"+type+" for "+f);
	        		stdout().println(f);
	        		return;
	    			}
	    		}
	    	catch(Exception err)
	    		{
	    		stdout().println(f);
	    		return;
	    		}
    		
    		//TEST BGZF termination
    		return;
    		}
    	
    	
    	GZIPInputStream in2=null;
    	try
    		{
    		in2=new GZIPInputStream(new FileInputStream(f));
    		testVcf(f,in2);
    		LOG.warning("gzip but not bgzip :"+f);
    		return;
    		}
    	catch(Exception err)
    		{
    		err.printStackTrace();
    		}
    	finally
    		{
    		if(in2!=null) try{ in2.close();} catch(IOException e) {}
    		}
    	LOG.fine("Not a BGZIP file / Error in VCF: "+f);
    	stdout().println(f);
		}

    private void testTbi(File f)
    	{
		String filename=f.getName();
		//remove .tbi suffix
		filename=filename.substring(0,
				filename.length()-TabixUtils.STANDARD_INDEX_EXTENSION.length());

    	File baseFile=new File(f.getParentFile(),filename);
    	if(baseFile.exists() && baseFile.isFile()) return;
    	LOG.fine("Missing associated file for : "+f);
    	stdout().println(f);		
		}
    
    private void testBai(File f)
		{
		String filename=f.getName();
		//remove .bai suffix
		filename=filename.substring(0,
				filename.length()-BAMIndex.BAMIndexSuffix.length());

		//ends with bam.bai
		if(filename.endsWith(BamFileIoUtils.BAM_FILE_EXTENSION))
			{
			//nothing
			}
		else // file.bai and file.bam
			{
			filename+=BamFileIoUtils.BAM_FILE_EXTENSION;
			}
    	File bam=new File(f.getParentFile(),filename);
    	if(bam.exists() && bam.isFile()) return;
    	LOG.fine("Missing associated BAM file for : "+f);
    	stdout().println(f);		
		}
    
	private void analyze(File f)
		{
		LOG.fine("Scanning "+f);
		if(f.isDirectory())
			{
			LOG.info("skipping "+f+" (directory)");
			}
		else
			{	
			String name=f.getName();
			if(name.endsWith(".bam"))
				{
				testBam(f);
				}
			else if(name.endsWith(".vcf.gz"))
				{
				testVcf(f);
				}
			else if(name.endsWith(".fastq") || name.endsWith(".fastq.gz") || name.endsWith(".fq") || name.endsWith(".fq.gz"))
				{
				testFastq(f);
				}
			else if(name.endsWith(".bed") || name.endsWith(".bed.gz"))
				{
				testBed(f);
				}
			else if(name.endsWith(TabixUtils.STANDARD_INDEX_EXTENSION))
				{
				testTbi(f);
				}
			else if(name.endsWith(BAMIndex.BAMIndexSuffix))
				{
				testBai(f);
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(!args.isEmpty())
			{
			LOG.fatal("Too many arguments");
			return -1;
			}
		LOG.info("reading from stdin");
		BufferedReader in=new BufferedReader(new InputStreamReader(stdin()));
		String line;
		try
			{
			while((line=in.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				File f=new File(line);
				analyze(f);
				}
			}
		catch(final IOException err)
			{
			err.printStackTrace();
			LOG.severe("I/O Error "+err.getMessage());
			return -1;
			}
		
		return 0;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FindCorruptedFiles().instanceMainWithExit(args);

	}

}
