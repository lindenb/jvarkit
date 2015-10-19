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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Collection;
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

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class FindCorruptedFiles extends AbstractFindCorruptedFiles
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(FindCorruptedFiles.class);



	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractFindCorruptedFiles.AbstractFindCorruptedFilesCommand
	 	{		
		private PrintStream outputstream=null;
		private ValidationStringency validationStringency=ValidationStringency.LENIENT;
    
    
    private void emptyFile(File f)
    	{
    	if(super.emptyIsError)
    		{
    		this.outputstream.println(f);
    		}
    	else
    		{
    		LOG.warn("Empty content:"+f);
    		}
    	}
	
    private void testBam(File f)
    	{
    	LOG.info("Test BAM for "+f);
    	
    	//Test BGZ-EOF
    	try
    		{
    		BlockCompressedInputStream.FileTermination type=BlockCompressedInputStream.checkTermination(f);
    		if(type!=BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK)
    			{
    			LOG.warn("bgz:"+type+" for "+f);
        		this.outputstream.println(f);
        		return;
    			}
    		}
    	catch(IOException err)
    		{
    		LOG.warn("Error in "+f);
    		this.outputstream.println(f);
    		return;
    		}
    	
    	
    	long n=0L;
    	SamReader r=null;
    	SAMRecordIterator iter=null;
    	try {
			r=SamFileReaderFactory.mewInstance().stringency(this.validationStringency).open(f);
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
    	catch (Exception e)
    		{
    		LOG.warn( "Error in "+f);
    		this.outputstream.println(f);
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
					LOG.warn( "BED error. Line "+(n+1)+" not enough columns in "+f);
					this.outputstream.println(f);
					break;
					}
				if(tokens[0].trim().isEmpty())
					{
					LOG.warn( "BED error. Line "+(n+1)+" Bad chrom in "+f);
					this.outputstream.println(f);
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
					LOG.warn( "BED error. Line "+(n+1)+" Bad start in "+f);
					this.outputstream.println(f);
					break;
					}
				
				if(start<0)
					{
					LOG.warn( "BED error. Line "+(n+1)+" Bad start in "+f);
					this.outputstream.println(f);
					break;
					}
				
				try
					{
					end=Integer.parseInt(tokens[2]);
					}
				catch(NumberFormatException err)
					{
					LOG.warn( "BED error. Line "+(n+1)+" Bad end in "+f);
					this.outputstream.println(f);
					break;
					}
				if(end<start)
					{
					LOG.warn( "BED error. Line "+(n+1)+" end<start in "+f);
					this.outputstream.println(f);
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
    		LOG.warn( "Error in "+f);
			this.outputstream.println(f);
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}


    
    private void testVcf(File f,InputStream in) throws IOException,TribbleException
    	{
    	long n=0;
    	VcfIterator iter=new VcfIterator(in);
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
    	LOG.info("Test Fastq for "+f);
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
    		LOG.info("Cannot read "+f);
    		this.outputstream.println(f);
    		}
    	finally
    		{
    		CloserUtil.close(r);
    		}
		}
    
    private void testVcf(File f)
		{
    	LOG.info("Test VCF for "+f);
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
	    			LOG.warn("bgz:"+type+" for "+f);
	        		this.outputstream.println(f);
	        		return;
	    			}
	    		}
	    	catch(Exception err)
	    		{
	    		this.outputstream.println(f);
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
    		LOG.warn("gzip but not bgzip :"+f);
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
    	LOG.info("Not a BGZIP file / Error in VCF: "+f);
    	this.outputstream.println(f);
		}

    private void testTbi(File f)
    	{
		String filename=f.getName();
		//remove .tbi suffix
		filename=filename.substring(0,
				filename.length()-TabixUtils.STANDARD_INDEX_EXTENSION.length());

    	File baseFile=new File(f.getParentFile(),filename);
    	if(baseFile.exists() && baseFile.isFile()) return;
    	LOG.info("Missing associated file for : "+f);
    	this.outputstream.println(f);		
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
    	LOG.info("Missing associated BAM file for : "+f);
    	this.outputstream.println(f);		
		}
    
	private void analyze(File f)
		{
		LOG.info("Scanning "+f);
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
	public Collection<Throwable> call(String inputName) throws Exception
			{
			if(super.validationStringencyStr!=null)
				{
				this.validationStringency=ValidationStringency.valueOf(super.validationStringencyStr);

				}
			
			BufferedReader in=null;
			String line;
			this.outputstream =null;
			try
				{
				this.outputstream =  openFileOrStdoutAsPrintStream();
				if(inputName==null)
					{
					in = new BufferedReader(new InputStreamReader(stdin()));
					}
				else
					{
					in = IOUtils.openURIForBufferedReading(inputName);
					}
				LOG.info("reading from stdin");
				
				while((line=in.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					File f=new File(line);
					analyze(f);
					}
				this.outputstream.flush();
				return RETURN_OK;
				}
			catch(IOException err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(in);
				CloserUtil.close(this.outputstream);
				}
			}
	 	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FindCorruptedFiles().instanceMainWithExit(args);

	}

}
