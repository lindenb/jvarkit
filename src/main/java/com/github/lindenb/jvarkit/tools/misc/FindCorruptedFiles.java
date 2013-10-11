package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.broad.tribble.TribbleException;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.BlockCompressedInputStream;

import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class FindCorruptedFiles extends AbstractCommandLineProgram
	{
	private static Log LOG=Log.getInstance(FindCorruptedFiles.class);

    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + " Search for corrupted NGS files (VCF/BAM/FASTQ) ";

	
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,doc="File and/or directories",minElements=0)
    public Set<File> IN=new HashSet<File>();

    @Option(shortName="N",doc="number of features (samrecord, variant) to read. -1= read everything. ",optional=true)
    public long NUM=100L;

    
    private FileFilter ngsFileFilter=new FileFilter()
    	{
		@Override
		public boolean accept(File f)
			{
			if(f.isDirectory()) return true;
			if(!f.isFile())  return false;
			if(f.isHidden()) return false;
			if(!f.canRead()) return false;
			String name=f.getName();
			if(name.endsWith(".bam")) return true;
			if(name.endsWith(".vcf.gz")) return true;
			if(name.endsWith(".fastq")) return true;
			if(name.endsWith(".fastq.gz")) return true;
			if(name.endsWith(".fq")) return true;
			if(name.endsWith(".fq.gz")) return true;
			return false;
			}
    	};
	
    private void testBam(File f)
    	{
    	LOG.debug("Test BAM for "+f);
    	
    	//Test BGZ-EOF
    	try
    		{
    		BlockCompressedInputStream.FileTermination type=BlockCompressedInputStream.checkTermination(f);
    		if(type!=BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK)
    			{
    			LOG.warn("bgz:"+type+" for "+f);
        		System.out.println(f);
        		return;
    			}
    		}
    	catch(IOException err)
    		{
    		LOG.warn(err, "Error in "+f);
    		System.out.println(f);
    		return;
    		}
    	
    	
    	long n=0L;
    	SAMFileReader r=null;
    	SAMRecordIterator iter=null;
    	try {
			r=new SAMFileReader(f);
			r.setValidationStringency(super.VALIDATION_STRINGENCY);
			r.getFileHeader();
			iter=r.iterator();
			while(iter.hasNext() && (NUM<0 || n<NUM))
				{
				iter.next();
				++n;
				}
			if(n==0)
				{
				LOG.warn("No SAM record in "+f);
				}
			} 
    	catch (Exception e)
    		{
    		LOG.warn(e, "Error in "+f);
    		System.out.println(f);
			}
    	finally
    		{
    		if(iter!=null) iter.close();
    		if(r!=null) r.close();
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
			LOG.warn("No Variant in "+f);
			}
    	iter.close();
    	}
    
    private void testFastq(File f)
		{
    	long n=0;
    	LOG.debug("Test Fastq for "+f);
    	FastqReader r=null;
    	try
    		{
    		r=new FastqReader(f);
    		while(r.hasNext() && (NUM<0 || n<NUM))
    			{
    			r.next();
    			++n;
    			}
    		}
    	catch(Exception err)
    		{
    		LOG.info(err,"Cannot read "+f);
    		System.out.println(f);
    		}
    	finally
    		{
    		if(r!=null) r.close();
    		}
		}
    
    private void testVcf(File f)
		{
    	LOG.debug("Test VCF for "+f);
    	Exception error1=null;
    	Exception error2=null;
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
	        		System.out.println(f);
	        		return;
	    			}
	    		}
	    	catch(Exception err)
	    		{
	    		System.out.println(f);
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
    		error2=err;
    		}
    	finally
    		{
    		if(in2!=null) try{ in2.close();} catch(IOException e) {}
    		}
    	LOG.debug(error1,"Not a BGZIP file / Error in VCF: "+f);
    	LOG.debug(error2,"Not a GZIP file "+f);
    	System.out.println(f);
		}

    
    
	private void analyze(File f)
		{
		LOG.debug("Scanning "+f);
		if(f.isDirectory())
			{
			for(File c:f.listFiles(ngsFileFilter))
				{
				analyze(c);
				}
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
			}
		}
	@Override
	protected int doWork()
		{
		for(File in:IN) analyze(in);
		return 0;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FindCorruptedFiles().instanceMainWithExit(args);

	}

}
