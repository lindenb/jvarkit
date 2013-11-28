package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.zip.GZIPInputStream;

import org.broad.tribble.TribbleException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class FindCorruptedFiles extends AbstractCommandLineProgram
	{

    //public String USAGE = " Reads filename from stdin and search for corrupted NGS files (VCF/BAM/FASTQ). Prints file with problem on stdout.";


    //@Option(shortName="N",doc="number of features (samrecord, variant) to read. -1= read everything. ",optional=true)
	private long NUM=100L;
    
    private ValidationStringency validationStringency=ValidationStringency.LENIENT;
    
	
    private void testBam(File f)
    	{
    	getLogger().fine("Test BAM for "+f);
    	
    	//Test BGZ-EOF
    	try
    		{
    		BlockCompressedInputStream.FileTermination type=BlockCompressedInputStream.checkTermination(f);
    		if(type!=BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK)
    			{
    			getLogger().warning("bgz:"+type+" for "+f);
        		System.out.println(f);
        		return;
    			}
    		}
    	catch(IOException err)
    		{
    		getLogger().warning("Error in "+f);
    		System.out.println(f);
    		return;
    		}
    	
    	
    	long n=0L;
    	SAMFileReader r=null;
    	SAMRecordIterator iter=null;
    	try {
			r=new SAMFileReader(f);
			r.setValidationStringency(this.validationStringency);
			r.getFileHeader();
			iter=r.iterator();
			while(iter.hasNext() && (NUM<0 || n<NUM))
				{
				iter.next();
				++n;
				}
			if(n==0)
				{
				getLogger().warning("No SAM record in "+f);
				}
			} 
    	catch (Exception e)
    		{
    		getLogger().warning( "Error in "+f);
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
			getLogger().warning("No Variant in "+f);
			}
    	iter.close();
    	}
    
    private void testFastq(File f)
		{
    	long n=0;
    	getLogger().fine("Test Fastq for "+f);
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
    		}
    	catch(Exception err)
    		{
    		getLogger().fine("Cannot read "+f);
    		System.out.println(f);
    		}
    	finally
    		{
    		CloserUtil.close(r);
    		}
		}
    
    private void testVcf(File f)
		{
    	getLogger().fine("Test VCF for "+f);
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
	    			getLogger().warning("bgz:"+type+" for "+f);
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
    		getLogger().warning("gzip but not bgzip :"+f);
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
    	getLogger().fine("Not a BGZIP file / Error in VCF: "+f);
    	System.out.println(f);
		}

    
    
	private void analyze(File f)
		{
		getLogger().fine("Scanning "+f);
		if(f.isDirectory())
			{
			getLogger().info("skipping "+f+" (directory)");
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
	public String getProgramDescription() {
		return "Reads filename from stdin and prints corrupted NGS files (VCF/BAM/FASTQ)";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FindCorruptedFiles";
	}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -h get help (this screen)");
		
		out.println(" -V (validation stringency). Optional");
		out.println(" -N (number). Number of records to test in each file (BAM, FASTQ, VCF...). ptional");

		
		out.println(" -v print version and exit.");
		out.println(" -L (level) log level. One of java.util.logging.Level . currently:"+getLogger().getLevel());
		out.println(" -r (file)  fasta sequence file indexed with picard. Required.");
		}

	
	@Override
	public int doWork(String args[])
		{
		GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args, "vhN:V:N:L:"))!=-1)
			{
			switch(c)
				{
				case 'v':
					{
					System.out.println(getVersion());
					return 0;
					}
				case 'h':
					{
					printUsage();
					return 0;
					}
				case 'N':
					{
					this.NUM=Integer.parseInt(getopt.getOptArg());
					break;
					}
				case 'V':
					{
					this.validationStringency=ValidationStringency.valueOf(getopt.getOptArg());
					getLogger().info("setting validation stringency to "+this.validationStringency);
					break;
					}
				case 'L':
					{
					getLogger().setLevel(Level.parse(getopt.getOptArg()));
					break;
					}
				default:
					{
					System.err.println("Unknown Option:"+getopt.getOptOpt()+" or missing argument");
					return -1;
					}
				}
			}
		if(getopt.getOptInd()!=args.length)
			{
			System.err.println("Too many arguments");
			return -1;
			}
		getLogger().info("reading from stdin");
		BufferedReader in=new BufferedReader(new InputStreamReader(System.in));
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
		catch(IOException err)
			{
			err.printStackTrace();
			getLogger().severe("I/O Error "+err.getMessage());
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
