package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Date;

import org.broadinstitute.variant.vcf.VCFHeader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class NgsFilesSummary extends AbstractCommandLineProgram
	{

    private enum InfoType { BAM,FASTQ,VCF};
    
    private NgsFilesSummary()
    	{
    	
    	}		
    @Override
    protected String getOnlineDocUrl() {
    	return "https://github.com/lindenb/jvarkit/wiki/NgsFilesSummary";
    	}
    
    @Override
    public String getProgramDescription() {
    	return "Scan folders and generate a summary of the files (SAMPLE/BAM SAMPLE/VCF etc..)";
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
    	
    
	
    private void readBam(File f)
    	{
    	if(!f.canRead()) return;
    	SAMFileReader r=null;
    	try {
			r=new SAMFileReader(f);
			r.setValidationStringency(ValidationStringency.LENIENT);
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
    		warning(e, "Error in "+f);
			}
    	finally
    		{
    		if(r!=null) r.close();
    		}
    	}
    
    private void readFastq(File f)
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
    
    private void readVCF(File f)
		{
    	if(!f.canRead()) return;
    	debug("readVCF  "+f);
    	    	

    	VcfIterator r=null;
    	InputStream in=null;
    	try
    		{
    		in=IOUtils.openFileForReading(f);
    		
    		r=new VcfIterator(in);
        	VCFHeader header=r.getHeader();
        	for(String sample:header.getSampleNamesInOrder())
	        	{
	        	print(sample,InfoType.VCF, f);
	    		}
        	
    		}
    	catch(Exception err)
    		{
    		error(err,"Error in VCF "+f);
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
    
	private void analyze(File f)
		{
		if(f==null) return;
		if(!f.canRead() || !f.exists() || f.isDirectory()) return;
		debug("Scanning "+f);
		
		
		String name=f.getName();
		if(name.endsWith(".bam") || name.endsWith(".sam"))
			{
			readBam(f);
			}
		else if(name.endsWith(".vcf.gz") || name.endsWith(".vcf"))
			{
			readVCF(f);
			}
		else if(name.endsWith(".fastq") || name.endsWith(".fastq.gz") ||
				name.endsWith(".fq") || name.endsWith(".fq.gz"))
			{
			readFastq(f);
			}
			
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()))!=-1)
			{
			switch(c)
				{
				default:
					{
					switch(super.handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS:return 0;
						case OK:break;
						}
					}
				}
			}
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				scan(new BufferedReader(new InputStreamReader(System.in)));
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("Reading from "+filename);
					BufferedReader r=IOUtils.openURIForBufferedReading(filename);
					scan(r);
					r.close();
					}
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
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
