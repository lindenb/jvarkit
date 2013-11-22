package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.InputStream;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;

import org.broadinstitute.variant.vcf.VCFHeader;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class NgsFilesSummary extends AbstractCommandLineProgram
	{
	private static Log LOG=Log.getInstance(NgsFilesSummary.class);

    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + " Scan folders and generate a summary of the files (SAMPLE/BAM SAMPLE/VCF etc..). ";

	
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,doc="File and/or directories",minElements=0)
    public Set<File> IN=new HashSet<File>();

    private enum InfoType { BAM,FASTQ,VCF};
    
  
    
    
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
			r.setValidationStringency(super.VALIDATION_STRINGENCY);
			SAMFileHeader h=r.getFileHeader();
			for(SAMReadGroupRecord rg: h.getReadGroups())
				{
				String sample=rg.getSample();
				if(sample==null || sample.isEmpty()) continue;
				print(sample,InfoType.BAM, f);
				}
			} 
    	catch (Exception e)
    		{
    		LOG.warn(e, "Error in "+f);
			}
    	finally
    		{
    		if(r!=null) r.close();
    		}
    	}
    
    private final Pattern uscore=Pattern.compile("_");
    private void readFastq(File f)
		{
    	//File parent=f.getParentFile();
    	//if(parent==null || super.VERBOSITY==Log.LogLevel.) return;
    	
    	
    	String tokens[]=this.uscore.split(f.getName());
		
		if(tokens.length<5)
			{
			//bad name
			return;
			}
		else if(tokens.length>5)
			{
			String tokens2[]=new String[5];
			tokens2[0]=tokens[0];
			int name_count=(tokens.length-5);
			for(int i=1;i<= name_count;++i)
				{
				tokens2[0]+="_"+tokens[i];
				}
			for(int i=name_count+1; i< tokens.length;++i)
				{
				tokens2[i-name_count]=tokens[i];
				}
			tokens=tokens2;
			}

		if(tokens[0].equalsIgnoreCase("Undetermined") ||  tokens[1].equalsIgnoreCase("Undetermined"))
			{
			return;
			}
    	print(tokens[0],InfoType.FASTQ, f);
    	
		}
    
    private void readVCF(File f)
		{
    	if(!f.canRead()) return;
    	LOG.debug("readVCF  "+f);
    	    	

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
    		LOG.info(err,"Error in VCF "+f);
    		}
    	finally
    		{
    		CloserUtil.close(r);
    		CloserUtil.close(in);
    		}
    	
		}

    
    
	private void analyze(File f)
		{
		if(f==null) return;
		if(!f.canRead()) return;
		LOG.debug("Scanning "+f);
		if(f.isDirectory())
			{
			File children[]=f.listFiles();
			if(children==null) return;
			for(File c:f.listFiles())
				{
				analyze(c);
				}
			}
		else
			{
			
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
		new NgsFilesSummary().instanceMainWithExit(args);

	}

}
