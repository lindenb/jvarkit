package com.github.lindenb.jvarkit.tools.ngsfiles;

import java.io.File;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

public abstract class AbstractScanNgsFilesProgram  extends AbstractCommandLineProgram
	{
    protected enum InfoType { BAM,FASTQ,VCF};
    
    protected AbstractScanNgsFilesProgram()
    	{
    	
    	}		
 
    	
    protected boolean isAcceptBam()
    	{
    	return true;
    	}
    protected boolean isAcceptFastq()
		{
		return true;
		}
    
    protected boolean isAcceptVcf()
		{
		return true;
		}
	
    protected void readBam(File f)
    	{
    	}
    
    protected void readFastq(File f)
		{
		}
    
    protected void readVCF(File f)
		{
    	
		}

    
	protected void analyze(File f)
		{
		if(f==null) return;
		if(!f.canRead() || !f.exists() || f.isDirectory()) return;
		debug("Scanning "+f);
		
		
		String name=f.getName();
		if(isAcceptBam() && (name.endsWith(".bam") || name.endsWith(".sam")))
			{
			readBam(f);
			}
		else if(isAcceptVcf() && (name.endsWith(".vcf.gz") || name.endsWith(".vcf")))
			{
			readVCF(f);
			}
		else if(isAcceptFastq() && (name.endsWith(".fastq") || name.endsWith(".fastq.gz") ||
				name.endsWith(".fq") || name.endsWith(".fq.gz")))
			{
			readFastq(f);
			}
			
		}
	
	
	}
