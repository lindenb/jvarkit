package com.github.lindenb.jvarkit.tools.ngsfiles;

import java.io.File;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

public abstract class AbstractScanNgsFilesProgram  extends Launcher
	{
	private static final Logger LOG = Logger.build(AbstractScanNgsFilesProgram.class).make();

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
	
    protected void readBam(final File f)
    	{
    	}
    
    protected void readFastq(final File f)
		{
		}
    
    protected void readVCF(final File f)
		{
    	
		}

    
	protected void analyze(final File f)
		{
		if(f==null) return;
		if(!f.canRead() || !f.exists() || f.isDirectory()) return;
		LOG.debug("Scanning "+f);
		
		
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
