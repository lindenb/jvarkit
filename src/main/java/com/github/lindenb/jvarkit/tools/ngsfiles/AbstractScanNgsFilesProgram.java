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

import java.io.File;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;

import htsjdk.samtools.util.FileExtensions;

public abstract class AbstractScanNgsFilesProgram  extends Launcher
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
		
		final String name=f.getName();
		if(isAcceptBam() && (name.endsWith(FileExtensions.SAM) || name.endsWith(FileExtensions.BAM)))
			{
			readBam(f);
			}
		else if(isAcceptVcf() && FileExtensions.VCF_LIST.stream().anyMatch(E->name.endsWith(E)))
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
