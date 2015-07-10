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
package com.github.lindenb.jvarkit.util.picard;

import java.io.File;
import java.io.PrintStream;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.cli.GetOpt;

/**
 * @author lindenb
 * A program producing a BAM/SAM
 *
 */
public abstract class AbstractBamWriterProgram extends AbstractCommandLineProgram
	{
	private File outputFile=null;
	private SamFormat outputFormat=null;
	protected AbstractBamWriterProgram()
		{
		
		}
	public File getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(File outputFile)
		{
		this.outputFile = outputFile;
		}
	
	public SamFormat getOutputSamFormat()
		{
		if(outputFormat!=null) return outputFormat;
		if( getOutputFile()!=null)
			{
			SamFormat sf = SamFormat.getSamFormatFromFile(getOutputFile());
			if(sf!=null) return sf;
			}
		return SamFormat.sam;
		}
	
	public void setOutputSamFormat(SamFormat outputFormat) {
		this.outputFormat = outputFormat;
		}
	
	protected SAMFileWriterFactory createSAMFileWriterFactory()
		{
		SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
		sfwf.setCreateIndex(false);
		sfwf.setCreateMd5File(false);
		File o = getOutputFile();
		if(o!=null && o.getParentFile()!=null)
			{
			sfwf.setTempDirectory(o.getParentFile());
			}
		return sfwf;
		}
	
	protected SAMFileWriter openSAMFileWriter(SAMFileHeader header,boolean presorted)
		{
		SAMFileWriterFactory sfwf=createSAMFileWriterFactory();
		if(getOutputFile()==null)
			{
			return getOutputSamFormat().openSAMFileWriterToStdout(sfwf, header, presorted);
			}	
		else
			{
			return getOutputSamFormat().openSAMFileWriterToFile(sfwf, getOutputFile(),header, presorted);
			}
		}
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -o (file) output file. Default: stdout");
		out.print(  " -T (format) output format. optional. One of "); 
			for(SamFormat f:SamFormat.values()) out.print(f.name()+" ");
			out.println();
		super.printOptions(out);
		}
	@Override
	protected String getGetOptDefault() {
		return super.getGetOptDefault()+"o:T:";
		}
	@Override
	protected GetOptStatus handleOtherOptions(int c, GetOpt opt, String[] args) {
		switch(c)
			{
			case 'o':
				{
				File f= new File(opt.getOptArg());
				SamFormat sf = SamFormat.getSamFormatFromFile(f);
				if(sf==null)
					{
					error("The extension of "+f+" is not a valid extension.");
					return GetOptStatus.EXIT_FAILURE;
					}
				this.setOutputFile(f);
				this.setOutputSamFormat(sf);
				return GetOptStatus.OK;
				}
			case 'T':
				{
				try {
					SamFormat sf = SamFormat.valueOf(opt.getOptArg());
					this.setOutputSamFormat(sf);
					} 
				catch (Exception e) {
					error("No a valid format.");
					error(e);
					return GetOptStatus.EXIT_FAILURE;
					}
				return GetOptStatus.OK;
				}
			default: return super.handleOtherOptions(c, opt, args);
			}
		}
	}
