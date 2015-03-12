/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.vcf;


import java.io.IOException;
import java.io.PrintStream;
import java.util.LinkedHashSet;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;


/**
 * extends com.github.lindenb.jvarkit.util.AbstractCommandLineProgram
 * while AbstractVCFFilter extends com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram
 * @author lindenb
 *
 */
public abstract class AbstractVCFFilter3
	extends AbstractKnimeApplication
	{
	/** number of variants filterer */
	private int countFilteredVariants=0;

	
	protected AbstractVCFFilter3()
		{
		}
	
	/** return the number of variants in the output vcf */
	public int getVariantCount()
		{
		return this.countFilteredVariants;
		}
	
	/** set the number of variants in the output vcf */
	public void setVariantCount(int countFilteredVariants)
		{
		this.countFilteredVariants = countFilteredVariants;
		}
	
	/** increase the number of variants count produced at the end */
	protected void incrVariantCount()
		{
		setVariantCount(1+getVariantCount());
		}
	
	@Override
	public String getProgramDescription() {
		return "Another VCF filter.";
		}
	
	/** open VariantContextWriter */
	protected VariantContextWriter createVariantContextWriter()
		throws IOException
		{
		if(getOutputFile()==null)
			{
			return VCFUtils.createVariantContextWriterToStdout();
			}
		else
			{
			info("opening vcf writer to "+getOutputFile());
			return VCFUtils.createVariantContextWriter(getOutputFile());
			}

		}
	
	/** called filterVcfIterator(in) */
	protected void doWork( String inputSource,VcfIterator vcfIn,VariantContextWriter out )
		throws IOException
		{
		throw new IOException(
				"This is the default method. It hasn't been implemented.");
		}
	
	/** called by executeKnime . Default behaviour: open VariantCountextWriter
	 *  and call filterVcfIterator(in,out) */
	protected void filterVcfIterator( String inputSource, VcfIterator vcfIn )
		throws IOException
		{
		VariantContextWriter out=null;
		try
			{
			out = createVariantContextWriter();
			setVariantCount(0);
			doWork(inputSource,vcfIn,out);
			info("Number of Variants:"+getVariantCount());
			}
		catch(IOException err)
			{
			throw err;
			}
		catch(Exception err)
			{
			throw err;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	
	/* give a chance to unroll files ending with '.list' */
	protected LinkedHashSet<String> unrollInputFiles(List<String> args)
		{
		return IOUtils.unrollFiles(args); 
		}
	
	
	@Override
	public int executeKnime(List<String> args)
		{
		VcfIterator vcfIn=null;
		try
			{
			LinkedHashSet<String> unrolledVcfFiles = unrollInputFiles(args); 
			if(unrolledVcfFiles.isEmpty() && args.isEmpty())
				{
				vcfIn = VCFUtils.createVcfIteratorStdin();
				this.filterVcfIterator(INPUT_SOURCE_STDIN,vcfIn);
				}
			else if(unrolledVcfFiles.size()==1)
				{
				String filename= unrolledVcfFiles.iterator().next();
				vcfIn = VCFUtils.createVcfIterator(filename);
				this.filterVcfIterator(filename,vcfIn);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfIn);
			}
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -o (fileout). Optional. Default: stdout");
		super.printOptions(out);
		}

	}
