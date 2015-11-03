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


import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.command.CommandFactory;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;


/**
 * extends com.github.lindenb.jvarkit.util.AbstractCommandLineProgram
 * while AbstractVCFFilter extends com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram
 * @author lindenb
 *
 */
@Deprecated
public abstract class AbstractVCFFilter4
	extends CommandFactory
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(AbstractVCFFilter4.class);


	protected static abstract class AbstractVCFCommand extends 
		Command
		{
		public abstract File getOutputFile();
		
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
		
		protected VCFHeader addMetaData(final VCFHeader header)
			{
			header.addMetaDataLine(new VCFHeaderLine(getName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header.addMetaDataLine(new VCFHeaderLine(getName()+"Version",String.valueOf(getVersion())));
			header.addMetaDataLine(new VCFHeaderLine(getName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			header.addMetaDataLine(new VCFHeaderLine(getName()+"HtsJdkHome",HtsjdkVersion.getHome()));
			return header;
			}

		
		/** increase the number of variants count produced at the end */
		protected void incrVariantCount()
			{
			setVariantCount(1+getVariantCount());
			}
			
		/** open VariantContextWriter */
		protected VariantContextWriter createVariantContextWriter()
			throws IOException
			{
			if(getOutputFile()==null)
				{
				return VCFUtils.createVariantContextWriterToOutputStream(this.stdout());
				}
			else
				{
				LOG.info("opening vcf writer to "+getOutputFile());
				return VCFUtils.createVariantContextWriter(getOutputFile());
				}

			}
		

		
		/** number of variants filterer */
		private int countFilteredVariants=0;
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
				LOG.info("Number of Variants:"+getVariantCount());
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
		LinkedHashSet<String> unrollInputFiles(List<String> args)
			{
			return IOUtils.unrollFiles(args); 
			}
		
		protected Throwable validateOptions()
			{
			return null;
			}
		
		@Override
		public Collection<Throwable> call() throws Exception {
			List<String> args = this.getInputFiles();
			VcfIterator vcfIn=null;
			try
				{
				Throwable errValidation = validateOptions();
				if(errValidation!=null)
					{
					return wrapException(errValidation);
					}
				LinkedHashSet<String> unrolledVcfFiles = unrollInputFiles(args); 
				if(unrolledVcfFiles.isEmpty() && args.isEmpty())
					{
					LOG.info("Reading from stdin");
					vcfIn = new VcfIterator(this.stdin());
					this.filterVcfIterator("<STDIN>",vcfIn);
					}
				else if(unrolledVcfFiles.size()==1)
					{
					String filename= unrolledVcfFiles.iterator().next();
					LOG.info("Reading from "+filename);
					vcfIn = VCFUtils.createVcfIterator(filename);
					this.filterVcfIterator(filename,vcfIn);
					}
				else
					{
					return wrapException(getMessageBundle("illegal.number.of.arguments"));
					}
				
				return Collections.emptyList();
				}
			catch(Exception err)
				{
				LOG.error(err);
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(vcfIn);
				}
			}

		
		}
	
	protected AbstractVCFFilter4()
		{
		}
	
	
	/* give a chance to unroll files ending with '.list' */
	protected LinkedHashSet<String> unrollInputFiles(List<String> args)
		{
		return IOUtils.unrollFiles(args); 
		}
	
	
	
	
	}
