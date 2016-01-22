package com.github.lindenb.jvarkit.util.vcf;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.IOException;

import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.StandardOptionDefinitions;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;


public abstract class AbstractVCFFilter
	extends AbstractCommandLineProgram
	{
	private static final Log LOG=Log.getInstance(AbstractVCFFilter.class);
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF file/URL to process. Default stdin. ",optional=true)
	public String IN=null;
	@Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF file to generate. Default stdout. ",optional=true)
	public File OUT=null;
	
	protected abstract void doWork(
			VcfIterator in,
			VariantContextWriter out
			) throws IOException;
	
	protected  VcfIterator createVcfIterator() throws IOException
		{
		if(IN==null)
			{
			LOG.info("reading from stdin");
			return VCFUtils.createVcfIteratorStdin();
			}
		else
			{
			LOG.info("reading from "+IN);
			return VCFUtils.createVcfIterator(IN);
			}
		}
	
	protected VariantContextWriter createVariantContextWriter() throws IOException
		{
		if(OUT==null)
			{
			LOG.info("writing to stdout");
			return VCFUtils.createVariantContextWriterToStdout();
			}
		else
			{
			LOG.info("writing to "+OUT);
			return VCFUtils.createVariantContextWriter(OUT);
			}
		}
	
	@Override
	protected int doWork()
		{
		VcfIterator r=null;
		VariantContextWriter w=null;
		try
			{
			r=this.createVcfIterator();
			w=this.createVariantContextWriter();
			doWork(r,w);
			}
		catch (Exception e)
			{
			e.printStackTrace();
			LOG.error(e);
			LOG.error("Commandline was : "+getCommandLine());
			testRemoteGit();
			return -1;
			}
		finally
			{
			if(w!=null) w.close();
			CloserUtil.close(r);
			}	
		return 0;
		}
	}
