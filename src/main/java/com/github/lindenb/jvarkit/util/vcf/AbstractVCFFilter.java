package com.github.lindenb.jvarkit.util.vcf;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.util.Log;
import net.sf.samtools.util.BlockCompressedOutputStream;

import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

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
			return new VcfIterator(System.in);
			}
		else
			{
			LOG.info("reading from "+IN);
			return new VcfIterator(IOUtils.openURIForReading(IN));
			}
		}
	
	protected VariantContextWriter createVariantContextWriter() throws IOException
		{
		if(OUT==null)
			{
			LOG.info("writing to stdout");
			return VariantContextWriterFactory.create(System.out,null,EnumSet.noneOf(Options.class));
			}
		else if(OUT.getName().endsWith(".gz"))
			{
			LOG.info("writing to "+OUT+" as bgz file.");
			BlockCompressedOutputStream bcos=new BlockCompressedOutputStream(OUT);
			return VariantContextWriterFactory.create(bcos,null,EnumSet.noneOf(Options.class));
			}
		else
			{
			LOG.info("writing to "+OUT);
			return  VariantContextWriterFactory.create(OUT,null,EnumSet.noneOf(Options.class));
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
			LOG.error(e);
			testRemoteGit();
			return -1;
			}
		finally
			{
			if(w!=null) w.close();
			if(r!=null) r.close();
			}	
		return 0;
		}
	}
