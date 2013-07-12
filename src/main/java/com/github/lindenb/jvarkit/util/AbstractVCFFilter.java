package com.github.lindenb.jvarkit.util;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.util.Log;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import com.github.lindenb.jvarkit.util.picard.IOUtils;


public abstract class AbstractVCFFilter
	extends CommandLineProgram
	{
	private final Log LOG=Log.getInstance(AbstractVCFFilter.class);
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF file/URL to process. Default stdin. ",optional=true)
	public String IN=null;
	@Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="VCF file to generate. Default stdout. ",optional=true)
	public File OUT=null;
	
	protected abstract void doWork(
			LineReader in,
			VariantContextWriter out
			) throws IOException;
	
	@Override
	protected int doWork()
		{
		AsciiLineReader r=null;
		VariantContextWriter w=null;
		try
			{
			if(IN==null)
				{
				LOG.info("reading from stdin");
				r=new AsciiLineReader(System.in);
				}
			else
				{
				LOG.info("reading from "+IN);
				r=new AsciiLineReader(IOUtils.openURIForReading(IN));
				}
			if(OUT==null)
				{
				LOG.info("writing to stdout");
				w=VariantContextWriterFactory.create(System.out,null,EnumSet.noneOf(Options.class));
				}
			else
				{
				w=VariantContextWriterFactory.create(OUT,null,EnumSet.noneOf(Options.class));
				}
			doWork(r,w);
			}
		catch (Exception e)
			{
			LOG.error(e);
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
