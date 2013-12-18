package com.github.lindenb.jvarkit.util.vcf;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;

import net.sf.samtools.util.BlockCompressedOutputStream;

import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;


/**
 * extends com.github.lindenb.jvarkit.util.AbstractCommandLineProgram
 * while AbstractVCFFilter extends com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram
 * @author lindenb
 *
 */
public abstract class AbstractVCFFilter2
	extends AbstractCommandLineProgram
	{
	protected AbstractVCFFilter2()
		{
		
		}
	
	@Override
	public String getProgramDescription() {
		return "Another VCF filter.";
		}
	
	protected abstract void doWork(
			VcfIterator in,
			VariantContextWriter out
			) throws IOException;
	
	protected  VcfIterator createVcfIterator(String IN) throws IOException
		{
		return VCFUtils.createVcfIterator(IN);
		}
	
	protected VariantContextWriter createVariantContextWriter(File OUT) throws IOException
		{
		return VCFUtils.createVariantContextWriter(OUT);
		}
	
	protected int doWork(String IN,File OUT)
		{
		VcfIterator r=null;
		VariantContextWriter w=null;
		try
			{
			r=this.createVcfIterator(IN);
			w=this.createVariantContextWriter(OUT);
			doWork(r,w);
			}
		catch (Exception e)
			{
			this.error(e);
			return -1;
			}
		finally
			{
			if(w!=null) w.close();
			if(r!=null) r.close();
			}	
		return 0;
		}
	
	protected int doWork(int optind,String args[])
		{

		try
			{
			if(optind==args.length)
				{
				info("reading from stdin.");
				return doWork((String)null, null);
				}
			else if(optind+1==args.length)
				{
				info("reading from "+args[optind]);
				return doWork(args[optind],null);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		}
	
	}
