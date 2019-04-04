package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.ContigPos;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

@Program(name="bluntgvcf",description="blunt a gvcf")
public class BluntGVcf extends Launcher
	{
	private static final Logger LOG = Logger.build(BluntGVcf.class).make();
	@Parameter(names = { "-o", "--out" }, description = OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names = { "-E", "--stop" }, description = "Blunt end this position",converter=ContigPos.Converter.class,required=true)
	private ContigPos pos = null;

	
	@Override
	protected int doVcfToVcf(
			final String inputName, 
			final VCFIterator in,
			final  VariantContextWriter out
			)
		{
		try
			{
			final LinkedList<String>  buffer= new LinkedList<>();
			final VCFHeader header= in.getHeader();
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if( dict == null ) {
				throw new IOException("Dict missing in input vcf.");
				}
			final int end_tid = dict.getSequenceIndex(this.pos.getContig());
			
			if( end_tid == -1 ) {
				throw new IOException("No such contig in dictionary "+this.pos);
			}

			out.writeHeader(header);
			while(in.hasNext())
				{
				final VariantContext ctx = in.next();
				int tid = dict.getSequenceIndex(ctx.getContig());
				if( tid == -1 ) {
					throw new IOException("No such contig in dictionary "+ctx);
					}
				if( tid > end_tid) break;
				if( tid < end_tid) {
					out.add(ctx);
					continue;
					}
				
				}
			
			return 0;
			}
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
		
			}
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		try
			{
			return doVcfToVcf(args,outputFile);
			}
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		}
	public static void main(String[] args)
		{
		new BluntGVcf().instanceMainWithExit(args);
		}
	}
