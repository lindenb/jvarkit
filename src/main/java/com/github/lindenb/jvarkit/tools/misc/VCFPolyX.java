package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.StandardOptionDefinitions;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VCFPolyX extends AbstractVCFFilter
	{
	private static Log LOG=Log.getInstance(VCFPolyX.class);

    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + " Number of repeated REF bases around POS. ";

    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME,doc="Reference",optional=false)
    public File REF=null;
    
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;

	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		GenomicSequence genomicSequence=null;
		LOG.info("opening reference "+REF);
		this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);

		final String TAG="POLYX";
		VCFHeader header=r.getHeader();
		
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG,1,
				VCFHeaderLineType.Integer,
				"number of repeated bases around REF")
				);
		
		w.writeHeader(h2);
	
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			if(genomicSequence==null || !ctx.getChr().equals(genomicSequence.getChrom()))
				{
				LOG.info("loading chromosome "+ctx.getChr());
				genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, ctx.getChr());
				}
			int pos0=ctx.getStart()-1;
			int count=1;
			char c0=Character.toUpperCase(genomicSequence.charAt(pos0));
			//go left
			pos0--;
			while(pos0>=0 && c0==Character.toUpperCase(genomicSequence.charAt(pos0)))
				{
				++count;
				pos0--;
				}
			//go right
			pos0=ctx.getStart()-1;
			pos0++;
			while(pos0< genomicSequence.getSAMSequenceRecord().getSequenceLength()
				&& c0==Character.toUpperCase(genomicSequence.charAt(pos0)))
				{
				++count;
				++pos0;
				}
			b.attribute(TAG,count);
			w.add(b.make());
			}		
		CloserUtil.close(indexedFastaSequenceFile);
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VCFPolyX().instanceMainWithExit(args);
		}

}
