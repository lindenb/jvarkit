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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.util.Collection;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VCFPolyX extends AbstractVCFPolyX
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VCFPolyX.class);
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	
	public VCFPolyX()
		{
		}
	
	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	@Override
	protected Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator r,
			final VariantContextWriter w
			) throws IOException
		{
		GenomicSequence genomicSequence=null;

		final VCFHeader header=r.getHeader();
		final VCFHeader h2=new VCFHeader(header);
		addMetaData(h2);
		
		final VCFInfoHeaderLine infoHeaderLine = new VCFInfoHeaderLine(
				super.polyXtag.trim(),
				1,
				VCFHeaderLineType.Integer,
				"Number of repeated bases around REF")
				;
		
		final VCFFilterHeaderLine filterHeaderLine = new VCFFilterHeaderLine(
				infoHeaderLine.getID()+"_ge_"+this.filterTrehsold,
				"Number of repeated bases around REF is greater or equal to "+this.filterTrehsold)
				;
		
		
		h2.addMetaDataLine(infoHeaderLine);
		
		if( this.filterTrehsold>-1) {
			h2.addMetaDataLine(filterHeaderLine);
		}

		w.writeHeader(h2);
		final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
		while(r.hasNext())
			{
			final VariantContext ctx=progress.watch(r.next());
			final VariantContextBuilder b=new VariantContextBuilder(ctx);
			if(genomicSequence==null || !ctx.getContig().equals(genomicSequence.getChrom()))
				{
				LOG.info("loading chromosome "+ctx.getContig());
				genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, ctx.getContig());
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
			pos0=ctx.getEnd()-1;
			c0=Character.toUpperCase(genomicSequence.charAt(pos0));
			pos0++;
			while(pos0< genomicSequence.getSAMSequenceRecord().getSequenceLength()
				&& c0==Character.toUpperCase(genomicSequence.charAt(pos0)))
				{
				++count;
				++pos0;
				}
			b.attribute(infoHeaderLine.getID(),count);
			
			/* filter */
			if(super.filterTrehsold>-1 && count>=super.filterTrehsold)
				{
				b.filter(filterHeaderLine.getID());
				}
			
			w.add(b.make());
			if(w.checkError()) break;
			}		
		progress.finish();
		genomicSequence=null;
		return RETURN_OK;
		}


	
	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.polyXtag.trim().isEmpty()) {
			return wrapException("Empty tag in option -"+OPTION_POLYXTAG);
		}
		try {
			if(this.REF==null)
				{
				return wrapException("Undefined Reference");
				}
			LOG.info("opening reference "+REF);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);
			}
		catch (final Exception e) {
			return wrapException(e);
			}
		return super.initializeKnime();
		}
	@Override
	public void disposeKnime() {
		CloserUtil.close(this.indexedFastaSequenceFile);
		this.indexedFastaSequenceFile=null;
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VCFPolyX().instanceMainWithExit(args);
		}

}
