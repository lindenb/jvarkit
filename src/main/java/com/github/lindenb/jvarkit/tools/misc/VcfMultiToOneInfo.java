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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfMultiToOneInfo
	extends AbstractVcfMultiToOneInfo
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfMultiToOneInfo.class);

	 public VcfMultiToOneInfo()
		{
		}
	 
	@Override
	/* public to be called by knime */
	public Collection<Throwable> doVcfToVcf(final String inputName,
			final VcfIterator in, final VariantContextWriter out) throws IOException {
		
		final VCFHeader srcHeader=in.getHeader();
		final VCFInfoHeaderLine srcInfo = srcHeader.getInfoHeaderLine(this.infoTag);
		if( srcInfo == null )
			{
			return wrapException("Cannot find INFO FIELD '"+ this.infoTag+"'");
			}
		switch( srcInfo.getCountType() )
			{
			case INTEGER:break;
			case UNBOUNDED:break;
			default: return wrapException("CountType is not supported '"+ srcInfo.getCountType() +"'");
			}
		switch( srcInfo.getType())
			{
			case Flag: return wrapException("Type is not supported '"+ srcInfo.getType() +"'");
			default:break;
			}
		
		final VCFHeader destHeader= new VCFHeader(srcHeader);
		super.addMetaData(destHeader);
		
		final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(srcHeader);
		out.writeHeader(destHeader);
		while(in.hasNext())
			{
			final VariantContext ctx=progess.watch(in.next());
			final List<Object> L=ctx.getAttributeAsList(srcInfo.getID());
			if( L.isEmpty() || L.size()==1)
				{
				out.add(ctx);
				continue;
				}
			for(final Object o: L)
				{
				final VariantContextBuilder vcb = super.getVariantContextBuilderFactory().newVariantContextBuilder(ctx);
				vcb.attribute(srcInfo.getID(), o);
				out.add(vcb.make());
				}
			}
		progess.finish();
		LOG.info("done");
		return RETURN_OK;
		}
	
	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.infoTag==null || super.infoTag.isEmpty())
			{
			return wrapException("No info tag defined");
			}
		return super.initializeKnime();
		}

	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		
		return doVcfToVcf(inputName);
		}
	
	public static void main(String[] args)
		{
		new VcfMultiToOneInfo().instanceMainWithExit(args);
		}
	
	}
