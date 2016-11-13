/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 * VcfMoveFiltersToInfo
 *
 */
public class VcfMoveFiltersToInfo
	extends AbstractVcfMoveFiltersToInfo
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfMoveFiltersToInfo.class);
	
	
	public VcfMoveFiltersToInfo()
		{
		}
	 
	
	/* public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator in,
			final VariantContextWriter out
			) throws IOException {
		if(StringUtil.isBlank(super.infoName)) {
			return wrapException("undefined option -"+OPTION_INFONAME);
			}
		final VCFHeader header = in.getHeader();
		final Set<String> limitToThoseFilters = new HashSet<>();
		if(!StringUtil.isBlank(super.onlyThoseFiltersTagStr)) {
			for(final String f:super.onlyThoseFiltersTagStr.split("[, ]"))
				{
				if(StringUtil.isBlank(f)) continue;
				limitToThoseFilters.add(f.trim());
				}
		}
		
		try {
			final VCFInfoHeaderLine infoHeaderLine = new VCFInfoHeaderLine(
					super.infoName.trim(),
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Variant was previously FILTERed with the given values."
					);
			if(header.getInfoHeaderLine(infoHeaderLine.getID())!=null)
				{
				return wrapException("INFO["+infoHeaderLine.getID()+"] already exists in input VCF.");
				}
			
			
			final VCFHeader h2= new VCFHeader(header);
			
			h2.addMetaDataLine(infoHeaderLine);
			
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			out.writeHeader(h2);
			while(in.hasNext() &&  !out.checkError())
				{
				final VariantContext ctx = progess.watch(in.next());
				if(ctx.isNotFiltered())
					{
					out.add(ctx);
					}
				else
					{
					final Set<String> INFOfilters = new HashSet<>();
					final Set<String> FILTERfilters = new HashSet<>();
					for(final String filter : ctx.getFilters()) {
						if( filter.equals(VCFConstants.UNFILTERED) ||
							filter.equals(VCFConstants.PASSES_FILTERS_v3) ||
							filter.equals(VCFConstants.PASSES_FILTERS_v4)
							)
							{
							continue;
							}
						if(!limitToThoseFilters.isEmpty() && !limitToThoseFilters.contains(filter)) {
							FILTERfilters.add(filter);
							}
						else
							{
							INFOfilters.add(filter);
							}
						}
					
					
					final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
					if(FILTERfilters.isEmpty()) {
						vcb.unfiltered();
						}
					else
						{
						vcb.filters(FILTERfilters);
						}
					if(!INFOfilters.isEmpty()) {
						vcb.attribute(infoHeaderLine.getID(), new ArrayList<>(INFOfilters));
						}
					out.add(vcb.make());
					}
				}
			progess.finish();
			LOG.info("done.");
			return RETURN_OK;
			} catch(Exception err) {
				return wrapException(err);
			} finally {
				CloserUtil.close(in);
				CloserUtil.close(out);
			}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfMoveFiltersToInfo().instanceMainWithExit(args);
		}
	}
