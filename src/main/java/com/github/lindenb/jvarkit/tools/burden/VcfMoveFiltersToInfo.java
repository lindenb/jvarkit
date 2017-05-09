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

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
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
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
/**

BEGIN_DOC


For Matilde: move the information in FILTER to the INFO column to keep a trace of the FILTERs.


END_DOC
*/

@Program(name="vcfmovefilterstoinfo",
		description="Move any FILTER to the INFO column. reset FILTER to PASS",
		keywords={"vcf","burden","format","info"}
		)
public class VcfMoveFiltersToInfo
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfMoveFiltersToInfo.class).make();


	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;


	@Parameter(names={"-f","--filter"},description="INFO name. This tag will be used to store the previous filters")
	private String infoName = "PREVIOUSLY_FILTERED_AS";

	@Parameter(names={"-t","--limitto"},description="If not empty, limit to those FILTERS. Multiple separated by comma/space.")
	private String onlyThoseFiltersTagStr = null;
	
	
	public VcfMoveFiltersToInfo()
		{
		}
	 
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		if(StringUtil.isBlank(this.infoName)) {
			LOG.error("undefined option for infoName");
			return -1;
			}
		final VCFHeader header = in.getHeader();
		final Set<String> limitToThoseFilters = new HashSet<>();
		if(!StringUtil.isBlank(this.onlyThoseFiltersTagStr)) {
			for(final String f:this.onlyThoseFiltersTagStr.split("[, ]"))
				{
				if(StringUtil.isBlank(f)) continue;
				limitToThoseFilters.add(f.trim());
				}
		}
		
		try {
			final VCFInfoHeaderLine infoHeaderLine = new VCFInfoHeaderLine(
					this.infoName.trim(),
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Variant was previously FILTERed with the given values."
					);
			if(header.getInfoHeaderLine(infoHeaderLine.getID())!=null)
				{
				LOG.error("INFO["+infoHeaderLine.getID()+"] already exists in input VCF.");
				return -1;
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
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(in);
				CloserUtil.close(out);
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		return doVcfToVcf(args,outputFile);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfMoveFiltersToInfo().instanceMainWithExit(args);
		}
	}
