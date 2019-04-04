/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.PostponedVariantContextWriter;
import htsjdk.variant.vcf.VCFIterator;

@Program(
		name="vcfmulti2oneinfo",
		description="'one INFO with N values' to 'N variants with one INFO'",
		keywords={"vcf"}
		)
public class VcfMultiToOneInfo
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfMultiToOneInfo.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-i","--info"},description="The INFO tag",required=true)
	private String infoTag = null;
	@ParametersDelegate
	private PostponedVariantContextWriter.WritingVcfConfig writingVcfArgs = new PostponedVariantContextWriter.WritingVcfConfig();

	 public VcfMultiToOneInfo()
		{
		}
	 
	 @Override
	protected int doVcfToVcf(final String inputName,final VCFIterator in,final VariantContextWriter out) {
		final VCFHeader srcHeader=in.getHeader();
		final VCFInfoHeaderLine srcInfo = srcHeader.getInfoHeaderLine(this.infoTag);
		if( srcInfo == null )
			{
			LOG.error("Cannot find INFO FIELD '"+ this.infoTag+"'");
			return -1;
			}
		switch( srcInfo.getCountType() )
			{
			case INTEGER:break;
			case UNBOUNDED:break;
			default: {
				LOG.error("CountType is not supported '"+ srcInfo.getCountType() +"'");
				return -1;
				}
			}
		switch( srcInfo.getType())
			{
			case Flag: 
				{
				LOG.error("Type is not supported '"+ srcInfo.getType() +"'");
				return -1;
				}
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
				final VariantContextBuilder vcb =new VariantContextBuilder(ctx);
				vcb.attribute(srcInfo.getID(), o);
				out.add(vcb.make());
				}
			}
		progess.finish();
		LOG.info("done");
		return RETURN_OK;
		}
	 
	@Override
	protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
		return new PostponedVariantContextWriter(this.writingVcfArgs,stdout(),this.outputFile);
		}

	@Override
	public int doWork(final List<String> args) {
		if(this.infoTag==null || this.infoTag.isEmpty())
			{
			LOG.error("No info tag defined");
			return -1;
			}
		return doVcfToVcf(args,outputFile);
		}
	
	public static void main(String[] args)
		{
		new VcfMultiToOneInfo().instanceMainWithExit(args);
		}
	
	}
