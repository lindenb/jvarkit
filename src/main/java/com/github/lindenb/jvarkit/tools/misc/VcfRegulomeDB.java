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

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 * VcfRegulomeDB
 * @author lindenb
 *
 */
@Program(name="ccfregulomedb",description="Annotate a VCF with the Regulome data (http://regulome.stanford.edu/")
public class VcfRegulomeDB extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfRegulomeDB.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;
	@Parameter(names="-b",description=" bed indexed with tabix. Format: chrom(tab)start(tab)end(tab)rank",required=true)
	private String bedFile=null;
	@Parameter(names="-T",description="tag in vcf INFO.")
	private String infoTag="REGULOMEDB";
	@Parameter(names="-x",description="(int) base pairs. look.for data around the variation +/- 'x' ")
	private int extend=5;
	@Parameter(names="-r",description="if defined, only accept the rank matching the regular expression ")
	private String acceptRegexStr=null;	
	private Pattern acceptRegex=null;
	
	private RegDataTabixFileReader regDataTabixFileReader=null;
	
	
	private static class RegData
		{
		@SuppressWarnings("unused")
		String chrom;
		int chromSart;
		@SuppressWarnings("unused")
		int chromEnd;
		String rank;
		RegData(String tokens[])
			{
			this.chrom=tokens[0];
			this.chromSart=Integer.parseInt(tokens[1]);
			this.chromEnd=Integer.parseInt(tokens[2]);
			this.rank=tokens[3];
			}
		}
	
	private static class RegDataTabixFileReader
		extends AbstractTabixObjectReader<RegData>
		{
		RegDataTabixFileReader(String uri) throws IOException
			{
			super(uri);
			}
		@Override
		protected Iterator<RegData> iterator(Iterator<String> delegate) {
			return new MyIterator(delegate);
			}
		private class MyIterator
	    	extends AbstractMyIterator
	    	{
	    	private Pattern tab=Pattern.compile("[\t]");
	    	MyIterator(Iterator<String> delegate)
	    		{
	    		super(delegate);
	    		}
	    	
	    	@Override
	    	public RegData next() {
	    		return new RegData(this.tab.split(delegate.next(),5));
	    		}
	    	}

		}
	
	private VcfRegulomeDB()
		{
		
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in,
			VariantContextWriter out)
		{
		
		VCFHeader header=in.getHeader();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		header.addMetaDataLine(new VCFInfoHeaderLine(
				this.infoTag,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				"Format: Position|Distance|Rank"
				));
		
		
		out.writeHeader(header);
		
		while(in.hasNext())
			{
			List<String>  regDataList=new ArrayList<String>();
			VariantContext ctx=in.next();
			
			progress.watch(ctx.getContig(),ctx.getStart());
			
			int start=Math.max(0,ctx.getStart()-this.extend);
			int end=ctx.getEnd()+this.extend;
			
			for(Iterator<RegData> iter=this.regDataTabixFileReader.iterator(ctx.getContig(), start, end);
					iter.hasNext();
					)
				{
				RegData curr=iter.next();
				if(this.acceptRegex!=null && 
				   !this.acceptRegex.matcher(curr.rank).matches()
				   )
					{
					continue;
					}
				String str=
						String.valueOf(curr.chromSart)+"|"+
						String.valueOf(Math.abs(curr.chromSart-(ctx.getStart()-1)))+"|"+
						curr.rank;
				regDataList.add(str);
				}
			if(regDataList.isEmpty())
				{
				out.add(ctx);
				continue;
				}
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.attribute(this.infoTag, regDataList.toArray());
			out.add(vcb.make());
			}
		progress.finish();
		return 0;
		}
	
		

	@Override
	public int doWork(List<String> args)
		{
		
		
		if(bedFile==null)
			{
			LOG.error("Bed file indexed with tabix is missing");
			return -1;
			}
		
		try
			{
			if(this.acceptRegexStr!=null)
				{
				this.acceptRegex=Pattern.compile(this.acceptRegexStr);
				}
			
			LOG.info("Opening "+bedFile);
			this.regDataTabixFileReader=new RegDataTabixFileReader(bedFile);
			return doVcfToVcf(args, outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.regDataTabixFileReader);
			}
		}
	public static void main(String[] args)
		{
		new VcfRegulomeDB().instanceMainWithExit(args);
		}
	}
