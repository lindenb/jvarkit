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
package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * 
 * VCFBedSetFilter
 *
 */
public class VCFBedSetFilter extends AbstractVCFFilter3
	{
	private IndexedBedReader bedReader =null;
	private String filterName = null; 
	private String filterDescription = null; 
	private boolean invertSelection = false;
	private File TABIX;
	
	public VCFBedSetFilter()
		{
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VCFBedSetFilter";
		}
	
	@Override
	public String getProgramDescription() {
		return "Set FILTER for VCF if it intersects with BED.";
		}
	
	public void setFilterDescription(String filterDescription) {
		this.filterDescription = filterDescription;
	}
	
	public void setFilterName(String filterName) {
		this.filterName = filterName;
		}

	
	public void setInvertSelection(boolean invertSelection) {
		this.invertSelection = invertSelection;
		}
	@Override
	protected void doWork(String inputSource,VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		VCFHeader header=r.getHeader();

		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

		
		h2.addMetaDataLine(new VCFFilterHeaderLine(this.filterName, this.filterDescription));

		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
		w.writeHeader(h2);
		while(r.hasNext())
			{
			VariantContext ctx= progress.watch(r.next());
			boolean set_filter=false;
			
			CloseableIterator<IndexedBedReader.BedLine> iter = this.bedReader.iterator(
					ctx.getContig(),
					ctx.getStart()-1,
					ctx.getEnd()+1
					);
			while(iter.hasNext())
				{
				IndexedBedReader.BedLine bed=iter.next();
				if(!ctx.getContig().equals(bed.getContig())) continue;
				if(ctx.getStart()-1 >= bed.getEnd() ) continue;
				if(ctx.getEnd()-1 < bed.getStart() ) continue;
				set_filter=true;
				break;
				}
			CloserUtil.close(iter);
			
			if(this.invertSelection) set_filter=!set_filter;
			
			if(checkOutputError()) break;
			
			if(!set_filter)
				{
				incrVariantCount();
				w.add(ctx);
				continue;
				}
			
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.filter(this.filterName);
			w.add(vcb.make());
			incrVariantCount();
			
			}
		progress.finish();
		}
	
	
	
	public void setBedFile(File tabix) {
		this.TABIX = tabix;
		}
	

	
	@Override
	public int initializeKnime() {
		try
			{
			if(this.TABIX==null)
				{
				error("Undefined tabix file");
				return -1;
				}
			if(this.filterName==null || this.filterName.trim().isEmpty())
				{
				error("Undefined filter name.");
				return -1;
				}
			if(this.filterDescription==null)
				{
				this.filterDescription="Undefined";
				}
			this.info("opening Bed "+this.TABIX);
			this.bedReader= new IndexedBedReader(this.TABIX);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime()
		{
		CloserUtil.close(this.bedReader);
		this.bedReader = null;
		super.disposeKnime();
		}
	
	public void printOptions(PrintStream out) {
		out.println(" -f (filter name).FILTER Name. Required.");
		out.println(" -d (filter description). Optional.");
		out.println(" -B (path). BED file indexed with tribble/tabix");
		out.println(" -V  inverse selection.");
		out.println(" -o (file) output file. default: stdout");
		super.printOptions(out);
		}

	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:d:B:Vf:"))!=-1)
			{
			switch(c)
				{
				case 'f': this.setFilterName(opt.getOptArg());break;
				case 'd': this.setFilterDescription(opt.getOptArg());break;
				case 'B': this.setBedFile(new File(opt.getOptArg()));break;
				case 'o': this.setOutputFile(new File(opt.getOptArg()));break;
				case 'V': this.setInvertSelection(true);break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
		}

	
	
	public static void main(String[] args) throws Exception
		{
		new VCFBedSetFilter().instanceMainWithExit(args);
		}
}
