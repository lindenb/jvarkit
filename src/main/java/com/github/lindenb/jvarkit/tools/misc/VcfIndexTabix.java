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
* 2015 : moving to knime
* 2014 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFRecordCodec;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfIndexTabix
	extends AbstractVCFFilter3
	{
	private boolean sort=false;
	private int MAX_RECORDS_IN_RAM=10000;
	
	private class SortingVCFWriter implements VariantContextWriter
		{
		VariantContextWriter delegate;
		SortingCollection<VariantContext> sorter=null;
		SortingVCFWriter(VariantContextWriter delegate)
			{
			this.delegate=delegate;
			}
		@Override
		public void writeHeader(VCFHeader header) {
			this.delegate.writeHeader(header);			
			this.sorter =
	                SortingCollection.newInstance(
	                        VariantContext.class,
	                        new VCFRecordCodec(header),
	                        header.getVCFRecordComparator(),
	                        VcfIndexTabix.this.MAX_RECORDS_IN_RAM,
	                        VcfIndexTabix.this.getTmpDirectories().get(0)
	                        );
			
			}
		@Override
		public void add(VariantContext vc)
			{
			this.sorter.add(vc);
			}
		
		@Override
		public void close()
			{
			CloseableIterator<VariantContext> iter=null;
			try
				{
				this.sorter.doneAdding();
				this.sorter.setDestructiveIteration(true);
				iter= this.sorter.iterator();
				while(iter.hasNext())
					{
					this.delegate.add(iter.next());
					}
				this.sorter.cleanup();
				}
			catch(Exception err)
				{
				throw new RuntimeException(err);
				}
			finally
				{
				CloserUtil.close(iter);
				CloserUtil.close(this.delegate);
				}
			}
		}
	
	
	public VcfIndexTabix()
		{
		}
	
	
	@Override
	public String getProgramDescription() {
		return "Index and sort a VCF on the fly with Tabix.";
		}
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VcfIndexTabix";
		}
	
	public void setSort(boolean sort) {
		this.sort = sort;
		}
	
	public void setMaxRecordsInRam(int n) {
		this.MAX_RECORDS_IN_RAM = n;
		}
	
	
	@Override
	protected void doWork(String inputSource, VcfIterator vcfIn,
			VariantContextWriter out) throws IOException {
		throw new IllegalStateException("should not be invoked");
		}
	
	@Override
	protected void filterVcfIterator(String inputSource, VcfIterator vcfIn)
			throws IOException
		{
		File out=getOutputFile();
		if(out==null)
			{
			throw new IOException(getProgramName()+" output file undefined.");
			}
		if(!out.getName().endsWith(".vcf.gz"))
			{
			throw new IOException("output file should en with .vcf.gz but got "+out);
			}

		SortingVCFWriter sortingVCW=null;
		VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
		VariantContextWriter w=null;
		try {
			SAMSequenceDictionary dict = vcfIn.getHeader().getSequenceDictionary();
			if(dict!=null) vcwb.setReferenceDictionary(dict);
			
			vcwb.setOutputFile(out);
			vcwb.setOutputFileType(OutputType.BLOCK_COMPRESSED_VCF);
			w = vcwb.build();
			
			if(this.sort)
				{
				info("Creating a sorting writer");
				sortingVCW= new SortingVCFWriter(w);
				w=sortingVCW;
				}
			
			w.writeHeader(vcfIn.getHeader());
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(vcfIn.getHeader());
			while(vcfIn.hasNext())
				{
				w.add(progress.watch(vcfIn.next()));
				incrVariantCount();
				}
			progress.finish();
			
			w.close();
			w=null;
			} 
		catch (Exception e) {
			throw new IOException(e);
			}
		finally
			{
			CloserUtil.close(w);
			}
		}
	
	
	@Override
	public int initializeKnime() {
		return super.initializeKnime();
		}
	@Override
	public void disposeKnime() {
		super.disposeKnime();
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -o (file) filename.vcf.gz out. Required.");
		out.println(" -s sort VCF prior to saving.");
		out.println(" -n (max-record-in-ram). If defined, vcf will be sorted using this cache size. Default:"+MAX_RECORDS_IN_RAM);
		out.println(" -T (tmp).add tmp dir. optional)");
		super.printOptions(out);
		}
		
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:n:sT:"))!=-1)
			{
			switch(c)
				{
				case 'T': this.addTmpDirectory(new File(opt.getOptArg()));break;
				case 's': this.setSort(true);break;
				case 'n': this.setMaxRecordsInRam(Integer.parseInt(opt.getOptArg()));break;
				case 'o': this.setOutputFile(new File(opt.getOptArg()));break;
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
		}
		
		
		
	public static void main(String[] args)
		{
		new VcfIndexTabix().instanceMainWithExit(args);
		}
	}
