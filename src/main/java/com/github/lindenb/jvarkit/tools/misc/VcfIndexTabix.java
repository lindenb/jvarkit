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
import java.util.Collection;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.util.TabixUtils;
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
	extends AbstractVcfIndexTabix
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfIndexTabix.class);

	private class SortingVCFWriter implements VariantContextWriter
		{
		VariantContextWriter delegate;
		SortingCollection<VariantContext> sorter=null;
		SortingVCFWriter(VariantContextWriter delegate)
			{
			this.delegate=delegate;
			}
		
		@Override
		public boolean checkError()
			{
			return  delegate.checkError();
			}

		
		@Override
		public void writeHeader(VCFHeader header) {
			this.delegate.writeHeader(header);			
			this.sorter =
	                SortingCollection.newInstance(
	                        VariantContext.class,
	                        new VCFRecordCodec(header),
	                        header.getVCFRecordComparator(),
	                        VcfIndexTabix.this.getMaxRecordsInRam(),
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
	
	/* public for knime */
	public Collection<Throwable> doVcfToVcf(String inputName, final VcfIterator vcfIn,final File outFile) throws IOException {
		if(outFile==null)
			{
			throw new IOException(getName()+" output file undefined.");
			}
		if(!outFile.getName().endsWith(".vcf.gz"))
			{
			return wrapException("output file should en with .vcf.gz but got "+outFile);
			}

		SortingVCFWriter sortingVCW=null;
		VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
		VariantContextWriter w=null;
		try {
			SAMSequenceDictionary dict = vcfIn.getHeader().getSequenceDictionary();
			if(dict!=null) vcwb.setReferenceDictionary(dict);
			
			vcwb.setOutputFile(outFile);
			vcwb.setOutputFileType(OutputType.BLOCK_COMPRESSED_VCF);
			w = vcwb.build();
			
			if(this.sort)
				{
				LOG.info("Creating a sorting writer");
				sortingVCW= new SortingVCFWriter(w);
				w=sortingVCW;
				}
			
			w.writeHeader(vcfIn.getHeader());
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(vcfIn.getHeader());
			while(vcfIn.hasNext())
				{
				w.add(progress.watch(vcfIn.next()));
				}
			progress.finish();
			
			w.close();
			w=null;
			return RETURN_OK;
			} 
		catch (Exception e)
			{
			if(getOutputFile().exists() && getOutputFile().isFile())
				{
				LOG.warn("Deleting "+outFile);
				getOutputFile().delete();
				File tbi = new File(getOutputFile().getPath()+TabixUtils.STANDARD_INDEX_EXTENSION);
				if(tbi.exists() && tbi.isFile()) tbi.delete();
				}
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(w);
			}
		}
	
	@Override
	public Collection<Throwable> initializeKnime() {
		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime() {
		super.disposeKnime();
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(getOutputFile()==null)
			{
			return wrapException("output file is required");
			}
		VcfIterator iter = null;
		try {
			iter =  super.openVcfIterator(inputName);
			return doVcfToVcf(inputName==null?"STDIN":inputName, iter, getOutputFile());
		} catch (Exception e) {
			return wrapException(e);
		} finally {
			CloserUtil.close(iter);
		}
		
		
		}
		
		
		
	public static void main(String[] args)
		{
		new VcfIndexTabix().instanceMainWithExit(args);
		}
	}
