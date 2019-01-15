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


*/
package com.github.lindenb.jvarkit.tools.springbatch;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.batch.item.ExecutionContext;
import org.springframework.batch.item.ItemReader;
import org.springframework.batch.item.ItemStream;
import org.springframework.batch.item.ItemStreamException;
import org.springframework.batch.item.NonTransientResourceException;
import org.springframework.batch.item.ParseException;
import org.springframework.batch.item.UnexpectedInputException;
import org.springframework.batch.item.file.ResourceAwareItemReaderItemStream;
import org.springframework.core.io.Resource;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.FilteringVariantContextIterator;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextBatchReader implements 
	ResourceAwareItemReaderItemStream<List<VariantContext>>,
	ItemReader<List<VariantContext>>,
	ItemStream
	{
	private static final Log LOG = LogFactory.getLog(VariantContextBatchReader.class);
    private static final String CURRENT_INDEX = "current.vcf.index";

	private Resource rsrc = null;
	private CloseableIterator<VariantContext> iter = null;
	private VCFFileReader vcfFileReader = null;
	private Interval interval = null;
	private VariantContextFilter filter = V->true;
    private long currentIndex=0L;
	
	public void setInterval(final Interval interval) {
		this.interval = interval;
		}
	
	public void setFilter(final VariantContextFilter filter) {
		this.filter = filter;
		}
	
	@Override
	public void setResource(final Resource rsrc) {
		this.rsrc = rsrc;
		}
	
	@Override
	public void open(final ExecutionContext executionContext) throws ItemStreamException {
		if(this.rsrc==null) throw new ItemStreamException("resource is not defined");
		final File vcfFile;
		try {
			if(LOG.isInfoEnabled()) LOG.info("Opening "+this.rsrc);
			vcfFile = this.rsrc.getFile();
			IOUtil.assertFileIsReadable(vcfFile);
			this.vcfFileReader = new VCFFileReader(
				vcfFile,
				this.interval!=null
				);
			final VCFHeader header = this.vcfFileReader.getFileHeader();
			
			executionContext.put(SpringBatchUtils.VCF_HEADER_KEY, header);
			
			if(this.interval == null) {
				this.iter = this.vcfFileReader.iterator();
				}
			else
				{
				this.iter = this.vcfFileReader.query(
					this.interval.getContig(),
					this.interval.getStart(),
					this.interval.getEnd()
					);
				}
			if(this.filter!=null) {
				this.iter = new FilteringVariantContextIterator(this.iter, this.filter);
				}
			 if(!executionContext.containsKey(CURRENT_INDEX)){
				 this.currentIndex = 0L;
			 	}
			 else {
				 this.currentIndex = executionContext.getLong(CURRENT_INDEX);
				 if(LOG.isInfoEnabled()) LOG.info("skipping "+this.currentIndex+" variants");
				 for(long n=0L;n< this.currentIndex ;n++)
					{
					if(!this.iter.hasNext()) {
						throw new IllegalStateException("no more variants");
						}
					this.iter.next();
					}
			 	}
			
			}
		catch(final IOException err)
			{
			if(LOG.isErrorEnabled()) LOG.error("Cannot open Vcf",err);
			priv_close();
			throw new ItemStreamException(err);
			}
		}
	@Override
	public List<VariantContext> read()
			throws Exception, UnexpectedInputException, ParseException, NonTransientResourceException {
		if(iter!=null && iter.hasNext())
				{
				this.currentIndex ++;
				return Collections.singletonList(this.iter.next());
				}
		else
				{
				return null;
				}
		}
	
	@Override
	public void update(final ExecutionContext executionContext) throws ItemStreamException {
		 executionContext.putLong(CURRENT_INDEX, this.currentIndex);
		}
	
	private void priv_close()
		{
		CloserUtil.close(this.iter);
		this.iter=null;
		CloserUtil.close(this.vcfFileReader);
		this.vcfFileReader=null;
		}
	
	@Override
	public void close() throws ItemStreamException {
		priv_close();
		}
	}
