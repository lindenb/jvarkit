/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.batch.item.ExecutionContext;
import org.springframework.batch.item.ItemStreamException;
import org.springframework.batch.item.file.ResourceAwareItemWriterItemStream;
import org.springframework.core.io.Resource;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextBatchWriter implements 
	ResourceAwareItemWriterItemStream<List<VariantContext>> {
	private static final Log LOG = LogFactory.getLog(VariantContextBatchWriter.class);

	private Resource resource = null;
	private VariantContextWriter vcw = null;
	private boolean createMD5 = false;
	
@Override
public void open(final ExecutionContext executionContext) throws ItemStreamException {
	if(this.resource==null) throw new ItemStreamException("resource is not defined");
	try {
		if(LOG.isInfoEnabled()) LOG.info("Opening "+this.resource);
		final File vcfFile = this.resource.getFile();
		if(!Arrays.stream(IOUtil.VCF_EXTENSIONS).anyMatch(SUFF->vcfFile.getName().endsWith(SUFF)))
			{
			throw new ItemStreamException("Bad extension for a VCF file:" + vcfFile);
			}
		final VCFHeader header= SpringBatchUtils.getVcfHeader(executionContext);
		final VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
		vcwb.setOutputFile(vcfFile);
		vcwb.setReferenceDictionary(header.getSequenceDictionary());
		vcwb.setCreateMD5(this.createMD5);
		this.vcw = vcwb.build();
		this.vcw.writeHeader(header);
		}
	catch(final Exception err)
		{
		priv_close();
		throw new ItemStreamException(err);
		}
	}

@Override
public void close() throws ItemStreamException {
	priv_close();
	}
@Override
public void update(final ExecutionContext executionContext) throws ItemStreamException {
	// TODO Auto-generated method stub
	
	}
@Override
public void write(List<? extends List<VariantContext>> variants) throws Exception {
	if(this.vcw==null) return;
	variants.stream().flatMap(L->L.stream()).forEach(ctx->this.vcw.add(ctx));
	}
@Override
public void setResource(final Resource resource) {
	this.resource = resource;
	}

public void setCreateMD5(boolean createMD5) {
	this.createMD5 = createMD5;
	}

private void priv_close() {
	CloserUtil.close(this.vcw);
	this.vcw = null;
	}
}
