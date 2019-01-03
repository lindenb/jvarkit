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
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.batch.item.ExecutionContext;
import org.springframework.batch.item.ItemStreamException;
import org.springframework.batch.item.file.ResourceAwareItemWriterItemStream;
import org.springframework.core.io.Resource;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextBatchWriter implements 
	ResourceAwareItemWriterItemStream<List<VariantContext>> {
	private static final Log LOG = LogFactory.getLog(VariantContextBatchWriter.class);

	private Function<ExecutionContext,String> filenameFactory = null;
	private VariantContextWriter vcw = null;
	private boolean createMD5 = false;
	private File reference = null;
	
@Override
public void open(final ExecutionContext executionContext) throws ItemStreamException {
	if(this.filenameFactory==null) throw new ItemStreamException("resource is not defined");
	try {
		final String filename = this.filenameFactory.apply(executionContext);
		if(StringUtil.isBlank(filename)) throw new ItemStreamException("No output file defined.");
		if(LOG.isInfoEnabled()) LOG.info("Opening "+filename+" for writing");
		
		if(!Arrays.stream(IOUtil.VCF_EXTENSIONS).anyMatch(SUFF->filename.endsWith(SUFF)))
			{
			throw new ItemStreamException("Bad extension for a VCF file:" + filename);
			}
		final File vcfFile = new File(filename);
		VCFHeader header= SpringBatchUtils.getVcfHeader(executionContext);
		final VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
		vcwb.setOutputFile(vcfFile);
		if(this.reference!=null)
			{
			final SAMSequenceDictionary dic=SAMSequenceDictionaryExtractor.extractDictionary(this.reference);
			vcwb.setReferenceDictionary(dic);
			if(header.getSequenceDictionary()==null) {
				header = new VCFHeader(header);
				header.setSequenceDictionary(dic);
				}
			}
		else
			{
			vcwb.setReferenceDictionary(header.getSequenceDictionary());
			}
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

public void setFilenameFactory(final Function<ExecutionContext, String> resourceFactory) {
	this.filenameFactory = resourceFactory;
}

@Override
public void setResource(final Resource resource) {
	this.setFilenameFactory( CTX-> resource.getFilename() );
	}

public void setCreateMD5(boolean createMD5) {
	this.createMD5 = createMD5;
	}
public void setReference(final File reference) {
	this.reference = reference;
	}

private void priv_close() {
	CloserUtil.close(this.vcw);
	this.vcw = null;
	}
}
