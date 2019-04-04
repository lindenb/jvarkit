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
* 2016 creation

*/
package com.github.lindenb.jvarkit.util.vcf;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
 * used to store a VCF in memory. If there are
 * too many variant, the content is flushed to disk.
 * 
 * @author lindenb
 *
 */
public class VCFBuffer implements VariantContextWriter{
private static final Logger LOG = Logger.build(VCFBuffer.class).make();

/** vcf header */
private VCFHeader header=null;
/** memory buffer */
private final List<VariantContext> buffer = new ArrayList<>();
/** file buffer */
private File tmpFile = null;
/** directory where to create tmpFile */
private final File tmpDir;
/** variant context writer. Null if buffered in memory */
private VariantContextWriter variantContextWriter= null;
/** max number of variants in memory */
private final int maxRecordsInRam;
/** done adding , we can always call 'iterator()' */
private boolean done_adding=false;
/** we cannot use it anymore, tmp File have been deleted */
private boolean disposed=false;
public VCFBuffer(final int maxRecordsInRam,final File tmpDir) {
	this.maxRecordsInRam = maxRecordsInRam;
	this.tmpDir = tmpDir;
	}

public VCFBuffer() {
	this(1000,null);
}

/** close and delete temporary files */
public void dispose() {
	close();
	if(this.tmpFile!=null) this.tmpFile.delete();
	if(this.buffer!=null) this.buffer.clear();
	this.tmpFile=null;
	disposed=true;
	}

@Override
protected void finalize() throws Throwable {
	dispose();
	super.finalize();
	}

public Stream<VariantContext> stream()
	{
	close();
	if(this.disposed) throw new IllegalStateException("buffer was disposed");
	if(this.header==null) {
		throw new IllegalStateException("header was not set");	
		}
	if(this.tmpFile==null) 
		{
		return this.buffer.stream();
		}
	else
		{
		final VCFFileReader reader= new VCFFileReader(this.tmpFile,false);
		final CloseableIterator<VariantContext> iter = reader.iterator();
		return  StreamSupport.stream(new IterableAdapter<VariantContext>(iter).spliterator(), false).onClose(
				()->{CloserUtil.close(iter);CloserUtil.close(reader);}
				);
		}
	}

public VCFIterator iterator() {
close();
if(this.disposed) throw new IllegalStateException("buffer was disposed");	
if(this.header==null) {
	throw new IllegalStateException("header was not set");	
	}
if(this.tmpFile==null) 
	{
	return new ArrayIterator();
	}
else
	{
	try {
		return VCFUtils.createVCFIteratorFromFile(this.tmpFile);
	} catch (final IOException e) {
		throw new RuntimeIOException(e);
	}
	}
}

@Override
public void setHeader(final VCFHeader header) {
	throw new IllegalArgumentException("setHeader shouldn't be called");	
	}

@Override
public void writeHeader(final VCFHeader header) {
	if(this.done_adding) throw new IllegalArgumentException("iterator() already called");
	if(this.header!=null) throw new IllegalArgumentException("Header already set");
	this.header = header;
	}

public VCFHeader getHeader()
	{
	return this.header;
	}

@Override
public void close() {
	CloserUtil.close(this.variantContextWriter);
	this.variantContextWriter=null;
	this.done_adding=true;
	}

@Override
public boolean checkError() {
	return false;
}
@Override
public void add(final VariantContext vc) {
	if(this.done_adding) throw new IllegalArgumentException("iterator() already called");
	if(this.header==null) throw new IllegalArgumentException("Header wasn't set");
	if(this.variantContextWriter!=null) {
		this.variantContextWriter.add(vc);
		}
	else if(this.buffer.size()+1>= this.maxRecordsInRam )
		{
		try {
			this.tmpFile = File.createTempFile("buffer.", ".vcf.gz",this.tmpDir);
			LOG.debug("Flushing to disk "+this.tmpFile);
			this.tmpFile.deleteOnExit();
			this.variantContextWriter = VCFUtils.createVariantContextWriter(this.tmpFile);
			this.variantContextWriter.writeHeader(this.header);
			for(final VariantContext bvc:this.buffer) {
				this.variantContextWriter.add(bvc);
			}
			this.buffer.clear();
			this.variantContextWriter.add(vc);
		} catch (IOException e) {
			throw new RuntimeIOException(e);
			}
		}
	else
		{
		this.buffer.add(vc);
		}
	}

private class ArrayIterator implements VCFIterator {
	int index=-1;
	
	@Override
	public VCFHeader getHeader() {
		return VCFBuffer.this.header;
		}
	
	//@Override
	//public AbstractVCFCodec getCodec() { return VCFUtils.createDefaultVCFCodec(); }
	@Override
	public boolean hasNext() {
		return index+1< VCFBuffer.this.buffer.size();
		}
	@Override
	public VariantContext next() {
		index++;
		return VCFBuffer.this.buffer.get(index);
		}
	@Override
	public VariantContext peek() {
		return  VCFBuffer.this.buffer.get(index+1);
		}
	@Override
	public void close() {
		
	}
}


}




