/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.variant.vcf;

import java.io.ByteArrayInputStream;
import java.io.IOException;

import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;


public class TabixVcfReader implements VCFReader
	{
	private final TabixReader tabix;
    private final VCFCodec vcfCodec = new VCFCodec();
    private final VCFHeader vcfHeader;
    
    public TabixVcfReader(final String uri) {
    	try {
    		this.tabix = new TabixReader(uri);
    		}
    	catch(final IOException err) {
    		throw new RuntimeIOException(err);
    		}
    	
    	try {
	    	final StringBuilder sb = new StringBuilder();
	    	String line;
	    	while((line=this.tabix.readLine())!=null && line.startsWith(VCFHeader.HEADER_INDICATOR))
	    		{
	    		sb.append(line).append('\n');
	    		if(line.startsWith("#CHROM\t")) break;
	    		}
	    	try(ByteArrayInputStream bais= new ByteArrayInputStream(sb.toString().getBytes())) {
	    		this.vcfHeader = (VCFHeader) this.vcfCodec.readActualHeader(this.vcfCodec.makeSourceFromStream(bais));
	    		}
    		}
    	catch(final IOException err) {
    		this.tabix.close();
    		throw new RuntimeIOException(err);
    		}
    	
    	}
    
    @Override	
	public VCFHeader getHeader()
		{
	    return this.vcfHeader;
		}
    
    @Override
    public CloseableIterator<VariantContext> query(String chrom, int start, int end) {
    	final TabixReader.Iterator tbxIter = this.tabix.query(chrom, start, end);
    	return new AbstractCtxIterator() {
    		@Override
    		protected String nextLine() throws IOException {
    			return tbxIter.next();
    			}
    		};
    	}
    
    @Override
    public CloseableIterator<VariantContext> iterator() {
    	return new AbstractCtxIterator() {
    		@Override
    		protected String nextLine() throws IOException {
    			return tabix.readLine();
    			}
    		};
    	}
    
    @Override
    public boolean isQueryable() {
    	return true;
    	}
    
    @Override
    public void close() {
    	this.tabix.close();
    	}
    
    @Override
    public String toString() {
    	return this.tabix.getSource();
    	}
    
    
    private abstract class AbstractCtxIterator extends AbstractCloseableIterator<VariantContext>
		{
	    protected abstract String nextLine() throws IOException;
		@Override
		protected VariantContext advance() {
			try{
				final String line = this.nextLine();
	    		if(line==null) return null;
	    		return TabixVcfReader.this.vcfCodec.decode(line);
	    		}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() {
			//nothing
			}
		}	
    
	}
