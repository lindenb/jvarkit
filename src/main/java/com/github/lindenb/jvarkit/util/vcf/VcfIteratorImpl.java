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
* 2014 creation

*/
package com.github.lindenb.jvarkit.util.vcf;

import java.io.BufferedReader;
import java.io.InputStream;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

public class VcfIteratorImpl implements htsjdk.variant.vcf.VCFIterator
	{
	/** associated VCF codec */
    private final AbstractVCFCodec vcfCodec;
	/** associated VCF header */
    private final VCFHeader vcfHeader;
	/** associated line iterator */
    private final LineIterator lineIterator;
    
	public VcfIteratorImpl(final InputStream vcfStream)
		{
		this( new LineIteratorImpl(new SynchronousLineReader(vcfStream)));
		}
	
	public VcfIteratorImpl(final BufferedReader vcfStream)
		{
		this(IOUtils.toLineIterator(vcfStream));
		}
	
	public VcfIteratorImpl(final LineIterator r)
		{
		this.lineIterator = r;
	   // this.vcfHeader = (VCFHeader) vcfCodec.readActualHeader(lineIterator);
	    
	    final VCFUtils.CodecAndHeader cah=VCFUtils.parseHeader(r);
	    this.vcfHeader=cah.header;
	    this.vcfCodec=cah.codec;
		}
	
	//@Override
	public AbstractVCFCodec getCodec()
		{
		return this.vcfCodec;
		}
	
	@Override
    public VCFHeader getHeader()
    	{
        return this.vcfHeader;
    	}

	@Override
    public VariantContext peek()
    	{
    	return vcfCodec.decode(lineIterator.peek());
		}
    
	@Override
	public boolean hasNext() {
        return lineIterator.hasNext();
		}
	@Override
	public VariantContext next() {
		return vcfCodec.decode(lineIterator.next());
		}
	
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
		}
	
	@Override
	public void close()
		{
		CloserUtil.close(lineIterator);
		}
	@Override
	public String toString() {
		return "VCF Iterator. Codec: "+vcfCodec;
		}
	}
