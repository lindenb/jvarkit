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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;


public class TabixVcfFileReader extends AbstractTabixObjectReader<VariantContext>
	//,Iterable<VariantContext> NO, not a true iterator
	{
    private final AbstractVCFCodec vcfCodec;
    private final VCFHeader vcfHeader;
   
    
    public TabixVcfFileReader(final String uri) throws IOException
    	{
    	super(uri);
    	
    	final List<String> stack=new ArrayList<String>();
    	String line;
    	while((line=super.tabix.readLine())!=null && line.startsWith("#"))
    		{
    		stack.add(line);
    		if(line.startsWith("#CHROM\t")) break;
    		}
    	final VCFUtils.CodecAndHeader cah=VCFUtils.parseHeader(stack);
    	this.vcfHeader=cah.header;
    	this.vcfCodec=cah.codec;
    	}
    
    
	public AbstractVCFCodec getCodec()
		{
		return this.vcfCodec;
		}
	
	public VCFHeader getHeader()
		{
	    return this.vcfHeader;
		}
    
    @Override
    protected  Iterator<VariantContext> iterator(final Iterator<String> delegate)
		{
		return new MyIterator(delegate);
		}
    
    
    private class MyIterator
    	extends AbstractMyIterator
    	{
    	MyIterator(Iterator<String> delegate)
    		{
    		super(delegate);
    		}
    	@Override
    	public VariantContext next() {
    		return getCodec().decode(delegate.next());
    		}
    	}	
    
	}
