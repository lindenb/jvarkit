/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;
import java.util.function.Predicate;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

/** net.sf.picard.vcf.VCFIterator deleted from from picard 1.100 */
public interface VcfIterator extends Iterator<VariantContext>,Closeable
	{
	/** get the codec for this vcf */
	public AbstractVCFCodec getCodec();
	/** get the VCF header */ 
	public VCFHeader getHeader();
	/** peek the next Variant */
    public VariantContext peek();
    
    /** wrap as a filtering iterator */
    public static VcfIterator filter(final VcfIterator delegate,final Predicate<VariantContext> predicate)
    	{
         class FilteringVcfIterator
         	extends AbstractIterator<VariantContext>
    		implements VcfIterator
    		{
        	private final VcfIterator delegate;
        	private final Predicate<VariantContext> predicate;
        	FilteringVcfIterator(final VcfIterator delegate,final Predicate<VariantContext> predicate) {
    			this.delegate = delegate;
    			this.predicate=predicate;
    			}
    		@Override
    		public AbstractVCFCodec getCodec() { return this.delegate.getCodec();}
    		@Override
    		public VCFHeader getHeader() { return this.delegate.getHeader();}
    		@Override
    		public void close() throws IOException {this.delegate.close();}
    		@Override
    		protected VariantContext advance() {
    			while(this.delegate.hasNext())
    				{
    				final VariantContext ctx=this.delegate.next();
    				if(this.predicate.test(ctx)) return ctx;;
    				}
    			return null;
    			}
    		}
    	return new FilteringVcfIterator(delegate,predicate);
    	}
    }
