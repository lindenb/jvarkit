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
package com.github.lindenb.jvarkit.tools.vcfcmp;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.Closeable;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import htsjdk.variant.vcf.VCFIterator;

public class EqualRangeVcfIterator 
	implements Closeable
	{	
	private VCFIterator in;
	private final List<VariantContext> buffer=new ArrayList<>();
	private VariantContext prev_ctx=null;
	private final Comparator<VariantContext> ctxComparator;
	public EqualRangeVcfIterator(
			final VCFIterator in,
			final Comparator<VariantContext> ctxComparator)
		{
		this.in = in;
		this.ctxComparator=ctxComparator;
		}
	
	public List<VariantContext> next(final VariantContext ctx) throws IOException
		{
		if(this.in==null) return Collections.emptyList();
		if(!this.buffer.isEmpty())
			{
			final int d= this.ctxComparator.compare(buffer.get(0), ctx);
			if( d == 0) return this.buffer;
			else if( d >0 ) return Collections.emptyList();
			//all variants in buffer was before ctx, continue
			this.buffer.clear();
			}
		for(;;)
			{
			if(!this.in.hasNext())
				{
				close();
				break;
				}
			final VariantContext ctx2= this.in.peek();
			if(this.prev_ctx!=null && this.ctxComparator.compare(this.prev_ctx, ctx2)>0)
				{
				throw new IOException(
						"Bad order : got\n "+this.prev_ctx+"\nbefore\n "+ctx2+"\n"
						);
				}
			this.prev_ctx =  ctx2;
			int d = this.ctxComparator.compare(ctx2, ctx);
			if( d < 0)
				{
				in.next();//consumme
				continue;
				}
			else if(d>0)
				{
				break;
				}
			else //equal
				{
				this.buffer.add(in.next());
				}
			}
		return this.buffer;
		}
	
	@Override
	public void close() {
		CloserUtil.close(this.in);
		in=null;
		}
	}
