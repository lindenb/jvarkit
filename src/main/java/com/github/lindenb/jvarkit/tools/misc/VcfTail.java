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

import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfTail
	extends AbstractVcfTail
	{
	public VcfTail()
		{
		}
	@Override
	public Collection<Throwable> initializeKnime()
		{
		if(this.count<0) return wrapException("bad value for count "+this.count);
		return super.initializeKnime();
	 	}

	@Override
	protected Collection<Throwable> doVcfToVcf(String inputName,
			VcfIterator in, VariantContextWriter out) throws IOException
			{
			try {
				final VCFHeader header=in.getHeader();
				final VCFHeader h2= addMetaData(new VCFHeader(header));
				final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header);
				out.writeHeader(h2);
				final LinkedList<VariantContext> L=new LinkedList<VariantContext>();
				while(in.hasNext() && L.size()< this.count && !out.checkError())
					{	
					L.add(progess.watch(in.next()));
					}
				while(in.hasNext())
					{
					L.add(progess.watch(in.next()));
					L.removeFirst();
					}
				for(VariantContext ctx:L)
					{
					out.add(ctx);
					}
				progess.finish();
				return RETURN_OK;
				}
			finally
				{
				CloserUtil.close(out);
				out=null;
				}
			}
	
	 @Override
	protected Collection<Throwable> call(String inputName) throws Exception
		{
		return doVcfToVcf(inputName);
		}
		
		
	public static void main(String[] args)
		{
		new VcfTail().instanceMainWithExit(args);
		}
	}
