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
* 2014 creation
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.util.Collection;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfHead
	extends AbstractVcfHead
	{
	public VcfHead()
		{
		}
	 
		 @Override
		public Collection<Throwable> initializeKnime() {
			if(this.count<0) return wrapException("bad value found count "+this.count);
			return super.initializeKnime();
		 	}
		/* public for knime */
		@Override
		public Collection<Throwable> doVcfToVcf(
				String inputName,
				VcfIterator in,
				VariantContextWriter out
				) throws IOException {
			final VCFHeader header=in.getHeader();
			VCFHeader h2=addMetaData(new VCFHeader(header));
			SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			out.writeHeader(h2);
			int n=0;
			while(in.hasNext() && n< super.count  && !out.checkError())
				{
				out.add(progess.watch(in.next()));
				++n;
				}
			progess.finish();
			return RETURN_OK;
			}
		
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			return doVcfToVcf(inputName);
			}
	 	
	
	public static void main(String[] args)
		{
		new VcfHead().instanceMainWithExit(args);
		}
	}
