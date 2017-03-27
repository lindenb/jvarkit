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
import java.util.LinkedList;
import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**

```
 BEGIN_DOC

## Example

```
$  curl -s "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
java -jar dist/vcftail.jar -n2 |\
grep -v "##"| cut -f 1,2,4,5

#CHROM  POS REF ALT
chr1    935492  G   T
chr1    1334052 CTAGAG  C
```

END_DOC

**/
public class VcfTail
	{
	@Parameter(names={"-n","--count"},description="number of variants")
	private int _count=10;
	@Parameter(names={"-c","--bycontig"},descriptionKey="Print last variant for each contig; Implies VCF is sorted",order=1,description="number of variants")
	private boolean _by_contig=false;

	
	public VcfTail()
		{
		}

	public VcfTail setByContig(boolean _by_contig) {
		this._by_contig = _by_contig;
		return this;
	}
	
	public VcfTail setCount(int _count) {
		this._count = _count;
		return this;
		}
	
	public VariantContextWriter open(VariantContextWriter delegate)
		{
		if(this._count<0) throw new JvarkitException.CommandLineError("bad value for count "+this._count);
		return new MyWriter(delegate,_count,_by_contig);
		}
	
	private static class MyWriter extends DelegateVariantContextWriter
		{
		private final int count;
		private final boolean by_contig;
		private String prev_contig=null;
		final LinkedList<VariantContext> buffer=new LinkedList<VariantContext>();
		MyWriter(VariantContextWriter w,int count,boolean by_contig) {
			super(w);
			this.count=count;
			this.by_contig=by_contig;
			}
		private void dump()
			{
			for(VariantContext ctx:buffer) 
				{
				getDelegate().add(ctx);
				}
			buffer.clear();
			}
		@Override
		public void add(final VariantContext ctx) {
			if(isClosed()) {
				return;
				}
			if(by_contig && (this.prev_contig!=null ||
					!this.prev_contig.equals(ctx.getContig())))
				{
				dump();
				prev_contig=ctx.getContig();
				}
			buffer.add(ctx);
			if(buffer.size()>this.count)
				{
				buffer.removeFirst();
				}
			}
		@Override
		public void close() {
			dump();
			super.close();
			}
		}

	
	@Program(name="vcftail",description="print the last variants of a vcf")
	public static class Launcher extends com.github.lindenb.jvarkit.util.jcommander.Launcher
		{
		@ParametersDelegate
		private VcfTail instance=new VcfTail();
		@Parameter(names={"-o","--out"},required=false,description="Output vcf , ot stdin")
		private VariantContextWriter out=new VcfWriterOnDemand();
		
		@Override
		public int doWork(final List<String> args) {
			try {
				final VcfIterator in = VCFUtils.createVcfIterator(super.oneFileOrNull(args));
				final VariantContextWriter out=instance.open(this.out);
				out.writeHeader(in.getHeader());
				while(in.hasNext()) {
					out.add(in.next());
					}
				in.close();
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			return 0;
			}	
		}
		
	public static void main(String[] args)
		{
		new VcfTail.Launcher().instanceMainWithExit(args);
		}
	}
