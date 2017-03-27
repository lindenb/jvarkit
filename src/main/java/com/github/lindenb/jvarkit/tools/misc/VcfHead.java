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
import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.validators.FilesValidators;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 BEGIN_DOC
 
 Hello world
 
 END_DOC
 */

public class VcfHead
	{
	@Parameter(names={"-n","--count"},description="number of variants",validateValueWith=FilesValidators.X.class)
	private long _count=10;
	@Parameter(names={"-c","--bycontig"},descriptionKey="Print first variant for each contig; Implies VCF is sorted",order=1,description="number of variants")
	private boolean _by_contig=false;;

	public VcfHead()
		{
		}
	 
	public VcfHead count(long n) { this._count=n;return this;}
	public VcfHead byContig() { this._by_contig=true;return this;}
	
	public VariantContextWriter open(final VariantContextWriter w) {
		if(this._count<0L) throw new JvarkitException.CommandLineError("bad value found count "+this._count);
		return new MyWriter(w,_count,this._by_contig);
		}
	
	private static class MyWriter extends DelegateVariantContextWriter
		{
		private final long count;
		private final boolean by_contig;
		private String prev_contig=null;
		private long n=0L;
		MyWriter(VariantContextWriter w,long n,boolean by_contig) {
			super(w);
			this.count=n;
			this.by_contig=by_contig;
			}
		@Override
		public void add(final VariantContext ctx) {
			if(isClosed()) {
				return;
				}
			if(by_contig && (this.prev_contig==null ||
					!this.prev_contig.equals(ctx.getContig())))
				{
				this.prev_contig = ctx.getContig();
				this.n=0L;
				};
				
			if(n<count) {
				super.add(ctx);
				++n;
				}
			else if(!this.by_contig)
				{
				close();
				}
			}
		
		}
		 
		 
		
	 	
	public static  class Launcher extends com.github.lindenb.jvarkit.util.jcommander.Launcher
		{
		@ParametersDelegate
		private VcfHead instance=new VcfHead();
		@Parameter(names={"-o","--out"},required=false,description="Output vcf , ot stdin")
		private VariantContextWriter out=new VcfWriterOnDemand();
		
		Launcher() {
			
			}
		@Override
		public int doWork(List<String> args) {
			try {
				VcfIterator in = VCFUtils.createVcfIterator(super.oneFileOrNull(args));
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
		new VcfHead.Launcher().instanceMainWithExit(args);
		}
	}
