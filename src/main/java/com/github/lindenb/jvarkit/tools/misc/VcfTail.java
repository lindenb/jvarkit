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
* 2015 : moving to knime
* 2014 : creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.PostponedVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import htsjdk.variant.vcf.VCFIterator;

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
@Program(
	name="vcftail",
	description="print the last variants of a vcf",keywords={"vcf"})
public class VcfTail extends Launcher
	{
	private static final Logger LOG=Logger.build(VcfTail.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT,required=false)
	private File output=null;
	@ParametersDelegate
	private PostponedVariantContextWriter.WritingVcfConfig writingVcfArgs = new PostponedVariantContextWriter.WritingVcfConfig();
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	
	@XmlType(name="vcftail")
	@XmlRootElement(name="vcftail")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
			{
			@XmlElement(name="count")
			@Parameter(names={"-n","--count"},description="number of variants")
			private long count=10;
			@XmlElement(name="by-contig")
			@Parameter(names={"-c","--bycontig"},descriptionKey="Print last variant for each contig; Implies VCF is sorted",order=1,description="number of variants")
			private boolean by_contig=false;
			
			public void setCount(long count) {
				this.count = count;
				}
			
			public void setByContig(boolean by_contig) {
				this.by_contig = by_contig;
				}
			
			private class CtxWriter extends DelegateVariantContextWriter
				{
				private String prev_contig=null;
				private final LinkedList<VariantContext> buffer=new LinkedList<VariantContext>();

				private void dump()
					{
					for(final VariantContext ctx:this.buffer) 
						{
						super.add(ctx);
						}
					this.buffer.clear();
					}
				
				CtxWriter(final VariantContextWriter delegate) {
					super(delegate);
					}
				@Override
				public void writeHeader(final VCFHeader header) {
					super.writeHeader(header);
					}
				
				@Override
				public void add(final VariantContext ctx) {
					if(CtxWriterFactory.this.by_contig &&
							(this.prev_contig==null || !this.prev_contig.equals(ctx.getContig())))
						{
						dump();
						this.prev_contig = ctx.getContig();
						}
					this.buffer.add(ctx);
					if(buffer.size() > CtxWriterFactory.this.count)
						{
						this.buffer.removeFirst();
						}
					}
				
				@Override
				public void close() {
					dump();
					super.close();
					}
				}

			@Override
			public VariantContextWriter open(final VariantContextWriter delegate) {
				return new CtxWriter(delegate);
				}
			
			
		}
	
		
		
		public VcfTail()
			{
			}
		
		@Override
		protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
			return new PostponedVariantContextWriter(writingVcfArgs,stdout(),outorNull);
			}
		
		@Override
		protected int doVcfToVcf(
				final String inputName, 
				final VCFIterator in,
				final VariantContextWriter delegate) {
	
			try
				{
				final VariantContextWriter out  = this.component.open(delegate);
				out.writeHeader(in.getHeader());
				
				final ProgressFactory.Watcher<VariantContext> progress= ProgressFactory.newInstance().dictionary(in.getHeader()).logger(LOG).build();
				while(in.hasNext())
					{
					out.add(progress.apply(in.next()));
					}
				progress.close();
				out.close();
				return 0;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			}
		
		@Override
		public int doWork(final List<String> args) {
			return doVcfToVcf(args,output);
			}
		
	
		public static void main(final String[] args)
			{
			new VcfTail().instanceMainWithExit(args);
			}
		}
