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

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.PostponedVariantContextWriter;
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
@Program(
	name="vcftail",
	description="print the last variants of a vcf",keywords={"vcf"})
public class VcfTail extends com.github.lindenb.jvarkit.util.jcommander.Launcher
	{
	private static final Logger LOG=Logger.build(VcfTail.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT,required=false)
	private File output=null;
	@Parameter(names={"-n","--count"},description="number of variants")
	private int count=10;
	@Parameter(names={"-c","--bycontig"},descriptionKey="Print last variant for each contig; Implies VCF is sorted",order=1,description="number of variants")
	private boolean by_contig=false;
	@ParametersDelegate
	private PostponedVariantContextWriter.WritingVcfConfig writingVcfArgs = new PostponedVariantContextWriter.WritingVcfConfig();

	
	public VcfTail()
		{
		}
	
	@Override
	protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
		return new PostponedVariantContextWriter(writingVcfArgs,stdout(),outorNull);
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		String prev_contig=null;
		final LinkedList<VariantContext> buffer=new LinkedList<VariantContext>();

		try
			{
			final Runnable dump=()->{
				for(final VariantContext ctx:buffer) 
					{
					out.add(ctx);
					}
				buffer.clear();
				};
			
			final VCFHeader header2=  new VCFHeader(in.getHeader());
			out.writeHeader(header2);
			while(in.hasNext())
				{
				final VariantContext ctx = in.next();
				if(by_contig && (prev_contig==null ||
						!prev_contig.equals(ctx.getContig())))
					{
					dump.run();;
					prev_contig=ctx.getContig();
					}
				buffer.add(ctx);
				if(buffer.size()>this.count)
					{
					buffer.removeFirst();
					}

				}
			dump.run();
			return 0;
			}
		catch(Exception err)
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
