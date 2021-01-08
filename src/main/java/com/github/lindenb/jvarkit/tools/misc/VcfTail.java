/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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

import java.util.LinkedList;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
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
	description="print the last variants of a vcf",
	keywords={"vcf"},
	modificationDate="20200518",
	creationDate="20131210"
	)
public class VcfTail extends OnePassVcfLauncher
	{
	private static final Logger LOG=Logger.build(VcfTail.class).make();
	@Parameter(names={"-n","-N","--count"},description="number of variants")
	private long count=10;
	@Parameter(names={"-c","--bycontig"},description="Print first variant for each contig; Implies VCF is sorted",order=1)
	private boolean by_contig=false;

	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter w) {
		if(count<0) {
			LOG.error("count < 0");
			return -1;
			}
		try {
			final VCFHeader header = in.getHeader();			
			final LinkedList<VariantContext> buffer=new LinkedList<VariantContext>();
			
			JVarkitVersion.getInstance().addMetaData(this, header);
			String prev_contig=null;
			w.writeHeader(header);
			while(in.hasNext())
				{
				final VariantContext ctx = in.next();
				if(this.by_contig &&
						(prev_contig==null || !prev_contig.equals(ctx.getContig())))
					{
					for(final VariantContext v2:buffer) w.add(v2);
					buffer.clear();
					prev_contig = ctx.getContig();
					}
				buffer.add(ctx);
				if(buffer.size() > this.count)
					{
					buffer.removeFirst();
					}
				}
			for(final VariantContext v2:buffer) w.add(v2);
			return 0;
		} catch (final Throwable err) {
			LOG.error(err);
			return -1;
		} finally
			{
			CloserUtil.close(in);
			CloserUtil.close(w);
			}
		}
		
	
		public static void main(final String[] args)
			{
			new VcfTail().instanceMainWithExit(args);
			}
		}
