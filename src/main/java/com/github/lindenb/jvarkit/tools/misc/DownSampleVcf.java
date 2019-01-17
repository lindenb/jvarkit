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

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC


## Example

```bash
$ curl -skL "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz" |\
  gunzip -c |\
java -jar dist/downsamplevcf.jar -n 100 > output.vcf
```

END_DOC
 */

@Program(
	name="downsamplevcf",
	description="DownSample a VCF. Will keep 'n' random variants in a vcf.",
	keywords={"vcf"}
	)
public class DownSampleVcf extends Launcher
	{
	private final Logger LOG=Logger.build(DownSampleVcf.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected File outputFile = null;
	@Parameter(names="-n",description="output size. keep 'n' random variants in the input vcf")
	private int reservoir_size=10;
	@Parameter(names="-N",description="random seed. -1==use current time")
	private long seed=-1L;
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator in, final VariantContextWriter out) {
		final Random rand=new Random(this.seed==-1L?System.currentTimeMillis():this.seed);
		final List<VariantContext>  buffer=new ArrayList<>(this.reservoir_size);
		
		final VCFHeader h1=in.getHeader();
		final SAMSequenceDictionary dict=h1.getSequenceDictionary();
		final VCFHeader h2=new VCFHeader(h1);
		super.addMetaData(h2);
		
		final ProgressFactory.Watcher<VariantContext> progess=ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
		out.writeHeader(h2);
		if(this.reservoir_size!=0)
			{
			while(in.hasNext())
				{	
				final VariantContext ctx = progess.apply(in.next());
				if(buffer.size() < this.reservoir_size)
					{
					buffer.add(ctx);
					}
				else
					{
					buffer.set(rand.nextInt(buffer.size()), ctx);
					}
				}
			}
		
		buffer.
			stream().
			sorted(dict==null?VCFUtils.createChromPosRefComparator():VCFUtils.createTidPosRefComparator(dict)).
			forEach(V->out.add(V));
		progess.close();
		return 0;
		}
	@Override
	public int doWork(final List<String> args) {
		return doVcfToVcf(args, outputFile);
		}

	public static void main(final String[] args)
		{
		new DownSampleVcf().instanceMainWithExit(args);
		}
	}
