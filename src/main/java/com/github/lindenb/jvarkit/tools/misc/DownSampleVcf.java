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

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

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
	description="DownSample a VCF",
	keywords={"vcf"}
	)
public class DownSampleVcf extends Launcher
	{
	private final Logger LOG=Logger.build(DownSampleVcf.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected File outputFile = null;

	@Parameter(names="-n",description="output size")
	private int reservoir_size=10;
	@Parameter(names="-N",description=" random seed")

	private long seed=System.currentTimeMillis();
	
	private DownSampleVcf()
		{
		}
	

	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		final Random rand=new Random(this.seed);
		final List<VariantContext>  buffer=new ArrayList<VariantContext>(this.reservoir_size);
		final VCFHeader h2=new VCFHeader(in.getHeader());
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(in.getHeader());
		out.writeHeader(h2);
		if(this.reservoir_size!=0)
			{
			while(in.hasNext())
				{	
				if(buffer.size() < this.reservoir_size)
					{
					buffer.add(progess.watch(in.next()));
					}
				else
					{
					buffer.set(rand.nextInt(buffer.size()), progess.watch(in.next()));
					}
				}
			
			}
		for(VariantContext ctx:buffer)
			{
			out.add(ctx);
			}
		progess.finish();
		LOG.info("done");
		return 0;
		}
	@Override
	public int doWork(List<String> args) {
		
		return doVcfToVcf(args, outputFile);
		}

	public static void main(String[] args)
		{
		new DownSampleVcf().instanceMainWithExit(args);
		}
	}
