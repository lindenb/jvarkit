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

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.Writer;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
/**
BEGIN_DOC


## Example

```
$ java -jar dist/vcfmakedict.jar test.vcf  mutations.vcf

@HD	VN:1.5
@SQ	SN:3	LN:124490595
@SQ	SN:19	LN:57230811
@SQ	SN:rotavirus	LN:1064
```


END_DOC

*/
@Program(name="vcfmakedict",
	description="Create a SAM Sequence Dictionary from a set of VCF files.",
	keywords={"vcf","dict","fai"},
	modificationDate="2090528"
	)
public class VcfCreateDictionary extends Launcher
	{
	private static final Logger LOG=Logger.build(VcfCreateDictionary.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	
	public VcfCreateDictionary()
		{
		}


	@Override
	public int doWork(final List<String> args) {
		VCFIterator in = null;
		Writer out= null;
		
		if (this.outputFile != null && !outputFile.getFileName().endsWith(".dict")) {
				LOG.error("output file should end with .dict :" + this.outputFile);
				return -1;
			}
		
		try {
			final LinkedHashMap<String, Integer> buildNewDictionary =new LinkedHashMap<String, Integer>();
			
			int optind=0;
			do
				{
				in = super.openVCFIterator(args.isEmpty()?null:args.get(optind));
				final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(in.getHeader()).logger(LOG).build();
	
				while(in.hasNext())
					{
					final VariantContext ctx = progress.apply(in.next());
					Integer length= buildNewDictionary.getOrDefault(ctx.getContig(),0);
					if(ctx.getEnd()>length)
						{
						buildNewDictionary.put(ctx.getContig(),ctx.getEnd());
						}
					}
				progress.close();
				in.close();
				in=null;
				++optind;
				} while(optind < args.size());
			

			final List<SAMSequenceRecord> list = new ArrayList<SAMSequenceRecord>(buildNewDictionary.size());
			for (final String k : buildNewDictionary.keySet()) {
				list.add(new SAMSequenceRecord(k, buildNewDictionary.get(k)));
				}
			final SAMFileHeader sfh = new SAMFileHeader();
			sfh.setSequenceDictionary(new SAMSequenceDictionary(list));
			final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
			codec.setValidationStringency(htsjdk.samtools.ValidationStringency.SILENT);
			out = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			codec.encode(out, sfh);
			out.flush();
			out.close();
			out=null;
			
			return 0;
			} 
		catch (final Exception err2) {
			LOG.error(err2);
			return -1;
			} 
		finally {
			CloserUtil.close(out);
			CloserUtil.close(in);
			}
		
	}

	public static void main(final String[] args)
	{
		new VcfCreateDictionary().instanceMainWithExit(args);
	}
}
