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
package com.github.lindenb.jvarkit.tools.vcfconcat;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Example

### From stdin

```bash
$ find ./ -name "*.vcf" | grep Sample1 | java -jar dist/vcfconcat.jar > out.vcf
```

### From files

```bash
$ java -jar dist/vcfconcat.jar Sample1.samtools.vcf Sample1.gatk.vcf > out.vcf
```


END_DOC
*/
@Program(name="vcfconcat",
	keywords={"vcf"},
		description="Concatenante sorted VCF with same sample, does NOT merge genotypes"
		)
public class VcfConcat extends Launcher
	{
	private static final Logger LOG =Logger.build(VcfConcat.class).make();
	
	public static final String VARIANTSOURCE="VARIANTSOURCE";
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputfile=null;
	
	private Set<String> inputFiles=new HashSet<>();
	
	public VcfConcat()
		{
		}
	
	
	

	
	private int fromFiles(final VariantContextWriter out) throws IOException
		{
		List<VCFIterator> inputs=new ArrayList<VCFIterator>(this.inputFiles.size());
		List<String> inputFiles=new ArrayList<>(this.inputFiles.size());
		List<String> samples=new ArrayList<>();
		SAMSequenceDictionary dict=null;
		try
			{
			final Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
			
			/* open each vcf file */
			for(String vcfFile:this.inputFiles)
				{
				LOG.info("Opening "+vcfFile);
				VCFIterator r=VCFUtils.createVCFIterator(vcfFile);
				
				/* check VCF dict */
				VCFHeader header = r.getHeader();
				if(dict==null && inputs.isEmpty())
					{
					dict = header.getSequenceDictionary();
					}
				else if(!inputs.isEmpty() &&
					(
					(dict==null && header.getSequenceDictionary()!=null) ||
					(dict!=null && header.getSequenceDictionary()==null))
					)
					{
					LOG.error("not.the.same.sequence.dictionaries");
					return -1;
					}
				else if(!inputs.isEmpty() && dict!=null && 
					!SequenceUtil.areSequenceDictionariesEqual(dict, header.getSequenceDictionary()))
					{
					LOG.error("not.the.same.sequence.dictionaries");
					return -1;
					}
				/* check samples */
				if(inputs.isEmpty())
					{
					samples = header.getSampleNamesInOrder();
					}
				else if(!header.getSampleNamesInOrder().equals(samples))
					{
					LOG.error("No same samples");
					return -1;
					}
				
				metaData.addAll(header.getMetaDataInInputOrder());
				inputs.add(r);
				inputFiles.add(VCFUtils.escapeInfoField(vcfFile));
				}
			/* create comparator according to dict*/
			final Comparator<VariantContext> comparator = (dict==null?
					VCFUtils.createChromPosRefComparator():
					VCFUtils.createTidPosRefComparator(dict)
					);
			metaData.addAll(JVarkitVersion.getInstance().getMetaData(getClass().getSimpleName()));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			metaData.add(new VCFInfoHeaderLine(
					VARIANTSOURCE,1,VCFHeaderLineType.String,
					"Origin File of Varant"
					));
			VCFHeader h2 = new VCFHeader(
					metaData,
					samples
					);
			out.writeHeader(h2);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			for(;;)
				{
				/* get 'smallest' variant */
				VariantContext smallest = null;
				int idx=0;
				int best_idx = -1;
				
				while(idx < inputs.size())
					{
					VCFIterator in= inputs.get(idx);
					if(!in.hasNext())
						{
						CloserUtil.close(in);
						inputs.remove(idx);
						inputFiles.remove(idx);
						}
					else
						{
						VariantContext ctx = in.peek();
						if( smallest==null ||
							comparator.compare(smallest,ctx)>0)
							{
							smallest = ctx;
							best_idx = idx;
							}
						++idx;
						}
					}
				
				if(smallest==null ) break;
				
				final VariantContext ctx = progress.watch(inputs.get(best_idx).next());
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(VARIANTSOURCE, inputFiles.get(best_idx));
				out.add(vcb.make());
				}
			progress.finish();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(inputs);
			}

		}
	
		

	@Override
	public int doWork(final List<String> args) {
		VariantContextWriter w=null;
		BufferedReader r=null;
		try
			{
			if(args.isEmpty())
				{
				LOG.info("Reading filenames from stdin");
				r= IOUtils.openStdinForBufferedReader();
				this.inputFiles.addAll(r.lines().
					filter(line->!(line.trim().isEmpty() || line.startsWith("#"))).
					collect(Collectors.toSet())
					);
				r.close();
				r=null;
				}
			else
				{
				for(final String filename:args)
					{
					this.inputFiles.addAll(IOUtils.unrollFiles(Arrays.asList(filename)));	
					}
				}
			if(this.inputFiles.isEmpty())
				{
				LOG.error("No input");
				return -1;
				}
			w= super.openVariantContextWriter(this.outputfile);
			return fromFiles(w);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(r);
			}
		}

	public static void main(String[] args)	{
	new VcfConcat().instanceMainWithExit(args);
	}

	
	}
