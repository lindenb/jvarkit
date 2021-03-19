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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/biostar497922.jar -n 10 -o TMP src/test/resources/rotavirus_rf.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000001.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000002.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000003.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000004.vcf.gz
[INFO][Biostar497922]Writing TMP/split.000005.vcf.gz
[INFO][Biostar497922]. Completed. N=45. That took:0 second


t$ find TMP/ -type f -name "*.vcf.gz" | sort | while read F; do echo -n "$F " && gunzip -c $F | grep -v "#" | wc -l  ; done
TMP/split.000001.vcf.gz 10
TMP/split.000002.vcf.gz 10
TMP/split.000003.vcf.gz 10
TMP/split.000004.vcf.gz 10
TMP/split.000005.vcf.gz 5

```


END_DOC
*/
@Program(
		name="biostar497922",
		description="Split VCF into separate VCFs by SNP count",
		keywords={"vcf"},
		biostars=497922,
		modificationDate="20210319",
		creationDate="20210319"
		)
public class Biostar497922 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar497922.class).make();

	@Parameter(names={"-o","--output"},description="Output directory",required=true)
	private Path outputDirectory = null;
	@Parameter(names={"--count","-n"},description="number of variants per vcf")
	private long count=-1L;

	@Parameter(names={"--index"},description="start numering file index from 'x'")
	private int start_index=1;
	@Parameter(names={"--prefix"},description="file prefix")
	private String prefix="split";
	@Parameter(names={"-C","--chrom"},description="by chromosomes when using '-n'.")
	private boolean consider_chromosome = false;
	@Parameter(names={"-m","--manifest"},description="output manifest to this file.")
	private Path manifestPath = null;
	@Parameter(names={"-D","--length"},description="max distance beween first and last variant (ignore if <=0) . " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int max_interval_length = -1;
	@Parameter(names={"-d","--distance"},description="max distance beween consecutive variants (ignore if <=0) . " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int max_distance = -1;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();

	
	
	@Override
	public int doWork(final List<String> args) {
			if(this.count<1L && max_interval_length<0 && max_distance<0)
				{
				LOG.error("no criteria was specified");
				return -1;
				}
			if(StringUtils.isBlank(this.prefix))
				{
				LOG.error("empty prefix");
				return -1;
				}
			IOUtil.assertDirectoryIsWritable(this.outputDirectory);
			PrintWriter manifest = null;
			VCFIterator in=null;
			final String inputName= oneFileOrNull(args);
			try
				{
				in = super.openVCFIterator(inputName);
				final VCFHeader header=in.getHeader();
				manifest = this.manifestPath==null?new PrintWriter(new NullOuputStream()):IOUtils.openPathForPrintWriter(this.manifestPath);
				VariantContextWriter w = null;
				int nFiles = this.start_index;
				VariantContext prev = null;
				VariantContext first = null;
				long n = 0L;
				final ProgressFactory.Watcher<VariantContext> progress=ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
				for(;;)
					{
					final VariantContext ctx= in.hasNext()?progress.apply(in.next()):null;
					if(ctx==null ||
						(this.count>0 && n>=this.count) || 
						(this.consider_chromosome &&  prev!=null && !ctx.contigsMatch(prev)) ||
						(this.max_distance>0 && prev!=null && !ctx.withinDistanceOf(prev, this.max_distance)) ||
						(this.max_interval_length>0 && first!=null && !ctx.withinDistanceOf(first, this.max_interval_length))
						) {
						if(w!=null) {
							manifest.print(prev.getContig());
							manifest.print("\t");
							manifest.print(prev.getStart());
							manifest.print("\t");
							manifest.print(prev.getReference().getDisplayString());
							manifest.print("\t");
							manifest.print(n);
							manifest.println();
							w.close();
							}
						if(ctx==null) {
							break;
							}
						n=0L;
						prev=null;
						first=null;
						w=null;
						}
					if(w==null) {
						final Path outputFile = this.outputDirectory.resolve(Paths.get(String.format("%s.%06d.vcf.gz",this.prefix,nFiles)));
						LOG.info("Writing "+ outputFile);
						final VCFHeader copy = new VCFHeader(header);
						header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+".Nsplit", String.valueOf(nFiles)));
						JVarkitVersion.getInstance().addMetaData(this, copy);
						w = this.writingVariantsDelegate.dictionary(header).open(outputFile);
						w.writeHeader(copy);
						first = ctx;
						manifest.print(outputFile);	
						manifest.print("\t");
						manifest.print(first.getContig());
						manifest.print("\t");
						manifest.print(first.getStart());
						manifest.print("\t");
						manifest.print(first.getReference().getDisplayString());
						manifest.print("\t");
						prev = null;
						n=0L;
						++nFiles;
						}
					w.add(ctx);
					prev=ctx;
					n++;
					}
				
				progress.close();
				manifest.flush();
				manifest.close();
				manifest=null;
				return 0;
				}
			catch (final Throwable e)
				{
				LOG.error(e);
				return -1;
				}
			finally
				{
				CloserUtil.close(in);
				}
			}

	public static void main(final String[] args)
		{
		new Biostar497922().instanceMainWithExit(args);
		}
	}
