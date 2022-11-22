/*
The MIT License (MIT)

Copyright (c) 2022 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfsplit;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.UnaryOperator;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;


/**
BEGIN_DOC

# motivation

Split VCF to 'N' VCF files

# example

```
$ java -jar dist/vcfsplitnvariants.jar -N 10 src/test/resources/rotavirus_rf.vcf.gz -o SPLITVCF --manifest jeter.mf --index -R src/test/resources/rotavirus_rf.fa
[INFO][VcfSplitNVariants]open SPLITVCF.00001.vcf.gz
[INFO][VcfSplitNVariants]open SPLITVCF.00002.vcf.gz
[INFO][VcfSplitNVariants]open SPLITVCF.00003.vcf.gz
[INFO][VcfSplitNVariants]open SPLITVCF.00004.vcf.gz
[INFO][VcfSplitNVariants]open SPLITVCF.00005.vcf.gz
[INFO][VcfSplitNVariants]open SPLITVCF.00006.vcf.gz
[INFO][VcfSplitNVariants]open SPLITVCF.00007.vcf.gz
[INFO][VcfSplitNVariants]open SPLITVCF.00008.vcf.gz
[INFO][VcfSplitNVariants]open SPLITVCF.00009.vcf.gz
[INFO][VcfSplitNVariants]open SPLITVCF.00010.vcf.gz

$ cat jeter.mf
vcf	count	contigs
/home/lindenb/src/jvarkit-git/SPLITVCF.00001.vcf.gz	5	RF01,RF03,RF04,RF06,RF09
/home/lindenb/src/jvarkit-git/SPLITVCF.00002.vcf.gz	5	RF02,RF03,RF05,RF06,RF10
/home/lindenb/src/jvarkit-git/SPLITVCF.00003.vcf.gz	5	RF02,RF03,RF05,RF07,RF10
/home/lindenb/src/jvarkit-git/SPLITVCF.00004.vcf.gz	5	RF02,RF03,RF05,RF07,RF10
/home/lindenb/src/jvarkit-git/SPLITVCF.00005.vcf.gz	5	RF02,RF04,RF05,RF07,RF11
/home/lindenb/src/jvarkit-git/SPLITVCF.00006.vcf.gz	4	RF02,RF04,RF05,RF07
/home/lindenb/src/jvarkit-git/SPLITVCF.00007.vcf.gz	4	RF03,RF04,RF05,RF08
/home/lindenb/src/jvarkit-git/SPLITVCF.00008.vcf.gz	4	RF03,RF04,RF06,RF08
/home/lindenb/src/jvarkit-git/SPLITVCF.00009.vcf.gz	4	RF03,RF04,RF06,RF09
/home/lindenb/src/jvarkit-git/SPLITVCF.00010.vcf.gz	4	RF03,RF04,RF06,RF09
```


END_DOC
*/

@Program(
		name="vcfsplitnvariants",
		description="Split VCF to 'N' VCF files ",
		creationDate = "202221122",
		modificationDate="202221122",
		keywords= {"vcf"}
		)
public class VcfSplitNVariants 	extends Launcher {
private static final Logger LOG = Logger.build(VcfSplitNVariants.class).make();	

@Parameter(names={"-o","--output","--prefix"},description="files prefix",required = true)
private String outputFile = null;
@Parameter(names={"--vcf-count"},description="number of vcf files. Or use --variant-count")
private int split_n_files = -1;
@Parameter(names={"--vcf-variants"},description="number of vcf files. Or use --variant-count")
private int split_n_variants = -1;
@Parameter(names={"--tbi","--index"},description="index on the fly as tbi files")
private boolean index_on_fly = false;
@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
private Path faidx=null;
@Parameter(names={"-m","--manifest"},description="write optional manifest")
private Path manifest=null;


private static class OutputVCF {
	Path path;
	VariantContextWriter w;
	long n_variants = 0L;
	final Set<String> contigs = new LinkedHashSet<>();
}

private OutputVCF newOutputVCF(final VCFHeader srcHeader, SAMSequenceDictionary dict,final int writer_id) {
	OutputVCF outvcf = new OutputVCF();
	final VariantContextWriterBuilder vcb = new VariantContextWriterBuilder();
	outvcf.path = Paths.get(this.outputFile + "." + String.format("%05d",(writer_id+1)) + ".vcf.gz");
	final VCFHeader header2 = new VCFHeader(srcHeader);
	JVarkitVersion.getInstance().addMetaData(this, header2);
	
	if(dict!=null) {
		vcb.setReferenceDictionary(dict);
		}
	if(index_on_fly) {
		final TabixIndexCreator tabixIndexCreator = new TabixIndexCreator(dict, TabixFormat.VCF);
		vcb.setIndexCreator(tabixIndexCreator);
		vcb.setOption(Options.INDEX_ON_THE_FLY);
		}
	
	outvcf.w = vcb. setOutputPath(outvcf.path).
			setCreateMD5(false).
			build();
	outvcf.w.writeHeader(header2);
	return outvcf;
	}

@Override
public int doWork(final List<String> args) {
	if(split_n_files<=0 && split_n_variants<=0) {
		LOG.error("--vcf-count or --vcf-variants must be defined");
		return -1;
		}
	if(split_n_files>0 && split_n_variants>0) {
		LOG.error("--vcf-count *OR* --vcf-variants must be defined but not both");
		return -1;
		}
	try {
		Path reference = null;
		int N_FILES = 10;
		
		final SAMSequenceDictionary dict;
		if(this.index_on_fly) {
			if(faidx==null) {
				LOG.error("reference is required if you need indexing.");
				return -1;
				}
			dict = SequenceDictionaryUtils.extractRequired(faidx);
			if(dict==null) {
				LOG.error("cannot load dict for " + reference);
				return -1;
				}
			}
		else
			{
			dict = null;
			}
		final String input = super.oneFileOrNull(args);
		/* get number of variants */
		try(VCFIterator iter = (input==null?
				new VCFIteratorBuilder().open(stdin()):
				new VCFIteratorBuilder().open(Paths.get(input)))) {
			final VCFHeader header = iter.getHeader();
			final List<OutputVCF> writers = new ArrayList<>(N_FILES);
			
			final UnaryOperator<VariantContext> checker;
			if(this.index_on_fly) {
				checker = new OrderChecker<VariantContext>(dict,false);
				}
			else
				{
				checker = T->T;
				}
			long n_variants = 0L;
			while(iter.hasNext()) {
				final VariantContext ctx = checker.apply(iter.next());
				final OutputVCF outvcf;
				
				if(this.split_n_files>0) {
					final int writer_id = (int)(n_variants%this.split_n_files);
					if(writer_id >= writers.size()) {
						outvcf =newOutputVCF(header,dict,writer_id);
						writers.add(outvcf);
						LOG.info("open "+outvcf.path);
						}
					else {
						outvcf = writers.get(writer_id);	
						}
					}
				else if(this.split_n_variants>0) {
					if(writers.isEmpty() || writers.)
					}
				else
					{
					throw new IllegalStateException();
					}
				outvcf.contigs.add(ctx.getContig());
				outvcf.w.add(ctx);
				outvcf.n_variants++;
				n_variants++;
				if(this.split_n_variants>0) {
					
					}
				}
				

			for(OutputVCF out: writers) {
				out.w.close();
				}
			if(manifest!=null) {
				try(final PrintWriter pw = IOUtils.openPathForPrintWriter(this.manifest)) {
					pw.println("vcf\tcount\tcontigs");
					for(OutputVCF x:writers) {
						pw.print(x.path.toAbsolutePath().toString());
						pw.print("\t");
						pw.print(x.n_variants);
						pw.print("\t");
						pw.print(String.join(",", x.contigs));
						pw.println();
						}
					pw.flush();
					}
				}
			}//end iter
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	}
public static void main(String[] args) {
	new VcfSplitNVariants().instanceMainWithExit(args);
	}
}
