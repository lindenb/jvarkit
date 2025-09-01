/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.UnaryOperator;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.variant.vcf.BcfIteratorBuilder;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIterator;


/**
BEGIN_DOC

# motivation

Split VCF to 'N' VCF files

# example

```
$ java -jar dist/vcfsplitnvariants.jar --vcf-count 10 src/test/resources/rotavirus_rf.vcf.gz -o SPLITVCF --manifest jeter.mf --index -R src/test/resources/rotavirus_rf.fa
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
		description="Split VCF to 'N' VCF files, or by number fo variant of by distance between variants ",
		creationDate = "20221122",
		modificationDate="20250901",
		biostars={9548193},
		keywords= {"vcf"},
		jvarkit_amalgamion =  true,
		menu="VCF Manipulation"
		)
public class VcfSplitNVariants 	extends Launcher {
private static final Logger LOG = Logger.of(VcfSplitNVariants.class);	

@Parameter(names={"-o","--output","--prefix"},description="files prefix",required = true)
private String outputFile = null;
@Parameter(names={"--distance"},description="distance between variants. Or use --variants-count or --vcf-count. "+DistanceParser.OPT_DESCRIPTION,splitter = NoSplitter.class,converter = DistanceParser.StringConverter.class)
private int distance_between_variants = -1;
@Parameter(names={"--vcf-count","--n-vcfs"},description="number of output vcf files. Or use --variants-count or --distance")
private int split_n_files = -1;
@Parameter(names={"--variants-count","--vc-count","--n-variants"},description="number of variants. Or use --vcf-count or --distance")
private int split_n_variants = -1;
@Parameter(names={"--tbi","--index"},description="index on the fly as tbi files")
private boolean index_on_fly = false;
@Parameter(names={"-m","--manifest"},description="write optional manifest to that file.")
private Path manifest=null;
@Parameter(names={"-f","--force"},description="overwrite existing files")
private boolean overwrite_existing_files = false;


private static class OutputVCF {
	Path path;
	VariantContextWriter w;
	long n_variants = 0L;
	final Set<String> contigs = new LinkedHashSet<>();
	void add(final VariantContext ctx) {
		this.n_variants++;
		this.contigs.add(ctx.getContig());
		this.w.add(ctx);
		}
	
	void close(final PrintWriter pw) {
		pw.print(this.path.toAbsolutePath().toString());
		pw.print("\t");
		pw.print(this.n_variants);
		pw.print("\t");
		pw.print(String.join(",", this.contigs));
		pw.println();
		w.close();
		w=null;
	}
}

private OutputVCF newOutputVCF(final VCFHeader srcHeader, final SAMSequenceDictionary dict,final int writer_id) throws IOException {
	final OutputVCF outvcf = new OutputVCF();
	final VariantContextWriterBuilder vcb = new VariantContextWriterBuilder();
	outvcf.path = Paths.get(this.outputFile + "." + String.format("%05d",(writer_id+1)) + ".vcf.gz");
	if(!overwrite_existing_files && Files.exists(outvcf.path)) {
		throw new IOException("File already exists "+outvcf.path);
		}
	LOG.info("open "+outvcf.path);
	final VCFHeader header2 = new VCFHeader(srcHeader);
	header2.addMetaDataLine(new VCFHeaderLine("vcfsplitnvariants.id", String.valueOf(writer_id+1)));
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
	int check = 0;
	if(split_n_files>0) check++;
	if(split_n_variants>0) check++;
	if(distance_between_variants>0) check++;
	if(check==0) {
		LOG.error("--vcf-count or --distance or --vcf-variants must be defined");
		return -1;
		}
	if(check!=1) {
		LOG.error("--vcf-count *OR* --vcf-variants **OR** --distance must be defined but not both");
		return -1;
		}
	try {
		final String input = super.oneFileOrNull(args);

		try(final PrintWriter pw = this.manifest==null?
				new PrintWriter(new NullOuputStream()):
				IOUtils.openPathForPrintWriter(this.manifest)) {
		pw.println("vcf\tcount\tcontigs");
		/* get number of variants */
		try(VCFIterator iter = (input==null?
				new  BcfIteratorBuilder().open(stdin()):
					new BcfIteratorBuilder().open(Paths.get(input)))) {
				final VCFHeader header = iter.getHeader();
				
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
				final UnaryOperator<VariantContext> checker;
				if(this.index_on_fly) {
					checker = new OrderChecker<VariantContext>(dict,false);
					}
				else
					{
					checker = T->T;
					}
				if(this.distance_between_variants>0) {
					VariantContext prev=null;
					OutputVCF outvcf=null;
					int writer_id=0;
					for(;;) {
						final VariantContext ctx = iter.hasNext()?checker.apply(iter.next()):null;
						if(outvcf==null || ctx==null ||
							(prev!=null && (!prev.contigsMatch(ctx) || CoordMath.getLength(prev.getEnd(), ctx.getStart()) > this.distance_between_variants)) 
							)
							{
							if(outvcf!=null) {
								outvcf.close(pw);
								outvcf=null;
								}
							if(ctx==null) break;
							outvcf = newOutputVCF(header,dict,writer_id++);
							}
						outvcf.add(ctx);
						prev=ctx;
						}
					}
				else if(this.split_n_files>0) {
					final List<OutputVCF> writers = new ArrayList<>(this.split_n_files);
					long n_variants = 0L;
					while(iter.hasNext()) {
						final VariantContext ctx = checker.apply(iter.next());
						final OutputVCF outvcf;
						
						final int writer_id = (int)(n_variants%this.split_n_files);
						if(writer_id >= writers.size()) {
							outvcf = newOutputVCF(header,dict,writer_id);
							writers.add(outvcf);
							}
						else {
							outvcf = writers.get(writer_id);	
							}
						outvcf.add(ctx);
						n_variants++;
						}
					for(OutputVCF out: writers) {
						out.close(pw);
						}
					}
				else if(this.split_n_variants>0) {
					OutputVCF outvcf=null;
					int writer_id=0;
					for(;;) {
						final VariantContext ctx = iter.hasNext()?checker.apply(iter.next()):null;
						if(outvcf==null || ctx==null ||  outvcf.n_variants>=this.split_n_variants) {
							if(outvcf!=null) {
								outvcf.close(pw);
								outvcf=null;
								}
							if(ctx==null) break;
							outvcf = newOutputVCF(header,dict,writer_id++);
							}
						outvcf.add(ctx);
						}
					}
				else
					{
					throw new IllegalStateException();
					}
				}
			pw.flush();
			} /* end write pw */
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	}
public static void main(final String[] args) {
	new VcfSplitNVariants().instanceMainWithExit(args);
	}
}
