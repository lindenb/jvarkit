/*
 
BEGIN_DOC

## Motivation

Extract the first and last variants of a tabix-index-VCF file for each chromosome.
Output is a BED file contig/start/end/vcf.
Input can be:

  - stdin (the path to the vcfs, one per line)
  - the vcfs
  - a file with the suffix '.list' containing the path to the vcfs, one per line.

## Example:


END_DOC

*/
package com.github.lindenb.jvarkit.tools.tbi2bed;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.BinningIndexContent;
import htsjdk.samtools.LinearIndex;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

/**
# Motivation

find the first and last position of each chromosome from a set of tabix-indexed VCF files.

Input is either:

  - stdin . One path to the VCF per line
  - the path to the VCFs
  - a file with the '.list' suffix containing the path to the VCFs, one per line

# warning

the END position of the last variant of the contig is not always the MAX(END)  for this contig. 
For example if the penultimate variant is a very large variant.

# Example
 
```
$ find . -type f -name "S[12].vcf.gz" | java -jar dist/jvarkit.jar  vcftbi2bed

RF02	876	1965	./src/test/resources/S1.vcf.gz
RF03	2149	2150	./src/test/resources/S1.vcf.gz
RF04	886	1920	./src/test/resources/S1.vcf.gz
RF06	516	517	./src/test/resources/S1.vcf.gz
RF09	293	294	./src/test/resources/S1.vcf.gz
RF10	45	46	./src/test/resources/S1.vcf.gz
RF11	655	656	./src/test/resources/S1.vcf.gz
RF03	1220	2573	./src/test/resources/S2.vcf.gz
RF05	498	1297	./src/test/resources/S2.vcf.gz
RF06	542	543	./src/test/resources/S2.vcf.gz
RF07	224	225	./src/test/resources/S2.vcf.gz
RF08	991	992	./src/test/resources/S2.vcf.gz
```

```
$ java -jar dist/jvarkit.jar  vcftbi2bed src/test/vcftbi2bed/S[45].vcf.gz

RF02	577	2333	src/test/resources/S4.vcf.gz
RF03	1241	1242	src/test/resources/S4.vcf.gz
RF04	990	1262	src/test/resources/S4.vcf.gz
RF05	40	795	src/test/resources/S4.vcf.gz
RF06	694	1132	src/test/resources/S4.vcf.gz
RF07	683	684	src/test/resources/S4.vcf.gz
RF08	1013	1014	src/test/resources/S4.vcf.gz
RF10	138	175	src/test/resources/S4.vcf.gz
RF11	73	79	src/test/resources/S4.vcf.gz
RF01	969	3246	src/test/resources/S5.vcf.gz
RF02	2661	2662	src/test/resources/S5.vcf.gz
RF03	1687	1708	src/test/resources/S5.vcf.gz
RF04	1899	1900	src/test/resources/S5.vcf.gz
RF05	328	329	src/test/resources/S5.vcf.gz
RF06	667	668	src/test/resources/S5.vcf.gz
RF07	951	952	src/test/resources/S5.vcf.gz
RF08	925	926	src/test/resources/S5.vcf.gz
RF09	316	317	src/test/resources/S5.vcf.gz
```

*/

@Program(name="vcftbi2bed",
description="extracts BED for each contig in a tabix-indexed VCF peeking first of last variant for each chromosome.",
keywords={"bed","vcf","tabix"},
creationDate="20230214",
modificationDate="20230214",
jvarkit_amalgamion = true
)
public class VcfTbiToBed extends Launcher {
	private static final Logger LOG = Logger.of(VcfTbiToBed.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--POS","-POS"},description="for the end position, default is to use the END position (think about SV, INDEL...) of the variant. This option just use the POS.")
	private boolean use_POS = false;

	
	private void oneVcf(final PrintWriter out,final String vcf) throws IOException {
		LOG.info(vcf);
		final TabixIndex tbi = new TabixIndex(Paths.get( vcf + FileExtensions.TABIX_INDEX));
		
		final List<String> contigs = tbi.getSequenceNames();
		final BinningIndexContent[] binIndexContents = tbi.getIndices();
		try(final SeekableStream seekable = SeekableStreamFactory.getInstance().getStreamFor(vcf)) {
			try(BlockCompressedInputStream bci = new BlockCompressedInputStream(seekable)) {
				final VCFCodec vcfCodec = new VCFCodec();
				LineIterator li = vcfCodec.makeSourceFromStream(bci);
				vcfCodec.readActualHeader(li);
			
				for(int tid = 0; tid < contigs.size();tid++) {
					final String contig = contigs.get(tid);
					final LinearIndex linearIndex = binIndexContents[tid].getLinearIndex();
				
					long offset=linearIndex.get(0);
					//System.err.println("offset 0 is "+offset);
					bci.seek(offset);
					li = vcfCodec.makeSourceFromStream(bci);
					if(!li.hasNext()) continue;
					final VariantContext ctx0 = vcfCodec.decode(li);
					if(ctx0==null) {
						//System.err.println("No variant for "+contig);
						continue;
						}
					if(!ctx0.getContig().equals(contig)) throw new IllegalStateException("expected contig "+ contig + " but got "+ctx0.getContig());
					
					offset=linearIndex.get(linearIndex.size()-1);
					//System.err.println("offset 1 is "+offset);
					bci.seek(offset);
					li = vcfCodec.makeSourceFromStream(bci);
					VariantContext ctx1 = ctx0;
					while(li.hasNext()) {
						final VariantContext ctx2 = vcfCodec.decode(li);
						if(ctx2==null) break;
						if(!ctx2.getContig().equals(contig)) break;
						ctx1 = ctx2;
						}
					out.print(contig);
					out.print("\t");
					out.print(ctx0.getStart()-1);
					out.print("\t");
					out.print(Math.max(
							(this.use_POS?ctx0.getStart():ctx0.getEnd()),
							(this.use_POS?ctx1.getStart():ctx1.getEnd())
							));
					out.print("\t");
					out.print(vcf);
					out.println();
					}
				}
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try  {
			final List<Path> vcfs = IOUtils.unrollPaths(args);
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				if(vcfs.isEmpty()) {
					try(BufferedReader br = IOUtils.openStreamForBufferedReader(stdin())) {
						String line;
						while((line=br.readLine())!=null) {
							if(line.startsWith("#")) continue;
							if(StringUtils.isBlank(line)) continue;
							oneVcf(out,line);
							}
						}
					}
				else
					{
					for(final Path p: vcfs) {
						oneVcf(out,p.toString());
						}
					}
				out.flush();
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return 1;
			}
		}
	
	public static void main(final String[] args) {
		new VcfTbiToBed().instanceMainWithExit(args);
	}

}
