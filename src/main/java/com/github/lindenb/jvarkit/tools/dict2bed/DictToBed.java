/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.dict2bed;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
BEGIN_DOC

## Motivation

extract SAM Sequence dictionaries from SAM/BAM/FASTA/VCF files and convert them to bed

## Example

```
$ java -jar dist/jvarkit.jar dict2bed ~/src/jvarkit-git/src/test/resources/*.bam |  head | column -t
chrom  start  end        path                                                                         tid  buildName  AS  M5  SP  UR
chr1   0      248956422  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  0    GRCh38     .   .   .   .
chr2   0      242193529  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  1    GRCh38     .   .   .   .
chr3   0      198295559  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  2    GRCh38     .   .   .   .
chr4   0      190214555  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  3    GRCh38     .   .   .   .
chr5   0      181538259  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  4    GRCh38     .   .   .   .
chr6   0      170805979  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  5    GRCh38     .   .   .   .
chr7   0      159345973  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  6    GRCh38     .   .   .   .
chr8   0      145138636  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  7    GRCh38     .   .   .   .
chr9   0      138394717  /home/lindenb/src/jvarkit-git/src/test/resources/ENCFF331CGL.rnaseq.b38.bam  8    GRCh38     .   .   .   .
```

```
$ find src/test/resources/ -type f -name "*.vcf.gz" | java -jar dist/jvarkit.jar dict2bed | head | column -t
chrom  start  end        path                                                     tid  buildName
1      0      249250621  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  0    GRCh37
2      0      243199373  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  1    GRCh37
3      0      198022430  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  2    GRCh37
4      0      191154276  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  3    GRCh37
5      0      180915260  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  4    GRCh37
6      0      171115067  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  5    GRCh37
7      0      159138663  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  6    GRCh37
8      0      146364022  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  7    GRCh37
9      0      141213431  src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz  8    GRCh37
```

END_DOC
 */
@Program(name="dict2bed",
	description="convert a SAM dictionary from vcf,sam,bam,dict, etc.. to bed.",
	keywords={"dict","bed","sam","bam","vcf"},
	creationDate="20240603",
	modificationDate="20240603",
	nfcore = "https://nf-co.re/modules/jvarkit_dict2bed.html",
	jvarkit_amalgamion =  true
	)
public class DictToBed extends Launcher {
	private static Logger LOG=Logger.of(DictToBed.class);

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"--no-header"},description="disable print header")
	private boolean disable_header = false;
	@Parameter(names={"--skip-attributes"},description="skip attributes output")
	private boolean skip_attributes = false;
	@Parameter(names={"--ignore-errors"},description="ignore errors, skip files that don't have a dictionary")
	private boolean skip_error = false;

	
	private SAMSequenceDictionary extract(final Path path) {
		if(!skip_error) {
			return SequenceDictionaryUtils.extractRequired(path);
			}
		else
			{
			return SequenceDictionaryUtils.extractDictionary(path).orElse(null);
			}
	}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Set<String> attributes = new TreeSet<>();
			final List<Path> paths;
			
			if(args.isEmpty()) {
				try(BufferedReader br = IOUtils.openStdinForBufferedReader()) {
					paths = br.lines().map(S->Paths.get(S)).collect(Collectors.toList());
					}
				}
			else
				{
				paths = IOUtils.unrollPaths(args);
				}
			
			if(!skip_attributes) {
				for(Path path: paths) {
					final SAMSequenceDictionary dict = extract(path);
					if(dict==null) continue;
					for(SAMSequenceRecord ssr: dict.getSequences()) {
						ssr.getAttributes().stream().forEach(KV->attributes.add(KV.getKey()));
						}
					}
				}
			try (PrintWriter pw = super.openPathOrStdoutAsPrintWriter(outputFile)) {
				if(!disable_header) {
					pw.print("chrom");
					pw.print("\t");
					pw.print("start");
					pw.print("\t");
					pw.print("end");
					pw.print("\t");
					pw.print("path");
					pw.print("\t");
					pw.print("tid");
					pw.print("\t");
					pw.print("buildName");
					for(final String key:attributes) {
						pw.print("\t");
						pw.print(key);
						}
					pw.println();
					}
				
				for(Path path: paths) {
					final SAMSequenceDictionary dict = extract(path);
					if(dict==null) continue;
					for(SAMSequenceRecord ssr: dict.getSequences()) {
						pw.print(ssr.getContig());
						pw.print("\t0\t");
						pw.print(ssr.getSequenceLength());
						pw.print("\t");
						pw.print(path.toString());
						pw.print("\t");
						pw.print(ssr.getSequenceIndex());
						pw.print("\t");
						pw.print(SequenceDictionaryUtils.getBuildName(dict).orElse("."));
						for(final String key:attributes) {
							pw.print("\t");
							pw.print(StringUtils.ifBlank(ssr.getAttribute(key),"."));
							}
						pw.println();
						}
					}
				pw.flush();
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new DictToBed().instanceMainWithExit(args);
	}

}
