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
package com.github.lindenb.jvarkit.tools.dict2R;


import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;

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
@Program(name="dict2r",
	description="convert a SAM dictionary from vcf,sam,bam,dict, etc.. to R code functions.",
	keywords={"dict","bed","sam","bam","vcf","R"},
	creationDate="20241212",
	modificationDate="20241212",
	jvarkit_amalgamion =  true,
	jvarkit_hidden = true,
	generate_doc = false
	)
public class DictToR extends Launcher {
	private static Logger LOG=Logger.of(DictToR.class);

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"--prefix"},description="prefix name. If empty try to use build name, if it is available")
	private String prefix="";
	@Parameter(names={"--regex"},description="keep chromosomes matching that regular expression")
	private String contig_regex="(chr)?[0-9XY]+";	
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Path path = Paths.get(oneAndOnlyOneFile(args));
			final SAMSequenceDictionary dict0 = SequenceDictionaryUtils.extractRequired(path);
			if(StringUtils.isBlank(prefix)) {
				prefix = SequenceDictionaryUtils.getBuildName(dict0).orElse(null);
				}
			if(StringUtils.isBlank(prefix)) {
				LOG.error("prefix is empty");
				return -1;
				}
			final String dictName = prefix+".dict";
			
			final SAMSequenceDictionary dict = new SAMSequenceDictionary(
					dict0.getSequences().stream().
					filter(SSR->SSR.getContig().matches(contig_regex)).
					collect(Collectors.toList())
					);
			
			
			try (PrintWriter pw = super.openPathOrStdoutAsPrintWriter(outputFile)) {
				
				
				pw.append(prefix).append(".genomeLength <- "+dict.getSequences().stream().mapToLong(R->R.getLengthOnReference()).sum());
				
				pw.append(dictName).append(" <- ");
				pw.println("data.frame(");
				pw.println("    name=c("+dict.getSequences().stream().map(R->StringUtils.doubleQuote(R.getContig())).collect(Collectors.joining(","))+"),");
				pw.println("    tid=c("+dict.getSequences().stream().map(R->String.valueOf(R.getSequenceIndex())).collect(Collectors.joining(","))+"),");
				pw.println("    length=c("+dict.getSequences().stream().map(R->String.valueOf(R.getLengthOnReference())).collect(Collectors.joining(","))+"),");
				pw.print(  "    genomicIndex=c(");
				
				long x=0;
				for(int i=0;i< dict.size();i++) {
					if(i>0) pw.print(",");
					pw.print(String.valueOf(x));
					x+= dict.getSequence(i).getLengthOnReference();
					}
				pw.print(")\n)");
			
				pw.println("row.names("+prefix+".dict)<- "+dictName+"$name;");
				
				pw.append(prefix).println(".pos2idx <- function(chrom,pos,width) {");
				pw.println("    if(!any(row.names("+dictName+") == chrom))) { return -1;}");
				pw.println("    return ("+dictName+"[chrom]$genomicIndex+pos);");
				pw.println("}");
				
				pw.append(prefix).println(".pos2x <- function(chrom,pos,width) {");
				pw.println("    if(!any(row.names("+prefix+".dict) == chrom))) { return -1;}");
				pw.println("    idx <- "+prefix+".pos2idx(chrom,pos);");
				pw.println("    if(idx<0) return -1;");
				pw.println("    (idx/(1.0 * "+prefix+".genomeLength))*width;");
				pw.println("}");
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
		new DictToR().instanceMainWithExit(args);
		}
}
