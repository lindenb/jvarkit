# GatkHaplotypeCaller

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Wrapper for GATK HaplotypeCaller


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar gatkhc  [options] Files

Usage: gatkhc [options] Files
  Options:
  * -L, -bed, --bed
      restrict to bed
    -dbsnp, --dbsnp
      path to dbsnp
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      mapping quality
      Default: 10
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --references
      Other references. If a reference is different from the main reference, 
      the contigs of a  GVCF file will be converted (e.g: 1 -> chr1) to the 
      main reference dictionary.
      Default: []
    --tmp, --tmp-dir
      temporary directory
      Default: /tmp
    --version
      print version and exit

```


## Keywords

 * gatk
 * vcf
 * bam



## Creation Date

20240625

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gatk/GatkHaplotypeCaller.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gatk/GatkHaplotypeCaller.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gatkhc** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


This tool is a wrapper over gatk4 HaplotypeCaller

Not All versions of jvarkit will work for this tool. Because it needs to
be compiled against the java code from gatk using a command like:

```
./gradlew jvarkit -Dgatk4.local.jar=/path/to/gatk/gatk-package-4.*-local.jar
```

Furthermore, don't trust the automatic generated documentation,  the tool must be invoked using the following syntax:

```
java -cp /path/to/jvarkit.jar:/path/to/gatk-package-4.xxx-local.jar \
	com.github.lindenb.jvarkit.tools.gatk.GatkHaplotypeCaller \
	-R src/test/resources/rotavirus_rf.fa \
	--bed input.bed \
	src/test/resources/S*.bam
```

This tool:

* call each BAM into a .vcf.gz using haplotype caller
  * convert the dictionary/chromosomes if needed
  * convert the sample name if the same sample if present more than once
* group g.vcf.gz files by the sqrt(number-of-gvcf-files)
* combine each group of gvcf files
* genotypegvcfs for the final gvcf file


