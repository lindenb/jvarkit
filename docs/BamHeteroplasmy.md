# BamHeteroplasmy

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Call a VCF for the Mitochondria (Experimental)


## Usage

```
Usage: bamheteroplasmy [options] Files
  Options:
    --all-alleles
      when calculating the depth, use all alleles, not just the major/minor
      Default: false
    -B, --baq
      min base quality
      Default: 0
    --discordant
      accept discordant alignments (mate mapping another contig)
      Default: false
    --fisher-treshold
      set filter if fisher test for strand bias is lower than 'x'.
      Default: 0.01
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -Q, --mapq
      min mapping quality
      Default: 30
    --mate-unmapped
      accept mate unmapped
      Default: false
    --max-clipping
      max clip per read
      Default: 2
    --min-allele-dp
      min-dp per alt allele
      Default: 20
    --organelle
      Organelle name; if empty, the program tries to find the best name in the 
      REF 
    -o, --output
      Output file. Optional . Default: stdout
    -partition, --partition
      Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -sa, --sa
      accept read having 'SA:' supplementary alignments mapping another contig
      Default: false
    --secondary
      accept secondary alignments.
      Default: false
    --supplementary
      accept supplementary alignments.
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * sam
 * bam
 * mitochondria


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bamheteroplasmy
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190910

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/mito/BamHeteroplasmy.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/mito/BamHeteroplasmy.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bamheteroplasmy** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Experimental.

Input is one or more indexed BAM/CRAM files or a file with the suffix '.list' containing the path to the BAMs.


