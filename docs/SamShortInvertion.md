# SamShortInvertion

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Scan short inversions in SAM using supplementary reads.


## Usage

```
Usage: samshortinvert [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -B, --bed, -r, --rgn
      A source of intervals. The following suffixes are recognized: vcf, 
      vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise it could be an 
      empty string (no interval) or a list of plain interval separated by '[ 
      \t\n;,]' 
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      min mapping quality
      Default: 1
    -m, --maxsize
      max size of inversion.
      Default: 10000
    -o, --output
      Output file. Optional . Default: stdout
    -partition, --partition
      Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -F, --ratio
      Two intervals are the same if they both have more or equals of this 
      fraction of length in common. A decimal number between 0.0 and 1.0. If 
      the value ends with '%' it is interpretted as a percentage eg. '1%' => 
      '0.01'. A slash '/' is interpretted as a ratio. e.g: '1/100' => '0.01'.
      Default: 0.75
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit
    -s, -supporting
      Don't print the variant if INFO/DP <= 's'
      Default: 1

```


## Keywords

 * sam
 * bam
 * sv
 * inversion


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samshortinvert
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20140228

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamShortInvertion.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamShortInvertion.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samshortinvert** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

finds the regions having some short inversions.

input is a set of BAM files. One file ending with '.list' is interpreted as a file containing some path to the bams.

output is a VCF file

## Example:

```
$ find DIR -type f -name "*.bam" > bams.list
$ java -jar ${JVARKIT_DIST}/samshortinvert.jar -R ref.fasta bams.list |\
 	bcftools view -i 'INFO/DPMAX>10' > out.vcf
```

## Screenshot

* https://twitter.com/yokofakun/status/1222848286048112641
![https://twitter.com/yokofakun/status/1222848286048112641](https://pbs.twimg.com/media/EPhuCJnX4AA3Brc?format=png&name=medium)

* https://twitter.com/yokofakun/status/1222832425518141442
![https://twitter.com/yokofakun/status/1222832425518141442](https://pbs.twimg.com/media/EPhfm8EW4AAiaBq?format=png&name=medium)

* https://twitter.com/yokofakun/status/1222853635656364032
![https://twitter.com/yokofakun/status/1222853635656364032](https://pbs.twimg.com/media/EPhy5fAUYAAMunf?format=png&name=medium)

