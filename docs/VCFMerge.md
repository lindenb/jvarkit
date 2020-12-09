# VCFMerge

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Merge a large number of VCF Files


## Usage

```
Usage: vcfmerge [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --fields
      print the following INFO/FORMAT Fields.
      Default: AC,AN,AF,DP,GQ,AD,PL
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -homref, --homref, -hr
      Use HomRef 0/0 for unknown variant
      Default: false
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    --ploidy
      Ploidy
      Default: 2
    -region, --region, -r
      Merge in that region: An interval as the following syntax : 
      "chrom:start-end" or "chrom:middle+extend"  or "chrom:start-end+extend" 
      or "chrom:start-end+extend-percent%".A program might use a Reference 
      sequence to fix the chromosome name (e.g: 1->chr1)
      Default: <empty string>
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * sort
 * merge


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfmerge
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20130916

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMerge.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMerge.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMergeTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfmerge/VCFMergeTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfmerge** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


The motivation for this is to merge a large number of VCF files without opening a bunch of temporary files.

For a regular normal number of files you should use  GATK combineVariants or bcftools merge
 
## Example


```bash
$ java -jar dist/vcfmerge.jar -hr src/test/resources/S*.vcf.gz | more
[INFO][VCFMerge]merging...5 vcfs
##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF01	970	.	A	C	.	.	AC=2;AF=0.200;AN=10	GT	0/0	0/0	0/0	0/0	1/1
RF01	3246	.	A	G	.	.	AC=2;AF=0.200;AN=10	GT	0/0	0/0	0/0	0/0	1/1
RF02	578	.	G	A	.	.	AC=2;AF=0.200;AN=10	GT	0/0	0/0	0/0	1/1	0/0
RF02	877	.	T	A	.	.	AC=1;AF=0.100;AN=10	GT	0/1	0/0	0/0	0/0	0/0
RF02	1962	.	TACA	TA	.	.	AC=1;AF=0.100;AN=10	GT	0/1	0/0	0/0	0/0	0/0
RF02	2332	.	AT	A	.	.	AC=1;AF=0.100;AN=10	GT	0/0	0/0	0/0	0/1	0/0
RF02	2662	.	G	C	.	.	AC=2;AF=0.200;AN=10	GT	0/0	0/0	0/0	0/0	1/1

```

