# IndexCovToVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert indexcov data to vcf


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar indexcov2vcf  [options] Files

Usage: indexcov2vcf [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --cases
      File or comma-separated list of control samples
    --controls
      File or comma-separated list of control samples
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    --no-merge
      disable adjacent block merging for the same variant (keep the original 
      bed structure)
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    -t, --treshold
      DUP if 1.5-x<=depth<=1.5+x . HET_DEL if 0.5-x<=depth<=0.5+x HOM_DEL if 
      0.0-x<=depth<=0.0+x.. A decimal number between 0.0 and 1.0. If the value 
      ends with '%' it is interpretted as a percentage eg. '1%' => '0.01'. A 
      slash '/' is interpretted as a ratio. e.g: '1/100' => '0.01'.
      Default: 0.05
    --version
      print version and exit

```


## Keywords

 * cnv
 * duplication
 * deletion
 * sv
 * indexcov



## Creation Date

20200528

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/IndexCovToVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/IndexCovToVcf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/IndexCovToVcfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/IndexCovToVcfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **indexcov2vcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a tab-delimited file created by e.g: indexcov (https://github.com/brentp/goleft/tree/master/indexcov)

```
#chrom  start  end     SampleBB  SampleBC  SampleBD  SampleBE  SampleBF  SampleBG  SampleBH
chr1    23778  40778   1.59      1.31      1.67      1.61      1.83      1.52      1.48
chr1    29106  46106   1.9       1.54      1.72      1.97      1.88      1.53      1.95
chr1    84581  101581  0.764     0.841     1.2       1.16      1.18      1.13      1.23
chr1    15220  32220   0.355     0.704     1.09      0.784     0.81      1.37      0.954
chr1    58553  75553   0.353     0.436     0.912     0.836     1.16      1.09      0.611
chr1    19347  36347   0.381     0.411     0.811     0.795     1.16      1.22      0.495
chr1    81062  98062   1.09      0.972     1.35      1.22      1.66      1.76      1.1
chr1    17353  34353   1.06      1.06      1.23      1.26      1.44      1.43      1.03
chr1    48498  65498   1.08      0.996     1.28      1.44      1.52      1.57      1.05
```

output:

```
##fileformat=VCFv4.2
##FILTER=<ID=ALL_DEL,Description="number of samples >1 and all are deletions">
##FILTER=<ID=ALL_DUP,Description="number of samples >1 and all are duplication">
##FILTER=<ID=NO_SV,Description="There is no DUP or DEL in this variant">
##FORMAT=<ID=DEL,Number=1,Type=Integer,Description="set to 1 if relative number of copy <= 0.6">
##FORMAT=<ID=DUP,Number=1,Type=Integer,Description="set to 1 if relative number of copy >= 1.9">
##FORMAT=<ID=F,Number=1,Type=Float,Description="Relative number of copy: 0.5 deletion 1 normal 2.0 duplication">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=NDEL,Number=1,Type=Integer,Description="Number of samples being deleted">
##INFO=<ID=NDUP,Number=1,Type=Integer,Description="Number of samples being duplicated">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SampleBB	SampleBC	SampleBD	SampleBE	SampleBF	(...)
chr1	0	.	N	<DUP>	.	.	END=16384;NDEL=0;NDUP=8	GT:DUP:F	0:0:1.59	0:0:1.31	0:0:1.67	0:0:1.61	0:0:1.83 (...)
```


