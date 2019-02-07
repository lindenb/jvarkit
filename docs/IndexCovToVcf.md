# IndexCovToVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert indexcov data to vcf


## Usage

```
Usage: indexcov2vcf [options] Files
  Options:
    -del, --deletion
      Deletion treshold
      Default: 0.6
    -dup, --duplication
      Duplication treshold
      Default: 1.9
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -ped, --pedigree
      Optional pedigree. A pedigree is a text file delimited with tabs. No 
      header. Columns are (1) Family (2) Individual-ID (3) Father Id or '0' 
      (4) Mother Id or '0' (5) Sex : 1 male/2 female / 0 unknown (6) Status : 
      0 unaffected, 1 affected,-9 unknown
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```


## Keywords

 * cnv
 * jfx
 * duplication
 * deletion
 * sv


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew indexcov2vcf
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/IndexCovToVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/IndexCovToVcf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/IndexCovToVcfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/IndexCovToVcfTest.java)


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

## history

  * 20191112 : add pedigree

