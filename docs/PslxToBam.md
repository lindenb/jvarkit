# PslxToBam

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert PSL to SAM/BAM


## Usage

```
Usage: psl2bam [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -b, --baq
      default base symbol
      Default: 2
    --disable-pslx
      disable use of bases provided by the pslX format.
      Default: false
    -D, --disable-secondary
      disable 'consecutive reads with same name will be flagged as secondary'
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --insert-base
      default base for insertion
      Default: N
    --intron
      use 'N' operator instead of 'D' when deletion are larger than 'x'
      Default: 50
    -q, --mapq
      default mapping quality
      Default: 255
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -g, --rg
      read group header line such as 'RG\tID:foo\tSM:bar' . '@' prefix should 
      be ignore because of this bug/feature: '@ to refer to contents in a 
      file': https://github.com/cucumber/cucumber-jvm/issues/266
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```


## Keywords

 * blat
 * sam
 * bam
 * psl
 * pslx


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew psl2bam
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190918

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/PslxToBam.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/PslxToBam.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/PslxToBamTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/PslxToBamTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **psl2bam** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## motivation

convert the psl  (or pslx) format to sam/bam.

When the input format is pslx. The program will use the bases to set the SEQ column. Nevertheless, there are sometimes fewers bases in the pslx format than expected (I don't understand why), so I sometimes fill the SEQ with 'N'.

## Example

### remapping reads with blat

```
samtools fasta src/test/resources/S1.bam |\
	blat -out=pslx  src/test/resources/rotavirus_rf.fa stdin stdout |\
	java -jar dist/psl2bam.jar -R src/test/resources/rotavirus_rf.fa 
```

### example
 
```
java -jar dist/psl2bam.jar -R src/test/resources/rotavirus_rf.fa   input.psl

@HD	VN:1.6	SO:unsorted
@SQ	SN:RF01	LN:3302	M5:59dccb944425dd61f895a564ad7b56a7	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
(...)
@SQ	SN:RF11	LN:666	M5:7a7cf2c7813f2e8bd74be383014202ca	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@PG	ID:0	CL:-R src/test/resources/rotavirus_rf.fa	PN:psl2bam	VN:a08d9c1
@CO	psl2bam. compilation:20190918074721 githash:a08d9c1 htsjdk:2.20.1 date:20190918075100. cmd:-R src/test/resources/rotavirus_rf.fa
RF01:100-200	0	RF01	100	255	101M	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:0
RF01:100-200/rc	16	RF01	100	255	101M	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:0
RF01:100-200+N	0	RF01	100	255	5H101M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:0
RF01:100-200/rc+N	16	RF01	100	255	5H101M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:0
RF01:100-200+N+mismatch	0	RF01	100	255	5H101M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:1
RF01:100-200/rc+N+mismatch	16	RF01	100	255	5H101M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:1
RF01:100-200+N+DEL	0	RF01	100	255	5H33M29D39M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTCTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:29
RF01:100-200/rc+N+DEL	16	RF01	100	255	5H33M29D39M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTCTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:29
RF01:100-200+N+INS	0	RF01	100	255	5H31M29I70M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:29
RF01:100-200/rc+N+INS	16	RF01	100	255	5H67M29I34M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:29
```

## see also

 * https://github.com/samtools/samtools/blob/develop/misc/psl2sam.pl
 * https://github.com/bsipos/uncle_psl

