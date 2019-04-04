# Biostar352930

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Fills the empty SEQ(*) and QUAL(*) in a bam file using the the reads with the same name carrying this information.


## Usage

```
Usage: biostar352930 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam



## See also in Biostars

 * [https://www.biostars.org/p/352930](https://www.biostars.org/p/352930)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar352930
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar352930.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar352930.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar352930Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar352930Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar352930** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Motivation

Fills the empty SEQ(*) and QUAL(*) in a bam file using the the reads with the same name carrying this information.

## Input

a sam or a bam file. Input **MUST** be sorted on query name using picard. ( see see https://github.com/samtools/hts-specs/issues/5 )

## Example

A read `A00223:8:H7YG3DMXX:1:1101:1163:35383` in the remote bam file below is missing SEQ/QUAL while another has this information

```
$ wget -q -O - "https://gist.githubusercontent.com/toddknutson/90430a0dd736898037ed18bcd044df7f/raw/87c1ea5a548ac71c628d2b72f5bd6ee6415efbcd/gistfile1.txt" | grep -F "A00223:8:H7YG3DMXX:1:1101:1163:35383"
A00223:8:H7YG3DMXX:1:1101:1163:35383	16	chr7	6027025	9	151M	*	0	0	GCTAGAAGACAGCAGACCCCTTGTCTGTCCTAGAGGGCTCCTTCTTGGTTCTGGAGTCTTTGGGCTGTGAGGCTTGTTCTCTGTTGTGTGACGAAGAGAAAAGGCCTCTCGCAGTCTGGAAATGGACACGTCTTTTTTTTCTTCTCCAGTC	FFFFFFFFFFFFFF,:FFFF,F,FFFFFFFF:FFFFFFFFF:FF::FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFF	NM:i:2	MD:Z:14T7T128	AS:i:2978	XS:i:2894	RG:Z:INV_4SMALL_Aug_17_2017_1
A00223:8:H7YG3DMXX:1:1101:1163:35383	256	chr7	6776783	0	151M	*	0	0	*	*	NM:i:6	MD:Z:17G0G109A7A2T0C10	AS:i:2894	RG:Z:INV_4SMALL_Aug_17_2017_1

```

get the sam/bam file, sort it on queryname using picard

```
$ wget -q -O - "https://gist.githubusercontent.com/toddknutson/90430a0dd736898037ed18bcd044df7f/raw/87c1ea5a548ac71c628d2b72f5bd6ee6415efbcd/gistfile1.txt" |\
 	java -jar /path/to/picard.jar SortSam I=/dev/stdin O=/dev/stdout SO=queryname VALIDATION_STRINGENCY=LENIENT |\
 	java -jar dist/biostar352930.jar
```

output:

```
@HD	VN:1.5	SO:queryname
(...)
@RG	ID:INV_4SMALL_Aug_17_2017_1	SM:0X	PL:ILLUMINA	LB:libeX
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -t 12 -R @RG\tID:sample1\tSM:0X\tPL:ILLUMINA\tLB:libeX -Y -a -k 11 -O 2 -B 1 -A 20 hg19_canonical_PARmaskedonchrY.fasta R1.fastq
@CO	biostar352930. compilation: 2018-12-05:10-12-50 githash: 03f79caf6a62458937b2d55d0661737617334733 htsjdk: 2.15.0. cmd:
A00223:8:H7YG3DMXX:1:1101:1163:35383	256	chr7	6776783	0	151M	*	0	0	GACTGGAGAAGAAAAAAAAGACGTGTCCATTTCCAGACTGCGAGAGGCCTTTTCTCTTCGTCACACAACAGAGAACAAGCCTCACAGCCCAAAGACTCCAGAACCAAGAAGGAGCCCTCTAGGACAGACAAGGGGTCTGCTGTCTTCTAGC	FFFFFF,FF:FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF::FF:FFFFFFFFF:FFFFFFFF,F,FFFF:,FFFFFFFFFFFFFF	MD:Z:17G0G109A7A2T0C10	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:6	AS:i:2894
A00223:8:H7YG3DMXX:1:1101:1163:35383	16	chr7	6027025	9	151M	*	0	0	GCTAGAAGACAGCAGACCCCTTGTCTGTCCTAGAGGGCTCCTTCTTGGTTCTGGAGTCTTTGGGCTGTGAGGCTTGTTCTCTGTTGTGTGACGAAGAGAAAAGGCCTCTCGCAGTCTGGAAATGGACACGTCTTTTTTTTCTTCTCCAGTC	FFFFFFFFFFFFFF,:FFFF,F,FFFFFFFF:FFFFFFFFF:FF::FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFF	MD:Z:14T7T128	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:2	AS:i:2978	XS:i:2894
A00223:8:H7YG3DMXX:1:1101:15338:10113	256	chr7	6777077	0	151M	*	0	0	GTTCAGCATCCCAGACACGGGCAGTCACTGCAGCAGCGAGTATGCGGCCAGCTCCCCAGGGGACAGGGGCTCGCAGGAACATGTGGACTCTCAGGAGAAAGCGCCTGAAACTGACGACTCTTTTTCAGATGTGGACTGCCATTCAAACCAG	,FFFFF,F,FFFFFFFFFFFF:F,F:F::F,FFFFFFFFFFFFFFFFFFFFFFF,F:FFF,:FFFFFFFFFF,FFF:FFFFF,FF,F:FFFFFF:F,FFFF:,::,,FFFFF:FFFF::F,FFF,:,,F:,FF,F,FF,FFF,FFFF,FF,	MD:Z:41G2T7A98	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:3	AS:i:2957
A00223:8:H7YG3DMXX:1:1101:15338:10113	16	chr7	6026731	4	151M	*	0	0	CTGGTTTGAATGGCAGTCCACATCTGAAAAAGAGTCGTCAGTTTCAGGCGCTTTCTCCTGAGAGTCCACATGTTCCTGCGAGCCCCTGTCCCCTGGGGAGCTGGCCGCATACTCGCTGCTGCAGTGACTGCCCGTGTCTGGGATGCTGAAC	,FF,FFFF,FFF,FF,F,FF,:F,,:,FFF,F::FFFF:FFFFF,,::,:FFFF,F:FFFFFF:F,FF,FFFFF:FFF,FFFFFFFFFF:,FFF:F,FFFFFFFFFFFFFFFFFFFFFFF,F::F:F,F:FFFFFFFFFFFF,F,FFFFF,	MD:Z:44T106	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:1	AS:i:2999	XS:i:2957
A00223:8:H7YG3DMXX:1:1101:15881:12242	0	chr7	6026908	11	151M	*	0	0	GTGCCCCGAGTCCTTCTCCACCTCCGCTCTGTCCGTAGGGTCACTGGGTCCGTGACTGGAACTCACTGCCTCTTTCTGAGATCTCAGGACGCCTTTGTCAGAGATGGCACCTGAAGTGCTAGAAGACAGCATACCCCTTTTCTGTCCTAGA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F:FFFFFFFFFFF:F,FFFFFFFFFFFFFFF:FFFFFFFFF:FF:FFF:FFFF:FFFFFFFF,FFFFFFFFF	MD:Z:80G70	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:1	AS:i:2999	XS:i:2895
A00223:8:H7YG3DMXX:1:1101:15881:12242	272	chr7	6776900	0	151M	*	0	0	TCTAGGACAGAAAAGGGGTATGCTGTCTTCTAGCACTTCAGGTGCCATCTCTGACAAAGGCGTCCTGAGATCTCAGAAAGAGGCAGTGAGTTCCAGTCACGGACCCAGTGACCCTACGGACAGAGCGGAGGTGGAGAAGGACTCGGGGCAC	FFFFFFFFF,FFFFFFFF:FFFF:FFF:FF:FFFFFFFFF:FFFFFFFFFFFFFFF,F:FFFFFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	MD:Z:22T0C17A28C28G50T0	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:6	AS:i:2895
```

the read `A00223:8:H7YG3DMXX:1:1101:1163:35383` is now fixed.


