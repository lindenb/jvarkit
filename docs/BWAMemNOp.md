# BWAMemNOp

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Merge the other BWA-MEM alignements with its SA:Z:* attributes to an alignment containing a cigar string with 'N' (  Skipped region from the reference.)


## Usage

```
Usage: bwamemnop [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit
    -S
       only print converted reads
      Default: false
    -m
      min soft clip length
      Default: 10
    -o
      Output file. Optional . Default: stdout

```


## Keywords

 * bwa
 * sam
 * bam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bwamemnop
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/mem/BWAMemNOp.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/mem/BWAMemNOp.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bwamemnop** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

The original SAM output:

```bash
$ samtools view file.bam | grep -F "HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811"
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	163	8	833200	60	62M39S	=	833209	62	ATCCCAGCAAGTATTGTACTCAGGCGTCTAGAAGCTCCAATCGCAGAATCCACAGAAAGAGCGCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACA	D@CCFFFFFGHGFHJGIHGGGIIIJGIGGHGHGDHIJJJJIJIIIJJJIGIGCDEHIJEHEDDEFFEECACDCDDDDD@CCDDDDDDDDDDDBDDDDDD>A	NM:i:0	MD:Z:62AS:i:62	XS:i:19	SA:Z:8,836189,+,62S34M1D5M,49,1;
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	83	8	833209	60	53M48S	=	833200	-62	AGTATTGTACTCAGGCGTCTAGAAGCTCCAATCGCAGAATCCACAGAAAGAGCGCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACAGGGCGCCCC	DBDBBB?BCCDDDCC=5CCCFFDFB?ACHECGGGIHCGGCCHGBGGJHGHGGGJJIFBJJIIHFGGIIJIGIGIIIGHGCIIGIHGIIHGHHHFFFFFCCC	NM:i:0	MD:Z:53AS:i:53	XS:i:19	SA:Z:8,836189,-,53S34M1D14M,60,1;
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	419	8	836189	49	62H34M1D5M	=	833209	-2929	GCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACA	DEFFEECACDCDDDDD@CCDDDDDDDDDDDBDDDDDD>A	NM:i:1	MD:Z:34^T5	AS:i:35	XS:i:26	SA:Z:8,833200,+,62M39S,60,0;
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	339	8	836189	60	53H34M1D14M	=	833200	-3038	GCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACAGGGCGCCCC	JJIFBJJIIHFGGIIJIGIGIIIGHGCIIGIHGIIHGHHHFFFFFCCC	NM:i:1	MD:Z:34^A14	AS:i:41	XS:i:26	SA:Z:8,833209,-,53M48S,60,0;
```

Using BWAMemNOp:

```
$ java -jar dist/bwamemnop.jar -S  file.bam | grep -F "HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811"
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	163	8	833200	60	62M2927N34M1D5M	=	833209	62	ATCCCAGCAAGTATTGTACTCAGGCGTCTAGAAGCTCCAATCGCAGAATCCACAGAAAGAGCGCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACA	D@CCFFFFFGHGFHJGIHGGGIIIJGIGGHGHGDHIJJJJIJIIIJJJIGIGCDEHIJEHEDDEFFEECACDCDDDDD@CCDDDDDDDDDDDBDDDDDD>A	MD:Z:62NM:i:0	AS:i:62	XS:i:19
HWI-1KL149:61:D2C11ACXX:4:1205:18030:70811	83	8	833209	60	53M2927N34M1D14M	=	833200	-62	AGTATTGTACTCAGGCGTCTAGAAGCTCCAATCGCAGAATCCACAGAAAGAGCGCCTGTAGCAGAAGGACTTGATTGATGTTGAATGCAACAGGGCGCCCC	DBDBBB?BCCDDDCC=5CCCFFDFB?ACHECGGGIHCGGCCHGBGGJHGHGGGJJIFBJJIIHFGGIIJIGIGIIIGHGCIIGIHGIIHGHHHFFFFFCCC	MD:Z:53	NM:i:0	AS:i:53	XS:i:19

```

converted to BED12 with bedtools/**bamToBed** and visualized in the UCSC genome browser... :


![http://i.imgur.com/NccxWdT.jpg](http://i.imgur.com/NccxWdT.jpg)


