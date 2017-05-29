# BWAMemNOp

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
      What kind of help
      Possible Values: [usage, markdown, xml]
    --samoutputformat
      Sam output format.
      Default: TypeImpl{name='SAM', fileExtension='sam', indexExtension='null'}
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

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make bwamemnop
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/mem/BWAMemNOp.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/mem/BWAMemNOp.java
)
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



