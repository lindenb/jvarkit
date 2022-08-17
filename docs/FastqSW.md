# FastqSW

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

align fasta sequences vs fastq


## Usage

```
Usage: java -jar dist/fastqsw.jar  [options] Files
Usage: fastqsw [options] Files
  Options:
    --disable-count-gap
      Do not count gaps when calculating percentage of identity.
      Default: false
    --disable-reverse-complement
      disable search in reverse-complement.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --inverse
      Inverse logic, discard matching reads.
      Default: false
    --max-gap
      Maximum number of gaps. -1 to ignore
      Default: -1
    -md5, --md5
      write md5 file
      Default: false
    --min-distance
      Minimum distance
      Default: 0.0
    --min-identity
      Minimum identity. A decimal number between 0.0 and 1.0. If the value 
      ends with '%' it is interpretted as a percentage eg. '1%' => '0.01'. A 
      slash '/' is interpretted as a ratio. e.g: '1/100' => '0.01'.
      Default: 9.0
    --min-length
      Minimum length
      Default: 0
    --min-similarity
      Minimum similarity. A decimal number between 0.0 and 1.0. If the value 
      ends with '%' it is interpretted as a percentage eg. '1%' => '0.01'. A 
      slash '/' is interpretted as a ratio. e.g: '1/100' => '0.01'.
      Default: 0.0
    -o, --out, -R1
      Output file for R1 fastq record or interleaved output.Output file. 
      Optional . Default: stdout
    --pair-logical
      Logical choice for pair of reads.
      Default: OR
      Possible Values: [OR, AND, XOR, XNOR, NAND]
    --paired
      assume input is paired end: we expect two files, or the input is assumed 
      interleaved fastq.
      Default: false
    --pairwise-aligner-type
      One item from 
      org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType 
      Default: LOCAL
      Possible Values: [GLOBAL, GLOBAL_LINEAR_SPACE, LOCAL, LOCAL_LINEAR_SPACE]
    -f, --queries
      Fasta file containing queries
    --save-align
      Save alignments in this file (mostly for debugging)
    -s, --string-queries
      Query as a DNA string
      Default: []
    --version
      print version and exit
    -R2
      Output file for R2 fastq record. If input is paired but R2 is omitted, 
      output will be interleaved.
    -extension-penalty
      Gap extension penalty.
      Default: 1
    -open-penalty
      Gap open penalty.
      Default: 10

```


## Keywords

 * fastq
 * align
 * sw



## See also in Biostars

 * [https://www.biostars.org/p/9534472](https://www.biostars.org/p/9534472)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew fastqsw
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20220207

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastqsw/FastqSW.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastqsw/FastqSW.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastqsw** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example:

```
$ java -jar dist/fastqsw.jar -s 'TCGTGACCTCTCAGGTTAGACCTACAAAATCTTCGATCAAT' \
	--min-identity 0.9 --save-align jeter.txt --min-length 40 --pair-logical OR \
	src/test/resources/S1.R1.fq.gz src/test/resources/S1.R2.fq.gz


(...)

@RF11_137_644_0:0:0_2:0:0_23/1
GCTCCCTCGTGACCTCTCAGGTTAGACCTACAAATCTTCGATCAATAGCATTGCGACTTGTTTCATTCTC
+
2222222222222222222222222222222222222222222222222222222222222222222222
@RF11_137_644_0:0:0_2:0:0_23/2
GTACATTTCACCAGATGCAGAAGCATTCAGTAAATACATGCTGTCAAAGTCTCCAGAAGATATTGGACCA
+
2222222222222222222222222222222222222222222222222222222222222222222222

$ tail -6 jeter.txt 
#RF11_137_644_0:0:0_2:0:0_23/1 distance:0.46 score:189.0 similarity:0.54 length:41 identicals:40 similars:40 identity:0.975609756097561
    1                                      41
  7 TCGTGACCTCTCAGGTTAGACCTAC-AAATCTTCGATCAAT 46
    ||||||||||||||||||||||||| |||||||||||||||
  1 TCGTGACCTCTCAGGTTAGACCTACAAAATCTTCGATCAAT 41
```



