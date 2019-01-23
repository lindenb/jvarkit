# FastqToFasta

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

fastq -> fasta


## DEPRECATED

use awk, samtools...

## Usage

```
Usage: fastq2fasta [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -N
      fasta line length
      Default: 50
    -b
      trim fasta header after space
      Default: false

```


## Keywords

 * fastq
 * fasta


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew fastq2fasta
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FastqToFasta.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FastqToFasta.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastq2fasta** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example

```bash
$ java -jar dist/fastq2fasta.jar  -N 60 -b file.fastq.gz

>HWI-1KL149:61:D2C11TCXX:2:1213:4591:29626
GAGTTGCTTTGTTTGAATATAGGTTGACTATACGAAGTGTGCGAGGACCTGCACCACGCA
GTAGGCCAAGATCAACTGAAACAGTGCTATCTGCACGACAA
>HWI-1KL149:61:D2C11TCXX:2:1213:4525:29650
CCTAGTAGTTCGTGGCCCCGGGCCCCTACTTAAACTCCTAGAACCACTCCTAGAAAGGGG
TGTTGCAGTTCGGCTCAGTCCCCGTGGTCGACTACTGTTTC
>HWI-1KL149:61:D2C11TCXX:2:1213:4569:29706
GCGCAGAGTTGTTTTAGCTATGCTGTGTTTGCATGGTTAGGTGGTGTACCTAGTGGTTTT
CTGAGACTTCTCTGAGGTTCTTGAGTAGATTAATACATCCC
>HWI-1KL149:61:D2C11TCXX:2:1213:4594:29713
```

