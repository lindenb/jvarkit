# BedMergeCnv

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Merge Bed records if they overlap a fraction of their lengths.


## Usage

```
Usage: bedmergecnv [options] Files
  Options:
    -f, --fraction, --overlap
      Intervals will be merged if they both overlap this fraction of their 
      lengths. A decimal number between 0.0 and 1.0. If the value ends with 
      '%' it is interpretted as a percentage eg. '1%' => '0.01'. A slash '/' 
      is interpretted as a ratio. e.g: '1/100' => '0.01'.
      Default: 0.9
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * bed
 * chromosome
 * contig


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bedmergecnv
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200330

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedtools/BedMergeCnv.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedtools/BedMergeCnv.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedmergecnv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


input is bed on standard input or it's a set of interval files (.bed, .interval_list, .gtf, etc... )

output is a BED file:

* 4th column: number of items merged
* 5th column: fraction overlap with previous record
* 6th column an optional label. Concatenation of the distinct label. For a bed record the label is the 4th column. For a VCF record, the label is the name of the samples

## Example

```
$ java -jar dist/bedmergecnv.jar src/test/resources/manta.B00*.vcf.gz | more

chr21	9653162	9653235	1	0.0
chr21	9653169	9653243	1	0.89
chr21	9653357	9653502	1	0.00
chr21	9660834	9660835	3	0.00
chr21	9841718	9842058	2	0.00
chr21	9846433	9846531	1	0.00
chr21	9846456	9846557	1	0.74
chr21	9846481	9846580	1	0.75
chr21	9851247	9854008	1	0.00
chr21	9862002	9862003	1	0.00
chr21	9881949	9882060	1	0.00
chr21	9936419	9936560	1	0.00
chr21	10459126	10459235	3	0.00
chr21	10475193	10475513	1	0.00
chr21	10492877	10493042	3	0.00
chr21	10493815	10493874	1	0.00
chr21	10618665	10618720	1	0.00
chr21	10633657	10633709	1	0.00
chr21	10699831	10700200	2	0.00
(...)
```


