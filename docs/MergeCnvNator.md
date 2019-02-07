# MergeCnvNator

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Merge CNVNator results


## Usage

```
Usage: mergecnvnator [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    -r, --ratio
      two intervals are the same if they both have more or equals of this 
      fraction of length in common
      Default: 0.75
    --version
      print version and exit
    -R, -reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary

```


## Keywords

 * cnv
 * indel
 * cnvnator


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew mergecnvnator
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/MergeCnvNator.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/MergeCnvNator.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/MergeCnvNatorTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/MergeCnvNatorTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **mergecnvnator** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## input

Input is the tabular output of CNVNator

```
(...)
deletion     chr2:1-2649000        2.649e+06  0        6.01633e-14  0            6.02087e-14  0           -1
duplication  chr2:3712001-3721000  9000       1.89036  0.0274362    5.17568e-42  0.137369     3.6838e-64  0.00821444
(...)
```

The name of each sample is the `basename` of the file, before the first `.`

a list of paths can be specified if the only input file ends with '.list' 


## Example

```
find DIR1 DIR2 -type f -name "*.tsv" > in.list
java -jar dist/mergecnvnator.jar -R ref.fasta in.list > out.vcf

(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample7	Sample8	Sample9
chr1	1	.	N	<DEL>	.	.	END=10000;IMPRECISE;SVLEN=10000;SVTYPE=DEL	GT:CN:P1:P2:P3:P4:Q0:RD	1/1:0:1.594e-11:0.00:1.992e-11:0.00:-1.000e+00:0.00	1/1:0:1.594e-11:0.00:1.992e-11:0.00:-1.000e+00:0.00	1/1:0:1.594e-11:0.00:1.992e-11:0.00:-1.000e+00:0.00
chr1	40001	.	N	<DEL>	.	.	END=48000;IMPRECISE;SVLEN=8000;SVTYPE=DEL	GT:CN:P1:P2:P3:P4:Q0:RD	./.	./.	0/1:1:34.95:5.850e-06:323.37:2.143e-03:0.889:0.603
(...)
```

## See also

  * https://github.com/abyzovlab/CNVnator

