# FindAVariation

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Finds a specific mutation in a list of VCF files


## Usage

```
Usage: findavariation [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -homref, --homref
      Hide HOM_REF genotypes
      Default: false
    -indexed, --indexed
      [20171020] Search only in indexed vcf
      Default: false
    -nocall, --nocall
      Hide NO_CALL genotypes
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -f, --posfile
      Add this file containing chrom:position
      Default: []
    -p, --position
      A list of 'chrom/position'
      Default: []
    -snp, --snp
      Search only variant have the very same position (ignore overlapping 
      variants) 
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * variation
 * search


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew findavariation
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindAVariation.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindAVariation.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/FindAVariationTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/FindAVariationTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **findavariation** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Example

```
$ find ./ -name "*.vcf" -o -name "*.vcf.gz" |\
   java -jar dist/findavariation.jar -p "chr1:1234" 


htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12878	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12891	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12892	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12878	HOM_REF	C C
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12891	HET	C T
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12892	HET	C T
```

## History

  * 20180914 : replace DP4 with AD
 
