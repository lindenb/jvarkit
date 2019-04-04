# VcfGapFrequent

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filter VCF annotated with external (AF or AC/AN) frequency information like vcfgnomad


## Usage

```
Usage: vcfgapfrequent [options] Files
  Options:
    -D, --database
      VCF database(s) used as a reference of frequent variants. One file 
      ending with '.list' is interpretted as a list of path.
      Default: []
    -f, --fields
      Where to peek the frequencies from the database.How to extract the 
      AlleleFrequencies from a variant. Multiple separated with comma or 
      semicolon. e.g: "AC/AN;exome_CEU_*;genome_NFE_AF;another_AC/another/AN". 
      Input is a set of AC/AN field pairs or/and AF field separated by 
      semicolon. 'x/y' means AC/AN fields. '*' will be replaced with AC and 
      AN, hence, 'exome_CEU_*' will be interpreted as 
      exome_CEU_AC/exome_CEU_AN. Other field will be interpreted as an AF 
      field. 
      Default: AC/AN;AF
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -M, --max-length
      Max segment length of NO_CALL
      Default: 1000000
    -m, --min-length
      Min segment length of NO_CALL
      Default: 1000
    -o, --output
      Output file. Optional . Default: stdout
    -C, -skip, --skip
      Skip missing whole chromosome
      Default: false
    -t, --max, --treshold
      Allele Frequency Treshold. Only Variant in database(s) with extracted AF 
      greater than this value are considered.
      Default: 0.4
    --version
      print version and exit

```


## Keywords

 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfgapfrequent
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/VcfGapFrequent.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/VcfGapFrequent.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/VcfGapFrequentTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/VcfGapFrequentTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgapfrequent** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ find /commun/data/pubdb/broadinstitute.org/gnomad/release-170228/ -name "*.vcf.gz" > gnomad.list

$ java -jar dist/vcfgapfrequent.jar --fields "AF_POPMAX" -C -D gnomad.list genotyped.22.vcf.gz  > out.bed

```

