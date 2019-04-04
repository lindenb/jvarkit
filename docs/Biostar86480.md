# Biostar86480

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Genomic restriction finder


## Usage

```
Usage: biostar86480 [options] Files
  Options:
    -E, --enzyme
      restrict to that enzyme.
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * rebase
 * genome
 * enzyme
 * restricion
 * genome



## See also in Biostars

 * [https://www.biostars.org/p/86480](https://www.biostars.org/p/86480)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar86480
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar86480.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar86480.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar86480Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar86480Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar86480** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example
```bash
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz" |\
gunzip -c  |\
java -jar dist/biostar86480.jar -E AarI -E EcoRI  

chr3	60645	60651	GAATTC	1000	+	EcoRI	G^AATTC
chr3	60953	60959	GAATTC	1000	+	EcoRI	G^AATTC
chr3	68165	68172	GCAGGTG	1000	-	AarI	CACCTGC(4/8)
chr3	70263	70269	GAATTC	1000	+	EcoRI	G^AATTC
chr3	70945	70952	GCAGGTG	1000	-	AarI	CACCTGC(4/8)
chr3	71140	71146	GAATTC	1000	+	EcoRI	G^AATTC
chr3	72264	72270	GAATTC	1000	+	EcoRI	G^AATTC
chr3	74150	74156	GAATTC	1000	+	EcoRI	G^AATTC
chr3	75063	75069	GAATTC	1000	+	EcoRI	G^AATTC
chr3	78438	78444	GAATTC	1000	+	EcoRI	G^AATTC
chr3	81052	81059	CACCTGC	1000	+	AarI	CACCTGC(4/8)
chr3	84498	84504	GAATTC	1000	+	EcoRI	G^AATTC
chr3	84546	84552	GAATTC	1000	+	EcoRI	G^AATTC
chr3	84780	84787	CACCTGC	1000	+	AarI	CACCTGC(4/8)
chr3	87771	87777	GAATTC	1000	+	EcoRI	G^AATTC
chr3	95344	95351	GCAGGTG	1000	-	AarI	CACCTGC(4/8)
chr3	96358	96364	GAATTC	1000	+	EcoRI	G^AATTC
chr3	96734	96740	GAATTC	1000	+	EcoRI	G^AATTC
chr3	105956	105962	GAATTC	1000	+	EcoRI	G^AATTC
chr3	107451	107457	GAATTC	1000	+	EcoRI	G^AATTC
(...)
```
