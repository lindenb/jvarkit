# KnownGenesToBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

converts UCSC knownGenes file to BED.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar kg2bed  [options] Files

Usage: kg2bed [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -hide, --hide
      don't show the following items (comma separated, one of 
      'INTRON,UTR,CDS,EXON,TRANSCRIPT,NON_CODING,CODING'). Empty don't hide 
      anything 
      Default: <empty string>
    -o, --output
      Output file. Optional . Default: stdout
    -s, --select
      JEXL select expression. Object 'kg' is an instance of KnownGene (https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/ucsc/KnownGene.java).JEXL 
      stands for Java EXpression Language.  See 
      https://commons.apache.org/proper/commons-jexl/reference/syntax.html 
      Default: <empty string>
    -sql, --sql
      SQL Schema URI. Each instance of transcript can be associated to a .sql 
      schema to help the software to decode the semantics of the columns. Eg.: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV20.sql
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * ucsc
 * bed
 * knownGenes



## See also in Biostars

 * [https://www.biostars.org/p/151628](https://www.biostars.org/p/151628)
 * [https://www.biostars.org/p/9557497](https://www.biostars.org/p/9557497)



## Creation Date

20140311

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/kg2bed/KnownGenesToBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/kg2bed/KnownGenesToBed.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/kg2bed/KnownGenesToBedTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/kg2bed/KnownGenesToBedTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **kg2bed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


##  See also

**ucsc** tools :  genePredToBed

### Example

```
$ curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" |\
  gunzip -c |\
  java -jar dist/jvarkit.jar kg2bed
chr1	11873	14409	+	uc001aaa.3	TRANSCRIPT	uc001aaa.3
chr1	11873	12227	+	uc001aaa.3	EXON	Exon 1
chr1	12227	12612	+	uc001aaa.3	INTRON	Intron 1
chr1	11873	12227	+	uc001aaa.3	UTR	UTR3
chr1	12612	12721	+	uc001aaa.3	EXON	Exon 2
chr1	12721	13220	+	uc001aaa.3	INTRON	Intron 2
chr1	12612	12721	+	uc001aaa.3	UTR	UTR3
chr1	13220	14409	+	uc001aaa.3	EXON	Exon 3
chr1	13220	14409	+	uc001aaa.3	UTR	UTR3
chr1	11873	14409	+	uc010nxr.1	TRANSCRIPT	uc010nxr.1
chr1	11873	12227	+	uc010nxr.1	EXON	Exon 1
chr1	12227	12645	+	uc010nxr.1	INTRON	Intron 1
chr1	11873	12227	+	uc010nxr.1	UTR	UTR3
chr1	12645	12697	+	uc010nxr.1	EXON	Exon 2
chr1	12697	13220	+	uc010nxr.1	INTRON	Intron 2
```

