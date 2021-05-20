# GroupByGene

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff


## Usage

```
Usage: groupbygene [options] Files
  Options:
    -e, -E, --extractors
      [20190626]Gene Extractors Name. Space/semicolon/Comma separated
      Default: ANN/GeneId VEP/GeneId BCSQ/gene SMOOVE
    --filtered
      ignore FILTERED variants
      Default: false
    --fisher
      [20170726] Print fisher for case/control (experimental, need to work on 
      this) 
      Default: false
    --gtFiltered
      [20170725] ignore FILTERED genotypes
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -l, --list
      [20190626]list all available gene extractors
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    -p, --ped, --pedigree
      [20170725] A pedigree file. tab delimited. Columns: 
      family,id,father,mother, sex:(0:unknown;1:male;2:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit
    -positions
      include variants positions in the output table.
      Default: false

```


## Keywords

 * vcf
 * gene



## See also in Biostars

 * [https://www.biostars.org/p/342790](https://www.biostars.org/p/342790)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew groupbygene
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20131209

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/groupbygene/GroupByGene.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/groupbygene/GroupByGene.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/groupbygene/GroupByGeneTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/groupbygene/GroupByGeneTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **groupbygene** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Group VCF data by gene/transcript. By default it tries to use data from VEP and SnpEff

## Example

### Delimited output

```
$ curl -s -k "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" |\
java -jar dist/groupbygene.jar |\
head | column  -t

#chrom  min.POS    max.POS    gene.name  gene.type         samples.affected  count.variations  M10475  M10478  M10500  M128215
chr10   52004315   52004315   ASAH2      snpeff-gene-name  2                 1                 0       0       1       1
chr10   52004315   52004315   ASAH2      vep-gene-name     2                 1                 0       0       1       1
chr10   52497529   52497529   ASAH2B     snpeff-gene-name  2                 1                 0       1       1       0
chr10   52497529   52497529   ASAH2B     vep-gene-name     2                 1                 0       1       1       0
chr10   48003992   48003992   ASAH2C     snpeff-gene-name  3                 1                 1       1       1       0
chr10   48003992   48003992   ASAH2C     vep-gene-name     3                 1                 1       1       1       0
chr10   126678092  126678092  CTBP2      snpeff-gene-name  1                 1                 0       0       0       1
chr10   126678092  126678092  CTBP2      vep-gene-name     1                 1                 0       0       0       1
chr10   135336656  135369532  CYP2E1     snpeff-gene-name  3                 2                 0       2       1       1
```

```

## History

* 201707: added pedigree, removed XML output


