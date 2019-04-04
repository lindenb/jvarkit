# VCFComposite

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

(in developpement) Finds Variants involved in a Het Composite Disease


## Usage

```
Usage: vcfcomposite [options] Files
  Options:
    --filter
      [20180718] set FILTER for the variants that are not part of a composite 
      mutation. 
      Default: NOT_COMPOSITE
    -gf, --genotype-filter
      A Java EXpression Language (JEXL) expressions to filter a genotye in a 
      VCF. Empty string will accept all genotypes. Expression returning a TRUE 
      will accept the genotypes. See 
      https://gatkforums.broadinstitute.org/gatk/discussion/1255 
      Default: <empty string> (ACCEPT ALL)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -max, --max, --max-variants
      [20180718] Max variants per gene. If different than -1, used to set a 
      maximum number of variants per gene; The idea is to filter out the gene 
      having a large number of variants.
      Default: -1
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --out
      Output file. Optional . Default: stdout
    -p, -ped, --pedigree
      A pedigree is a text file delimited with tabs. No header. Columns are 
      (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother Id or '0' 
      (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 unaffected, 1 
      affected,-9 unknown
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    -vf, --variant-filter
      A Java EXpression Language (JEXL) expressions to filter the variants 
      from a VCF. Empty string will accept all variants. Expression returning 
      a TRUE will accept the variant. See 
      https://gatkforums.broadinstitute.org/gatk/discussion/1255 
      Default: <empty string> (ACCEPT ALL)
    --version
      print version and exit

```


## Keywords

 * vcf
 * disease
 * annotation
 * pedigree


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfcomposite
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcomposite/VCFComposite.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcomposite/VCFComposite.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcomposite/VCFCompositeTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcomposite/VCFCompositeTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcomposite** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a VCF file annotated with SNPEff or VEP.



