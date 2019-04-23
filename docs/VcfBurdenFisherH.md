# VcfBurdenFisherH

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Fisher Case /Controls per Variant


## Usage

```
Usage: vcfburdenfisherh [options] Files
  Options:
    --attribute
      [20190418] Name of the attribue used as FILTER and INFO
      Default: BurdenHFisher
    -gtf, --gtf, --gtFiltered
      [20180115] Ignore FILTERed **Genotype**
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -ignoreFiltered, --ignoreFiltered
      [20171031] Don't try to calculate things why variants already FILTERed 
      (faster) 
      Default: false
    -lumpy-su-min, --lumpy-su-min
      [20180115] if variant identified as LUMPy-SV variant. This is the 
      minimal number of 'SU' to consider the genotype as a variant.
      Default: 1
    -fisher, --minFisherPValue
      if p-value fisher(case/control vs have alt/have not alt) lower than 
      'fisher' the FILTER Column is Filled
      Default: 0.05
    -o, --output
      Output file. Optional . Default: stdout
    -p, --pedigree
      [20180115] Pedigree file. Default: use the pedigree data in the VCF 
      header.A pedigree is a text file delimited with tabs. No header. Columns 
      are (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother Id or 
      '0' (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 unaffected, 1 
      affected,-9 unknown
    --report
      [20190418] save report as bed file
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * fisher


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfburdenfisherh
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFisherH.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFisherH.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFisherHTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFisherHTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdenfisherh** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

Variants in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.

VCF header must contain a pedigree ( see VCFinjectPedigree ) or a pedigree must be defined.

## Lumpy-SV

 * 20180115: this tools recognize lumpy-sv genotypes


### see also

 *  VcfBurdenMAF
 *  VcfBurdenFilterExac

