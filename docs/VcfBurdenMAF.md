# VcfBurdenMAF

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Burden : MAF for Cases / Controls 


## Usage

```
Usage: vcfburdenmaf [options] Files
  Options:
    -gtf, --gtf, --gtFiltered
      [20180117] Ignore FILTERed **Genotype**
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -c, --homref
      Treat No Call './.' genotypes as HomRef
      Default: false
    -ignoreFiltered, --ignoreFiltered
      [20171031] Don't try to calculate things why variants already FILTERed 
      (faster) 
      Default: false
    -lumpy-su-min, --lumpy-su-min
      [20180117] if variant identified as LUMPy-SV variant. This is the 
      minimal number of 'SU' to consider the genotype as a variant.
      Default: 1
    -maxMAF, --maxMAF
      if MAF of cases OR MAF of control is greater than maxMAF, the the FILTER 
      Column is Filled
      Default: 0.05
    -o, --output
      Output file. Optional . Default: stdout
    -p, --pedigree
      [20180117] Pedigree file. Default: use the pedigree data in the VCF 
      header.A pedigree is a text file delimited with tabs. No header. Columns 
      are (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother Id or 
      '0' (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 unaffected, 1 
      affected,-9 unknown
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * maf
 * case
 * control


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfburdenmaf
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenMAF.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenMAF.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdenmaf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



Variant in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.


### Output




#### INFO column


 *  BurdenMAFCas : MAF cases
 *  BurdenMAFControls : MAF controls





#### FILTER column


 *  BurdenMAFCas : MAF for cases  doesn't meet  user's requirements
 *  BurdenMAFControls : MAF for controls  doesn't meet  user's requirements
 *  BurdenMAFCaseOrControls : MAF for controls or cases  doesn't meet  user's requirements





### see also


 *  VcfBurdenFilterExac
 *  VcfBurdenFisherH



Variant in that VCF should have one and **only one** ALT allele. Use [https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele](https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele) if needed.

### Output


#### INFO column

  * **BurdenMAFCas** : MAF cases
  * **BurdenMAFControls** : MAF controls

#### FILTER column

  * **BurdenMAFCas** : MAF for cases  doesn't meet  user's requirements
  * **BurdenMAFControls** : MAF for controls  doesn't meet  user's requirements
  * **BurdenMAFCaseOrControls** : MAF for controls or cases  doesn't meet  user's requirements

