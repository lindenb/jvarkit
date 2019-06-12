# VcfBurdenRscriptV

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Fisher Case / Controls per Variant (Vertical)


## Usage

```
Usage: vcfburdenrscriptv [options] Files
  Options:
    --cadd, -cadd
      [20180831] Include CADD data, if available (INFO/CADD_PHRED 
      INFO/CADD_SCORE) 
      Default: false
    -cpm, --cadd-phred-missing
      [20180831] value for CADD / phred missing data
      Default: NA
    -csm, --cadd-score-missing
      [20180921] value for CADD / score missing data
      Default: NA
    -f, --function
      User defined R function to be called after each VCF
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -if, --ignorefilter
      accept variants having a FILTER column. Default is ignore variants with 
      a FILTER column
      Default: false
    -minusnineiszero, --minusnineiszero
      No Call is '0' (default is -9)
      Default: false
    --nfe
      [20180910] INCLUDE gnomad genome NFE AC
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -p, --pedigree
      [20180831] pedigree file (or I will try to extract the pedigree from the 
      vcf header. A pedigree is a text file delimited with tabs. No header. 
      Columns are (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother 
      Id or '0' (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 
      unaffected, 1 affected,-9 unknown
    -t, --title
      Try to find ##(TITLE)=abcdefghijk in the VCF header and use it as the 
      name of the VCF chunk
      Default: <empty string>
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
$ ./gradlew vcfburdenrscriptv
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenRscriptV.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenRscriptV.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdenrscriptv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Variant in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.


### Output


#### INFO column

 *  BurdenF1Fisher : Fisher test

#### FILTER column

 *  BurdenF1Fisher :Fisher test doesn't meet  user's requirements

### see also

 *  VcfBurdenFilter3




### Output

#### INFO column


 *  BurdenF1Fisher : Fisher test

#### FILTER column

 *  BurdenF1Fisher :Fisher test doesn't meet  user's requirements

### see also

 *  VcfBurdenFilter3

### History 

  * [20180910] add NFE https://www.youtube.com/watch?v=Fi5dLGAH8R0
  * [20180831] add CADD values from VcfCadd

