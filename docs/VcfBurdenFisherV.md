# VcfBurdenFisherV

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Fisher Case / Controls per Variant (Vertical)


## Usage

```
Usage: vcfburdenfisherv [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -if, --ignorefilter
      accept variants having a FILTER column. Default is ignore variants with 
      a FILTER column
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
  * -p, --pedigree
      A pedigree file. tab delimited. Columns: family,id,father,mother, 
      sex:(0:unknown;1:male;2:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    -table, --table
      Write statistics into that file instead of the VCF (faster)
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
$ ./gradlew vcfburdenfisherv
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20160418

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFisherV.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFisherV.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdenfisherv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

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


