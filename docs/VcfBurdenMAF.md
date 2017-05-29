# VcfBurdenMAF

Burden : MAF for Cases / Controls 


## Usage

```
Usage: vcfburdenmaf [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -c, --homref
      Treat No Call './.' genotypes as HomRef
      Default: false
    -maxMAF, --maxMAF
      if MAF of cases OR MAF of control is greater than maxMAF, the the FILTER 
      Column is Filled
      Default: 0.05
    -o, --output
      Output file. Optional . Default: stdout
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

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcfburdenmaf
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenMAF.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenMAF.java
)
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


