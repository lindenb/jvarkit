# VcfGeneOntology


## Usage

```
Usage: vcfgo [options] Files
  Options:
    -h, --help
      print help and exits
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exits
  * -A
      (goa input url)
      Default: http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD
    -C
      (Go:Term) Add children to the list of go term to be filtered. Can be 
      used multiple times.
      Default: []
    -F
       if -C is used, don't remove the variant but set the filter
  * -G
      (go  url)
      Default: http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz
    -T
      INFO tag.
      Default: GOA
    -r
      remove variant if no GO term is found associated to variant
      Default: false
    -v
      inverse filter if -C is used
      Default: false

```


## Description

Find the GO terms for VCF annotated with SNPEFF or VEP


## Keywords

 * vcf
 * go


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
$ make vcfgo
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

https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfgo/VcfGeneOntology.java

## Contribute

- Issue Tracker: http://github.com/lindenb/jvarkit/issues
- Source Code: http://github.com/lindenb/jvarkit

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgo** ? https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md

The current reference is:

http://dx.doi.org/10.6084/m9.figshare.1425030

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> http://dx.doi.org/10.6084/m9.figshare.1425030


