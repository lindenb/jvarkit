# VcfPeekVcf


## Usage

```
Usage: vcfpeekvcf [options] Files
  Options:
    -a, --alt
      **ALL** alt allele must be found in indexed file.
      Default: false
    -h, --help
      print help and exits
    -o, --output
      Output file. Optional . Default: stdout
    -p, --prefix
      prefix all database tags with this prefix to avoid collisions
      Default: <empty string>
    -i, --replaceid
      Replace the ID field if it exists
      Default: false
  * -f, --tabix
      The VCF file indexed with TABIX or tribble. Source of the annotations
    -t, --tags
      tag1,tag2,tag... the INFO keys to peek from the indexed file
    --version
      print version and exits

```


## DEPRECATED

Use GATK or vcftools

## Description

Get the INFO from a VCF and use it for another VCF


## Keywords

 * vcf


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
$ make vcfpeekvcf
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

https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekVcf.java

## Contribute

- Issue Tracker: http://github.com/lindenb/jvarkit/issues
- Source Code: http://github.com/lindenb/jvarkit

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfpeekvcf** ? https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md

The current reference is:

http://dx.doi.org/10.6084/m9.figshare.1425030

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> http://dx.doi.org/10.6084/m9.figshare.1425030


