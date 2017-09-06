# VcfFilterNotInPedigree

Adds a FILTER NotInPedigree if the only not(homref) genotypes are not in a pedigree


## Usage

```
Usage: vcffilternotinpedigree [options] Files
  Options:
    -f, --filter
      FILTER name. Will be set for variant where the only genotypes non-homref 
      are NOT in the pedigree
      Default: NoGenotypeInPedigree
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -r, --remove
      remove the variant instead of setting the FILTER
      Default: false
    -sf, --sfilter
      FILTER name for option singleton
      Default: SingletonAlt
    -s, --singleton
      Variant is flagged/FILTERed as SingletonAlt if the ALT if found in less 
      or equal times 'singleton-times' in the genotypes. -1:ignore
      Default: 1
    --version
      print version and exit

```


## Keywords

 * burden
 * vcf
 * pedigree


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
$ make vcffilternotinpedigree
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfFilterNotInPedigree.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfFilterNotInPedigree.java)

Git History for this file:
```
Fri Aug 4 16:40:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/57f08e720a97f952bab81961431d83accdefeae3
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Wed May 10 20:57:52 2017 +0200 ; YC tag ; https://github.com/lindenb/jvarkit/commit/a9515d969d27c76ccd0814a093e886d71904b0f2
Tue May 24 16:49:08 2016 +0200 ; error in fisherburdenv ; https://github.com/lindenb/jvarkit/commit/ad8d8c2252786a71d063854fddcf5c3e276a052e
Tue May 3 17:34:10 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/7d668372271a7ecc28da6051c0ef251f70bbece9
Thu Apr 21 10:39:25 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/7adf87adc987efbe89def5c530f5a84be0c841d4
Tue Apr 19 14:25:22 2016 +0200 ; vcfnotinped ; https://github.com/lindenb/jvarkit/commit/265b0d11a280ad1458038fbd838a7a866952facf
```

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcffilternotinpedigree** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




