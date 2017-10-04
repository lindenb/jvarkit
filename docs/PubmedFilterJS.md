# PubmedFilterJS

Filters Pubmed XML with a javascript  (java rhino) expression. Context contain 'article' a  PubmedBookArticle or a PubmedArticle and 'index', the index in the XML file.


## Usage

```
Usage: pubmedfilterjs [options] Files
  Options:
    -e, --expression
      Javascript expression
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -f, --scriptfile
      Javascript file
    --version
      print version and exit

```


## Keywords

 * pubmed
 * javascript
 * xml
 * ncbi


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
$ make pubmedfilterjs
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedFilterJS.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/pubmed/PubmedFilterJS.java)


<details>
<summary>Git History</summary>

```
Fri May 19 17:10:13 2017 +0200 ; cont doc ; https://github.com/lindenb/jvarkit/commit/d2aea1eaa554d0498b197fb8fac01893b10ceb83
Fri Apr 14 17:55:48 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/298923b61cb33a9a3d187baf115e13cff4b74251
Wed Apr 5 20:37:28 2017 +0200 ; fix pubmed dtd ; https://github.com/lindenb/jvarkit/commit/b9cdddc89997662bd29971e88da307d48de84b86
Tue Apr 4 17:09:36 2017 +0200 ; vcfgnomad ; https://github.com/lindenb/jvarkit/commit/eac33a01731eaffbdc401ec5fd917fe345b4a181
Mon Dec 14 12:37:25 2015 +0100 ; pubmed filter js ; https://github.com/lindenb/jvarkit/commit/4a801c730cafb01176646ea1969f8314bb1bbdf4
Mon Dec 14 12:34:42 2015 +0100 ; pubmed filter js ; https://github.com/lindenb/jvarkit/commit/b94e9618802392f54fdf44233072c4cb716bab6c
Mon Jun 15 17:24:26 2015 +0200 ; vep as xml base 64 ; https://github.com/lindenb/jvarkit/commit/d629603576c14a970208c9599b39ecb2a8b39994
Mon Sep 1 15:47:33 2014 +0200 ; fix error in shuffle. pubmed filter js ; https://github.com/lindenb/jvarkit/commit/1c690742a61ba809e342eacd5d2a214134bfab72
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pubmedfilterjs** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


