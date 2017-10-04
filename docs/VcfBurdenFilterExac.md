# VcfBurdenFilterExac

Burden filter 3 - Exac


## Usage

```
Usage: vcfburdenexac [options] Files
  Options:
    -d, --discardNotInExac
      if variant was not found in Exac, set the FILTER. Default: don't set the 
      FILTER. 
      Default: false
  * -exac, --exac
      Path to Exac VCF file. At the time of writing, you'd better use a 
      normalized version of Exac (see 
      https://github.com/lindenb/jvarkit/wiki/VCFFixIndels )
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -maxFreq, --maxFreq
      set FILTER if max(exac frequency in any pop) is greater than this value)
      Default: 0.001
    -o, --output
      Output file. Optional . Default: stdout
    -pop, --population
      comma separated populations in exac
      Default: AFR,AMR,EAS,FIN,NFE,SAS
    -tabix, --tabix
      use tabix index for Exac it is present. Might speed up things if the 
      number of variant is low.
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * exac


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
$ make vcfburdenexac
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFilterExac.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFilterExac.java)


<details>
<summary>Git History</summary>

```
Wed Sep 20 15:52:53 2017 +0200 ; moving to amalgamation ; https://github.com/lindenb/jvarkit/commit/fca74f53afa062f238c8a899ee0ee6e7cd15136c
Fri Sep 8 11:29:13 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/0c5035b59124ca18aba0405d0c616b565a32d10e
Fri Aug 4 16:40:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/57f08e720a97f952bab81961431d83accdefeae3
Mon Jun 26 17:29:03 2017 +0200 ; burden ; https://github.com/lindenb/jvarkit/commit/a3b7abf21d07f0366e81816ebbb2cce26b2341e7
Fri Jun 23 15:26:55 2017 +0200 ; updated vcf2multiallele ; https://github.com/lindenb/jvarkit/commit/775e8ddcc38a3e283cf49d9287b06510d7634e31
Mon May 22 14:34:36 2017 +0200 ; canvas burden ; https://github.com/lindenb/jvarkit/commit/ce6092de9cc4bea4d35d848410f5559d3d76e235
Wed May 10 20:57:52 2017 +0200 ; YC tag ; https://github.com/lindenb/jvarkit/commit/a9515d969d27c76ccd0814a093e886d71904b0f2
Thu Apr 21 10:39:25 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/7adf87adc987efbe89def5c530f5a84be0c841d4
Mon Apr 18 17:34:40 2016 +0200 ; cnot burden ; https://github.com/lindenb/jvarkit/commit/e0403a175b479d9e8bec1ced1e3f35715f404ad8
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdenexac** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


20170626: this tool now supports multiple ALT in the user VCF, however it's not been tested for choosing when to set the FILTER or the min value

### Output

#### INFO column


 *  FreqExac : Exac frequency.
 *  AC_* and AN_*: Transpose original population data from original Exac file


#### FILTER column

 *  BurdenExac : if FreqExac doesn't fit the criteria maxFreq


### see also


 *  VcfBurdenMAF



