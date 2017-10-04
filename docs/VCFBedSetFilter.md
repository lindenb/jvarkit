# VCFBedSetFilter

Set FILTER for VCF if it doesn't intersects with BED.


## DEPRECATED

use GATK FilterVariants

## Usage

```
Usage: vcfbedsetfilter [options] Files
  Options:
    -B, --bed
      Tribble or Tabix bed file
    -d, --discard
      Discard filtered variants
      Default: false
    -f, --filter
      FILTER name
      Default: VCFBED
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -i, --inverse
      inverse selection
      Default: false
    -m, --map
      unindexed bed file, will be loaded in memory (faster than tribble/tabix 
      but memory consumming)
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * bed
 * filter


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
$ make vcfbedsetfilter
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfbed/VCFBedSetFilter.java)


<details>
<summary>Git History</summary>

```
Tue Jun 6 18:06:17 2017 +0200 ; postponed vcf ; https://github.com/lindenb/jvarkit/commit/bcd52318caf3cd76ce8662485ffaacaabde97caf
Fri May 19 17:10:13 2017 +0200 ; cont doc ; https://github.com/lindenb/jvarkit/commit/d2aea1eaa554d0498b197fb8fac01893b10ceb83
Tue May 9 10:40:20 2017 +0200 ; moving to jcommander ; https://github.com/lindenb/jvarkit/commit/88cfdecb60c1f193ae8b3176ad86181c4a15256b
Wed Mar 23 11:45:49 2016 +0100 ; vcfbed cont ; https://github.com/lindenb/jvarkit/commit/8e7370663fa5f3824ebab9de304a2ee4d90619c1
Sun Feb 28 14:54:44 2016 +0100 ; burden week-end ; https://github.com/lindenb/jvarkit/commit/2f49ec436743934639d0adf51b55a577db7ee3d2
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Mon Jun 1 15:27:11 2015 +0200 ; change getChrom() to getContig() ; https://github.com/lindenb/jvarkit/commit/5abd60afcdc2d5160164ae6e18087abf66d8fcfe
Fri May 15 18:07:33 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/346bd7d374bfe8f2c969de98ed176060b234f0e1
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbedsetfilter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Examples

```
$java -jar dist/vcfbedsetfilter.jar -f MYFILTER - -B in.bed in.vcf 

```


