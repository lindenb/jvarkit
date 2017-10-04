# VcfBurdenRscriptV

Fisher Case / Controls per Variant (Vertical)


## Usage

```
Usage: vcfburdenrscriptv [options] Files
  Options:
    -f, --function
      User defined R function to be called after each VCF
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -if, --ignorefilter
      accept variants having a FILTER column. Default is ignore variants with 
      a FILTER column
      Default: false
    -minusnineiszero, --minusnineiszero
      No Call is '0' (default is -9)
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
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

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcfburdenrscriptv
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenRscriptV.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenRscriptV.java)


<details>
<summary>Git History</summary>

```
Fri Aug 4 16:40:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/57f08e720a97f952bab81961431d83accdefeae3
Tue Jun 27 17:36:29 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/278970358111f7e3eca02e77d9a238321668a2dd
Mon Jun 26 17:29:03 2017 +0200 ; burden ; https://github.com/lindenb/jvarkit/commit/a3b7abf21d07f0366e81816ebbb2cce26b2341e7
Tue May 23 18:35:23 2017 +0200 ; lowres bam2raster ; https://github.com/lindenb/jvarkit/commit/e39a8f964b4bb11b28700c37ce1f2a7ba16b4653
Mon May 15 20:23:58 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/bb4421e107f53c95efdcad8fb54f022f9642312c
Fri Jun 24 16:32:07 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/36ff02783ed74a02e0802609a4397507660aff50
Thu Jun 9 10:17:58 2016 +0200 ; fix bug in maf/filter ; https://github.com/lindenb/jvarkit/commit/fb7a97b7d46904fe8f7d2e6b517d6ca4eb0f1117
Fri Jun 3 16:35:43 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/8990c064a122e300faddcab38bd476a9f8b9758e
Thu May 26 16:43:07 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/60ada53779722d3b5f4bff4d31b08cb518a38541
Wed May 25 10:49:01 2016 +0200 ; vcf burden R ; https://github.com/lindenb/jvarkit/commit/9e18b322b437bd9269c36406e35cdafce2ba71e5
Tue Apr 26 17:21:33 2016 +0200 ; vcfbuffer ; https://github.com/lindenb/jvarkit/commit/3300512769fd3bb2ee4430c9474367b06f2edc7c
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



