# VCFReplaceTag

Replace the key for INFO/FORMAT/FILTER


## Usage

```
Usage: vcfreplacetag [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -k, --tag
      tag to replace . Format FROM/TO
      Default: []
  * -t, --type
      replace type: one of FORMAT,FILTER,INFO
      Default: INFO
    --version
      print version and exit

```


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
$ make vcfreplacetag
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstripannot/VCFReplaceTag.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstripannot/VCFReplaceTag.java)


<details>
<summary>Git History</summary>

```
Mon May 29 16:53:42 2017 +0200 ; moved to docs ; https://github.com/lindenb/jvarkit/commit/6c0535d7add884e75b424af89a4f00aff6fae75f
Mon May 22 17:20:59 2017 +0200 ; moving to jcommaner ; https://github.com/lindenb/jvarkit/commit/60cbfa764f7f5bacfdb78e48caf8f9b66e53a6a0
Tue Apr 4 17:45:44 2017 +0200 ; fix bug in replace tag ; https://github.com/lindenb/jvarkit/commit/8e2a57c24f82df69c730b944f4aed99f082fb3c5
Tue Apr 4 17:09:36 2017 +0200 ; vcfgnomad ; https://github.com/lindenb/jvarkit/commit/eac33a01731eaffbdc401ec5fd917fe345b4a181
Mon Jan 18 16:58:08 2016 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/83f80fdbe8d6be71539cfdbf60d61ce7ead9c0fd
Thu May 21 15:13:01 2015 +0200 ; vcf peek info ; https://github.com/lindenb/jvarkit/commit/b451191e148ad57f498032e2f81adc1051336c07
Tue May 19 17:40:33 2015 +0200 ; vcf index tabix ; https://github.com/lindenb/jvarkit/commit/9c0d532b734ce8f5d9e2a2b67c4b2e3a764ca1d6
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfreplacetag** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$   java -jar dist/vcfreplacetag.jar -t INFO -k VDB/NEWNAME ~/jeter.vcf 
```




