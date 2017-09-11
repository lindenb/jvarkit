# SortVcfOnInfo

Sort a VCF a field in the INFO column


## Usage

```
Usage: sortvcfoninfo [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
  * -T, --tag, -t
      INFO tag
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * sort
 * annotation


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
$ make sortvcfoninfo
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sortvcfonref/SortVcfOnInfo.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sortvcfonref/SortVcfOnInfo.java)

Git History for this file:
```
Tue Jun 6 18:06:17 2017 +0200 ; postponed vcf ; https://github.com/lindenb/jvarkit/commit/bcd52318caf3cd76ce8662485ffaacaabde97caf
Sun Jun 4 21:53:22 2017 +0200 ; writing bcf ; https://github.com/lindenb/jvarkit/commit/784fdac37cd7e6eca04e35d0a3ddad8637826b4a
Mon May 15 12:10:21 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/b4895dd40d1c34f345cd2807f7a81395ba27e8ee
Wed Apr 26 17:26:23 2017 +0200 ; cont jcommander ; https://github.com/lindenb/jvarkit/commit/ab6c7b760cd5376e08da24426cede7f84a6b3ae2
Fri May 15 18:07:33 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/346bd7d374bfe8f2c969de98ed176060b234f0e1
Thu Sep 11 09:36:01 2014 +0200 ; problem with java dataInputSTream: writeUTF requires line.length < SHORt_MAX ; https://github.com/lindenb/jvarkit/commit/19eac4ee36909a730903546b50461de3c19a5c1f
Mon Jun 23 12:34:44 2014 +0200 ; find-a-variation + using abstractcodec instead of vcfcodec ; https://github.com/lindenb/jvarkit/commit/da621ba8326d56da8f6907c845c539e4ea785284
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Wed Feb 26 11:54:02 2014 +0100 ; sort on info ; https://github.com/lindenb/jvarkit/commit/88011f4fe0612a962c335ea0f92b35828501aec8
Tue Feb 18 18:05:31 2014 +0100 ; sortvcf on info, vcfutils parse header, pad call/format for fixvcfformat, vcfcadd ; https://github.com/lindenb/jvarkit/commit/6123910f68df940c1f3986d142f9b0414f76a43a
```

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sortvcfoninfo** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


##Example

```bash
]$ curl  "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" |\
   java -jar dist/sortvcfoninfo.jar -F BaseQRankSum | grep -vE "^#" 

chr10	1142208	.	T	C	3404.30	.	AC=8;AF=1.00;AN=8;
chr10	135336656	.	G	A	38.34	.	AC=4;AF=1.00;AN=4;
chr10	52004315	.	T	C	40.11	.	AC=4;AF=1.00;AN=4;
chr10	52497529	.	G	C	33.61	.	AC=4;AF=1.00;AN=4;
chr10	126678092	.	G	A	89.08	.	AC=1;AF=0.13;AN=8;BaseQRankSum=-3.120;
chr16	72057435	.	C	T	572.98	.	AC=1;AF=0.13;AN=8;BaseQRankSum=-2.270;
chr10	48003992	.	C	T	1047.87	.	AC=4;AF=0.50;AN=8;BaseQRankSum=-0.053;
chr10	135210791	.	T	C	65.41	.	AC=4;AF=0.50;AN=8;BaseQRankSum=2.054;
chr10	135369532	.	T	C	122.62	.	AC=2;AF=0.25;AN=8;BaseQRankSum=2.118;
```


