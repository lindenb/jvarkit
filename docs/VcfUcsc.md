# VcfUcsc

annotate an VCF with mysql UCSC data


## Usage

```
Usage: vcfucsc [options] Files
  Options:
    -D, --database
      database name
      Default: hg19
    -e, --expression
      What should be displayed in the INFO field. expression string using 
      column offsets. or names e.g: "${chrom}:${2}-${3}"
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
  * -T, -t, --table
      table name
    -tag, --tag
      tag prefix. Default = UCSC_${db}_${table}
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * ucsc
 * mysql
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
$ make vcfucsc
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfucsc/VcfUcsc.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfucsc/VcfUcsc.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfucsc** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example


```
java -jar dist/vcfucsc.jar --table snp142 -e '${name}' input.vcf
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr3	124290753	.	G	C	579.77	.	UCSC_HG19_SNP142=rs145115089
chr3	124290943	.	A	G	491.77	.	UCSC_HG19_SNP142=rs7372055
chr3	124291069	.	G	A	266.77	.	UCSC_HG19_SNP142=rs7373767
chr3	124291171	.	C	CA	240.73	.	.
chr3	124291245	.	A	G	563.77	.	UCSC_HG19_SNP142=rs12695439
chr3	124291351	.	A	G	194.77	.	UCSC_HG19_SNP142=rs7613600
chr3	124291416	.	G	T	308.77	.	UCSC_HG19_SNP142=rs73189597
chr3	124291579	.	T	C	375.77	.	UCSC_HG19_SNP142=rs7649882
```


