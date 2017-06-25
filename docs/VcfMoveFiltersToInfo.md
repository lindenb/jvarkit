# VcfMoveFiltersToInfo

Move any FILTER to the INFO column. reset FILTER to PASS


## Usage

```
Usage: vcfmovefilterstoinfo [options] Files
  Options:
    -f, --filter
      INFO name. This tag will be used to store the previous filters
      Default: PREVIOUSLY_FILTERED_AS
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -t, --limitto
      If not empty, limit to those FILTERS. Multiple separated by comma/space.
      Default: []
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * format
 * info


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
$ make vcfmovefilterstoinfo
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfMoveFiltersToInfo.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfMoveFiltersToInfo.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfmovefilterstoinfo** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


For Matilde K: move the information in FILTER to the INFO column to keep a trace of the FILTERs.

## Example

```
$ cat input.vcf | grep -v "##"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16057607	rs201535778	G	GAAAA	.	.	AAC=2;AAF=0.0005;BEACON=T|solvebio,T|bob,T|solvebio-133,T|altruist,T|prism,T|kaviar;NS=3690;RAC=3688;RAF=0.9995;VTYPE=SNV
22	16057608	rs201535778	G	T	.	.	AAC=2;AAF=0.0005;BEACON=T|solvebio,T|bob,T|solvebio-133,T|altruist,T|prism,T|kaviar;NS=3690;RAC=3688;RAF=0.9995;VTYPE=SNV


$ java -jar dist/vcfmovefilterstoinfo.jar input.vcf | grep -v "##"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16057608	rs201535778	G	T	.	.	AAC=2;AAF=0.0005;NS=3690;PREVIOUSLY_FILTERED_AS=FT1;RAC=3688;RAF=0.9995;VTYPE=SNV
22	16058492	.	G	A	.	.	AAC=2;AAF=0.0005;NS=3708;PREVIOUSLY_FILTERED_AS=FT2;RAC=3706;RAF=0.9995;VTYPE=SNV


$ java -jar dist/vcfmovefilterstoinfo.jar -f OLDFILTER -t FT2 input.vcf | grep -v "##"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16057608	rs201535778	G	T	.	FT1	AAC=2;AAF=0.0005;NS=3690;RAC=3688;RAF=0.9995;VTYPE=SNV
22	16058492	.	G	A	.	.	AAC=2;AAF=0.0005;NS=3708;OLDFILTER=FT2;RAC=3706;RAF=0.9995;VTYPE=SNV

```


