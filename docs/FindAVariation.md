# FindAVariation

Finds a specific mutation in a list of VCF files


## Usage

```
Usage: findavariation [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -homref, --homref
      Hide HOM_REF genotypes
      Default: false
    -nocall, --nocall
      Hide NO_CALL genotypes
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    -f, --posfile
      Add this file containing chrom:position
      Default: []
    -p, --position
      A list of 'chrom/position'
      Default: []
    -snp, --snp
      Search only variant have the very same position (ignore overlapping 
      variants) 
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * variation
 * search


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
$ make findavariation
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindAVariation.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FindAVariation.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **findavariation** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Example

```
$ find ./ -name "*.vcf" -o -name "*.vcf.gz" |\
   java -jar dist/findamutation.jar -p "chr1:1234" 


htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12878	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12891	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12892	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12878	HOM_REF	C C
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12891	HET	C T
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12892	HET	C T
```

 

