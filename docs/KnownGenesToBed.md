# KnownGenesToBed

converts UCSC knownGenes file to BED.


## Usage

```
Usage: kg2bed [options] Files
  Options:
    -c, --cds
      Hide CDSs
      Default: false
    -x, --exon
      Hide Exons
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -i, --intron
      Hide Introns
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -t, --transcript
      Hide Transcript
      Default: false
    -u, --utr
      Hide UTRs
      Default: false
    --version
      print version and exit

```


## Keywords

 * ucsc
 * bed
 * knownGenes



## See also in Biostars

 * [https://www.biostars.org/p/151628](https://www.biostars.org/p/151628)


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
$ make kg2bed
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/KnownGenesToBed.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/KnownGenesToBed.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **kg2bed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Example



```
$ curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" |\
  gunzip -c |\
  java -jar dist/kg2bed.jar
chr1	11873	14409	+	uc001aaa.3	TRANSCRIPT	uc001aaa.3
chr1	11873	12227	+	uc001aaa.3	EXON	Exon 1
chr1	12227	12612	+	uc001aaa.3	INTRON	Intron 1
chr1	11873	12227	+	uc001aaa.3	UTR	UTR3
chr1	12612	12721	+	uc001aaa.3	EXON	Exon 2
chr1	12721	13220	+	uc001aaa.3	INTRON	Intron 2
chr1	12612	12721	+	uc001aaa.3	UTR	UTR3
chr1	13220	14409	+	uc001aaa.3	EXON	Exon 3
chr1	13220	14409	+	uc001aaa.3	UTR	UTR3
chr1	11873	14409	+	uc010nxr.1	TRANSCRIPT	uc010nxr.1
chr1	11873	12227	+	uc010nxr.1	EXON	Exon 1
chr1	12227	12645	+	uc010nxr.1	INTRON	Intron 1
chr1	11873	12227	+	uc010nxr.1	UTR	UTR3
chr1	12645	12697	+	uc010nxr.1	EXON	Exon 2
chr1	12697	13220	+	uc010nxr.1	INTRON	Intron 2

```



