# Biostar90204

Bam version of linux split. See also http://www.biostars.org/p/90204/


## Usage

```
Usage: biostar90204 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --samoutputformat
      Sam output format.
      Default: TypeImpl{name='SAM', fileExtension='sam', indexExtension='null'}
    --version
      print version and exit
    -M
       manifest file. Optional
    -a
      suffix length
      Default: 2
    -n
      Records per file
      Default: -1
    -p
      (prefix) output file prefix.
      Default: _splitbam

```


## Keywords

 * sam
 * bam
 * split
 * util



## See also in Biostars

 * [https://www.biostars.org/p/90204](https://www.biostars.org/p/90204)


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
$ make biostar90204
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar90204.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar90204.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar90204** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


##Example

```bash
$ java -jar dist/biostar90204.jar -m bam.manifest -n 3 -a 5 samtools-0.1.18/examples/toy.sam

$ cat bam.manifest
_splitbam.00001.bam	1	3
_splitbam.00002.bam	4	6
_splitbam.00003.bam	7	9
_splitbam.00004.bam	10	12

$ samtools-0.1.18/samtools view -h _splitbam.00003.bam 
@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.biostar.Biostar90204	VN:7e17f8bd273cf081d4415bc4f579cd34e2c681d1	CL:-m bam.manifest -n 3 -a 5 samtools-0.1.
18/examples/toy.sam
@CO	SPLIT:3
@CO	SPLIT:Starting from Read7
x1	0	ref2	1	30	20M	*	0	0	AGGTTTTATAAAACAAATAA	????????????????????
x2	0	ref2	2	30	21M	*	0	0	GGTTTTATAAAACAAATAATT	?????????????????????
x3	0	ref2	6	30	9M4I13M	*	0	0	TTATAAAACAAATAATTAAGTCTACA	??????????????????????????
```


