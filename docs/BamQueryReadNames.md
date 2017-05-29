# BamQueryReadNames

Query a Bam file indexed with BamIndexReadNames


## Usage

```
Usage:  [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --out
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: TypeImpl{name='SAM', fileExtension='sam', indexExtension='null'}
    --version
      print version and exit
    -N
       save unmatched names here
    -s
      user list of read names is sorted
      Default: false

```

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
$ make 
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamindexnames/BamQueryReadNames.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamindexnames/BamQueryReadNames.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
 
## Example



```bash
$cat read.names
HWI-1KL149:18:C0RNBACXX:1:1101:12800:7069
HWI-1KL149:18:C0RNBACXX:1:1101:15533:71685
HWI-1KL149:18:C0RNBACXX:1:1103:6001:91243
HWI-1KL149:18:C0RNBACXX:1:1107:2088:3461
HWI-1KL149:18:C0RNBACXX:1:1108:2098:26795
HWI-1KL149:18:C0RNBACXX:1:1110:10318:73043
HWI-1KL149:18:C0RNBACXX:1:1112:18688:6422
HWI-1KL149:18:C0RNBACXX:1:1116:8824:38450/1
HWI-1KL149:18:C0RNBACXX:1:1202:13982:33444/2
ZZZZ:X



$ java -jar dist/bamqueryreadnames.jar -b -s -N list.notfound input.bam read.names | samtools view |head

HWI-1KL149:18:C0RNBACXX:1:1101:12800:7069	77	*	0	0	*	*	0	0	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	##################################################	RG:Z:p1294	AS:i:0	XS:i:0
HWI-1KL149:18:C0RNBACXX:1:1101:12800:7069	141	*	0	0	*	*	0	0	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	##################################################	RG:Z:p1294	AS:i:0	XS:i:0
HWI-1KL149:18:C0RNBACXX:1:1101:15533:71685	99	X	238583	60	100M	=	23858972	274	(...)
(...)

$ cat list.notfound
ZZZZ:X
```


 

