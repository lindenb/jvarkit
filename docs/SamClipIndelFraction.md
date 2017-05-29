# SamClipIndelFraction

Extract clipping/indel fraction from BAM


## Usage

```
Usage: Samclipindelfraction [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -t
      type
      Default: allclip
      Possible Values: [leftclip, rightclip, allclip, insert, deletion]

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
$ make Samclipindelfraction
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFraction.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamClipIndelFraction.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **Samclipindelfraction** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ samtools view -h -F3844 my.bam  | java -jar dist/samclipindelfraction.jar 

##UNMAPPED_READS=0
##MAPPED_READS=3028359
##CLIPPED_READS=1182730
##CLIPPED_READS_5_PRIME=597757
##CLIPPED_READS_3_PRIME=617399
##UNCLIPPED_READS=1845629
##COUNT_BASES=338644685
#CLIP	COUNT	FRACTION_OF_MAPPED_READS
0	1845629	0.5
1	7	1.8963724562195327E-6
2	6756	0.0018302703306027376
3	695	1.8828269386751074E-4
4	794	2.1510281860547272E-4
5	819	2.2187557737768533E-4
6	471	1.275987752684857E-4
7	447	1.210969268471616E-4
(...)
```

plotting:
```bash
$ java -jar dist/samclipindelfraction.jar |\
   grep -v "##" | cut -f1,2 | tr -d '#' > output.txt
```

then, in R:
```R
T<-read.table('output.txt',header=TRUE)
plot(T[T$CLIP>0,])
```


