# PadEmptyFastq

Pad empty fastq sequence/qual with N/#


## DEPRECATED

use awk

## Usage

```
Usage: pademptyfastq [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -N
      number of bases/qual to be added.  -1=length of the first read
      Default: -1

```


## Keywords

 * fastq


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
$ make pademptyfastq
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/PadEmptyFastq.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/PadEmptyFastq.java)


<details>
<summary>Git History</summary>

```
Mon May 15 17:17:02 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/fc77d9c9088e4bc4c0033948eafb0d8e592f13fe
Tue Apr 25 15:40:45 2017 +0200 ; cont jcommander ; https://github.com/lindenb/jvarkit/commit/16aeb209fda502b60dd75689b85d1304f469775b
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Tue Dec 10 22:12:54 2013 +0100 ; vcf head tail ; https://github.com/lindenb/jvarkit/commit/25adbeaf3de1c64fd793c5ddf5236322adb430d1
Thu Nov 28 08:16:28 2013 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/d41deb4c340967592eb53e98101077ccbd84a3dd
Fri Nov 22 21:13:00 2013 +0100 ; pad fastq ; https://github.com/lindenb/jvarkit/commit/0d712b792c3c3ede6c78312611597a6dba893fab
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **pademptyfastq** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### AWK

compiled version of awk:

```awk
{if(NR%4==2 && length($0)==0) { printf("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n");} else if(NR%4==0 && length($0)==0) { printf("##################################################\n");} else {print;}}
```


### Example

```bash
$ cutadapt -a AGATCGGAAGAGCGTCGT  2> /dev/null in.fastq.gz |\
  java -jar dist/pademptyfastq.jar -o pad.fastq.gz
```


