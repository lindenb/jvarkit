# FastqToFasta

fastq -> fasta


## DEPRECATED

use awk, samtools...

## Usage

```
Usage: fastq2fasta [options] Files
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
    -N
      fasta line length
      Default: 50
    -b
      trim fasta header after space
      Default: false

```


## Keywords

 * fastq
 * fasta


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
$ make fastq2fasta
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FastqToFasta.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FastqToFasta.java)


<details>
<summary>Git History</summary>

```
Mon May 15 10:41:51 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/c13a658b2ed3bc5dd6ade57190e1dab05bf70612
Mon Apr 24 17:49:35 2017 +0200 ; cont jcommander ; https://github.com/lindenb/jvarkit/commit/d822a90a1eaba26a4d874472ccd45e689e8ba063
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Thu Feb 27 17:10:54 2014 +0100 ; cont, fix bug in bam2fastq, shortread, starting change-ref bam, extract clipped seq ; https://github.com/lindenb/jvarkit/commit/d83138c95883cf87078565b54614b2aa7aa04740
Sun Feb 2 18:55:03 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/abd24b56ec986dada1e5162be5bbd0dac0c2d57c
Thu Nov 28 14:54:21 2013 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/6bd741fe898f5d735e5ada6b59222f8818c08baf
Wed Nov 27 20:00:16 2013 +0100 ; abstract bam filter ; https://github.com/lindenb/jvarkit/commit/6da95f7c2f27ea15634c8f3504cdc71495020248
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastq2fasta** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example

```bash
$ java -jar dist/fastq2fasta.jar  -N 60 -b file.fastq.gz

>HWI-1KL149:61:D2C11TCXX:2:1213:4591:29626
GAGTTGCTTTGTTTGAATATAGGTTGACTATACGAAGTGTGCGAGGACCTGCACCACGCA
GTAGGCCAAGATCAACTGAAACAGTGCTATCTGCACGACAA
>HWI-1KL149:61:D2C11TCXX:2:1213:4525:29650
CCTAGTAGTTCGTGGCCCCGGGCCCCTACTTAAACTCCTAGAACCACTCCTAGAAAGGGG
TGTTGCAGTTCGGCTCAGTCCCCGTGGTCGACTACTGTTTC
>HWI-1KL149:61:D2C11TCXX:2:1213:4569:29706
GCGCAGAGTTGTTTTAGCTATGCTGTGTTTGCATGGTTAGGTGGTGTACCTAGTGGTTTT
CTGAGACTTCTCTGAGGTTCTTGAGTAGATTAATACATCCC
>HWI-1KL149:61:D2C11TCXX:2:1213:4594:29713

```



