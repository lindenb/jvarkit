# IGVServer

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Web Server displaying SAM/BAM/VCF file and using the IGV.js API. 


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar igvserver  [options] Files

Usage: igvserver [options] Files
  Options:
    --bed, --roi
      Bed of regions of interest displated as a html 'select' item
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --igv-js-url
      Where to download IGV js
      Default: https://cdn.jsdelivr.net/npm/igv@3.8.10/dist/igv.min.js
    -P, --port, -port
      Server listening port
      Default: 8080
    --shutdown-after
      Stop the server after 'x' seconds.
      Default: -1
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * table
 * visualization
 * server
 * web



## Creation Date

20261104

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/igvserver/IGVServer.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/igvserver/IGVServer.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **igvserver** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a set of bam/vcf/cram/fasta file. A fasta file compatible with the bam/vcf is required.


## Example:

```
 $ java -jar dist/jvarkit.jar igvserver --bed roi.bed src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam src/test/resources/rotavirus_rf.vcf.gz
 $ java -jar dist/jvarkit.jar igvserver igv_sessions1.xml igv_sessions2.xml 
 
$ find /path/to -type  f -name "*.cram" -o -name "*.fa" > input.list
$ java -jar dist/jvarkit.jar igvserver input.list 
 
```


