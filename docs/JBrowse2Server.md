# JBrowse2Server

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

create a run a local instance of jbrowse2


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar jbrowse2  [options] Files

Usage: jbrowse2 [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --port
      server port.
      Default: 8080
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --version
      print version and exit
    --zip
      JBrowse2 archive source
      Default: https://github.com/GMOD/jbrowse-components/releases/download/v3.2.0/jbrowse-web-v3.2.0.zip

```


## Keywords

 * jbrowse
 * browser
 * genome



## Creation Date

20250404

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/jbrowse2/JBrowse2Server.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/jbrowse2/JBrowse2Server.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **jbrowse2** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Example

```
java -jar dist/jeter12.jar -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam
```


