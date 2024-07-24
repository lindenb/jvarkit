# SameDict

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

check if all HTS files share the same dictionary


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar samedict  [options] Files

Usage: samedict [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --verbose
      Be verbose
      Default: false
    --version
      print version and exit

```


## Keywords

 * dict
 * bed
 * sam
 * bam
 * vcf



## Creation Date

20240724

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samedict/SameDict.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/samedict/SameDict.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samedict** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Return value

return -1 if there is no argument on the command line, if a dictionary cannot
be extracted from BAM/SAM/CRAM/BCF/VCF/FASTA etc...

return 0 if all argument share the same dictionary


## Example

```
$ java -jar dist/jvarkit.jar samedict src/test/resources/S*.bam && echo "OK"
OK

$ java -jar dist/jvarkit.jar samedict || echo "ERROR"
[INFO][Launcher]samedict Exited with failure (-1)
ERROR

$ java -jar dist/jvarkit.jar samedict src/test/resources/S*.bam src/test/resources/ENCFF331CGL.rnaseq.b38.bam || echo "ERROR"
[INFO][Launcher]samedict Exited with failure (-1)
ERROR

$ java -jar dist/jvarkit.jar samedict src/test/resources/S*.bam src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz src/test/resources/rotavirus_rf.fa && echo OK
OK

```




