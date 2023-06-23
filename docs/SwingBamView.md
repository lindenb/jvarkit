# SwingBamView

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Read viewer using Java Swing UI


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar swingbamview  [options] Files

Usage: swingbamview [options] Files
  Options:
    --gff, --gff3
      Tabix-Indexed GFF3 File
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    -V, --variant
      Indexed VCF File
    --version
      print version and exit

```


## Keywords

 * bam
 * alignment
 * graphics
 * visualization
 * swing



## Creation Date

20220503

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfviewgui/SwingBamView.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfviewgui/SwingBamView.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **swingbamview** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example:

```
java -jar dist/jvarkit.jar swingbamview -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam
```

```
find dir -type f -name "*.bam" > out.list
java -jar dist/jvarkit.jar swingbamview -R src/test/resources/rotavirus_rf.fa out.list
```


