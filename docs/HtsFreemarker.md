# HtsFreemarker

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Apply Freemarker to VCF/BAM/JSON files.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar htsfreemarker  [options] Files

Usage: htsfreemarker [options] Files
  Options:
    --bam, --sam
      <name> <bam-file>. Add this Bam to the freemarker context. Object 
      created is [header:object,reads:list]
      Default: []
    --freemarker-version
      FreeMarker version.
      Default: 2.3.32
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --interval
      Restrict BAM/VCF/... file(s) to that interval. Files must be indexed.
    --json
      <name> <json-file> . Add this pair to the freemarker context.
      Default: []
    --json-string
      <name> <json-string> . Add this pair to the freemarker context.
      Default: []
    --output, -o
      Output file. Optional . Default: stdout
    --reference, -R
      For reading CRAM.Indexed fasta Reference file. This file must be indexed 
      with samtools faidx and with picard CreateSequenceDictionary
    --string
      <name> <string>. Add this pair to the freemarker context.
      Default: []
    --templates-dir, -T
      Templates base directory
      Default: .
    --vcf
      <name> <vcf-file>. Add this VCF to the freemarker context. Object 
      created is [header:object,variants:list]
      Default: []
    --version
      print version and exit

```


## Keywords

 * template
 * freemarker
 * vcf
 * bam



## Creation Date

20230616

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/htsfreemarker/HtsFreemarker.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/htsfreemarker/HtsFreemarker.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **htsfreemarker** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


todo


