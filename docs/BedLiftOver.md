# BedLiftOver

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

LiftOver a BED file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bedliftover  [options] Files

Usage: bedliftover [options] Files
  Options:
  * -f, --chain
      LiftOver chain file. Can be a local chain file, a URL 'https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToCriGri1.over.chain.gz', 
      or a chain identifier like 'hg19ToHg38'.
    --chainvalid
      Ignore LiftOver chain validation
      Default: false
    -x, --failed
        write bed failing the liftOver here. Optional.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --minmatch
      lift over min-match.
      Default: 0.95
    --original, --src
      Append original interval
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
  * -D, -R, -r, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --version
      print version and exit

```


## Keywords

 * bed
 * liftover



## Creation Date

20140311

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/BedLiftOver.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/liftover/BedLiftOver.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedliftover** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example


```
cat in.bed | java -jar jvarkit.jar bedliftover --chain x.chain -R ref.fa
```



