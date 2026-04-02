# AddLinearIndexToBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Use a Sequence dictionary to create a linear index for a BED file. Can be used as a X-Axis for a chart.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar addlinearindextobed  [options] Files

Usage: addlinearindextobed [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --ignore
      skip unknown chromosomes.
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference, --dict
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    --regex
      keep chromosomes matching that regular expression. (use --ignore too 
      prevent error the other chromosomes)
      Default: (chr)?[0-9XY]+
    --version
      print version and exit

```


## Keywords

 * bed
 * reference



## Creation Date

20140201

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/AddLinearIndexToBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/AddLinearIndexToBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **addlinearindextobed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## EXAMPLE

```
$ cat input.bed | java -jar dist/addlinearindextobed.jar -R  human_g1k_v37.fasta 

10000   1       10000   177417
227417  1       227417  267719
317719  1       317719  471368
(...)
3060255274      Y       23951428        28819361
3095123207      Y       58819361        58917656
3095271502      Y       58967656        59363566
```


