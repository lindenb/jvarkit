# BigwigMerge

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

merge several Bigwig files using different descriptive statistics (mean, median, etc..)


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bigwigmerge  [options] Files

Usage: bigwigmerge [options] Files
  Options:
    --header
      write track header
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --interval, -r
       process only this interval
    --method, -m
       how to compute merge of wiggle values ?
      Default: median
      Possible Values: [median, average, count, sum, min, max, random]
    --min-item-count
       skip output data if thre is less than 'x' bigwig file at a location
      Default: 0
    --min-value
       skip output data with values lower than 'x'
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    --version
      print version and exit

```


## Keywords

 * wig
 * bigwig



## Creation Date

20240417

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bigwigmerge/BigwigMerge.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bigwigmerge/BigwigMerge.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bigwigmerge** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# motivation

merge several Bigwig files using different descriptive statistics (mean, median, etc..)

Output is a BedGraph file.

Input is a set of bigwig file or a file with the '.list' suffix containing the path to the bigwig


## Example

```
find DIR -type f -name "*.bigWig" > tmp.list
java -jar jvarkit.jar bigwigmerge -R genome.fa tmp.list --interval "chr1:234-567" --header --method median  > bedGraph.bed
```



