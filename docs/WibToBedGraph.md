# WibToBedGraph

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extract Wib files to bedgraph or wig


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar wib2bedgraph  [options] Files

Usage: wib2bedgraph [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
  * -r, --regions
      An interval as the following syntax : "chrom:start-end" or 
      "chrom:middle+extend"  or "chrom:start-end+extend" or 
      "chrom:start-end+extend-percent%".A program might use a Reference 
      sequence to fix the chromosome name (e.g: 1->chr1)
      Default: <empty string>
  * --tabix
      A wib associated indexed tabix file.e.g: wget 
      "https://hgdownload.soe.ucsc.edu/gbdb/hg38/multiz100way/phastCons100way.wib" 
      && wget -O - "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/phastCons100way.txt.gz" 
      | gunzip -c | bgzip > phastCons100way.txt.gz && tabix --force -0 -s 2 -e 
      3 -e 4 phastCons100way.txt.gz

      Default: <empty string>
    --version
      print version and exit
  * --wib
      The wib file itself e.g: wget 
      'https://hgdownload.soe.ucsc.edu/gbdb/hg38/multiz100way/phastCons100way.wib' 
      . 
      Default: <empty string>
    --wig
      output WIG format instead of BED graph
      Default: false

```


## Keywords

 * wib
 * wig
 * bed



## Creation Date

20230819

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/wib/WibToBedGraph.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/wib/WibToBedGraph.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **wib2bedgraph** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Example:

```
$ java -jar dist/jvarkit.jar wib2bedgraph  --tabix ~/phastCons100way.txt.gz --wib ~/phastCons100way.wib -r "chr2:1-125168903" | head

chr2	11391	11392	0.14969291
chr2	11392	11393	0.14667717
chr2	11393	11394	0.14064567
chr2	11394	11395	0.13461417
chr2	11395	11396	0.12858267
chr2	11396	11397	0.13461417
chr2	11397	11398	0.14366141
chr2	11398	11399	0.14969291
chr2	11399	11400	0.1557244
chr2	11400	11401	0.15874015

$ java -jar  dist/jvarkit.jar wib2bedgraph  --tabix ~/phastCons100way.txt.gz --wib ~/phastCons100way.wib -r "chr2:1-125168903" --wig | head

variableStep chrom=chr2 span=1
11392	0.14969291
11393	0.14667717
11394	0.14064567
11395	0.13461417
11396	0.12858267
11397	0.13461417
11398	0.14366141
11399	0.14969291
11400	0.1557244
```
 

