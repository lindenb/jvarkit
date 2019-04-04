# FaidxSplitter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split bed of reference genome into overlapping parts


## Usage

```
Usage: faidxsplitter [options] Files
  Options:
    --exclude
      A regular expression to exclude some chromosomes from the dictionary.
      Default: (chr)?(M|MT)$
    -gap, --gaps, --gap
      gap file. A bed file containing the known gaps in the genome. E.g: 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz 
    -gene, --gene, --genes
      gene file. You shouldn't break a record in this file.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -x, --overlap
      Overlap  Size. The resulting BED region should overlap with 'x' bases.
      Default: 1000
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -w, --size
      BED Size. The genome should be split into BED of 'w' size.
      Default: 1000000
    -s, --small
      If it remains 's' bases in the BED split to the end of the chromosome, 
      extends the current BED.
      Default: 1000
    --version
      print version and exit

```


## Keywords

 * vcf
 * reference
 * bed


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew faidxsplitter
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FaidxSplitter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FaidxSplitter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/FaidxSplitterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/FaidxSplitterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **faidxsplitter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Motivation

For WGS sequencing, a tool I used to split the genome and parallelize things...

# Example

```
$ wget -q -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz" | gunzip -c | cut -f 2,3,4 > jeter.gaps.txt

$ wget -q -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz"  | gunzip -c | cut -f 3,5,6 > jeter.genes.txt

$ java -jar dist/faidxsplitter.jar -gap jeter.gaps.txt -gene jeter.genes.txt -R src/test/resources/human_b37.dict

1	10000	177417
1	227417	267719
1	317719	471368
1	521368	1521369
1	1520369	2634220
1	2684220	3689209
1	3688209	3845268
1	3995268	4995269
1	4994269	6241184
1	6240184	7830767
(...)

```

 

