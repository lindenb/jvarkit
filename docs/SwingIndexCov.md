# SwingIndexCov

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

indexcov visualization


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar swingindexcov  [options] Files

Usage: swingindexcov [options] Files
  Options:
    --gff, --gff3
      GFF3 file indexed with tabix to plot the genes.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --jexl
      JEXL code for highlighting. JEXL stands for Java EXpression Language.  
      See https://commons.apache.org/proper/commons-jexl/reference/syntax.html
    -R, --reference
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    --roi, --bed
      BED file with regions of interest, to plot
    --version
      print version and exit

```


## Keywords

 * cnv
 * duplication
 * deletion
 * sv



## Creation Date

2020511

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/SwingIndexCov.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/indexcov/SwingIndexCov.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **swingindexcov** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Examples:

```
java -jar dist/swingindexcov.jar indexcov.bed.gz
```
or 
```
java -jar dist/swingindexcov.jar -R reference.fa indexcov.bed.gz --gff indexed.gff3.gz --helper "/TMP/myhelper.jar MyHelper"
```


## Helper

helper is based on jexl, it is used to colorize the peaks

the following JEX scrip show dot for the samples below

```jexl
var cases = ["SAMPLE1","SAMPLE2","SAMPLE3"];

for(var c:calls) {
        var is_case= cases.indexOf(c.getSample())!=-1;
        var is_var = c.isVariant();
        if(is_case && !is_var) return false;
        if(!is_case && is_var) return false;
        }

return true;
```


