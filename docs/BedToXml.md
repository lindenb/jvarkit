# BedToXml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert BED to XML


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bed2xml  [options] Files

Usage: bed2xml [options] Files
  Options:
    --columns
      Optional Comma separated list names of column starting from the 3rd 
      column (after 'end'). Use '.' to ignore.
      Default: <empty string>
    --dict-out
      Save reduced/modified dict to this file
    -d, --distance
      if >=0, add a 'y' attribute that could be used to display the bed 
      records in a browser, 	this 'y' is the graphical row where the item 
      should be displayed. 	This distance is the distance between two item 
      where there is a collision. Memory consuming . A distance specified as a 
      positive integer.Commas are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb,g,gb
      Default: -1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --min-contig-length
      keep chromosomes which length is greater than 'x'
      Default: 0
    --omit-xml-desclaration
      Don't print XM declaration.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    --regex
      keep chromosomes matching that regular expression
    --type
      BedType. bedGraph will print the min/max of value. bed12 will display 
      each block
      Default: generic
      Possible Values: [generic, bedGraph, bed12]
    --version
      print version and exit

```


## Keywords

 * bed
 * xml



## Creation Date

20251128

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bed2xml/BedToXml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bed2xml/BedToXml.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bed2xml/BedToXmlTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bed2xml/BedToXmlTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bed2xml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


