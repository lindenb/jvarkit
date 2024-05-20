# BedRenameChromosomes

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert the names of the chromosomes in a Bed file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bedrenamechr  [options] Files

Usage: bedrenamechr [options] Files
  Options:
    -c, --column
      1-based chromosome column(s), multiple separated by commas
      Default: 1
    -convert, --convert
      What should I do when  a converstion is not found
      Default: RAISE_EXCEPTION
      Possible Values: [RAISE_EXCEPTION, SKIP, RETURN_ORIGINAL]
    -d, --delim
      field delimiter.
      Default: 	
    -s, --header
      Ignore lines starting with this java regular expression but print them 
      anyway 
      Default: (#|browser|track)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -f, --mapping, -m, -R
      load a custom name mapping.Chromosome mapping file. If the file looks 
      like a NGS file (fasta, vcf, bam...) the mapping is extracted from a 
      dictionary; Otherwise, it is interpreted as a mapping file ( See 
      https://github.com/dpryan79/ChromosomeMappings )
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * bed
 * chromosome
 * contig
 * convert



## Creation Date

20190503

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedrenamechr/BedRenameChromosomes.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedrenamechr/BedRenameChromosomes.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedrenamechr** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


