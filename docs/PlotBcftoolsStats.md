# PlotBcftoolsStats

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Plot bcftools stats output


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar plotbcftoolsstats  [options] Files

Usage: plotbcftoolsstats [options] Files
  Options:
    --categories, --phenotypes
      A tab delimited file (sample)<tab>(category)
    --format
      output format.
      Default: PDF
      Possible Values: [PDF, PNG, SVG]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --prefix
      output prefix.
      Default: <empty string>
    --sections
      Limit to those sections. eg. 'QUAL,PSC,PSI'. Comma-separated values. 
      Ignore if empty. Inverse the selection if it starts with '^'.
      Default: <empty string>
    --version
      print version and exit
    -D
      set some css style elements. '-Dkey=value'. Undocumented.
      Syntax: -Dkey=value
      Default: {}

```


## Keywords

 * bcftools
 * qc
 * vcf



## Creation Date

20210622

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bcftools/PlotBcftoolsStats.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bcftools/PlotBcftoolsStats.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **plotbcftoolsstats** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Warning

this tool is not stable

## Example

todo


