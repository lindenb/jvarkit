# HaploGroupCaseControl

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Run Fisher test for Haplogroup input.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar haplogroupcasectrl  [options] Files

Usage: haplogroupcasectrl [options] Files
  Options:
    --cases
      File or comma-separated list of control samples
    --controls
      File or comma-separated list of control samples
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --phylotree, --tree
      Path/URL to a valid XML formatted phylotree ( 
      https://github.com/genepi/phylotree-rcrs-17 ).
      Default: https://raw.githubusercontent.com/genepi/phylotree-rcrs-17/main/src/tree.xml
    --version
      print version and exit

```


## Keywords

 * haplogroup
 * burden
 * mitochondrial
 * 



## Creation Date

20240610

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/haplogroupcasectrl/HaploGroupCaseControl.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/haplogroupcasectrl/HaploGroupCaseControl.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **haplogroupcasectrl** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Examples

```
 java -jar jvarkit.jar haplogroupcasectrl --cases input.cases  --controls input.ctrls  --tree input.xml   input.tsv

```



### Examples




