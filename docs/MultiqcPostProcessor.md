# MultiqcPostProcessor

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Enhances multiqc output by reading the data folder and producing new plots (eg. boxplot per population.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar multiqcpostproc  [options] Files

Usage: multiqcpostproc [options] Files
  Options:
    --beeswarm
      use plot_type=beeswarm instead of boxplot.
      Default: false
    --custom
      custom mapping file (undocumented
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -o, --output
      output directory
    --sample2collection
      tab delimited file containing (sample-name)(TAB)(collection-name). Empty 
      lines or starting with '#' are skipped
    --version
      print version and exit

```


## Keywords

 * multiqc



## Creation Date

20240708

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/multiqc/MultiqcPostProcessor.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/multiqc/MultiqcPostProcessor.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **multiqcpostproc** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

Enhances multiqc output by reading the data folder and producing new plots (eg. boxplot per population

If no group is defined and the tool can find the file `dragen_ploidy.txt`, this file is used to create a group with male and females

Doesn't work for now ? https://github.com/MultiQC/MultiQC/issues/2689

## Example

```
# run first time
multiqc --force --file-list input.list
java -jar jeter.jar multiqc_data -o OUTDIR2
find OUTDIR2 --type f -name "*.json" >> input.list
# run second time
multiqc --force --file-list input.list
```

