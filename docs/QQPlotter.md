# QQPlotter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

plot QQplot


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar qqplotter  [options] Files

Usage: qqplotter [options] Files
  Options:
    --column-label
      column label for 'name of the point'
      Default: <empty string>
    -p, --column-pvalue
      column label for 'p-value'
      Default: P
    --column-url
      column label for 'url of the point'
      Default: <empty string>
    -d, --delim
      column separator as a java regexe
      Default: [ 	]+
    --disable-log10
      data are already -log10(x), so do not apply -log10
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -w, -width
      image width
      Default: 1000

```


## Keywords

 * qqplot
 * gwas
 * statistics



## Creation Date

20250324

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/qqplot/QQPlotter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/qqplot/QQPlotter.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **qqplotter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




# Example



