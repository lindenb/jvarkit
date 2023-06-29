# BarcodeGenerator

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Barcode generator for EricCharp


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar barcodegenerator  [options] Files

Usage: barcodegenerator [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --max-gc
      max gc
      Default: 0.8
    --max-generation
      max generations. -1 : forever
      Default: -1
    --min-differences
      min differences between barcodes
      Default: 4
    --min-gc
      min gc
      Default: 0.2
    --number
      number of barcode to generated
      Default: 96
  * -o, --out
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
    --polyx
      max polyx
      Default: 3
    --random
      random seed. -1 == current time
      Default: -1
    --size
      barcode size
      Default: 8
    --version
      print version and exit

```


## Keywords

 * barcode



## Creation Date

20230629

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/barcode/BarcodeGenerator.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/barcode/BarcodeGenerator.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **barcodegenerator** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


