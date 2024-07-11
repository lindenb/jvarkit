# BigWigTView

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

view bigwig file coverage in a terminal


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bigwigtview  [options] Files

Usage: bigwigtview [options] Files
  Options:
    --height
      screen height
      Default: 20
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * --interval, --region
      interval chr:start-end.
      Default: <empty string>
    --max
      when binning data, use max value instead of average
      Default: false
    -o, --out
      Output is a setfile. Output file. Optional . Default: stdout
    --plain
      disable ansi
      Default: false
    --version
      print version and exit
    --width
      screen width
      Default: 100

```


## Keywords

 * wig
 * bigwig



## Creation Date

20240704

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bigwigtview/BigWigTView.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bigwigtview/BigWigTView.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bigwigtview** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

plot set of bigwig files in terminal

## EXAMPLES

## Input

intut is one or more bigwigfile or a file with the '.list' suffix containing the path to the bigwigs
 
### Example

```
$ java -jar dist/jvarkit.jar bigwigtview --plain -region "chr1:1014238-1014330" src/test/resources/Uniqueness35bp.bigWig 


> src/test/resources/Uniqueness35bp.bigWig chr1:1014238-1014330
  1.050000 |                                                                                                     
  0.997500 |                                                               ######### #########             ##### 
  0.945000 |                                                               ######### #########             ##### 
  0.892500 |                                                               ######### #########             ##### 
  0.840000 |                                                               ######### #########             ##### 
  0.787500 |                                                               ######### #########             ##### 
  0.735000 |                                                               ######### #########             ##### 
  0.682500 |                                                               ######### #########             ##### 
  0.630000 |                                                               ######### #########             ##### 
  0.577500 |                                                               ######### #########             ##### 
  0.525000 |                                                               ######### #########             ##### 
  0.472500 |                                ###                            ######### ############# ############# 
  0.420000 |                                ###                            ######### ############# ############# 
  0.367500 |                                ###                            ######### ############# ############# 
  0.315000 |                                ###              ######### ############# ############# ############# 
  0.262500 |                                ###              ######### ############# ############# ############# 
  0.210000 |  ##                            ########### ############## ############# ############# ############# 
  0.157500 |  ##                            ########### ############## ############# ############# ############# 
  0.105000 |  ##                            ########### ############## ############# ############# ############# 
  0.052500 |  ##                            ########### ############## ############# ############# ############# 
  0.000000 | ####################################################################################################
      chr1 | 1014238 1014245 1014252 1014260 1014267 1014275 1014282 1014290 1014297 1014304 1014312 1014319 
```



