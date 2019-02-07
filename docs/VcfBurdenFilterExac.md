# VcfBurdenFilterExac

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Burden filter 3 - Exac


## Usage

```
Usage: vcfburdenexac [options] Files
  Options:
    -d, --discardNotInExac
      if variant was not found in Exac, set the FILTER. Default: don't set the 
      FILTER. 
      Default: false
  * -exac, --exac
      Path to Exac VCF file. At the time of writing, you'd better use a 
      normalized version of Exac (see 
      https://github.com/lindenb/jvarkit/wiki/VCFFixIndels )
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -maxFreq, --maxFreq
      set FILTER if max(exac frequency in any pop) is greater than this value)
      Default: 0.001
    -o, --output
      Output file. Optional . Default: stdout
    -pop, --population
      comma separated populations in exac
      Default: AFR,AMR,EAS,FIN,NFE,SAS
    -tabix, --tabix
      use tabix index for Exac it is present. Might speed up things if the 
      number of variant is low.
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * exac


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfburdenexac
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFilterExac.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/VcfBurdenFilterExac.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfburdenexac** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


20170626: this tool now supports multiple ALT in the user VCF, however it's not been tested for choosing when to set the FILTER or the min value

### Output

#### INFO column


 *  FreqExac : Exac frequency.
 *  AC_* and AN_*: Transpose original population data from original Exac file


#### FILTER column

 *  BurdenExac : if FreqExac doesn't fit the criteria maxFreq


### see also


 *  VcfBurdenMAF


