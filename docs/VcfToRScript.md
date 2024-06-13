# VcfToRScript

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert VCF to R so it can be used for burden testing


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcf2r  [options] Files

Usage: vcf2r [options] Files
  Options:
    --cases
      File or comma-separated list of control samples
    --controls
      File or comma-separated list of control samples
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --info
      comma separated of extra INFO fields included in the variant dataframe. 
      Only with --variants
      Default: <empty string>
    -o, --output
      Output file. Optional . Default: stdout
    --variants
      include a table with variants
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * fisher
 * R


## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf2r/VcfToRScript.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcf2r/VcfToRScript.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2r** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





