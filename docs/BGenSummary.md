# BGenSummary

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

bgen file summary


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bgensummary  [options] Files

Usage: bgensummary [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -g
      print the genotypes for 'g' genotypes
      Default: 3
    -n
      print the summary for 'n' variants
      Default: 10

```


## Keywords

 * bgen



## Creation Date

20250507

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bgen/summary/BGenSummary.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bgen/summary/BGenSummary.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bgensummary** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/jvarkit.jar bgensummary *.bgen | xmllint --format -
$ java -jar dist/jvarkit.jar bgensummary ~/src/jbgen/src/test/resources/htsjdk/tribble/bgen/*bgen | xmllint --format -
<?xml version="1.0" encoding="UTF-8"?>
<bgen-summary>
  <bgen-file filename=complex.bgen" compression="e_ZlibCompression" n-variants="10" n-samples="4" layout="e_Layout2" snps-offset="68" anonymous="false">
    <samples n-samples="4">
      <samples idx="0">sample_0</samples>
      <samples idx="1">sample_1</samples>
      <samples idx="2">sample_2</samples>
      <samples idx="3">sample_3</samples>
    </samples>
  </bgen-file>
  <bgen-file filename="example.v11.bgen" compression="e_ZlibCompression" n-variants="199" n-samples="500" layout="e_Layout1" snps-offset="5920" anonymous="true"/>
  <bgen-file filename="haplotypes.bgen" compression="e_ZlibCompression" n-variants="4" n-samples="4" layout="e_Layout2" snps-offset="68" anonymous="false">
    <samples n-samples="4">
      <samples idx="0">sample_0</samples>
      <samples idx="1">sample_1</samples>
      <samples idx="2">sample_2</samples>
      <samples idx="3">sample_3</samples>
    </samples>
  </bgen-file>
</bgen-summary>

```


