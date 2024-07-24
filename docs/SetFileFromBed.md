# SetFileFromBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert bed chrom/start/end/name sorted on 4th column to set file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar setfilefrombed  [options] Files

Usage: setfilefrombed [options] Files
  Options:
    --disable-interval-merge
      Do not merge overlapping intervals in a setFile record
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout. For action=cluster, output is: 
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    -t, --trim-chr
      Remove chr prefix in chromosome names on output.
      Default: false
    --version
      print version and exit

```


## Keywords

 * setfile
 * bed



## Creation Date

20210125

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/setfile/SetFileFromBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/setfile/SetFileFromBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **setfilefrombed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


```
$ echo -e "RF01\t150\t200\tA\nRF01\t190\t300\tA\nRF01\t350\t400\tA\nRF01\t150\t200\tB" |\
	java -jar dist/jvarkit.jar setfilefrombed -R src/test/resources/rotavirus_rf.fa frombed

A	RF01:151-300,RF01:351-400
B	RF01:151-200

```



