# VcfSplitBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split VCF by BED/Region Name


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfsplitbed  [options] Files

Usage: vcfsplitbed [options] Files
  Options:
  * -B, --bed
      User bed fourth column is used as the key to aggregate the variants. 
      Empty is replace with {chrom}_{start}_{end}
    --disable-hash-directory, --dhd
      disable default which is to save each file in a checksum-based 
      directory-a-la-nextflow to avoid a large number of files in the same 
      directory. 
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -o, --output
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
    --prefix
      prefix each output VCF file with this string
      Default: <empty string>
    --version
      print version and exit

```


## Keywords

 * vcf
 * sliding
 * window



## Creation Date

20190619

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfsplitbed/VcfSplitBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfsplitbed/VcfSplitBed.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfsplitbed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


### Example

```
$ java -jar dist/jvarkit.jar vcfsplitbed --bed input.ned src/test/resources/rotavirus_rf.vcf.gz 
```

# see also

 * vcfgenesplitter


