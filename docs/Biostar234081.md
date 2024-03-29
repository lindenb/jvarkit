# Biostar234081

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert extended CIGAR to regular CIGAR ('X','=' -> 'M')


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar biostar234081  [options] Files

Usage: biostar234081 [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      For reading/writing CRAM files. Indexed fasta Reference file. This file 
      must be indexed with samtools faidx and with picard/gatk 
      CreateSequenceDictionary or samtools dict
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * cigar



## See also in Biostars

 * [https://www.biostars.org/p/234081](https://www.biostars.org/p/234081)



## Creation Date

20170130

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar234081.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar234081.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar234081Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar234081Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar234081** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
 $ cat toy.sam 
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
r001    163 ref 7   30  1M2X5=4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG*XX:B:S,12561,2,20,112

 $ java -jar dist/biostar234081.jar toy.sam 
@HD VN:1.5  SO:unsorted
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG*XX:B:S,12561,2,20,112
```


