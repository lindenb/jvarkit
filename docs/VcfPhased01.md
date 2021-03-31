# VcfPhased01

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

X10 Phased SVG to Scalar Vector Graphics (SVG)


## Usage

```
Usage: vcfphased01 [options] Files
  Options:
    -xp, --xpos, --extra-highligth
      Extra Highligth positions that are not always in the vcfs. (comma 
      separated) 
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -p, --pos, --highligth
      Highligth positions in the VCFs. (comma separated)
      Default: <empty string>
  * -r, --interval, --region
      interval CHROM:START-END
    -k, --knownGenes
      UCSC knownGene File/URL. The knowGene format is a compact alternative to 
      GFF/GTF because one transcript is described using only one line.	Beware 
      chromosome names are formatted the same as your REFERENCE. A typical 
      KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
      Default: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * x10
 * phased
 * genotypes
 * svg



## See also in Biostars

 * [https://www.biostars.org/p/9462569](https://www.biostars.org/p/9462569)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfphased01
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190710

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/phased/VcfPhased01.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/phased/VcfPhased01.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfphased01** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Motivation

For @154sns , 10 July 2019

Displays phased genotypes from X10 genomics into SVG.

# Input

input is a list of indexed vcf files or a file with the suffix '.list' containing the full path to the vcfs

# Example

```
find . -type f -name "*.vcf.gz" > in.list
java -jar dist/vcfphased01.jar -r "chr1:1000-2000" -xp '1001,1010' in.list 
```

# Screenshots

https://twitter.com/yokofakun/status/1148964221482414080

![https://pbs.twimg.com/media/D_Hwd2dXoAAzx8g.jpg](https://pbs.twimg.com/media/D_Hwd2dXoAAzx8g.jpg)

