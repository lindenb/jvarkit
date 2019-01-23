# GenScan

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Paint a Genome Scan picture from a Tab delimited file (CHROM/POS/VALUE1/VALUE2/....).


## Usage

```
Usage: genscan [options] Files
  Options:
    -dbc, --distancebetweencontigs
      number of pixels between contigs
      Default: 1.0
    -H, --height
      Image height
      Default: 700
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -maxy, --maxy, --ymax, -ymax
      max Y value
      Default: 100.0
    -miny, --miny, --ymin, -ymin
      min Y value
      Default: 0.0
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -r, --region
      One or more Specific region to observe. empty string is whole genome. Or 
      the name of a chromosome. Or an interval.
      Default: []
    --removeContigsSmallerThan
      When displaying a whole reference, don't use the configs having a length 
      lower than this number. Useful to only display the main chromosomes.
      Default: 0
    -style, --style
      Default style
      Default: <empty string>
    -track, --track
      Add a track by specifying it's name
      Default: []
    --version
      print version and exit
    -W, --width
      Image width
      Default: 1000

```


## Keywords

 * chromosome
 * reference
 * chart
 * visualization


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew genscan
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/genscan/GenScan.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/genscan/GenScan.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **genscan** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## INPUT



Input consists in 3 fields delimited with a tabulation. An optional 4th column can be used to set a specific style for the point

```
(CHROM)<tab>(POS)<tab>(VALUE)(<tab>STYLE)?
```

### directives

* `#!push` push a copy of the current style on the stack
* `#!pop` pop the current style off the stack

## style




## Example

```
$ samtools depth in.bam |\
		java -jar dist/genscan.jar  --removeContigsSmallerThan 500000 --width 2000 -R  human_g1k_v37.fasta -o out.png
```	

```
(echo "#!track:S1;"; samtools depth S1.bam ;echo "#!track:S2;"; samtools depth S2.bam ) | java -jar dist/genscan.jar --track "S1,S2" --width 2000 --removeContigsSmallerThan 500000 -R  human_g1k_v37.fasta -o out.png
```


