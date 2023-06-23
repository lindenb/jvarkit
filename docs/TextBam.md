# TextBam

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Write text in a bam. Mostly for fun...


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar texbam  [options] Files

Usage: texbam [options] Files
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
    -p, --pos
      use this position instead of a random one. Syntax: CHROM:POS
      Default: <empty string>
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --snp
      use SNV instead of deletion
      Default: false
    --version
      print version and exit

```


## Keywords

 * fun
 * bam
 * sam
 * txt



## Creation Date

20220708

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/textbam/TextBam.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/textbam/TextBam.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **texbam** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


```
$ ./gradlew textbam && \
	java -jar dist/textbam.jar -R src/test/resources/rotavirus_rf.fa -p 'RF01:100' "Hello world" | samtools view -O BAM -o jeter.bam &&  \
	samtools index jeter.bam && \
	samtools tview ~/jeter.bam src/test/resources/rotavirus_rf.fa


samtools tview ~/jeter.bam src/test/resources/rotavirus_rf.fa -d T -p RF01:100
 101       111       121       131       141       151       161                
TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATC
................................................................................
................................................................................
 ...............................................................................
  ..............*....*......******......*...........*............*****..........
   .............*....*......*...........*...........*............*....*.........
    ............*....*......*...........*...........*...........*.....*.........
     ...........******......****........*...........*...........*.....*.........
      ..........*....*......*...........*...........*...........*.....*.........
       .........*....*......*...........*...........*...........*....*..........
        ........*....*......*...........*...........*............**..*..........
         .......*....*......******......******......******.........**...........
          ......................................................................
           .....................................................................
```


