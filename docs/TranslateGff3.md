# TranslateGff3

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

translates the output of bcftools consensus


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar translategff3  [options] Files

Usage: translategff3 [options] Files
  Options:
  * --gff, --gff3
      Gff3 file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * gff
 * gff3
 * fasta
 * genpred
 * peptide
 * protein



## Creation Date

20241003

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/translategff3/TranslateGff3.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/translategff3/TranslateGff3.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **translategff3** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Experimental.
Translates the output of bcftools consensus.
Output is a BED file containing the cDNA, the translated peptide and the mean codon usage .

## Example
```
samtools faidx ref.fa 8:11870-11890 |\
	bcftools consensus in.vcf.gz  -s SAMPLE -H A --include 'TYPE="snp"' |\
	java -jar dist/jvarkit.jar translategff3 --gff3 annot.gff3
```

