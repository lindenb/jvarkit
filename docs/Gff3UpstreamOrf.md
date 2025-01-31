# Gff3UpstreamOrf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Takes a ucsc genpred file, scan the 5' UTRs and generate a GFF3 containing upstream-ORF. Inspired from https://github.com/ImperialCardioGenetics/uORFs 


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar gff3upstreamorf  [options] Files

Usage: gff3upstreamorf [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
  * -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --strength
      only accept events that are greater or equal to this Kozak strength.
      Default: nil
      Possible Values: [Strong, Moderate, Weak, nil]
    --version
      print version and exit

```


## Keywords

 * gff
 * gff3
 * uorf
 * uorf



## Creation Date

20220724

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/upstreamorf/Gff3UpstreamOrf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/upstreamorf/Gff3UpstreamOrf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **gff3upstreamorf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## inspiration

Part of this code was inspired from: https://github.com/ImperialCardioGenetics/uORFs/blob/master/5primeUTRannotator/five_prime_UTR_annotator.pm

Wikipedia:

> An Upstream Open Reading Frame (uORF) is an open reading frame (ORF) within the 5' untranslated region (5'UTR) of an mRNA. uORFs can regulate eukaryotic gene expression.
> Translation of the uORF typically inhibits downstream expression of the primary ORF. In bacteria, uORFs are called leader peptides, and were originally discovered on the basis of their impact on the regulation of genes involved in the synthesis or transport of amino acids. 


Input is a UCSC "genpred/knowngene" file, but if you only have a gff/gff3 file, you can use `gff2kg` to create one.

## Examples


```
wget -O - "https://hgdownload.soe.ucsc.edu/goldenPath/hg38//database/wgEncodeGencodeBasicV47.txt.gz" |\
	gunzip -c |\
	java -jar dist/jvarkit.jar gff3upstreamorf  -R GRCh38.fa  > uorf.gff3

bcftools csq --ncsq 1000 -l -f GRCh38.fa -g uorf.gff3 in.bcf |\
	bcftools annotate --rename-annots <(echo -e 'INFO/BCSQ\tUTR_BCSQ')

	
```

note to self: test ENSG00000141736 https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003529


