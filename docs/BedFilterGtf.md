# BedFilterGtf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filter/Annotate BED using a GTF file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar bedgtf  [options] Files

Usage: bedgtf [options] Files
  Options:
    --bed
      only use GTF lines overlaping that bed file. A Bed file: 
      (CHROM)<tab>(START 0-based)<tab>(END)[<tab>otherfields...].
    --biotype
      What feature biotype 'x'  should be used.
      Default: all
      Possible Values: [all, protein_coding]
    --cut
      The GTF feature should not enclose the variant. This is interesting for 
      the SV like the DUP or the INV where  we want that the variant 'cut' a 
      genomic feature. Features where the whole Variant is contained will be 
      ignored. 
      Default: false
    --extend, --extends, -x, -X
      extends each gtf feature by 'x' bases. Multiple separated by commas. 
      Empty will be replaced by '0'. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 0
    --gene-id
      Optional comma delimited string or File with the '.txt' suffix 
      containing a list of gene IDs . Filter on attribute gene_id. Trailing 
      '.[0-9]+' for version will be removed
      Default: <empty string>
    --gene-name
      Optional comma delimited string or File with the '.txt' suffix 
      containing a list of gene names . Filter on attribute gene_name
      Default: <empty string>
  * -gtf, --gtf
      GTF file
    --header-type
      Input bed header type. 0 : auto detect ; 1 : the first line is always 
      the header; 2: there is no header.
      Default: 0
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hide-no-overlap, --hide
      Don't show variant that have NO intersection with any of the GTF 
      features 
      Default: false
    --no-biotype
      Do not insert gene biotypes in the output
      Default: false
    --no-gene-id
      Do not insert gene ids in the output
      Default: false
    --no-gene-name
      Do not insert gene names in the output
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference
      A SAM Sequence dictionary source: it can be a *.dict file, a fasta file 
      indexed with 'picard CreateSequenceDictionary' or 'samtools dict', or 
      any hts file containing a dictionary (VCF, BAM, CRAM, intervals...)
    --show-protein-id
      Do show protein IDs in the output
      Default: false
    --show-transcript-id
      Do show transcript IDs in the output
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * gtf



## Creation Date

20261113

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedgtf/BedFilterGtf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bedgtf/BedFilterGtf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bedgtf/BedFilterGtfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bedgtf/BedFilterGtfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bedgtf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



## Motivation

filter/annotate a BED with a GTF. Very same idea than  `vcfgtf`


## Example

Filter to protein-coding genes only

```
cat in.bed |\
java -jar jvarkit.jar bedgtf \
    --gtf annotations.gtf.gz \
    --biotype protein_coding \
    --hide-no-overlap \
    --extend 0,1kb \
    -o annotated.bed
```

Annotate structural variants with cut-mode

```
cat in.bed |\
java -jar jvarkit.jar bedgtf \
    --gtf gencode.v44.gtf \
    --cut \
    --extend 0,100 \
    --show-transcript-id \
    --show-protein-id \
    -o sv_annotated.vcf
```

Filter to a subset of genes via file

```
# genes.txt contains one gene name per line
cat in.bed |\
java -jar jvarkit.jar bedgtf \
    --gtf ref.gtf \
    --gene-name genes.txt \
    --hide-no-overlap \
    -o filtered.vcf
```




