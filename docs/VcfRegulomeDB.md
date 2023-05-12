# VcfRegulomeDB

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate a VCF with the Regulome2 data (https://regulomedb.org/)


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfregulomedb  [options] Files

Usage: vcfregulomedb [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -x, --extends
      (int) base pairs. look.for data around the variation +/- 'x'. A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 0
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    -r, --ranking-regex
      if defined, only accept the rank matching the regular expression. see 
      https://regulomedb.org/regulome-help/ . For example: 1a	eQTL/caQTL + TF 
      binding + matched TF motif + matched Footprint + chromatin accessibility 
      peak 
      Default: <empty string>
  * -b, --bed, --tabix, --regulomedb
      RegulomeDB bed sorted, bgzipped and indexed with tabix.
    --version
      print version and exit

```


## Keywords

 * vcf
 * regulomedb



## Creation Date

20140709

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfregulomedb/VcfRegulomeDB.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfregulomedb/VcfRegulomeDB.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfregulomedb** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


Build the database, for grch38


```
# here, we use `head` to get a short example
$ wget -q  -O - "https://encode-public.s3.amazonaws.com/2023/03/03/d38f3202-3364-415d-86ed-8690330cb7a2/ENCFF250UJY.tsv" |\
	head -n 1000 | sed 's/^chrom/#chrom/' > regulome.bed
$ bgzip -f  regulome.bed 
$ tabix -f -p bed regulome.bed.gz 
```

for grch37, the database is available at http://legacy.regulomedb.org/downloads/RegulomeDB.dbSNP141.txt.gz .

## Example

```bash
$ wget -q -O - "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz" |\
	bcftools view --types snps |\
	java -jar dist/jvarkit.jar vcfregulomedb --regulomedb regulome.bed.gz |\
	bcftools query -i 'REGULOMEDB>0' -f '%CHROM\t%POS\t%REF\t%ALT\t%REGULOMEDB\n' 
chr1	10177	A	C	0.829
chr1	10177	A	G	0.829
chr1	10181	A	C	0.342
chr1	10181	A	G	0.342
chr1	10181	A	T	0.342
chr1	10248	A	C	0.609
chr1	10248	A	G	0.609
chr1	10248	A	T	0.609
chr1	10250	A	C	0.609
chr1	10250	A	T	0.609
chr1	10255	A	T	0.609
chr1	10257	A	C	0.609
chr1	10327	T	A	0.796
chr1	10327	T	C	0.796
chr1	10327	T	G	0.796


```



