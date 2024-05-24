# VcfCadd

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate VCF with  Combined Annotation Dependent Depletion (CADD) (Kircher & al. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892.PubMed PMID: 24487276.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfcadd  [options] Files

Usage: vcfcadd [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    -d, --buffer-size
      Buffer size / processing window size. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 1000
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
  * -u, --uri, --tabix, --cadd
      Combined Annotation Dependent Depletion (CADD) Tabix file URI
      Default: <empty string>
    --version
      print version and exit
    -na
      value  used for 'allele-not-found'.
      Default: -999.0

```


## Keywords

 * vcf
 * prediction
 * cadd
 * annotation



## Creation Date

20220119

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cadd/VcfCadd.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cadd/VcfCadd.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/cadd/VcfCaddTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/cadd/VcfCaddTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcadd** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

```bash
$ java -Dhttp.proxyHost=my.proxy.host.fr -Dhttp.proxyPort=1234 -jar dist/jvarkit.jar vcfcadd \
	-u "http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz"  \
	src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz 2> /dev/null | ~/package/bcftools/bcftools annotate -x '^INFO/CADD_SCORE,INFO/CADD_PHRED'

##fileformat=VCFv4.2
(...)
##INFO=<ID=CADD_PHRED,Number=A,Type=Float,Description="PHRED expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD scores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc.  URI was http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz">
##INFO=<ID=CADD_SCORE,Number=A,Type=Float,Description="Score suggests that that variant is likely to be  observed (negative values) vs simulated(positive values).However, raw values do have relative meaning, with higher values indicating that a variant is more likely to be simulated (or -not observed-) and therefore more likely to have deleterious effects. URI was http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz">
##VcfCaddCmdLine=-u http://krishna.gs.washington.edu/download/CADD/v1.3/1000G_phase3.tsv.gz src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	905606	rs540662886	G	C,A	41743.9	PASS	CADD_PHRED=3.426,.;CADD_SCORE=0.082875,.
(...)
1	905621	rs368876607	G	A	14291.5	PASS	CADD_PHRED=6.025;CADD_SCORE=0.334762
(...)
1	905669	rs111483874	C	G,T	86574.3	PASS	CADD_PHRED=12.77,.;CADD_SCORE=1.39614,.
(...)
1	905723	rs150703609	G	A	15622.1	PASS	CADD_PHRED=23.7;CADD_SCORE=4.05532
1	905726	rs751084833	C	T,A	8733.36	PASS	.
1	905727	rs761609807	G	A	12936.9	PASS	.
(..)
```
## Note to self

I got problem with the certificate. Fixed with `-Dcom.sun.security.enableAIAcaIssuers=true -Dcom.sun.net.ssl.checkRevocation=false `


