# VcfCadd

Annotate VCF with  Combined Annotation Dependent Depletion (CADD) (Kircher & al. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892.PubMed PMID: 24487276.


## Usage

```
Usage: vcfcadd [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -d
      processing window size
      Default: 1000
    -u
      Combined Annotation Dependent Depletion (CADD) Tabix file URI
      Default: http://krishna.gs.washington.edu/download/CADD/v1.2/whole_genome_SNVs.tsv.gz

```

## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcfcadd
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfCadd.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfCadd.java
)
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
$  curl -s  "https://raw2.github.com/arq5x/gemini/master/test/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.snippet.vcf" | \
 java -jar dist/vcfcadd.jar \
      -u "http://krishna.gs.washington.edu/download/CADD/v1.2/whole_genome_SNVs.tsv.gz" |\
 grep -E '(CADD|#)'

(...)
##INFO=<ID=CADD,Number=.,Type=String,Description="(Allele|Score|Phred) Score suggests that that variant is likely to be  observed (negative values) vs simulated(positive 
values).However, raw values do have relative meaning, with higher values indicating that a variant is more likely to be simulated (or -not observed-) and therefore more l
ikely to have deleterious effects.PHRED expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD sc
ores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc">
(..)
##VcfCaddCmdLine=-u http://krishna.gs.washington.edu/download/CADD/v1.0/1000G.tsv.gz
##VcfCaddVersion=6123910f68df940c1f3986d142f9b0414f76a43a
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	1308871	.	A	T	6.20	.	AC1=2;AF1=1;CADD=T|-0.549120|0.089;...
1	1308982	.	A	G	6.98	.	AC1=2;AF1=1;CADD=G|-0.000088|3.329;...
1	1657021	.	T	C	3.02	.	AC1=2;AF1=1;CADD=C|-0.271229|0.740;...
(..)
```

