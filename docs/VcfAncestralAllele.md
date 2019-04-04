# VcfAncestralAllele

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate a VCF with it's ancestral allele. Data from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.README


## Usage

```
Usage: vcfancestralalleles [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -m, --manifest
      Manifest file containing the path to the fasta files. See doc. ALL fasta 
      files must be indexed with `samtools faidx`
    -o, --output
      Output file. Optional . Default: stdout
    -t, --tag
      Ancestral allele INFO attribute
      Default: AA
    --version
      print version and exit

```


## Keywords

 * vcf
 * sort


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfancestralalleles
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/onekgenomes/VcfAncestralAllele.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/onekgenomes/VcfAncestralAllele.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/onekgenomes/VcfAncestralAlleleTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/onekgenomes/VcfAncestralAlleleTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfancestralalleles** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

For @isamtalves : annotate a VCF for ancestral allele using data from 1000 genomes: 

  * http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.README

## Manifest

manifest is a tab delimited file.

empty lines are ignored

lines starting with '#' are ignored

three columns:

  * 1st column: REF chromosome name, multiple/aliases can be separated by a pipe '|'
  * 2nd column: ancestral contig name
  * 3td column: File PATH to the equivalent ancestral contig.

### note to self: building the manifest

```bash
 find ${PWD} -name "*.fa" | sort -V | while read F; do echo "$F" | sed 's%/commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_%%' | sed 's/.fa//' | tr -d "\n" && echo -ne "\t" &&  head -n 1 $F | cut -c 2- | tr -d '\n' && echo -en "\t" && echo $F  ; done | sed 's/^\([^\t]*\)/\1|chr\1/'
```

```
1|chr1    ANCESTOR_for_chromosome:GRCh37:1:1:249250621:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_1.fa
2|chr2    ANCESTOR_for_chromosome:GRCh37:2:1:243199373:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_2.fa
3|chr3    ANCESTOR_for_chromosome:GRCh37:3:1:198022430:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_3.fa
4|chr4    ANCESTOR_for_chromosome:GRCh37:4:1:191154276:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_4.fa
5|chr5    ANCESTOR_for_chromosome:GRCh37:5:1:180915260:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_5.fa
```

## Example:

```
$ java -jar dist/vcfancestralalleles.jar \
	-m /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/manifest.mf \
	src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz |\
	bcftools annotate -x '^INFO/AA'

(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	905606	rs540662886	G	C,A	41743.9	PASS	AA=G
1	905608	rs770396126	G	T,A	1006.45	PASS	AA=G
1	905609	.	G	A	694.04	PASS	AA=G
1	905610	rs775689041	TG	AG,T,TGGGGGGCCCAG	4327.1	PASS	AA=TG
1	905611	rs749351425	G	T	2434.35	PASS	AA=G
1	905616	.	G	C	1801.88	PASS	AA=G
1	905617	rs376988925	C	T,G	4379.6	PASS	AA=C
1	905619	rs774441222	C	T	19350.2	PASS	AA=C
1	905621	rs368876607	G	A	14291.5	PASS	AA=G
1	905623	rs770778738	G	A	8255.09	PASS	AA=g
1	905626	.	G	A	257.48	PASS	AA=G
1	905627	rs755913930	GC	G,GCC	75797.7	PASS	AA=GC
1	905628	rs776376856	C	A	76618.4	AC0;RF	AA=C
1	905629	rs759517260	C	T	13583.4	PASS	AA=C
1	905632	.	C	T	2152.62	PASS	AA=C
1	905634	rs765105613	C	T	4534.75	PASS	AA=C
1	905635	rs752637312	C	A	4947.89	PASS	AA=C
1	905636	.	C	T	888.44	PASS	AA=C
1	905637	rs762369764	C	T	313.28	PASS	AA=C
(...)
```

