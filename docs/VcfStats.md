# VcfStats


## Usage

```
Usage: vcfstats [options] Files
  Options:
    -h, --help
      print help and exit
    -K, -kg, --knownGenes
      UCSC knownGene URI. Beware chromosome names are formatted the same as 
      your REFERENCE. A typical KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
    -mafTag, --mafTag
      Do not calculate MAF for controls, but use this tag to get Controls' MAF
    -nchr, --nocallishomref
      treat no call as HomRef
      Default: false
  * -o, --output
      output Directory or zip file
    -ped, --pedigree
      A pedigree is a text file delimited with tabs. No header. Columns are 
      (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother Id or '0' 
      (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 unaffected, 1 
      affected,-9 unknown
    --prefix
      File/zip prefix
      Default: <empty string>
    -so, --soterms
      Sequence ontology Accession to observe. VCF must be annotated with 
      SNPEFF or VEP. e.g: "SO:0001818" (protein altering variant) "SO:0001819" 
      (synonymouse variant)
      Default: []
    -tee, --tee
      output the incoming vcf to stdout. Useful to get intermediary stats in a 
      pipeline 
      Default: false
    --trancheAffected
      tranches for the number of affected. A range of is a list of integers is 
      ascending order separated with semicolons.
      Default: [[-Inf/0[, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10/20[, [20/50[, [50/100[, [100/Inf[]
    --trancheAlts
      tranches for the number of ALTs. A range of is a list of integers is 
      ascending order separated with semicolons.
      Default: [[-Inf/0[, 0, 1, 2, 3, 4, 5, [6/8[, 8, 9, [10/Inf[]
    --trancheDP
      tranches for the DEPTH. A range of is a list of integers is ascending 
      order separated with semicolons.
      Default: [[-Inf/0[, [0/10[, [10/20[, [20/30[, [30/50[, [50/100[, [100/200[, [200/Inf[]
    --trancheDistance
      tranches for the distance between the variants. A range of is a list of 
      integers is ascending order separated with semicolons.
      Default: [[-Inf/0[, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10/20[, [20/100[, [100/200[, [200/300[, [300/400[, [400/500[, [500/1000[, [1000/Inf[]
    --trancheIndelSize
      tranches for the Indel size A range of is a list of integers is 
      ascending order separated with semicolons.
      Default: [[-Inf/0[, 0, 1, 2, 3, 4, 5, [6/8[, 8, 9, [10/15[, [15/20[, [20/Inf[]
    --vckey
      Variant Context Key. if defined, I will look at this key in the INFO 
      column and produce a CASE/CTRL graf for each item. If undefined, I will 
      produce a default graph with all variant
    --version
      print version and exit

```


## Description

Produce VCF statitics


## Keywords

 * vcf
 * stats
 * burden
 * gnuplot


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
$ make vcfstats
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStats.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStats.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfstats** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

http://dx.doi.org/10.6084/m9.figshare.1425030

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Tip: Adding a new key in the INFO field

Using vcffilterjs :



the script:

```
var ArrayList = Java.type("java.util.ArrayList");
var VariantContextBuilder = Java.type("htsjdk.variant.variantcontext.VariantContextBuilder");


function addInfo(v)
	{
	var vcb = new VariantContextBuilder(v);
	var atts = new ArrayList();
	atts.add(v.getType().name()+ (variant.isFiltered()?"_FILTERED":"_UNFILTERED"));
	atts.add(v.getType().name()+ (variant.hasID()?"_ID":"_NOID"));
	vcb.attribute("MYKEY",atts);
	return vcb.make();
	}


addInfo(variant);
```

run the program, but first use awk to insert the new INFO definition for 'MYKEY'

```
cat input.vcf |\
	awk '/^#CHROM/ {printf("##INFO=<ID=MYKEY,Number=.,Type=String,Description=\"My key\">\n");} {print;}' |\
	java -jar dist/vcffilterjs.jar -f script.js 
```


## Example

```
$ java -jar $< -o tmp --soterms SO:0001818 --soterms SO:0001819 -K ucsc/hg19/database/knownGene_noPrefix.txt.gz input.vcf
$ (cd tmp && make) 
$ ls tmp/
ALL.affectedSamples.png         ALL.countDistances.png  ALL.geneLoc.png  ALL.predictionsBySample.png  ALL.sample2gtype.png  ALL.variant2type.png
ALL.affectedSamples.tsv         ALL.countDistances.tsv  ALL.geneLoc.tsv  ALL.predictionsBySample.tsv  ALL.sample2gtype.tsv  ALL.variant2type.tsv
ALL.countDistancesBySample.png  ALL.countIndelSize.png  ALL.maf.png      ALL.predictions.png          ALL.transvers.png     Makefile
ALL.countDistancesBySample.tsv  ALL.countIndelSize.tsv  ALL.maf.tsv      ALL.predictions.tsv          ALL.transvers.tsv

==> tmp/ALL.affectedSamples.tsv <==
1/532	56
2/532	4
[10/20[/532	1
[100/Inf[/532	1

==> tmp/ALL.countDistancesBySample.tsv <==
Sample	[-Inf/0[	0	1	2	3	4	5	6	7	8	9	[10/20[	[20/100[	[100/200[	[200/300[	[300/400[	[400/500[	[500/1000[	[1000/Inf[
B00GWF0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
B00GWFL	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
B00GWFO	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
B00GWFS	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1

==> tmp/ALL.countDistances.tsv <==
1	1
2	1
3	1
6	1
9	1

==> tmp/ALL.countIndelSize.tsv <==
2	4
4	2

==> tmp/ALL.geneLoc.tsv <==
first_exon	4
internal_exon	46
last_exon	9

==> tmp/ALL.maf.tsv <==
0.006097560975609756	0.014705882352941176
0.001524390243902439	0.0
0.001524390243902439	0.0
0.001524390243902439	0.0
0.001524390243902439	0.0

==> tmp/ALL.predictionsBySample.tsv <==
Sample	synonymous_variant	protein_altering_variant
B00GWE0	0	1
B00GWE1	0	1
B00GWE2	0	1
B00GWE4	0	1

==> tmp/ALL.predictions.tsv <==
protein_altering_variant	59

==> tmp/ALL.sample2gtype.tsv <==
Sample	NO_CALL	HOM_REF	HET	HOM_VAR	UNAVAILABLE	MIXED
B00GWDY	0	62	0	0	0	0
B00GWDZ	0	62	0	0	0	0
B00GWE0	0	61	1	0	0	0
B00GWE1	0	61	1	0	0	0

==> tmp/ALL.transvers.tsv <==
TYPE	TRANSITION	TRANSVERSION
ALL	44	9
CDS	44	9

==> tmp/ALL.variant2type.tsv <==
Type	Count
INDEL	6
SNP	56

```



