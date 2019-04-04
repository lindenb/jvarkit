# VcfStats

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Produce VCF statitics


## Usage

```
Usage: vcfstats [options] Files
  Options:
    --binSize
      [20170718] When plotting data over a genome, divide it into 'N' bp.
      Default: 1000000
    --disableGTConcordance
      Disable Plot Sample vs Sample Genotypes (Faster...)
      Default: false
    --disableMAFPlot
      Disable MAF plot
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -K, -kg, --knownGenes
      UCSC knownGene File/URL. The knowGene format is a compact alternative to 
      GFF/GTF because one transcript is described using only one line.	Beware 
      chromosome names are formatted the same as your REFERENCE. A typical 
      KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
    -mafTag, --mafTag
      Do not calculate MAF for controls, but use this tag to get Controls' MAF
    -nchr, --nocallishomref
      treat no call as HomRef
      Default: false
  * -o, --output
      output Directory or zip file. The output contains the data files as well 
      as a Makefile to convert the data files to graphics using gnuplot.
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
      tranches for the number of affected. A 'range of integers' is a list of 
      integers in ascending order separated with semicolons.
      Default: [[-Inf/0[, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10/20[, [20/50[, [50/100[, [100/200[, [200/300[, [300/400[, [400/500[, [500/1000[, [1000/Inf[]
    --trancheAlts
      tranches for the number of ALTs. A 'range of integers' is a list of 
      integers in ascending order separated with semicolons.
      Default: [[-Inf/0[, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10/Inf[]
    --trancheDP
      tranches for the DEPTH. A 'range of integers' is a list of integers in 
      ascending order separated with semicolons.
      Default: [[-Inf/0[, [0/10[, [10/20[, [20/30[, [30/50[, [50/100[, [100/200[, [200/300[, [300/400[, [400/500[, [500/600[, [600/700[, [700/800[, [800/900[, [900/1000[, [1000/2000[, [2000/3000[, [3000/4000[, [4000/5000[, [5000/10000[, [10000/20000[, [20000/30000[, [30000/40000[, [40000/50000[, [50000/100000[, [100000/Inf[]
    --trancheDistance
      tranches for the distance between the variants. A 'range of integers' is 
      a list of integers in ascending order separated with semicolons.
      Default: [[-Inf/0[, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10/20[, [20/100[, [100/200[, [200/300[, [300/400[, [400/500[, [500/1000[, [1000/Inf[]
    --trancheIndelSize
      tranches for the Indel size A 'range of integers' is a list of integers 
      in ascending order separated with semicolons.
      Default: [[-Inf/0[, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10/15[, [15/20[, [20/50[, [50/100[, [100/Inf[]
    --vckey
      Variant Context Key. if defined, I will look at this key in the INFO 
      column and produce a CASE/CTRL graf for each item. If undefined, I will 
      produce a default graph with all variant
    --version
      print version and exit

```


## Keywords

 * vcf
 * stats
 * burden
 * gnuplot


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfstats
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStats.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStats.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStatsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfstats/VcfStatsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfstats** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

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
ALL.affectedSamples.png     ALL.countDepthBySample.tsv      ALL.countDistances.png  ALL.geneLoc.tsv  ALL.predictionsBySample.png  ALL.sample2gtype.tsv  Makefile
ALL.affectedSamples.tsv     ALL.countDepth.png              ALL.countDistances.tsv  ALL.maf.png      ALL.predictionsBySample.tsv  ALL.transvers.png
ALL.countAltAlleles.png     ALL.countDepth.tsv              ALL.countIndelSize.png  ALL.maf.tsv      ALL.predictions.png          ALL.transvers.tsv
ALL.countAltAlleles.tsv     ALL.countDistancesBySample.png  ALL.countIndelSize.tsv  ALL.mendel.png   ALL.predictions.tsv          ALL.variant2type.png
ALL.countDepthBySample.png  ALL.countDistancesBySample.tsv  ALL.geneLoc.png         ALL.mendel.tsv   ALL.sample2gtype.png         ALL.variant2type.tsv

$ head -n3 tmp/*.tsv
==> tmp/ALL.affectedSamples.tsv <==
1/57	182265
2/57	98512
3/57	67449

==> tmp/ALL.countAltAlleles.tsv <==
2	11699
3	2406
4	679

==> tmp/ALL.countDepthBySample.tsv <==
Sample	[-Inf/0[	[0/10[	[10/20[	[20/30[	[30/50[	[50/100[	[100/200[	[200/Inf[
10_T1245	0	305229	97739	61972	75791	72779	14008	686
11AG09	0	400147	62261	39411	55254	86883	86192	46749

==> tmp/ALL.countDepth.tsv <==
[0/10[	19435
[10/20[	95369
[20/30[	92378

==> tmp/ALL.countDistancesBySample.tsv <==
Sample	[-Inf/0[	0	1	2	3	4	5	6	7	8	9	[10/20[	[20/100[	[100/200[	[200/300[	[300/400[	[400/500[	[500/1000[	[1000/Inf[
S1	0	0	2643	1755	1491	1430	1156	1106	893	858	867	6826	26696	18252	12411	9348	7270	20975	85769
S2	0	0	3632	2397	1967	1855	1516	1425	1246	1179	1109	8917	36006	25455	18652	14951	11912	37828	103341

==> tmp/ALL.countDistances.tsv <==
1	16275
2	11184
3	8505

==> tmp/ALL.countIndelSize.tsv <==
2	59420
3	22727
4	10216

==> tmp/ALL.geneLoc.tsv <==
first_exon	51180
internal_exon	70074
last_exon	75341

==> tmp/ALL.maf.tsv <==
0.029411764705882353	0.027777777777777776
0.5	0.43333333333333335
0.05263157894736842	0.0

==> tmp/ALL.mendel.tsv <==
Sample	synonymous_variant	protein_altering_variant
10_T1245	0	0
11AG09	170	298

==> tmp/ALL.predictionsBySample.tsv <==
Sample	synonymous_variant	protein_altering_variant
S1	13453	14527
S2	12820	14077

==> tmp/ALL.predictions.tsv <==
protein_altering_variant	59516
synonymous_variant	47454

==> tmp/ALL.sample2gtype.tsv <==
Sample	NO_CALL	HOM_REF	HET	HOM_VAR	UNAVAILABLE	MIXED
S1	345776	429002	98825	100945	0	0
S2	196925	504211	132576	140836	0	0

==> tmp/ALL.transvers.tsv <==
TYPE	TRANSITION	TRANSVERSION
ALL	581537	133472
CDS	533340	120720

==> tmp/ALL.variant2type.tsv <==
Type	Count
INDEL	126219
SNP	846677

```


