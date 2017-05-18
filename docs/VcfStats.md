# VcfStats


## Usage

```
Usage: vcfstats [options] Files
  Options:
    -h, --help
      print help and exit
    -kg, --knownGenes
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
    -tee, --tee
      output the incoming vcf to stdout. Useful to get intermediary stats in a 
      pipeline 
      Default: false
    --trancheAffected
      tranches for the number of affected. A range of is a list of integers is 
      ascending order separated with semicolons.
      Default: [[-Inf / 0[, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10 / 20[, [20 / 50[, [50 / 100[, [100 / Inf[]
    --trancheAlts
      tranches for the number of ALTs. A range of is a list of integers is 
      ascending order separated with semicolons.
      Default: [[-Inf / 0[, 0, 1, 2, 3, 4, 5, [6 / 8[, 8, 9, [10 / Inf[]
    --trancheDP
      tranches for the DEPTH. A range of is a list of integers is ascending 
      order separated with semicolons.
      Default: [[-Inf / 0[, [0 / 10[, [10 / 20[, [20 / 30[, [30 / 50[, [50 / 100[, [100 / 200[, [200 / Inf[]
    --trancheDistance
      tranches for the distance between the variants. A range of is a list of 
      integers is ascending order separated with semicolons.
      Default: [[-Inf / 0[, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, [10 / 20[, [20 / 100[, [100 / 200[, [200 / 300[, [300 / 400[, [400 / 500[, [500 / 1000[, [1000 / Inf[]
    --trancheIndelSize
      tranches for the Indel size A range of is a list of integers is 
      ascending order separated with semicolons.
      Default: [[-Inf / 0[, 0, 1, 2, 3, 4, 5, [6 / 8[, 8, 9, [10 / 15[, [15 / 20[, [20 / Inf[]
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
$ java -jar dist/vcfstats.jar ~karaka/BURDEN_JVARKIT/MVP/20170227.pct001gBED.Q.mvp_frex.vcf.gz -o tmp
$ ls tmp/
ALL.sample2gtype.tsv  ALL.variant2type.tsv  Makefile

$ head tmp/ALL.sample2gtype.tsv  tmp/ALL.variant2type.tsv

==> tmp/ALL.sample2gtype.tsv <==
Type	NO_CALL	HOM_REF	HET	HOM_VAR	UNAVAILABLE	MIXED
X0G73J	2538	3440	218	132	0	0
Y00G3I	2543	3462	193	130	0	0
Z03K	2547	3417	252	112	0	0
A0G3N	2547	3424	209	148	0	0
B980	1909	4068	202	149	0	0
C003P	2559	3417	216	136	0	0
D0741	2557	3428	204	139	0	0
E073O	2566	3433	212	117	0	0
F00G7	2560	3444	191	133	0	0
(...)

==> tmp/ALL.variant2type.tsv <==
Type	Count
MIXED	7
SNP	6125
INDEL	196

```



