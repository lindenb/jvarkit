# VcfRegulomeDB

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate a VCF with the Regulome data (http://regulome.stanford.edu/


## Usage

```
Usage: vcfregulomedb [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -T
      tag in vcf INFO.
      Default: REGULOMEDB
  * -b
       bed indexed with tabix. Format: chrom(tab)start(tab)end(tab)rank
    -r
      if defined, only accept the rank matching the regular expression
    -x
      (int) base pairs. look.for data around the variation +/- 'x'
      Default: 5

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfregulomedb
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRegulomeDB.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfRegulomeDB.java)


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


## Building the Tabix File for regulomeDB:
```bash
(for S in 1 2 3 4 5 6 ; do curl -s "http://regulome.stanford.edu/downloads/RegulomeDB.dbSNP132.Category${S}.txt.gz"  | gunzip -c |  cut -f 1,2,5 | sed -e 's/^chrX/23/'  -e 's/^chr//'  | awk -F ' ' '{printf("%s\t%d\t%s\t%s\n",$1,int($2)-1,$2,$3);}' | uniq ; done)| LC_ALL=C sort -t ' ' -k1,1n -k2,2n -k3,3n |  sed 's/^23/X/' | bgzip -c > regulomeDB.bed.gz && tabix  -p bed -f regulomeDB.bed.gz 
```

## Example

```bash
$   curl -kLs "https://raw.githubusercontent.com/arq5x/gemini/master/test/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.snippet.snpEff.vcf" |\
   java -jar dist/vcfstripannot.jar -k '*' |\
   java -jar dist/vcfregulomedb.jar -b regulomeDB.bed.gz -x 10  |\
   grep -i REGU | head

##INFO=<ID=REGULOMEDB,Number=.,Type=String,Description="Format: Position|Distance|Rank">
##VcfRegulomeDBCmdLine=-b regulomeDB.bed.gz -x 10
##VcfRegulomeDBVersion=235a18b083ea15c4ad94060de512f1edc74cec42
1	10583	rs58108140	G	A	100	PASS	REGULOMEDB=10582|0|5
1	13327	rs144762171	G	C	100	PASS	REGULOMEDB=13326|0|6
1	13980	rs151276478	T	C	100	PASS	REGULOMEDB=13971|8|6,13979|0|6,13980|1|6
1	46402	.	C	CTGT	31	PASS	REGULOMEDB=46402|1|6
1	55164	rs3091274	C	A	100	PASS	REGULOMEDB=55163|0|6
1	55299	rs10399749	C	T	100	PASS	REGULOMEDB=55298|0|6
1	55313	rs182462964	A	T	100	PASS	REGULOMEDB=55321|9|6
1	55326	rs3107975	T	C	100	PASS	REGULOMEDB=55321|4|6,55325|0|6
1	55330	rs185215913	G	A	100	PASS	REGULOMEDB=55321|8|6,55325|4|6

```

Those SNPs can be seen at:

* http://regulome.stanford.edu/snp/chr1/10582
* http://regulome.stanford.edu/snp/chr1/13326
* http://regulome.stanford.edu/snp/chr1/13971
* http://regulome.stanford.edu/snp/chr1/46402
* etc...

