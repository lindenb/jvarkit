# VcfMultiToOne

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

>Convert VCF with multiple samples to a VCF with one SAMPLE, duplicating variant and adding the sample name in the INFO column


## Usage

```
Usage: vcfmulti2one [options] Files
  Options:
    --disable-vc-attribute-recalc
      When genotypes are removed/changed, Dd not recalculate variant 
      attributes like DP, AF, AC, AN...
      Default: false
    -r, --discard_hom_ref
      discard if variant is hom-ref
      Default: false
    -c, --discard_no_call
      discard if variant is no-call
      Default: false
    -a, --discard_non_available
      discard if variant is not available
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --vc-attribute-recalc-ignore-filtered
      When recalculating variant attributes like DP AF, AC, AN, ignore 
      FILTERed **Genotypes**
      Default: false
    --vc-attribute-recalc-ignore-missing
      Ignore missing VCF headers (DP, AF, AC, AN). Default behavior: adding 
      VCF header if they're missing
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * sample



## See also in Biostars

 * [https://www.biostars.org/p/130456](https://www.biostars.org/p/130456)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfmulti2one
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/onesamplevcf/VcfMultiToOne.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/onesamplevcf/VcfMultiToOne.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfmulti2one** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ curl -s "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz" |\
gunzip -c |\
java -jar dist/vcfmulti2one.jar  -c -r -a  |\
grep -v '##' |\
grep -E '(CHROM|SAMPLENAME)' | head | verticalize 


>>> 2
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00096;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 2

>>> 3
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00097;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 0|1
<<< 3

>>> 4
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00099;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 0|1
<<< 4

>>> 5
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00100;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 5

>>> 6
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00102;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 6

>>> 7
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00103;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 7

>>> 8
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00105;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 8

>>> 9
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00106;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 1|0
<<< 9

>>> 10
$1   #CHROM : 1
$2      POS : 10177
$3       ID : .
$4      REF : A
$5      ALT : AC
$6     QUAL : 100
$7   FILTER : PASS
$8     INFO : AA=|||unknown(NO_COVERAGE);AC=2130;AF=0.425319;AFR_AF=0.4909;AMR_AF=0.3602;AN=5008;DP=103152;EAS_AF=0.3363;EUR_AF=0.4056;NS=2504;SAMPLENAME=HG00114;SAS_AF=0
.4949
$9   FORMAT : GT
$10  SAMPLE : 0|1
<<< 10
```



