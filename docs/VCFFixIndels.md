# VCFFixIndels

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Fix samtools indels (for @SolenaLS)


## DEPRECATED

use `bcftools norm`

## Usage

```
Usage: vcffixindels [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -t, --tag
      INFO/tag
      Default: INDELFIXED
    --version
      print version and exit

```


## Keywords

 * vcf
 * indel


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcffixindels
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcffixindels/VCFFixIndels.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcffixindels/VCFFixIndels.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcffixindels/VCFFixIndelsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcffixindels/VCFFixIndelsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcffixindels** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### See also

 *  https://github.com/lindenb/jvarkit/wiki/VCFFixIndels
 *  "Unified Representation of Genetic Variants" http://bioinformatics.oxfordjournals.org/content/early/2015/02/19/bioinformatics.btv112.abstract (hey ! it was published after I wrote this tool !)
 *  https://github.com/quinlan-lab/vcftidy/blob/master/vcftidy.py
 *  http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/

### Example

```
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/input_callsets/si/ALL.wgs.samtools_pass_filter.20130502.snps_indels.low_coverage.sites.vcf.gz" |\
 gunzip -c | java -jar dist/vcfstripannot.jar -k '*' 2> /dev/null |\
 java -jar dist/vcffixindels.jar  2> /dev/null | grep FIX | head -n 15

##INFO=<ID=INDELFIXED,Number=1,Type=String,Description="Fix Indels for @SolenaLS (position|alleles...)">
1   2030133 .   T   TTTTGT,TTTTG    999 PASS    INDELFIXED=2030101|CGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGT*|CGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGT|CGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTG
1   3046430 .   C   CCCT,CCC    999 PASS    INDELFIXED=3046429|TC*|TCCCT|TCCC
1   4258325 rs137902679;rs61115653  A   AAT,AA  999 PASS    INDELFIXED=4258316|CAAAAAAAAA*|CAAAAAAAAAA|CAAAAAAAAAAT
1   5374885 rs59294415  C   CCCC,CCCCA  999 PASS    INDELFIXED=5374881|TCCCC*|TCCCCCCC|TCCCCCCCA
1   5669438 rs143435517 C   CACAT,CAC   999 PASS    INDELFIXED=5669414|TACACACACACACACACACACACAC*|TACACACACACACACACACACACACAC|TACACACACACACACACACACACACACAT
1   5702062 .   A   AA,AAC  999 PASS    INDELFIXED=5702060|TAA*|TAAAC|TAAA
1   5713682 rs70977965  A   AAAAA,AAAAAC    999 PASS    INDELFIXED=5713678|CAAAA*|CAAAAAAAA|CAAAAAAAAC
1   5911136 .   T   TGCCATT,TGCCATTCCAAAGAGGCACTCA  999 PASS    INDELFIXED=5911135|CT*|CTGCCATTCCAAAGAGGCACTCA|CTGCCATT
1   6067269 rs34064079;rs59468731   G   GG,GGC  999 PASS    INDELFIXED=6067261|TGGGGGGGG*|TGGGGGGGGG|TGGGGGGGGGC
1   6069948 .   TC  T,TTC   999 PASS    INDELFIXED=6069933|CTTTTTTTTTTTTTTTC*|CTTTTTTTTTTTTTTTTC|CTTTTTTTTTTTTTTT
1   6480784 .   C   CGGGCCCCAGGCTGCCCGCC,CGGGCCCCAGGCTGCCCGCCT  999 PASS    INDELFIXED=6480783|GC*|GCGGGCCCCAGGCTGCCCGCCT|GCGGGCCCCAGGCTGCCCGCC
1   6829081 rs34184977;rs5772255    A   AAC,AA  999 PASS    INDELFIXED=6829070|TAAAAAAAAAAA*|TAAAAAAAAAAAA|TAAAAAAAAAAAAC
1   7086193 .   AG  A,AAG   999 PASS    INDELFIXED=7086179|TAAAAAAAAAAAAAAG*|TAAAAAAAAAAAAAAAG|TAAAAAAAAAAAAAA
1   8096161 .   T   TATATATATAC,TAT 999 PASS    INDELFIXED=8096143|CATATATATATATATATAT*|CATATATATATATATATATAT|CATATATATATATATATATATATATATAC

```


