# VcfAlleleBalance

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Insert missing allele balance annotation using FORMAT:AD


## Usage

```
Usage: vcfallelebalance [options] Files
  Options:
    -f, --filtered
      ignore FILTER-ed **GENOTYPES**
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -p, -ped, --pedigree, --ped
      A pedigree is a text file delimited with tabs. No header. Columns are 
      (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother Id or '0' 
      (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 unaffected, 1 
      affected,-9 unknown
    -s, --snp
      consider only snps
      Default: false
    --version
      print version and exit

```


## Keywords

 * vcf
 * allele-balance
 * depth


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfallelebalance
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfAlleleBalance.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfAlleleBalance.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfAlleleBalanceTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfAlleleBalanceTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfallelebalance** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

August 2018, for @MKarakachoff

part of the code was inspired from GATK public code :  https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/annotator/AlleleBalance.java

## Example

```
$ java -jar dist/vcfallelebalance.jar  src/test/resources/test_vcf01.vcf
$ java -jar dist/vcfallelebalance.jar -p src/test/resources/test_vcf01.ped src/test/resources/test_vcf01.vcf
```

