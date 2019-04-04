# VcfCompareCallers

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Compare two VCFs and print common/exclusive information for each sample/genotype


## Usage

```
Usage: vcfcomparecallers [options] Files
  Options:
    -B, --bed
      Limit to variants in that BED region
    --collapseGenotypeType
      collapse Genotype Type. Just show Same Genotype or Discordant, don't 
      print the type of genotype.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --jexl1
      An optional list of GATK-like JEXL expressions to filter the variants 
      from VCF File 1
      Default: []
    --jexl2
      An optional list of GATK-like JEXL expressions to filter the variants 
      from VCF File 2
      Default: []
    -c, --noCallIsHomRef
      No Call is HomRef (created when comparing merged vcf with GATK: there is 
      no homref, everything is nocall)
      Default: false
  * -o, --output
      Directory or zip file to save results to be plotted with gnuplot
    -p, --prefix
      Archive prefix (for option -d)
      Default: <empty string>
    -vcf1, --vcf1
      short descriptive name for VCF1
      Default: VCF1
    -vcf2, --vcf2
      short descriptive name for VCF2
      Default: VCF2
    --version
      print version and exit

```


## Keywords

 * vcf
 * compare
 * genotype


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfcomparecallers
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfCompareCallers.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfCompareCallers.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfCompareCallersTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfcmp/VcfCompareCallersTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcomparecallers** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Synopsis

```
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) stdin 
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) file2.vcf(.gz) 

```


both vcf **must** share the same sequence dictionary and must be sorted

#### History

* 20170704 : rewritten from scratch

### Example


```
$ java -jar dist/vcfcomparecallers.jar -o tmp Proj1.samtools.vcf.gz  Proj1.varscan.vcf.gz
$ (cd tmp && make)
```


