# FixVcfMissingGenotypes

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

After a VCF-merge, read a VCF, look back at some BAMS to tells if the missing genotypes were homozygotes-ref or not-called. If the number of reads is greater than min.depth, then a missing genotype is said hom-ref.


## Usage

```
Usage: fixvcfmissinggenotypes [options] Files
  Options:
    -B, --bams
      path of indexed BAM path with read Groups. You can put those paths in a 
      text file having a *.list sufffix
      Default: []
    -d, --depth
      minimal depth before setting a genotype to HOM_REF
      Default: 10
    --disable-vc-attribute-recalc
      When genotypes are removed/changed, Dd not recalculate variant 
      attributes like DP, AF, AC, AN...
      Default: false
    -filter, --filter
      A JEXL Expression that will be used to filter out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    --filtered
      Mark fixed genotypes as FILTERED with this FILTER
    --force, -f
      [20181120] Update all fields like DP even if the Genotype is called.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --partition
      Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -T, --tag
      FORMAT 'Tag' for fixed genotype
      Default: FXG
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

 * sam
 * bam
 * vcf
 * sv
 * genotype



## See also in Biostars

 * [https://www.biostars.org/p/119007](https://www.biostars.org/p/119007)
 * [https://www.biostars.org/p/263309](https://www.biostars.org/p/263309)
 * [https://www.biostars.org/p/276811](https://www.biostars.org/p/276811)
 * [https://www.biostars.org/p/302581](https://www.biostars.org/p/302581)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew fixvcfmissinggenotypes
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FixVcfMissingGenotypes.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FixVcfMissingGenotypes.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/FixVcfMissingGenotypesTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/FixVcfMissingGenotypesTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fixvcfmissinggenotypes** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Slow

This tool remains slow because there is a random-access in the bam for each './.' genotype.

You can always try to speed-up things by breaking your VCF in multiple regions and process them in parallel.

## Examples


### Example

```
$ find ~/src/gatk-ui/testdata/ -name "*.bam" > input.list

$ tail -2 input.vcf
rotavirus	1064	.	G	A	21.5606	.	DP=250;VDB=2.70971e-16;SGB=8.40135;RPB=0.935144;MQB=1;BQB=0.683886;MQ0F=0;AF1=0.25;G3=0.75,2.37734e-17,0.25;HWE=0.033921;AC1=2;DP4=0,219,0,31;MQ=60;FQ=22.8019;PV4=1,1.22605e-06,1,1	GT:PL	0/0:0,244,70	0/0:0,199,65	0/0:0,217,68	1/1:69,84,0
rotavirus	1064	.	G	A	21.5606	.	DP=250;VDB=2.70971e-16;SGB=8.40135;RPB=0.935144;MQB=1;BQB=0.683886;MQ0F=0;AF1=0.25;G3=0.75,2.37734e-17,0.25;HWE=0.033921;AC1=2;DP4=0,219,0,31;MQ=60;FQ=22.8019;PV4=1,1.22605e-06,1,1	GT:PL	./.	./.	./.	./.

$ java -jar dist/fixvcfmissinggenotypes.jar -d 50 --fixDP --filtered zz -B input.list input.vcf | tail -2
rotavirus	1064	.	G	A	21.56	.	AC1=2;AF1=0.25;BQB=0.683886;DP=188;DP4=0,219,0,31;FQ=22.8019;G3=0.75,2.37734e-17,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1.22605e-06,1,1;RPB=0.935144;SGB=8.40135;VDB=2.70971e-16	GT:DP:PL	0/0:48:0,244,70	0/0:63:0,199,65	0/0:53:0,217,68	1/1:24:69,84,0
rotavirus	1064	.	G	A	21.56	.	AC1=2;AF1=0.25;BQB=0.683886;DP=72;DP4=0,219,0,31;FQ=22.8019;G3=0.75,2.37734e-17,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1.22605e-06,1,1;RPB=0.935144;SGB=8.40135;VDB=2.70971e-16	GT:DP:FT:FXG	./.:48:PASS	0/0:63:zz:1	0/0:53:zz:1	./.:24:PASS
```

### Example


```
$ yourtool-mergingvcf 1.vcf 2.vcf 3.vcf > merged.vcf
$ find ./ -name "*.bam" > bams.list
$  java -jar dist/fixvcfmissinggenotypes.jar -B bams.list < merged.vcf > out.vcf
```

```
$ find DIR1 -name "PREFIX_*final.bam"  | grep -E '(S1|S2|S3|S4)' ) > bams.list

$ find DIR1 -name "PREFIX_*_variations.gatk.annotations.vcf.gz" |\
grep -E '(S1|S2|S3|S4)' |\
xargs perl  vcftools_0.1.12b/perl vcftools_0.1.12b/bin/vcf-merge |\
java -jar dist/fixvcfmissinggenotypes.jar -d 10 -f  bams.list |\
gzip --best > out.vcf.gz

```

## Cited in 

 * "Exome sequencing in genomic regions related to racing performance of Quarter Horses" Pereira, G.L., Malheiros, J.M., Ospina, A.M.T. et al. J Appl Genetics (2019). https://doi.org/10.1007/s13353-019-00483-1
 

### History

 * 2018-11-20 : adding features for structural variants
 * 2017-07-24 : rewrite whole program 
 * 2014: Creation


