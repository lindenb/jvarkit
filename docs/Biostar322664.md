# Biostar322664

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Extract PE Reads (with their mates) supporting variants in vcf file


## Usage

```
Usage: biostar322664 [options] Files
  Options:
    -all, --all
      used in complement of option -X . Will output any SamRecord, but some 
      reads will carrying the information about the variant in a 'Xx' 
      attribute. 
      Default: false
    --bamcompression
      Compression Level.
      Default: 5
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -index, --index
      Use the VCF input to query the BAM using bai index. Faster for large bam 
      + small VCF. Require option `--no-mate` and the bam file to be indexed.
      Default: false
    -nm, --no-mate
      Disable the 'mate' function. BAM is not expected to be sorted with 
      picard (bam can be sorted on coordinate), but mate will not be written.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -pair, --pair
      pair mode: the paired read and it's pair muts BOTH carry at least one 
      variant 
      Default: false
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
  * -V, --variant
      Variant VCF file. This tool **doesn't work** with INDEL/SV.
    --version
      print version and exit
    -x, -X
      If defined, add the variant(s) information in a 'X' metadata. One 
      character only.

```


## Keywords

 * sam
 * bam
 * vcf



## See also in Biostars

 * [https://www.biostars.org/p/322664](https://www.biostars.org/p/322664)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar322664
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar322664.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar322664.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar322664Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar322664Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar322664** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

  * BAM : MUST be sorted using Picard SortSam (see https://github.com/samtools/hts-specs/issues/5 )
  * VCF : only SNP are considered. Genotypes are ignored (all ALT alleles are observed regardless of the sample/genotype )

##Example

```
$ java -jar picard.jar SortSam I=src/test/resources/S1.bam O=query.bam SO=queryname
$ java -jar dist/biostar322664.jar -V src/test/resources/S1.vcf.gz query.bam  


RF02_358_926_2:0:0_2:1:0_83	83	RF02	857	60	70M	=	358	-569	GACGTGAACTATATAATTAAAATGGACAGAAATCTGCCATCAACAGCTAGATATATAAGACCTAATTTAC	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0
RF02_362_917_2:0:0_2:1:0_6f	147	RF02	848	60	70M	=	362	-556	ATAAGGAATCACGTTAACTATATACTTAAAATGGACTGAAATCTGCCATCAACAGCTAGATATATAAGAC	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0
(...)
```

```
$ java -jar dist/biostar322664.jar -nm  -index -V src/test/resources/S1.vcf.gz src/test/resources/S1.bam 
(...)
RF02_827_1385_4:1:0_2:0:0_3f    163     RF02    827     60      70M     =       1316    559     ATCAATTACATTCCTGAAAGGATAAGGAATGAGGTTAACTATCTACTTAAAATGGACAGAAATCTGCCAA  2222222222222222222222222222222222222222222222222222222222222222222222  RG:Z:S1 NM:i:5  AS:i:53XS:i:0
RF02_827_1292_4:0:0_1:0:0_50    163     RF02    827     60      70M     =       1223    466     TTCAATTACATTCCTGCAAGGATAAGGAATGCCGTTAACTATATACTTAATAAGGACAGAAATCTGGCAT  2222222222222222222222222222222222222222222222222222222222222222222222  RG:Z:S1 NM:i:4  AS:i:51XS:i:0
(...)
```


