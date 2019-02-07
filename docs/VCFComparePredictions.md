# VCFComparePredictions

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Compare predictions (SNPEff, VEP) for several VCFs


## Usage

```
Usage: vcfcmppred [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --maxRecordsInRam
      When writing  files that need to be sorted, this will specify the number 
      of records stored in RAM before spilling to disk. Increasing this number 
      reduces the number of file  handles needed to sort a file, and increases 
      the amount of RAM needed
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    --tmpDir
      tmp working directory. Default: java.io.tmpDir
      Default: []
    --version
      print version and exit

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfcmppred
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VCFComparePredictions.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfcmp/VCFComparePredictions.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcmppred** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

### VEP

```bash
$  java -jar dist/vcfcmppred.jar  f1.vcf f2.vcf 
(...)
7	8566286	rs2139	A	VEP discordant SO:terms between f1.vcf and f2.vcf	[SO:0001619, SO:0001632]
```


in f1.vcf (VEP 75) CSQ contains:

* intron_variant
* downstream_gene_variant
* nc_transcript_variant

in f2.vcf (VEP 71) CSQ contains:

* intron_variant

### SNPEFF

```bash
$  java -jar dist/vcfcmppred.jar  f1.vcf f2.vcf 
(...)
 8	1394127	.	G	SNPEff discordant SO:terms between between f1.vcf and f2.vcf	[SO:0001630]
```

in f1.vcf  (snpEff_3_6) EFF contains:

* downstream_gene_variant
* intron_variant
* splice_region_variant

in f2.vcf  (snpEff_3_4) EFF contains:

* downstream_gene_variant
* intron_variant




