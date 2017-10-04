# BimToVcf

convert a .bim to a .vcf


## Usage

```
Usage: bim2vcf [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit

```

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
$ make bim2vcf
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BimToVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BimToVcf.java)


<details>
<summary>Git History</summary>

```
Tue Jun 6 18:06:17 2017 +0200 ; postponed vcf ; https://github.com/lindenb/jvarkit/commit/bcd52318caf3cd76ce8662485ffaacaabde97caf
Sun Jun 4 21:53:22 2017 +0200 ; writing bcf ; https://github.com/lindenb/jvarkit/commit/784fdac37cd7e6eca04e35d0a3ddad8637826b4a
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Tue Apr 25 17:33:00 2017 +0200 ; fix make ; https://github.com/lindenb/jvarkit/commit/621fdfbe039f474e187e8f5c67ff17b368b77289
Tue Jul 12 16:29:16 2016 +0200 ; bim2vcf ; https://github.com/lindenb/jvarkit/commit/e319a093f69d1ffba0f31ea95b8c305922c17320
Tue Jul 12 16:26:04 2016 +0200 ; bim2vcf ; https://github.com/lindenb/jvarkit/commit/6df744a561505df03f56e0b18d75587f032fcdbe
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bim2vcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Contig conversion

chromosome 23 is converted to X or chrX, chromosome 24 is converted to Y or chrY, chromosome 25 is ignored, chromosome 26 is converted to chrM or MT.


### Example



```
$ java -jar dist/bim2vcf.jar -R human_g1k_v37.fasta input.bim 

##fileformat=VCFv4.2
##INFO=<ID=MORGAN,Number=1,Type=Float,Description="Centimorgan">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variation type">
##contig=<ID=1,length=249250621,assembly=human_g1k_v37>
(...)
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
1       12      rs73422 C       .       .       .       MORGAN=0.5224;SVTYPE=NOVARIATION
1       13      rs30315 G       A       .       .       MORGAN=0.530874;SVTYPE=SNV
1       14      rs14325 C       T       .       .       MORGAN=0.532596;SVTYPE=SNV
1       15      rs31319 A       G       .       .       MORGAN=0.532682;SVTYPE=SNV
1       16      rs954   C       T       .       .       MORGAN=0.537655;SVTYPE=SNV
1       17      rs62034 G       A       .       .       MORGAN=0.548645;SVTYPE=SNV
1       18      rs25996 A       G       .       .       MORGAN=0.575595;SVTYPE=SNV
1       19      rs12117 G       A       .       .       MORGAN=0.582608;SVTYPE=SNV
(...)

```









