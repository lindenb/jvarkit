# FixVcfMissingGenotypes

After a VCF-merge, read a VCF, look back at some BAMS to tells if the missing genotypes were homozygotes-ref or not-called. If the number of reads is greater than min.depth, then the missing genotypes is said hom-ref.


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
    -filter, --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    --filtered
      Mark fixed genotypes as FILTERED with this FILTER
    --fixDP
      Update/create DP field even if genotype is called but there is no DP
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -T, --tag
      FORMAT 'Tag' for fixed genotype
      Default: FXG
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * vcf



## See also in Biostars

 * [https://www.biostars.org/p/119007](https://www.biostars.org/p/119007)
 * [https://www.biostars.org/p/263309](https://www.biostars.org/p/263309)


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
$ make fixvcfmissinggenotypes
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FixVcfMissingGenotypes.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FixVcfMissingGenotypes.java
)
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
$  java -jar dist/fixvcfmissinggenotypes.jar -f bams.list < merged.vcf > out.vcf
```

```
$ find DIR1 -name "PREFIX_*final.bam"  | grep -E '(S1|S2|S3|S4)' ) > bams.list

$ find DIR1 -name "PREFIX_*_variations.gatk.annotations.vcf.gz" |\
grep -E '(S1|S2|S3|S4)' |\
xargs perl  vcftools_0.1.12b/perl vcftools_0.1.12b/bin/vcf-merge |\
java -jar dist/fixvcfmissinggenotypes.jar -d 10 -f  bams.list |\
gzip --best > out.vcf.gz

```



### History

 * 2017-07-24 : rewrite whole program 
 * 2014: Creation



