# VcfCalledWithAnotherMethod

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Compare one vcf with other , add a flag to tell if a variant was called with another method. Vcf must be sorted on the same Dict.


## Usage

```
Usage: vcfcalledwithanothermethod [options] Files
  Options:
    --filter
      FILTER name: the variant was NOT found in another VCF 
      (CONTIG/POS/REF/at-least-one-ALT). Empty: no filter
      Default: VariantNotFoundElseWhere
    --foundCount
      INFO name for the file identifiers where a variant was found
      Default: FOUND_COUNT
    --foundKey
      INFO name for the file identifiers where a variant was found
      Default: FOUND_KEY
    --gtDiscordant
      FORMAT name for the number of time we didn't find the same genotype
      Default: COUNT_DISCORDANT
    --gtSame
      FORMAT name for the number of time we found the same genotype
      Default: COUNT_SAME
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --nocallhomref
      NO_CALL is same as HOM_REF
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -f, --vcfs
      Add alternate VCF files. File ending with '.list' will be interpreted as 
      a list of path of vcf.
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * compare
 * concordance


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfcalledwithanothermethod
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfCalledWithAnotherMethod.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfCalledWithAnotherMethod.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfcalledwithanothermethod** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example


```make
SHELL=/bin/bash

define ff
dir1/dir2/sample_variations.$(1).annotations.vcf.gz
endef

all :
	java -jar  dist/vcfcalledwithanothermethod.jar \
		-f $(call ff,samtools) \
		-f $(call ff,varscan) \
		-f $(call ff,freebayes) \
			$(call ff,gatkHapCaller)
	
```

output:

```
(...)
1	12718	.	G	C	197.77	VariantNotFoundElseWhere	AC=1;AF=0.500;AN=2;BaseQRankSum=-1.418;ClippingRankSum=1.220;DP=22;QD=8.99;ReadPosRankSum=1.022;SEGDUP=1:10485-19844,1:10464-40733,1:10000-19844,1:10485-40733,1:10000-87112,1:10000-20818	GT:AD:COUNT_DISCORDANT:COUNT_SAME:DP:GQ:PL	0/1:12,10:0:0:22:99:226,0,286
1	23119	.	T	G	637.77	PASS	FOUND_COUNT=2;FOUND_KEY=sample_variations.varscan.annotations,sample_variations.samtools.annotations;FS=34.631;GERP_SCORE=-0.558;MLEAC=1;MLEAF=0.500;MQ=25.98;MQ0=0;MQRankSum=-2.888;POLYX=1;PRED=uc010nxq.1|||||intron_variant;QD=18.22;ReadPosRankSum=1.634;SEGDUP=1:10485-19844,1:10464-40733,1:10000-19844,1:10485-40733,1:10000-87112,1:10000-20818	GT:AD:COUNT_DISCORDANT:COUNT_SAME:DP:GQ:PL	0/1:17,18:0:2:35:99:666,0,727

```


