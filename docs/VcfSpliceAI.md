# VcfSpliceAI

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate VCF with spiceai web service


## Usage

```
Usage: vcfspliceai [options] Files
  Options:
    --base
      Base API
      Default: https://spliceailookup-api.broadinstitute.org/spliceai/
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --distance
      For each variant, SpliceAI looks within a window (+/- 50bp by default) 
      to see how the variant affects the probabilities of different positions 
      being splice acceptors or donors. The distance specified here controls 
      the size of this window. The maximum allowed value is 10,000bp
      Default: 50
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hg
      genome version. Must be 37 or 38. Otherwise, the dictionary is used to 
      detect the version.
      Default: -1
    -o, --out
      Output file. Optional . Default: stdout
    --raw
      Splicing changes corresponding to strengthening annotated splice sites 
      and weakening unannotated splice sites are typically much less 
      pathogenic than weakening annotated splice sites and strengthening 
      unannotated splice sites. Selecting 'masked' (default) will hide the 
      score for such splicing changes and show 0 instead. Selecting 'raw' will 
      show all scores. SpliceAI developers recommend using 'raw' scores for 
      alternative splicing analysis and 'masked' scores for variant 
      interpretation. 
      Default: false
    --tag
      INFO tag
      Default: SPLICEAI
    --version
      print version and exit

```


## Keywords

 * vcf
 * splice
 * splicing
 * spliceai


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfspliceai
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20201107

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/spliceai/VcfSpliceAI.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/spliceai/VcfSpliceAI.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfspliceai** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Example

```
java -jar dist/vcfspliceai.jar  src/test/resources/test_vcf01.vcf 

(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
(...)
1	866893	.	T	C	431	PASS	AA=t;AC=7;AF=0.7;AN=10;SPLICEAI=SAMD11|0.00|0.00|0.00|0.00|-13|24|-13|-48
1	870317	.	G	A	12	PASS	AC=11;AF=0.917;AN=12;SPLICEAI=SAMD11|0.00|0.00|0.00|0.00|2|17|16|-12
1	875770	.	A	G	338	PASS	AA=a;AC=8;AF=0.8;AN=10;SPLICEAI=SAMD11|0.00|0.00|0.01|0.00|-1|-45|-50|-46
1	903245	.	A	G	199	PASS	AA=a;AC=6;AF=0.6;AN=10;SPLICEAI=PLEKHN1|0.00|0.00|0.00|0.00|48|-37|-22|1
1	905130	.	ATG	A	487	PASS	AC=3;AF=0.5;AN=6;CIGAR=1M2D;IDREP=1;REFREP=2;RU=TG;SPLICEAI=PLEKHN1|0.00|0.00|0.00|0.00|-43|21|-33|-37
1	909238	.	G	C	229	PASS	AA=C;AC=8;AF=0.667;AN=12;SPLICEAI=PLEKHN1|0.00|0.01|0.00|0.00|-43|-50|39|-7
1	912049	.	T	C	400	PASS	AA=T;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.01|0.01|0.00|-28|-14|-27|-23
1	913889	.	G	A	372	PASS	AA=G;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.01|0.00|0.00|-46|9|2|-45
1	914333	.	C	G	556	PASS	AA=G;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.00|0.00|0.00|-3|27|-3|-38
1	914852	.	G	C	525	PASS	AA=C;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.00|0.00|0.00|22|-22|48|49
1	914940	.	T	C	488	PASS	AA=C;AC=5;AF=0.625;AN=8;SPLICEAI=PERM1|0.00|0.00|0.00|0.00|28|-30|-39|3
(...)
```

