# VcfNoCallToHomRef

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert the UNCALLED gentoypes in a VCF to HOM_REF. This tool can be used after using GATK CombineVariants.


## Usage

```
Usage: vcfnocall2homref [options] Files
  Options:
    -dp, --depth
      Default DEPTH. negative = don't set depth.
      Default: 10
    -f, --filter
      Set this Genotype FILTER for converted genotype
    -gq, --gq, --GT
      Default Genotype quality: negative : don't set GQ.
      Default: 1
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -s, --includeSamples
      only converts those samples. Default: all samples are converted.
      Default: []
    -sf, --includeSamplesFile
      only converts those samples. Default: all samples are converted. One 
      sample per line.
    --noRecount
      do not recount DP/AC/AN/AF atttributes
      Default: false
    -o, --out
      Output file. Optional . Default: stdout
    --outputbcf
      Output bcf (for streams)
      Default: false
    -p, --ploidy
      ploidy
      Default: 2
    --vcfcreateindex
      VCF, create tribble or tabix Index when writing a VCF/BCF to a file.
      Default: false
    --vcfmd5
      VCF, create MD5 checksum when writing a VCF/BCF to a file.
      Default: false
    --version
      print version and exit
    -x
      do not recount DP/AC/AN/AF atttributes
      Default: false

```


## Keywords

 * vcf



## See also in Biostars

 * [https://www.biostars.org/p/276811](https://www.biostars.org/p/276811)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfnocall2homref
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfNoCallToHomRef.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfNoCallToHomRef.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfnocall2homref** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

original VCF

```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	.	GT:PL	./.	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	.	GT:PL	0/0:0,255,133	0/1:40,0,31	./.	./.
```

default invocation

```
$ java -jar dist/vcfnocall2homref.jar input.vcf
##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC=2;AF=0.25;AN=8;DP=10	GT:DP:GQ:PL	0/0:10:1	0/0:.:.:0,255,127	0/0:.:.:0,255,137	1/1:.:.:70,255,0
rotavirus	91	.	A	T	5.45	.	AC=1;AF=0.125;AN=8;DP=20	GT:DP:GQ:PL	0/0:.:.:0,255,133	0/1:.:.:40,0,310/0:10:1	0/0:10:1
```

convert S3 and S4 only

```
$ java -jar dist/vcfnocall2homref.jar  -f CONVERTED -s S3 -s S4  ~/jeter.vcf 
##fileformat=VCFv4.2
##FILTER=<ID=CONVERTED,Description="NOCALL Genotypes converted to HOM_REF by com.github.lindenb.jvarkit.tools.misc.VcfNoCallToHomRef">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC=2;AF=0.33333334;AN=6;DP=0	GT:PL	./.	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	AC=1;AF=0.125;AN=8;DP=20	GT:DP:FT:GQ:PL	0/0:.:PASS:.:0,255,133	0/1:.:PASS:.:40,0,31	0/0:10:CONVERTED:1	0/0:10:CONVERTED:1
```

