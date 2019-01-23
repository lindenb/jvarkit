# Biostar251649

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

 Annotating the flanking bases of SNPs in a VCF file


## Usage

```
Usage: biostar251649 [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
  * -r, -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --version
      print version and exit
    -3
      Right tag
      Default: SEQ3_
    -5
      Left tag
      Default: SEQ5_
    -n
      number of bases
      Default: 1

```


## Keywords

 * vcf
 * annotation
 * sequence
 * reference



## See also in Biostars

 * [https://www.biostars.org/p/251649](https://www.biostars.org/p/251649)
 * [https://www.biostars.org/p/334253](https://www.biostars.org/p/334253)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar251649
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar251649.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar251649.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar251649Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar251649Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar251649** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Examples

### Example

```
$ java -jar dist/biostar251649.jar -n 10 -R tests/ref.fa tests/mutations.vcf
##INFO=<ID=SEQ3_10,Number=1,Type=String,Description="Sequence on the 3' of mutation">
##INFO=<ID=SEQ5_10,Number=1,Type=String,Description="Sequence on the 5' of mutation">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC1=2;AF1=0.25;BQB=1;DP=944;DP4=849,0,93,0;FQ=23.7972;G3=0.75,0,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.993129;SEQ3_10=GATGGTAAGC;SEQ5_10=TCTACTCAGC;SGB=-61.9012;VDB=3.53678e-05	GT:PL	0/0:0,255,134	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	AC1=1;AF1=0.124963;BQB=0.951201;DP=1359;DP4=1134,0,225,0;FQ=5.8713;MQ=60;MQ0F=0;MQB=1;PV4=1,4.80825e-05,1,1;RPB=0.0393173;SEQ3_10=GTTGTTGCTG;SEQ5_10=TTGAAGCTGC;SGB=-369.163;VDB=0.313337	GT:PL	0/0:0,255,133	0/1:40,0,31	0/0:0,255,134	0/0:0,255,82
```

### Example

convert to fasta:

```
java -jar dist/biostar251649.jar -R src/test/resources/rotavirus_rf.fa src/test/resources/rotavirus_rf.vcf.gz -n 20 |\
java -jar dist/bioalcidaejdk.jar -F VCF -e 'stream().forEach(V->println(">"+V.getContig()+":"+V.getStart()+"\n"+V.getAttribute("SEQ5_20")+"["+V.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/"))+"]"+V.getAttribute("SEQ3_20")));'
```


