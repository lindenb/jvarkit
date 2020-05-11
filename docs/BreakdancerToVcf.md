# BreakdancerToVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert output of breakdancer to VCF


## Usage

```
Usage: breakdancer2vcf [options] Files
  Options:
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -R, -reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary

```


## Keywords

 * cnv
 * sv
 * breakdancer
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew breakdancer2vcf
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20200511

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/breakdancer/BreakdancerToVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/breakdancer/BreakdancerToVcf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/breakdancer/BreakdancerToVcfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/breakdancer/BreakdancerToVcfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **breakdancer2vcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## input

Input is the tabular output of BreakDancer.

Genotypes are all set to 0/1.

## Example:

```
$ wget -O - -q "https://raw.githubusercontent.com/genome/breakdancer/master/test-data/expected_output" | java -jar dist/breakdancer2vcf.jar -R src/test/resources/human_b37.dict 

##fileformat=VCFv4.2
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimate allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=CHROM2,Number=1,Type=String,Description="Chromosome 2">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=ORIENT1,Number=1,Type=String,Description="Orientation 1">
##INFO=<ID=ORIENT2,Number=1,Type=String,Description="Orientation 2">
##INFO=<ID=POS2,Number=1,Type=String,Description="Position 2">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variation type">
##breakdancer.command=bdfast -o21 inv_del_bam_config
##breakdancer.version=1.4.1-unstable-10-fdfe9f2-dirty (commit fdfe9f2-dirty)
##breakdancer2vcf.meta=compilation:20200511120615 githash:c830e0b htsjdk:2.21.3 date:20200511120941 cmd:-R src/test/resources/human_b37.dict
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
(...)
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##contig=<ID=MT,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	H_IJ-NA19238-NA19238-extlibs	H_IJ-NA19240-NA19240-extlibs
21	29185056	.	N	<INS>	99	.	DP=2;END=29185377;SVLEN=226;SVTYPE=INS	GT:DP:GQ	0/1:1:99	0/1:1:99
21	29185462	.	N	<DEL>	99	.	DP=21;END=29186122;SVLEN=-545;SVTYPE=DEL	GT:AF:DP:GQ	0/1:176.58:21:99	0/1:167.89:.:99
21	34807694	.	N	<INS>	99	.	DP=3;END=34808852;SVLEN=304;SVTYPE=INS	GT:DP:GQ	0/1:1:99	0/1:2:99
21	34808937	.	N	<INV>	99	.	DP=2;END=34809799;SVLEN=-737;SVTYPE=INV	GT:AF:DP:GQ	0/1:847.39:.:99	0/1:878.83:2:99
```



