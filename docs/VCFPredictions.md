# VCFPredictions

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Basic Variant Effect prediction using ucsc-known gene


## Usage

```
Usage: vcfpredictions [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -k, --knownGene
      UCSC knownGene File/URL. The knowGene format is a compact alternative to 
      GFF/GTF because one transcript is described using only one line.	Beware 
      chromosome names are formatted the same as your REFERENCE. A typical 
      KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
      Default: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz
    -o, --output
      Output file. Optional . Default: stdout
    -os, --output-syntax, --syntax
      [20180122]output formatting syntax. SnpEff is still not complete.
      Default: Native
      Possible Values: [Native, Vep, SnpEff]
    -soacn, --printsoacn
      Print SO:term accession rather than label
      Default: false
  * -R, --reference
      [20180122](moved to faidx/DAS). Indexed Genome Reference. It can be a 
      the path to fasta file that must be indexed with samtools faidx and with 
      picard CreateSequenceDictionary. It can also be a BioDAS dsn url like 
      `http://genome.cse.ucsc.edu/cgi-bin/das/hg19/` . BioDAS references are 
      slower, but allow to work without a local reference file.
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation
 * prediction


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfpredictions
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannot/VCFPredictions.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfannot/VCFPredictions.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfpredictions** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



First described in http://plindenbaum.blogspot.fr/2011/01/my-tool-to-annotate-vcf-files.html.



### Examples




#### Example 1




```
$java -jar  dist/vcfpredictions.jar  \
	-R human_g1k_v37.fasta \
	-k <(curl  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c | awk -F '   ' '{if($2 ~ ".*_.*") next; OFS="        "; gsub(/chr/,"",$2);print;}'  ) \
	I=~/WORK/variations.gatk.annotations.vcf.gz | bgzip > result.vcf.gz

```



#### Annotate 1000 genomes




```
$  curl "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" |\
gunzip -c |\
java -jar dist/vcfpredictions.jar \
  -R  human_g1k_v37.fasta \
  -k <(curl  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c | awk -F '      ' '{if($2 ~ ".*_.*") next; OFS=" "; gsub(/chr/,"",$2);print;}'  )  |\
cut -d '    ' -f 1-8 | grep -A 10 CHROM


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10583	rs58108140	G	A	100	PASS	AA=.;AC=314;AF=0.14;AFR_AF=0.04;AMR_AF=0.17;AN=2184;ASN_AF=0.13;AVGPOST=0.7707;ERATE=0.0161;EUR_AF=0.21;LDAF=0.2327;PRED=|||||intergenic_variant;RSQ=0.4319;SNPSOURCE=LOWCOV;THETA=0.0046;VT=SNP
1	10611	rs189107123	C	G	100	PASS	AA=.;AC=41;AF=0.02;AFR_AF=0.01;AMR_AF=0.03;AN=2184;ASN_AF=0.01;AVGPOST=0.9330;ERATE=0.0048;EUR_AF=0.02;LDAF=0.0479;PRED=|||||intergenic_variant;RSQ=0.3475;SNPSOURCE=LOWCOV;THETA=0.0077;VT=SNP
1	13302	rs180734498	C	T	100	PASS	AA=.;AC=249;AF=0.11;AFR_AF=0.21;AMR_AF=0.08;AN=2184;ASN_AF=0.02;AVGPOST=0.8895;ERATE=0.0058;EUR_AF=0.14;LDAF=0.1573;PRED=uc010nxq.1|||||intron_variant;RSQ=0.6281;SNPSOURCE=LOWCOV;THETA=0.0048;VT=SNP
1	13327	rs144762171	G	C	100	PASS	AA=.;AC=59;AF=0.03;AFR_AF=0.02;AMR_AF=0.03;AN=2184;ASN_AF=0.02;AVGPOST=0.9698;ERATE=0.0012;EUR_AF=0.04;LDAF=0.0359;PRED=uc010nxq.1|||||intron_variant;RSQ=0.6482;SNPSOURCE=LOWCOV;THETA=0.0204;VT=SNP
1	13957	rs201747181	TC	T	28	PASS	AA=TC;AC=35;AF=0.02;AFR_AF=0.02;AMR_AF=0.02;AN=2184;ASN_AF=0.01;AVGPOST=0.8711;ERATE=0.0065;EUR_AF=0.02;LDAF=0.0788;PRED=uc010nxq.1|||||;RSQ=0.2501;THETA=0.0100;VT=INDEL
1	13980	rs151276478	T	C	100	PASS	AA=.;AC=45;AF=0.02;AFR_AF=0.01;AMR_AF=0.02;AN=2184;ASN_AF=0.02;AVGPOST=0.9221;ERATE=0.0034;EUR_AF=0.02;LDAF=0.0525;PRED=uc010nxq.1|||||3_prime_UTR_variant;RSQ=0.3603;SNPSOURCE=LOWCOV;THETA=0.0139;VT=SNP
1	30923	rs140337953	G	T	100	PASS	AA=T;AC=1584;AF=0.73;AFR_AF=0.48;AMR_AF=0.80;AN=2184;ASN_AF=0.89;AVGPOST=0.7335;ERATE=0.0183;EUR_AF=0.73;LDAF=0.6576;PRED=|||||intergenic_variant;RSQ=0.5481;SNPSOURCE=LOWCOV;THETA=0.0162;VT=SNP
1	46402	rs199681827	C	CTGT	31	PASS	AA=.;AC=8;AF=0.0037;AFR_AF=0.01;AN=2184;ASN_AF=0.0017;AVGPOST=0.8325;ERATE=0.0072;LDAF=0.0903;PRED=|||||intergenic_variant;RSQ=0.0960;THETA=0.0121;VT=INDEL
1	47190	rs200430748	G	GA	192	PASS	AA=G;AC=29;AF=0.01;AFR_AF=0.06;AMR_AF=0.0028;AN=2184;AVGPOST=0.9041;ERATE=0.0041;LDAF=0.0628;PRED=|||||intergenic_variant;RSQ=0.2883;THETA=0.0153;VT=INDEL
1	51476	rs187298206	T	C	100	PASS	AA=C;AC=18;AF=0.01;AFR_AF=0.01;AMR_AF=0.01;AN=2184;ASN_AF=0.01;AVGPOST=0.9819;ERATE=0.0021;EUR_AF=0.01;LDAF=0.0157;PRED=|||||intergenic_variant;RSQ=0.5258;SNPSOURCE=LOWCOV;THETA=0.0103;VT=SNP
```





### See also



 *  http://plindenbaum.blogspot.fr/2013/10/yes-choice-of-transcripts-and-software.html
 *  VcfFilterSequenceOntology
 *  GroupByGene



## History

 *  2013-Dec : moved to a standard arg/argv command line







