# SamFindClippedRegions

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Fins clipped position in one or more bam. 


## Usage

```
Usage: samfindclippedregions [options] Files
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
    --gtf
      Optional gtf file. Will be used to set a warning if the junction could 
      be a junction exon-exon of a retrogene. A GTF (General Transfer Format) 
      file. See https://www.ensembl.org/info/website/upload/gff.html . Please 
      note that CDS are only detected if a start and stop codons are defined.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --intron-distance
      when gtf is specified: max distance between breakend and the intron 
      bound 
      Default: 3
    --mapq
      min mapping quality
      Default: 1
    --min-clip-depth
      Ignore if number of clipped bases lower than 'x'
      Default: 10
    --min-depth
      Ignore if Depth lower than 'x'
      Default: 10
    --min-ratio
      Ignore genotypes where count(clip)/(count(clip)+DP) < x
      Default: 0.1
    -o, --out
      Output file. Optional . Default: stdout
    --groupby, --partition
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit
    -c
      consider only clip having length >= 'x'
      Default: 1

```


## Keywords

 * sam
 * bam
 * clip
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samfindclippedregions
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20140228

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamFindClippedRegions.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamFindClippedRegions.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamFindClippedRegionsTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/structvar/SamFindClippedRegionsTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samfindclippedregions** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


samfindclippedregions find 'blunt' regions where the reads are clipped. It can be used to find structural variations in exomes/genomes data.

input is a set of indexed BAM/CRAM files or a list with the '.list' suffix containing the path to the bam.

output is a VCF file


### Example

```
$ java -jar dist/samfindclippedregions.jar --min-depth 10 --min-ratio 0.2 src/test/resources/S*.bam


##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=CL,Number=1,Type=Integer,Description="Left Clip">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=RL,Number=1,Type=Integer,Description="Right Clip">
##FORMAT=<ID=TL,Number=1,Type=Integer,Description="Total Clip">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
(...)
##samfindclippedregions.meta=compilation:20191009163450 githash:b8d60cab htsjdk:2.20.1 date:20191009163608 cmd:--min-depth 10 --min-ratio 0.2 src/test/resources/S1.bam src/test/resources/S2.bam src/test/resources/S3.bam src/test/resources/S4.bam src/test/resources/S5.bam
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF01	996	.	N	<CLIP>	.	.	AC=1;AF=0.1;AN=10;DP=30	GT:AD:CL:DP:RL:TL	0/0:2,0:0:2:0:0	0/0:4,0:0:4:0:00/0:4,0:0:4:0:0	0/0:15,0:0:15:0:0	0/1:4,1:0:5:1:1
(...)
```


### Screenshot

https://twitter.com/yokofakun/status/1194921855875977216

![https://twitter.com/yokofakun/status/1194921855875977216](https://pbs.twimg.com/media/EJU3F9hWoAACgsd?format=png&name=large)

