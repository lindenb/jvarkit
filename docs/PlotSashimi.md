# PlotSashimi

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Print Sashimi plots from Bam


## Usage

```
Usage: plotsashimi [options] Files
  Options:
    --css
      Custom CSS stylesheet
    --force-max-coverage
      Force the maximum coverage to this value. ignored if <=0
      Default: 0
    -g, --gtf
      A GTF (General Transfer Format) file. See 
      https://www.ensembl.org/info/website/upload/gff.html . Please note that 
      CDS are only detected if a start and stop codons are defined.
    --gzip
      Generate gzipped compressed svg files.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -u, --url, --hyperlink
      creates a hyperlink an area is clicked. creates a hyperlink when 'click' 
      in an area. The URL must contains __CHROM__, __START__ and __END__ that 
      will be replaced by their values. Predefined values are 
      'hg19','hg38','igv'. IGV : 
      "http://localhost:60151/goto?locus=__CHROM__%3A__START__-__END__" , 
      UCSC: "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__"
      Default: <empty string>
  * -r, --region, --interval
      A source of intervals. The following suffixes are recognized: vcf, 
      vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise it could be an 
      empty string (no interval) or a list of plain interval separated by '[ 
      \t\n;,]' 
      Default: (empty)
    -m, --manifest
      Manifest Bed file output containing chrom/start/end of each gene
    --mapq
      Min mapping quality
      Default: 0
  * -o, --out
      An existing directory or a filename ending with the '.zip' or '.tar' or 
      '.tar.gz' suffix.
    --partition
      Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -R, --reference
      For Reading CRAM. Indexed fasta Reference file. This file must be 
      indexed with samtools faidx and with picard CreateSequenceDictionary
    --skip-empty
      Do not generate a SVG file if there is no read in the interval
      Default: false
    -D, --use-deletion
      also use the D operator in the cigar string (default is use only 'N').
      Default: false
    --version
      print version and exit
    -w, --width
      image width.
      Default: 1000

```


## Keywords

 * bam
 * visualization
 * svg
 * rna
 * exon
 * rnaseq



## See also in Biostars

 * [https://www.biostars.org/p/497894](https://www.biostars.org/p/497894)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew plotsashimi
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20191117

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sashimi/PlotSashimi.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sashimi/PlotSashimi.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sashimi/PlotSashimiTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sashimi/PlotSashimiTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **plotsashimi** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a set of indexed BAM files or a file with the '.list' suffix containing the path to the bams.

## Example

```
 java -jar dist/plotsashimi.jar -r "chr3:38597400-38599300" -m jeter.mf --gtf jeter.gtf ENCFF331CGL.rnaseq.bam -o TMP

$ find TMP/ -name "*.svg"
TMP/86/f9362065ddce1af4b31b47be80fff6/chr3_38599300_38599300.svg

$ cat jeter.mf 
#chrom	start	end	bam	Sample	Genes	svg
chr3	38597399	38599300	ENCFF331CGL.rnaseq.bam	SCN5A	.	86/f9362065ddce1af4b31b47be80fff6/chr3_38599300_38599300.svg

```


```
$ cat jeter.list
ENCFF331CGL.rnaseq.bam


$ cat jeter.bed 
chr3	38595150	38599347
chr3	38595350	38599500



$ java -jar dist/plotsashimi.jar -r jeter.bed -m jeter.mf --gtf jeter.gtf jeter.list -o jeter.zip

$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
    36566  2019-12-05 15:52   f6/d0bf04098aa24eb3111666ed171a11/chr3_38599347_38599347.svg
    36364  2019-12-05 15:52   82/fd0c7ad61ecc0fa22eda93e75a6943/chr3_38599500_38599500.svg
---------                     -------
    72930                     2 files


$ cat jeter.mf | column -t
#chrom  start     end       bam                     Genes   Samples  svg
chr3    38595150  38599347  ENCFF331CGL.rnaseq.bam  SCN5A   .        f6/d0bf04098aa24eb3111666ed171a11/chr3_38599347_38599347.svg
chr3    38595350  38599500  ENCFF331CGL.rnaseq.bam  SCN5A   .        82/fd0c7ad61ecc0fa22eda93e75a6943/chr3_38599500_38599500.svg

```


## Screenshot

* https://twitter.com/yokofakun/status/1202587424725127168

![https://twitter.com/yokofakun/status/1202587424725127168](https://pbs.twimg.com/media/ELBy4vAX0AABSzF?format=jpg&name=small)

* https://twitter.com/yokofakun/status/1202905778140712960

![https://twitter.com/yokofakun/status/1202905778140712960](https://pbs.twimg.com/media/ELGUbZ_W4AAMKxj?format=png&name=small)

* https://twitter.com/yokofakun/status/1207337424935936001

![https://twitter.com/yokofakun/status/1207337424935936001](https://pbs.twimg.com/media/EMFS-xGXsAAVR20?format=jpg&name=small)

