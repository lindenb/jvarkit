# BamToFastq

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Same as picard/SamToFastq but allow missing reads + shuffle reads using hash(name) so you can use them with bwa. 


## DEPRECATED

use picard

## Usage

```
Usage: bam2fastq [options] Files
  Options:
    -F, --forward
      Save fastq_R1 to file (default: stdout)
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -maxRecordsInRam, --maxRecordsInRam
      Max records in RAM
      Default: 50000
    -r, --repair
      repair: insert missing read
      Default: false
    -R, --reverse
      Save fastq_R2 to file (default: interlaced with forward)
    -T, --tmpDir
      tmp directory
      Default: /tmp
    --version
      print version and exit

```


## Keywords

 * fastq


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bam2fastq
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/BamToFastq.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/BamToFastq.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bam2fastq** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



# Deprecated: 

use picard please


# Warnings

Previous version was an Implementation of https://twitter.com/DNAntonie/status/402909852277932032

	illumina  read is filtered is always "n"
	illumina control number is always 0
	Illumina index sequence is lost.

## Example
piping bwa mem

```
$ bwa mem -M  human_g1k_v37.fasta  Sample1_L001_R1_001.fastq.gz Sample2_S5_L001_R2_001.fastq.gz |\
  java -jar dist/bam2fastq.jar  -F tmpR1.fastq.gz -R tmpR2.fastq.gz

```

before:

```
$ ls -lah Sample1_L001_R1_001.fastq.gz Sample2_S5_L001_R2_001.fastq.gz
-rw-r--r-- 1 lindenb lindenb 181M Jun 14 15:20 Sample1_L001_R1_001.fastq.gz
-rw-r--r-- 1 lindenb lindenb 190M Jun 14 15:20 Sample1_L001_R2_001.fastq.gz
```

after (these are Haloplex Data, with a lot of duplicates )

```
$ ls -lah tmpR1.fastq.gz  tmpR2.fastq.gz
-rw-rw-r-- 1 lindenb lindenb  96M Nov 20 17:10 tmpR1.fastq.gz
-rw-rw-r-- 1 lindenb lindenb 106M Nov 20 17:10 tmpR2.fastq.gz
```

using BZ2:

```
$  ls -lah *.bz2
-rw-rw-r-- 1 lindenb lindenb 77M Nov 20 17:55 tmpR1.fastq.bz2
-rw-rw-r-- 1 lindenb lindenb 87M Nov 20 17:55 tmpR2.fastq.bz2
```

check the number of reads

```
$ gunzip -c Sample1_L001_R1_001.fastq.gz | wc -l
5824676
$ gunzip -c tmpR1.fastq.gz | wc -l
5824676
```

verify one read

```
$ gunzip -c Sample1_L001_R1_001.fastq.gz | cat -n | head -n 4
     1	@M00491:25:000000000-A46H3:1:1101:11697:2045 1:N:0:5
     2	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACATTGGCAAATAGCATGCCGAGGTACGCTTAAAAAAAAAACGACGCGAGGCAGGGGGGGAGGAAGCAGGGGAGCAACAGGGGGAAGGGAAGGGAAGAGAAGAAGAACGAACGAAAG
     3	+
     4	AAAAAAAA1AC1FFGCGA0AFFBGAGHHFF2GBGHH0B2DBCF101111D211B////A11///B/1DE1E/>>E//?///</<><C////<?9-9-99A-;/---;---;-9--9=---------9:AF---9//:/9/:9---9-:-9-


$ gunzip -c tmpR1.fastq.gz | cat -n | grep  -A 3 -w "@M00491:25:000000000-A46H3:1:1101:11697:2045"
5771577	@M00491:25:000000000-A46H3:1:1101:11697:2045 1:N:0:1
5771578	AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACATTGGCAAATAGCATGCCGAGGTACGCTTAAAAAAAAAACGACGCGAGGCAGGGGGGGAGGAAGCAGGGGAGCAACAGGGGGAAGGGAAGGGAAGAGAAGAAGAACGAACGAAAG
5771579	+
5771580	AAAAAAAA1AC1FFGCGA0AFFBGAGHHFF2GBGHH0B2DBCF101111D211B////A11///B/1DE1E/>>E//?///</<><C////<?9-9-99A-;/---;---;-9--9=---------9:AF---9//:/9/:9---9-:-9-
```

## Example 2 from BAM

```
$ java -jar dist/bam2fastq.jar \
    -F tmpR1.fastq.gz -R tmpR2.fastq.gz file.bam
(...)
-rw-r--r-- 1 lindenb lindenb 565M Nov 18 10:44 Sample_S1_L001_R1_001.fastq.gz
-rw-r--r-- 1 lindenb lindenb 649M Nov 18 10:45 Sample_S1_L001_R2_001.fastq.gz
-rw-rw-r-- 1 lindenb lindenb 470M Nov 20 16:17 tmpR1.fastq.gz.fastq.gz
-rw-rw-r-- 1 lindenb lindenb 554M Nov 20 16:17 tmpR2.fastq.gz.fastq.gz
```

## Cited In:

  * "Plastomes of nine hornbeams and phylogenetic implications", Ying Li & al;  Ecology and Evolution, 2018; DOI: 10.1002/ece3.4414; https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.4414 

