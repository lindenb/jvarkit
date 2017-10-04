# BamToFastq

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
      What kind of help
      Possible Values: [usage, markdown, xml]
    -maxRecordsInRam, --maxRecordsInRam
      Max records in RAM
      Default: 50000
    -r, --repair
      repair: insert missing read
      Default: false
    -R, --reverse
      Save fastq_R2 to file (default: interlaced with forward)
    -T, --tmpDir
      mp directory
      Default: /tmp
    --version
      print version and exit

```


## Keywords

 * fastq


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
$ make bam2fastq
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/BamToFastq.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/fastq/BamToFastq.java)


<details>
<summary>Git History</summary>

```
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Fri May 12 18:07:46 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/ca96bce803826964a65de33455e5231ffa6ea9bd
Tue Apr 18 13:24:50 2017 +0200 ; cont-cleanup ; https://github.com/lindenb/jvarkit/commit/a86c8971fe5ebb3f8de175c75e78f2d0e5325cfd
Fri Mar 25 17:18:27 2016 +0100 ; sammask ; https://github.com/lindenb/jvarkit/commit/b9c834afec6c7c9904baecd2fb2b61e57261da0f
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Mon Oct 13 18:29:16 2014 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/c83f20cde867920870918ee6eb5e5406f554e2bb
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Thu Apr 3 17:24:00 2014 +0200 ; cont: doap, sam2fasta, sam2psl... ; https://github.com/lindenb/jvarkit/commit/bc9f11b1a0a1a7b0874a3b74d75b368e4de0bf98
Thu Feb 27 17:10:54 2014 +0100 ; cont, fix bug in bam2fastq, shortread, starting change-ref bam, extract clipped seq ; https://github.com/lindenb/jvarkit/commit/d83138c95883cf87078565b54614b2aa7aa04740
Wed Nov 20 17:26:42 2013 +0100 ; bam2fastq ; https://github.com/lindenb/jvarkit/commit/0fd5257cd14c23d833ee96e2a1e3c79a441a584d
```

</details>

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



Deprecated: use picard please

Previous version was an Implementation of https://twitter.com/DNAntonie/status/402909852277932032



Warnings

	illumina  read is filtered is always "n"
	illumina control number is always 0
	Illumina index sequence is lost.



Example
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



Example 2 from BAM


```

$ java -jar dist/bam2fastq.jar \
    -F tmpR1.fastq.gz -R tmpR2.fastq.gz file.bam

(...)
-rw-r--r-- 1 lindenb lindenb 565M Nov 18 10:44 Sample_S1_L001_R1_001.fastq.gz
-rw-r--r-- 1 lindenb lindenb 649M Nov 18 10:45 Sample_S1_L001_R2_001.fastq.gz
-rw-rw-r-- 1 lindenb lindenb 470M Nov 20 16:17 tmpR1.fastq.gz.fastq.gz
-rw-rw-r-- 1 lindenb lindenb 554M Nov 20 16:17 tmpR2.fastq.gz.fastq.gz

```



