# CompareBamAndBuild

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Compare two  BAM files mapped on two different builds. Requires a liftover chain file


## Usage

```
Usage: cmpbamsandbuild [options] Files
  Options:
  * -c, --chain
      ) Lift Over file from bam1 to bam2. REQUIRED
    -d, --distance
      distance tolerance between two alignments
      Default: 10
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -maxRecordsInRam, --maxRecordsInRam
      Max records in RAM
      Default: 50000
    -o, --output
      Output file. Optional . Default: stdout
    -r, --region
      restrict to that region chr:start-end
    -T, --tmpDir
      mp directory
      Default: /tmp
    --version
      print version and exit

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew cmpbamsandbuild
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cmpbams/CompareBamAndBuild.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cmpbams/CompareBamAndBuild.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **cmpbamsandbuild** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 

## Example

The following **Makefile** compares two BAMs mapped on hg19 and hg38.

```makefile
BWA=bwa
SAMTOOLS=samtools
.PHONY:all

all: file.diff

file.diff.gz : hg19ToHg38.over.chain hg19.bam hg38.bam
	java -jar jvarkit-git/dist/cmpbamsandbuild.jar -d 20 \
	    -c hg19ToHg38.over.chain hg19.bam hg38.bam | gzip --best > $@

hg38.bam:file1.fastq.gz file2.fastq.gz
	${BWA} mem  index-bwa-0.7.6a/hg38.fa $^ |\
		${SAMTOOLS} view -Sb - |\
		${SAMTOOLS}  sort - hg38 && \
		${SAMTOOLS} index hg38.bam
		
hg19.bam: file1.fastq.gz file2.fastq.gz
		${BWA} mem index-bwa-0.7.6a/hg19.fa $^ |\
		${SAMTOOLS} view -Sb - |\
		${SAMTOOLS}  sort - hg19 && \
		${SAMTOOLS} index hg19.bam


hg19ToHg38.over.chain:
	curl -o $@.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/$@.gz && gunzip -f $@.gz

```

file.diff :
```tsv
#READ-Name	COMPARE	hg19.bam	hg38.bam
HWI-1KL149:18:C0RNBACXX:3:1101:10375:80749/2	NE	chrMasked:581319279->chrMasked:581415967	chrMasked:581703411
HWI-1KL149:18:C0RNBACXX:3:1101:10394:60111/1	NE	chrMasked:581319428->chrMasked:581416116	chrMasked:581703560
HWI-1KL149:18:C0RNBACXX:3:1101:10394:60111/2	NE	chrMasked:581319428->chrMasked:581416116	chrMasked:581703560
HWI-1KL149:18:C0RNBACXX:3:1101:10413:69955/1	NE	chrMasked:581318674->chrMasked:581415362	chrMasked:581702806
HWI-1KL149:18:C0RNBACXX:3:1101:10413:69955/2	NE	chrMasked:581318707->chrMasked:581415395	chrMasked:581702839
HWI-1KL149:18:C0RNBACXX:3:1101:10484:89477/1	NE	chrMasked:581319428->chrMasked:581416116	chrMasked:581703560
HWI-1KL149:18:C0RNBACXX:3:1101:10484:89477/2	NE	chrMasked:581319428->chrMasked:581416116	chrMasked:581703560
HWI-1KL149:18:C0RNBACXX:3:1101:10527:3241/1	NE	chrMasked:581319279->chrMasked:581415967	chrMasked:581703411
HWI-1KL149:18:C0RNBACXX:3:1101:10527:3241/2	NE	chrMasked:581319279->chrMasked:581415967	chrMasked:581703411
HWI-1KL149:18:C0RNBACXX:3:1101:10580:13030/1	NE	chrMasked:581319331->chrMasked:581416019	chrMasked:581703463
HWI-1KL149:18:C0RNBACXX:3:1101:10580:13030/2	NE	chrMasked:581319331->chrMasked:581416019	chrMasked:581703463
HWI-1KL149:18:C0RNBACXX:3:1101:10618:51813/1	EQ	chrMasked:581318674->chrMasked:581415362	chrMasked:581415362
HWI-1KL149:18:C0RNBACXX:3:1101:10618:51813/2	NE	chrMasked:581318707->chrMasked:581415395	chrMasked:581702839
HWI-1KL149:18:C0RNBACXX:3:1101:10803:23593/1	NE	chrMasked:581319091->chrMasked:581415779	chrMasked:581703223
HWI-1KL149:18:C0RNBACXX:3:1101:10803:23593/2	NE	chrMasked:581319091->chrMasked:581415779	chrMasked:581703223
HWI-1KL149:18:C0RNBACXX:3:1101:11290:76217/1	EQ	chrMasked:581318674->chrMasked:581415362	chrMasked:581415362
HWI-1KL149:18:C0RNBACXX:3:1101:11290:76217/2	NE	chrMasked:581318707->chrMasked:581415395	chrMasked:581702839
HWI-1KL149:18:C0RNBACXX:3:1101:11307:71853/1	NE	chrMasked:581319254->chrMasked:581415942	chrMasked:581703386
HWI-1KL149:18:C0RNBACXX:3:1101:11307:71853/2	NE	chrMasked:581319258->chrMasked:581415946	chrMasked:581703390
HWI-1KL149:18:C0RNBACXX:3:1101:11359:100655/1	EQ	chrMasked:581318923->chrMasked:581415611	chrMasked:581415611
HWI-1KL149:18:C0RNBACXX:3:1101:11359:100655/2	EQ	chrMasked:581318923->chrMasked:581415611	chrMasked:581415611
HWI-1KL149:18:C0RNBACXX:3:1101:11467:8793/1	NE	chrMasked:581319428->chrMasked:581416116	chrMasked:581703560
HWI-1KL149:18:C0RNBACXX:3:1101:11467:8793/2	NE	chrMasked:581319428->chrMasked:581416116	chrMasked:581703560
HWI-1KL149:18:C0RNBACXX:3:1101:11560:69825/1	EQ	chrMasked:581319331->chrMasked:581416019	chrMasked:581416019
HWI-1KL149:18:C0RNBACXX:3:1101:11560:69825/2	EQ	chrMasked:581319335->chrMasked:581416023	chrMasked:581416023

```


## See also

* CmpBams

## History

* 2014: Creation
 
 
