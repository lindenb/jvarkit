# SamAddPI

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Add predicted median insert size 'PI' to SAM Read groups (RG).


## Usage

```
Usage: samaddpi [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    --filter
      A JEXL Expression that will be used to filter out some sam-records (see 
      https://software.broadinstitute.org/gatk/documentation/article.php?id=1255). 
      An expression should return a boolean value (true=exclude, false=keep 
      the read). An empty expression keeps everything. The variable 'record' 
      is the current observed read, an instance of SAMRecord (https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html).
      Default: record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -N, --num-reads
      Number of reads to test. Negative=all = memory consuming.
      Default: 1000000
    -o, --output
      Output file. Optional . Default: stdout
    -w, --overwrite
      Overwrite median insert size if it already exists
      Default: false
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew samaddpi
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamAddPI.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamAddPI.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamAddPITest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SamAddPITest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samaddpi** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example


```
$ samtools view -h ~/src/gatk-ui/testdata/S1.bam | head -n 10
@HD	VN:1.5	GO:none	SO:coordinate
@SQ	SN:rotavirus	LN:1074
@RG	ID:S1	SM:S1
@PG	ID:bwa	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_01_R1.fq.gz S1_01_R2.fq.gz
@PG	ID:bwa.1	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_02_R1.fq.gz S1_02_R2.fq.gz
@PG	ID:bwa.2	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_03_R1.fq.gz S1_03_R2.fq.gz
rotavirus_1_317_5:0:0_7:0:0_2de	99	rotavirus	1	60	70M	=	248	317	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:33G4A3T14A7T4	RG:Z:S1	NM:i:5	AS:i:45	XS:i:0
rotavirus_1_535_4:0:0_4:0:0_1a6	163	rotavirus	1	60	70M	=	466	535	GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:8A18G13C6G21	RG:Z:S1	NM:i:4	AS:i:50	XS:i:0
rotavirus_1_543_5:0:0_11:0:0_390	163	rotavirus	1	60	70M	=	487	530	GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:18G1G1T22T11A12	RG:Z:S1	NM:i:5	AS:i:45XS:i:0
rotavirus_1_578_3:0:0_7:0:0_7c	99	rotavirus	1	60	70M	=	509	578	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:43A2C5A17	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0



$ java -jar dist/samaddpi.jar S1.bam | head -n 10
@HD	VN:1.5	GO:none	SO:coordinate
@SQ	SN:rotavirus	LN:1074
@RG	ID:S1	SM:S1	PI:498
@PG	ID:bwa	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_01_R1.fq.gz S1_01_R2.fq.gz
@PG	ID:bwa.1	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_02_R1.fq.gz S1_02_R2.fq.gz
@PG	ID:bwa.2	PN:bwa	VN:0.7.12-r1044	CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_03_R1.fq.gz S1_03_R2.fq.gz
@CO	Processed with SamAddPI /home/lindenb/src/gatk-ui/testdata/S1.bam
rotavirus_1_317_5:0:0_7:0:0_2de	99	rotavirus	1	60	70M	=	248	317	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:33G4A3T14A7T4	RG:Z:S1	NM:i:5	AS:i:45	XS:i:0
rotavirus_1_535_4:0:0_4:0:0_1a6	163	rotavirus	1	60	70M	=	466	535	GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:8A18G13C6G21	RG:Z:S1	NM:i:4	AS:i:50	XS:i:0
rotavirus_1_543_5:0:0_11:0:0_390	163	rotavirus	1	60	70M	=	487	530	GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT	++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	MD:Z:18G1G1T22T11A12	RG:Z:S1	NM:i:5	AS:i:45XS:i:0

```
