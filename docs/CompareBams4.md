# CompareBams4

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Compare two BAM files. Print a tab-delimited report


## Usage

```
Usage: cmpbams4 [options] Files
  Options:
    -c, --chain
      Lift Over file from bam1 to bam2. Optional
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --mismatch
      Default Lift Over mismatch. negative=use default
      Default: 0.95
    -novalidchain, --novalidchain
      Disable Lift Over chain validation
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -sortmethod, --sortmethod
      [20171110]Method used to sort the read on query name. (samtools != 
      picard) see https://github.com/samtools/hts-specs/issues/5
      Default: picard
      Possible Values: [samtools, picard]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * compare


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew cmpbams4
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cmpbams/CompareBams4.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/cmpbams/CompareBams4.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **cmpbams4** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example
The following Makefile compare the bam for hg19 and hg38 on chr22 and 21

```Makefile

include ../../config/config.mk

CHROMS=21 22
OUTDIR=tmp

define run
${OUTDIR}/$(1).bam : ${OUTDIR}/$(1).fa.bwt R1.fastq.gz R1.fastq.gz
	${bwa.exe} mem -M   -R '@RG\tID:SAMPLE\tLB:SAMPLE\tSM:SAMPLE\tPL:illumina\tCN:Nantes' ${OUTDIR}/$(1).fa $$(word 2,$$^) $$(word 3,$$^) |  ${samtools.exe} view -b -u -S -F4 - | ${samtools.exe} sort -n -o $$@ -T ${OUTDIR}/$(1)_tmp -

${OUTDIR}/$(1).fa.bwt : ${OUTDIR}/$(1).fa
	${bwa.exe} index $$<

${OUTDIR}/$(1).dict : ${OUTDIR}/$(1).fa
	${java.exe} -jar $(picard.jar) CreateSequenceDictionary R=$$< O=$$@

${OUTDIR}/$(1).fa.fai : ${OUTDIR}/$(1).fa
	${samtools.exe} faidx $$<
	
${OUTDIR}/$(1).fa : 
	mkdir -p $$(dir $$@) && rm -f $$@
	$$(foreach C,${CHROMS}, curl "http://hgdownload.cse.ucsc.edu/goldenPath/$(1)/chromosomes/chr$${C}.fa.gz" | gunzip -c >> $$@;)

endef

all:  ${OUTDIR}/diff.txt 

${OUTDIR}/diff.txt : ${OUTDIR}/hg19.bam ${OUTDIR}/hg38.bam ${OUTDIR}/hg19ToHg38.over.chain ${OUTDIR}/hg19.dict ${OUTDIR}/hg38.dict
	mkdir -p $(dir $@) &&  $(call run_jvarkit,cmpbams4) --novalidchain -st -c $(word 3,$^) $(word 1,$^)  $(word 2,$^) > $@

${OUTDIR}/hg19ToHg38.over.chain :
	mkdir -p $(dir $@) && curl "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" | gunzip -c > $@

$(eval $(call run,hg19))
$(eval $(call run,hg38))

```

#### Output

```
onlyIn	liftover	compareContig	shift	diffCigarOperations	diffNM	diffFlags	diffChroms	Count
BOTH	SameChrom	DiscordantContig	.	-1	0	147/163	chr22/chr21	2
BOTH	SameChrom	SameContig	Gt100	3	15	83/83	chr22/chr22	1
BOTH	SameChrom	DiscordantContig	.	3	5	147/129	chr21/chr22	1
BOTH	SameChrom	SameContig	Gt100	-1	1	163/163	chr22/chr22	22
BOTH	SameChrom	DiscordantContig	.	0	1	83/99	chr21/chr22	32
BOTH	SameChrom	SameContig	Gt100	0	2	99/99	chr21/chr21	22
BOTH	SameChrom	DiscordantContig	.	2	6	81/65	chr22/chr21	1
BOTH	SameChrom	SameContig	Gt100	0	0	185/137	chr22/chr22	20
BOTH	SameChrom	SameContig	Zero	0	0	177/177	chr22/chr22	1417
(...)
```


