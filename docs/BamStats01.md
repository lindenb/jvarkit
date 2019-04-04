# BamStats01

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Statistics about the reads in a BAM.


## Usage

```
Usage: samstats01 [options] Files
  Options:
    -B, --bed
      capture bed file. Optional
    --groupby
      Group Reads by. Data partitioning using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -q, --qual
      min mapping quality
      Default: 30.0
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
$ ./gradlew samstats01
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats01.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bamstats01/BamStats01.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samstats01** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### History



* Dec 2013 Added PROPER_PAIR_HMQ for @SolenaLS
* Dec 2013 Added X and Y for @SolenaLS



### Output


See also: http://picard.sourceforge.net/explain-flags.html


#### Counts


* TOTAL : total number of reads (not PAIRS of reads)
* PAIRED: total number of reads paired (should be equals to ALL for a paired-end assay)
* UNMAPPED : count unmapped reads 
* MAPPED  : count mapped reads
* PROPER_PAIR  : count reads in proper pair (forward+reverse, same chromosome, distance is OK)
* PROPER_PAIR_HMQ  : proper pairs with mapping quality >= user qual
* PLUS : reads on plus strand
* MINUS : reads on minus strand
* PRIMARY_ALIGNMENT : alignment flagged as primary alignment (not alternative position)
* FAIL_MAPPING_QUALITY : MAQ < user qual
* DUPLICATE : the flag 'duplicate' was set
* FAIL_VENDOR_QUALITY : the flag "read fails platform/vendor quality checks" was set
* OK_FOR_PE_CALLING : reads ok for Paired-end mapping ( properly paired, not dup, not fails_vendor_qual,  not fails_mapping_qual, primary align )
* X and Y : number of reads mapping the chromosomes X/chrX and Y/chrY


#### Categories


* ALL: all reads
* IN_TARGET: reads overlapping user's BED (if provided)
* OFF_TARGET: reads with no overlap with user's BED (if provided)



### Example


```
$  java -jar dist/bamstats01.jar \
		-B capture.bed my.bam \
		

(...)
#Filename	Sample	ALL_TOTAL	ALL_PAIRED	ALL_UNMAPPED	ALL_MAPPED	ALL_PROPER_PAIR	ALL_PLUS_STRAND	ALL_MINUS_STRAND	ALL_PRIMARY_ALIGNMENT	ALL_FAIL_MAPPING_QUALITY	ALL_DUPLICATE	ALL_FAIL_VENDOR_QUALITY	IN_TARGET_TOTAL	IN_TARGET_PAIRED	IN_TARGET_UNMAPPED	IN_TARGET_MAPPED	IN_TARGET_PROPER_PAIR	IN_TARGET_PLUS_STRAND	IN_TARGET_MINUS_STRAND	IN_TARGET_PRIMARY_ALIGNMENT	IN_TARGET_FAIL_MAPPING_QUALITY	IN_TARGET_DUPLICATE	IN_TARGET_FAIL_VENDOR_QUALITY	OFF_TARGET_TOTAL	OFF_TARGET_PAIRED	OFF_TARGET_UNMAPPED	OFF_TARGET_MAPPED	OFF_TARGET_PROPER_PAIR	OFF_TARGET_PLUS_STRAND	OFF_TARGET_MINUS_STRAND	OFF_TARGET_PRIMARY_ALIGNMENT	OFF_TARGET_FAIL_MAPPING_QUALITY	OFF_TARGET_DUPLICATE	OFF_TARGET_FAIL_VENDOR_QUALITY
my.bam	Sample	1617984	1617984	3966	1614018	1407862	806964	807054	1614018	56980	0	0	1293922	1293922	0	1293922	1133808	644741	649181	1293922	14087	0	0	320096	320096	0	320096	274054	162223	157873	320096	42893	0	0
(...)

```


