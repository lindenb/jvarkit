# DepthOfCoverage

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

A custom 'Depth of Coverage'.


## Usage

```
Usage: depthofcoverage [options] Files
  Options:
    --auto-mask
      Use REFerence sequence to automatically mask bases that are not ATGC
      Default: false
    -B, --bed
      optional bed containing regions to be SCANNED (inverse of --mask)
    -ct, --ct
      summary Coverage Threshold. A 'range of integers' is a list of integers 
      in ascending order separated with semicolons.
      Default: [[-Inf/0[, [0/5[, [5/10[, [10/20[, [20/30[, [30/40[, [40/50[, [50/100[, [100/200[, [200/300[, [300/400[, [400/500[, [500/1000[, [1000/2000[, [2000/3000[, [3000/4000[, [4000/5000[, [5000/Inf[]
    --disable-paired-overlap
      Disable: Count overlapping bases with mate for paired-end
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
       min mapping quality.
      Default: 1
    -M, --mask
      optional bed containing regions to be MASKED
    --max-depth
      Ignore depth if it is bigger than this value.
      Default: 10000000
    --min-length
      Chromosomes to skip if their length is lower than this value.
      Default: 0
    -o, --out
      Output file. Optional . Default: stdout
    -R, --reference
      For reading CRAM. Indexed fasta Reference file. This file must be 
      indexed with samtools faidx and with picard CreateSequenceDictionary
    --skip
      Chromosomes to skip (regular expression)
      Default: (NC_007605|hs37d5)
    --use-index
      use bam index to query intervals if --bed is defined.
      Default: false
    --version
      print version and exit

```


## Keywords

 * depth
 * bam
 * sam
 * coverage


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew depthofcoverage
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20190927

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/coverage/DepthOfCoverage.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/coverage/DepthOfCoverage.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/coverage/DepthOfCoverageTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/coverage/DepthOfCoverageTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **depthofcoverage** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java  -jar dist/depthofcoverage.jar -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 2> /dev/null  | column -t 

#BAM                       Sample  Contig  Length  Count   Depth
src/test/resources/S1.bam  S1      RF01    3302    25037   7.582374318594791
src/test/resources/S1.bam  S1      RF02    2687    20275   7.545589877186453
src/test/resources/S1.bam  S1      RF03    2592    19583   7.55516975308642
src/test/resources/S1.bam  S1      RF04    2362    17898   7.577476714648603
src/test/resources/S1.bam  S1      RF05    1579    11887   7.528182393920202
src/test/resources/S1.bam  S1      RF06    1356    10201   7.522861356932153
src/test/resources/S1.bam  S1      RF07    1074    8115    7.555865921787709
src/test/resources/S1.bam  S1      RF08    1059    7980    7.5354107648725215
src/test/resources/S1.bam  S1      RF09    1062    7980    7.5141242937853105
src/test/resources/S1.bam  S1      RF10    751     5740    7.6431424766977365
src/test/resources/S1.bam  S1      RF11    666     5037    7.563063063063063
src/test/resources/S1.bam  S1      *       18490   139733  7.557220118983234
src/test/resources/S2.bam  S2      RF01    3302    25030   7.580254391278014
src/test/resources/S2.bam  S2      RF02    2687    20272   7.544473390398213
src/test/resources/S2.bam  S2      RF03    2592    19592   7.558641975308642
src/test/resources/S2.bam  S2      RF04    2362    17916   7.585097375105843
src/test/resources/S2.bam  S2      RF05    1579    11892   7.531348955034832
src/test/resources/S2.bam  S2      RF06    1356    10217   7.534660766961652
src/test/resources/S2.bam  S2      RF07    1074    8112    7.553072625698324

```

