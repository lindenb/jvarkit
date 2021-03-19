# CopyNumber01

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

experimental CNV detection.


## Usage

```
Usage: copynumber01 [options] Files
  Options:
    --bed, --capture
      Exome Capture as BED
    --gcDepthInterpolation
      Method to interpolate GC% and depth. See https://commons.apache.org/proper/commons-math/javadocs/api-3.0/org/apache/commons/math3/analysis/interpolation/UnivariateInterpolator.html
      Default: loess
      Possible Values: [loess, neville, difference, linear, spline, identity]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --mapq
      Min mapping quality
      Default: 1
    --max-depth
      Treat depth greater than this value as 'weird' and discard the sliding 
      windows at this place.
      Default: 500
    --max-gc
      Max GC%
      Default: 1.0
    --min-depth
      Treat depth lower than this value as 'weird' and discard the sliding 
      windows at this place.
      Default: 0
    --min-gc
      Min GC%
      Default: 0.0
    -o, --out
      Output file. Optional . Default: stdout
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --sex
      Sexual contigs, comma or space separated
      Default: chrX,chrY,X,Y
    --smooth
      Smooth normalized depth window. smooth normalized depth with the 'n' 
      neightbours 
      Default: 5
    --smooth-distance
      When using --smooth. Only merge if windows are within that distance.A 
      distance specified as a positive integer.Commas are removed. The 
      following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 1000
    --univariateDepth
      How to calculate depth in a BAM interval.
      Default: mean
      Possible Values: [mean, median]
    --univariateGC
      Loess needs only one GC value: we need to merge Depth with same GC%. How 
      do we merge ?
      Default: median
      Possible Values: [mean, median]
    --univariateMid
      Depth normalization. Used when we want to normalize the depths between 
      0.0 and 1.0
      Default: median
      Possible Values: [mean, median]
    --univariateSmooth
      How to smooth data with the --smooth option.
      Default: mean
      Possible Values: [mean, median]
    --version
      print version and exit
    --win-min
      Discard window where length on reference is lower than 'x'. A distance 
      specified as a positive integer.Commas are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 100
    -s, --win-shift
      window shift. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 500
    -w, --win-size
      window size. A distance specified as a positive integer.Commas are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 1000

```


## Keywords

 * cnv
 * bam
 * sam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew copynumber01
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20140201

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/redon/CopyNumber01.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/redon/CopyNumber01.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **copynumber01** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example:

```
$ java -jar dist/copynumber01.jar  -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam
[INFO][CopyNumber01]sorting...
[INFO][CopyNumber01]fill gc%
[INFO][CopyNumber01]remove high/low gc%
[INFO][CopyNumber01]Getting coverage for RF01 N=6
[INFO][CopyNumber01]Getting coverage for RF02 N=4
[INFO][CopyNumber01]Getting coverage for RF03 N=4
[INFO][CopyNumber01]Getting coverage for RF04 N=4
[INFO][CopyNumber01]Getting coverage for RF05 N=2
[INFO][CopyNumber01]Getting coverage for RF06 N=2
[INFO][CopyNumber01]Getting coverage for RF07 N=1
[INFO][CopyNumber01]Getting coverage for RF08 N=1
[INFO][CopyNumber01]Getting coverage for RF09 N=1
[INFO][CopyNumber01]removed 0. now N=25
[INFO][CopyNumber01]median norm depth : 8.331950991034539
#CHROM	START	END	Sample	IDX	GC	RAW-DEPTH	NORM-DEPTH
RF01	0	1001	S1	0	0.321	6.410	1.015
RF01	500	1501	S1	500	0.349	8.446	1.015
RF01	1000	2001	S1	1000	0.371	9.479	1.015
RF01	1500	2501	S1	1500	0.374	9.445	1.015
RF01	2000	3001	S1	2000	0.354	7.921	1.015
RF01	2500	3302	S1	2500	0.331	5.766	1.015
RF02	0	1001	S1	3302	0.347	7.380	0.998
RF02	500	1501	S1	3802	0.348	9.189	0.998
RF02	1000	2001	S1	4302	0.341	8.672	0.998
RF02	1500	2501	S1	4802	0.344	7.977	0.998
RF03	0	1001	S1	5989	0.314	7.060	0.931
RF03	500	1501	S1	6489	0.332	9.967	0.931
RF03	1000	2001	S1	6989	0.319	9.193	0.931
RF03	1500	2501	S1	7489	0.315	7.012	0.931
RF04	0	1001	S1	8581	0.352	6.119	1.021
RF04	500	1501	S1	9081	0.344	9.554	1.021
RF04	1000	2001	S1	9581	0.374	10.000	1.021
RF04	1500	2362	S1	10081	0.359	7.396	1.021
RF05	0	1001	S1	10943	0.326	8.827	0.980
RF05	500	1501	S1	11443	0.327	8.996	0.980
RF06	0	1001	S1	12522	0.363	8.571	1.035
RF06	500	1356	S1	13022	0.421	8.384	1.035
RF07	0	1001	S1	13878	0.329	7.900	0.996
RF08	0	1001	S1	14952	0.358	7.876	1.008
RF09	0	1001	S1	16011	0.374	7.919	1.089
```

