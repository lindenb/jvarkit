# BedNonOverlappingSet

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Split a Bed file into non-overlapping data set.


## Usage

```
Usage: bednonoverlappingset [options] Files
  Options:
    -x, --extend
      Extend intervals by 'x' bases
      Default: 0
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -m, --manifset
      Manifest file file containing the generated filenames/number of item.
  * -o, --out
      Output file. Filename *Must* contains the word __SETID__ and end with 
      '.bed' .
    -R, -r, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary If defined, will be used 
      to sort the bed record on chrom/pos before writing the bed records.
    --version
      print version and exit

```


## Keywords

 * bed


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bednonoverlappingset
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BedNonOverlappingSet.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/BedNonOverlappingSet.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/BedNonOverlappingSetTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/BedNonOverlappingSetTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bednonoverlappingset** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

GATK DepthOfCoverage merge overlapping segments (see https://gatkforums.broadinstitute.org/gatk/discussion/1865/ ). I wan't to get the coverage for a set of overlapping windows.

## EXAMPLES

### Example

```bash
$ awk '{printf("%s\t0\t%s\n",$1,$2);}' src/test/resources/rotavirus_rf.fa.fai |\
bedtools makewindows  -w 500 -s 100 -b - |\
java -jar dist/bednonoverlappingset.jar  -o tmp.__SETID__.bed -m tmp.manifest

[INFO][BedNonOverlappingSet]saving tmp.00001.bed 44
[INFO][BedNonOverlappingSet]saving tmp.00002.bed 39
[INFO][BedNonOverlappingSet]saving tmp.00003.bed 37
[INFO][BedNonOverlappingSet]saving tmp.00004.bed 36
[INFO][BedNonOverlappingSet]saving tmp.00005.bed 33

$ head -n 2 tmp.000*.bed
==> tmp.00001.bed <==
RF01	0	500
RF01	500	1000

==> tmp.00002.bed <==
RF01	100	600
RF01	600	1100

==> tmp.00003.bed <==
RF01	200	700
RF01	700	1200

==> tmp.00004.bed <==
RF01	300	800
RF01	800	1300

==> tmp.00005.bed <==
RF01	400	900
RF01	900	1400

$ cat tmp.manifest 
tmp.00001.bed	44
tmp.00002.bed	39
tmp.00003.bed	37
tmp.00004.bed	36
tmp.00005.bed	33
```

### Example

```bash
(...)
java -jar dist/bednonoverlappingset.jar -x 1 -R ref.fa -o "tmp.__SETID__.bed" -m tmp.manifest input.bed

cut -f 1 tmp.manifest | while read B
do
	${java_exe}   -Djava.io.tmpdir=.  -jar GenomeAnalysisTK.jar \
	   -T DepthOfCoverage -R "ref.fa" \
	   -o "SAMPLE"  -I input.bam  -L "${B}" --omitDepthOutputAtEachBase --omitLocusTable --omitPerSampleStats
 	   grep -v '^Target' "${sample}.sample_interval_summary" | awk -F '	' '{printf("%s\\t%s\\n",\$1,\$3);}' >> tmp.tsv	
	
	   rm "SAMPLE.sample_interval_summary"  "SAMPLE.sample_interval_statistics"
done
	
LC_ALL=C sort -t '	' -k1,1 tmp.tsv >> "SAMPLE.win.cov.tsv"

(...)
```

