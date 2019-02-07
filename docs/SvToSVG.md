# SvToSVG

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

BAM to SVG. Used to display the structural variations.


## Usage

```
Usage: sv2svg [options] Files
  Options:
    --coverage, --depth
      Coverage height. Don't print if cov '<=0'.
      Default: 70
    -d, --duration
      Animation duration, in secs
      Default: 10
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
  * -r, -i, --interval, --region
      interval CHROM:START-END
      Default: []
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary. Optional: if defined 
      will be used to display the mismatches.
    --repeat-count
      SVG animation repeat count
      Default: indefinite
    --variant, -V
      optional indexed VCF file.
    --version
      print version and exit
    -w, --width
      Page width
      Default: 1000

```


## Keywords

 * bam
 * alignment
 * graphics
 * visualization
 * svg


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew sv2svg
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2svg/SvToSVG.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2svg/SvToSVG.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2svg/SvToSVGTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2svg/SvToSVGTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sv2svg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ samtools view -b "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA20845/high_coverage_alignment/NA20845.wgs.ILLUMINA.bwa.GIH.high_cov_pcr_free.20140203.bam" "7:8352157-8356597" > jeter.bam && samtools index jeter.bam
$ java -jar dist/svg2svg.jar jeter.bam > jeter.svg
```

### A translocation


Translocation described in [PMC5932280](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5932280/)

> Identification of Balanced Chromosomal Rearrangements Previously Unknown Among Participants in the 1000 Genomes Project: Implications for Interpretation of Structural Variation in Genomes and the Future of Clinical Cytogenetics

```
$ samtools view -b "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02260/alignment/HG02260.mapped.ILLUMINA.bwa.PEL.low_coverage.20130415.bam" "9:137229907-137231907" > jeter1.bam
$ samtools view -b "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02260/alignment/HG02260.mapped.ILLUMINA.bwa.PEL.low_coverage.20130415.bam" "14:79838174-79840174" > jeter2.bam
$ samtools merge jeter3.bam jeter1.bam jeter2.bam
$ samtools index jeter3.bam
$ java -jar dist/sv2svg.jar -r "9:137229907-137231907" -r "14:79838174-79840174"  jeter3.bam > jeter.svg
```


## Gallery

[https://gist.github.com/lindenb/bf48989b8da31eeafdc2caa0694361eb](https://gist.github.com/lindenb/bf48989b8da31eeafdc2caa0694361eb)

[https://twitter.com/yokofakun/status/1063369955406725120](https://twitter.com/yokofakun/status/1063369955406725120)

![https://imgur.com/EVOrXuc](https://i.imgur.com/EVOrXuc.gif)

[https://twitter.com/yokofakun/status/1063511215832539136](https://twitter.com/yokofakun/status/1063511215832539136)

![https://pbs.twimg.com/media/DsJZKRrWsAE_QqA.jpg](https://pbs.twimg.com/media/DsJZKRrWsAE_QqA.jpg)

[https://gist.github.com/lindenb/877d1d00d9f19c618f2d8505a2fe5614](https://gist.github.com/lindenb/877d1d00d9f19c618f2d8505a2fe5614)

[https://twitter.com/yokofakun/status/1064484537059684355](https://twitter.com/yokofakun/status/1064484537059684355)

![https://twitter.com/yokofakun/status/1064484537059684355](https://pbs.twimg.com/media/DsXObwuXoAAepSw.jpg)

[https://twitter.com/yokofakun/status/1064503996285681666](https://twitter.com/yokofakun/status/1064503996285681666)

![https://pbs.twimg.com/media/DsXgnKrWwAAYiBS.jpg](https://pbs.twimg.com/media/DsXgnKrWwAAYiBS.jpg)
