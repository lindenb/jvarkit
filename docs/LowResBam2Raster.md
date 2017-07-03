# LowResBam2Raster

Low Resolution BAM to raster graphics


## Usage

```
Usage: lowresbam2raster [options] Files
  Options:
    -clip, --clip
      Show clipping
      Default: false
    -depth, --depth
      Depth track height.
      Default: 100
    -gcPercent, --gcPercent
      GC% track height.
      Default: 100
    -gcwin, --gcWindowSize
      GC% Window size
      Default: 10
    --groupby
      Group Reads by
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --highlight
      hightligth those positions.
      Default: []
    -kg, --knownGene
      UCSC knownGene URI. Beware chromosome names are formatted the same as 
      your REFERENCE. A typical KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)
    --mapqopacity
      How to handle the MAPQ/ opacity of the reads. all_opaque: no opacity, 
      handler 1: transparency under MAPQ=60
      Default: handler1
      Possible Values: [all_opaque, handler1]
    --limit, --maxrows
      Limit number of rows to 'N' lines. negative: no limit.
      Default: -1
    -minh, --minh
      Min. distance between two reads.
      Default: 10
    -noSuppl, --noSuppl
      Hide arcs of Supplementary alignments.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -printNames, --printNames
      Print Read Names (for debugging)
      Default: false
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
  * -r, --region
      Restrict to that region. An interval as the following syntax : 
      "chrom:start-end" or "chrom:middle+extend"  or 
      "chrom:start-end+extend".A program might use a Reference sequence to fix 
      the chromosome name (e.g: 1->chr1)
    -srf, --samRecordFilter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    --spaceyfeature
      number of pixels between features
      Default: 1
    -V, --variants, --vcf
      VCF files used to fill the position to hightlight with POS
      Default: []
    --version
      print version and exit
    -w, --width
      Image width
      Default: 1000

```


## Keywords

 * bam
 * alignment
 * graphics
 * visualization
 * png
 * knowngene


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
$ make lowresbam2raster
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/LowResBam2Raster.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/LowResBam2Raster.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **lowresbam2raster** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
java -jar dist/lowresbam2raster.jar \
	-o out.png -r "22:38999+10000" in.bam \
	 -clip -srf "" -R ref.fasta  -kg knownGene.txt.gz
```

## Screenshot

![https://pbs.twimg.com/media/DAldDxvXkAAGMoJ.jpg](https://pbs.twimg.com/media/DAldDxvXkAAGMoJ.jpg)


