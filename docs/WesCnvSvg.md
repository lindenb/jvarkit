# WesCnvSvg

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

SVG visualization of bam DEPTH for multiple regions


## Usage

```
Usage: wescnvsvg [options] Files
  Options:
    -cap, --cap
      Cap coverage to this value. Negative=don't set any limit
      Default: -1
    -B, --bed, -b, --capture
      BED Capture. BED file containing the Regions to be observed.
    -css, --css
      custom svg css stylesheet
    -x, --extend
      Extend each region in the bed by 'x' bases. If the argument ends with 
      '%' it is interpreted as a percentage.
    --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: Accept All/ Filter out nothing
    -height, --height
      Sample Track height
      Default: 100
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -u, --url, --hyperlink
      creates a hyperlink when 'click' in an area. The URL must contains 
      __CHROM__, __START__ and __END__ that will be replaced by their values. 
      IGV : "http://localhost:60151/goto?locus=__CHROM__%3A__START__-__END__" 
      , UCSC: "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__"
      Default: none
    -r, -rgn, --region, --interval
      Interval regions: 'CHR:START-END'. multiple separated with spaces or 
      semicolon 
    -o, --output
      Output file. Optional . Default: stdout
    -p, -percentile, --percentile
      How to compute the percentil of a region
      Default: AVERAGE
      Possible Values: [MIN, MAX, MEDIAN, AVERAGE, RANDOM, SUM]
    -R, --ref, --reference
      The parameter is the path to an Indexed fasta Reference file. This fasta 
      file must be indexed with samtools faidx and with picard 
      CreateSequenceDictionary. The parameter can also be a 'key' (matching 
      the regular expression `[A-Za-z][A-Za-z0-9_\\-]*`) in a catalog file. A 
      'catalog' file is a java property file ( 
      https://docs.oracle.com/javase/tutorial/essential/environment/properties.html 
      ) where the values are the path to the fasta file.  Catalogs are 
      searched in that order : `${PWD}/fasta-ref.properties`, 
      `${HOME}/.fasta-ref.properties`, `/etc/jvarkit/fasta-ref.properties`.  
      If the key or the path are not defined by the user, they will be 
      searched in that order 1) the java property 
      -Djvarkit.fasta.reference=pathTofastaOrCatalogKey . 2) the linux 
      environement variable $FASTA_REFERENCE=pathTofastaOrCatalogKey 3) The 
      catalogs. 
      Default: <<Default Fasta Reference Supplier>>
    -smooth, --smooth
      Smoothing pixel window size. Negative=don't smooth
      Default: 100
    --title
      document title
      Default: WesCnvSvg
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
 * wes
 * bed
 * capture
 * exome


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew wescnvsvg
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2svg/WesCnvSvg.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2svg/WesCnvSvg.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2svg/WesCnvSvgTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2svg/WesCnvSvgTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **wescnvsvg** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

input is a set of bam file or a file with the '*.list' suffix containing the path to the bam files.

## Example

```
$ find dir -name "*.bam"  > bam.list
$ java -jar dist/wescnvsvg.jar -R ref.fasta -B cnv.bed bam.list > cnv.svg 
```

## Screenshots

https://twitter.com/yokofakun/status/1022503372669300738 : 

![ScreenShot](https://pbs.twimg.com/media/DjCpKcYXgAAq4fw.jpg:large)

https://twitter.com/yokofakun/status/1022805656905150464

![ScreenShot](https://pbs.twimg.com/media/DjG8Do0XsAA4U46.jpg:large)

https://twitter.com/yokofakun/status/1023953180927963137

![ScreenShot](https://pbs.twimg.com/media/DjXP8_4X0AAtSQZ.jpg)

https://twitter.com/yokofakun/status/1024315948847849472

![ScreenShot](https://pbs.twimg.com/media/DjcZ6vNXcAEYjEt.jpg)

https://twitter.com/yokofakun/status/1024315948847849472

![ScreenShot](https://pbs.twimg.com/media/DjcZ6vNXcAEYjEt.jpg)

https://twitter.com/yokofakun/status/1025330308193779712

![ScreenShot](https://pbs.twimg.com/media/Djq0Se-W0AEAbyR.jpg)

https://twitter.com/yokofakun/status/1040592885786263554

![ScreenShot](https://pbs.twimg.com/media/DnDttNgX4AAtxax.jpg)

https://twitter.com/yokofakun/status/1040577235856580608

![ScreenShot](https://pbs.twimg.com/media/DnDfaGLXcAArg0P.jpg)

https://twitter.com/yokofakun/status/1057625407913111557

![ScreenShot](https://pbs.twimg.com/media/Dq1whOTX0AAzkZc.jpg)


