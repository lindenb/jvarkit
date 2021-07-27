# WesCnvSvg

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

SVG visualization of bam DEPTH for multiple regions


## Usage

```
Usage: wescnvsvg [options] Files
  Options:
    -css, --css
      custom svg css stylesheet
    -x, --extend
      Extending interval. The following syntaxes are supported: 1000; 1kb; 
      1,000; 30%(shrink); 150% (extend); 0.5 (shrink); 1.5 (extend)
      Default: 0
    -height, --height
      Sample Track height
      Default: 100
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -u, --url, --hyperlink
      creates a hyperlink an area is 'clicked'. creates a hyperlink when 
      'click' in an area. The URL must contains __CHROM__, __START__ and 
      __END__ that will be replaced by their values. Predefined values are 
      'hg19','hg38','igv'. IGV : 
      "http://localhost:60151/goto?locus=__CHROM__%3A__START__-__END__" , 
      UCSC: "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__"
      Default: <empty string>
    -B, --bed, -b, --capture, -r, -rgn, --region, --interval
      A source of intervals. The following suffixes are recognized: vcf, 
      vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise it could be an 
      empty string (no interval) or a list of plain interval separated by '[ 
      \t\n;,]' 
      Default: (unspecified)
    -Q, --mapq
      Min mapping quality
      Default: 1
    -o, --output
      Output file. Optional . Default: stdout
  * -R, --ref, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -smooth, --smooth
      how to smooth data
      Default: AVERAGE
      Possible Values: [AVERAGE, MEDIAN]
    --title
      document title
      Default: WesCnvSvg
    --vcf
      plot VCF data
    --version
      print version and exit
    -w, --width
      Page width
      Default: 1000
    -D
      other parameters. '-Dkey=value'. Undocumented.
      Syntax: -Dkey=value
      Default: {}

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


## Creation Date

20180726

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

https://twitter.com/yokofakun/status/1180046139502059521

![ScreenShot](https://pbs.twimg.com/media/EGBdtO7WoAEHjEE?format=jpg&name=small)


