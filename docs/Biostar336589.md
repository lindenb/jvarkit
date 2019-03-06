# Biostar336589

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

displays circular map as SVG from BED and REF file


## Usage

```
Usage: biostar336589 [options] Files
  Options:
    -css, --css
      custom svg css file
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -hist, --histogram
      histogram mode: score of each bed item must be defined. Items must not 
      overlap 
      Default: false
    -u, --url, --hyperlink
      creates a hyperlink when 'click' in an area. The URL must contains 
      __CHROM__, __START__ and __END__ that will be replaced by their values. 
      IGV : "http://localhost:60151/goto?locus=__CHROM__%3A__START__-__END__" 
      , UCSC: "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__"
      Default: none
    -o, --out
      Output file. Optional . Default: stdout
    --score-end-color
      When coloring with bed/score this is the color end.  A named color 
      ('red', 'blue'...) use the syntax 'rgb(int,int,int)'.
      Default: black
    --score-start-color
      When coloring with bed/score this is the color start.  A named color 
      ('red', 'blue'...) use the syntax 'rgb(int,int,int)'.
      Default: white
    --title
      document title
      Default: Biostar336589
    --version
      print version and exit
    --width
      Linear SVG size. When defined, this option switches to a 
      linear(!=circular) view
      Default: -1
  * -R
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -a
      rotate for 'x' seconds. ignore if <=0
      Default: -1
    -as
      arrow size
      Default: 10.0
    -ch
      contig height
      Default: 10.0
    -da
      distance between arcs
      Default: 10.0
    -fh
      arc height
      Default: 10.0
    -md
      min distance in bp between two features on the same arc.
      Default: 100
    -mr
      min internal radius
      Default: 100.0
    -ms
      skip chromosome reference length lower than this value. ignore if <=0
      Default: -1

```


## Keywords

 * genome
 * browser
 * circular
 * bed
 * svg



## See also in Biostars

 * [https://www.biostars.org/p/336589](https://www.biostars.org/p/336589)
 * [https://www.biostars.org/p/367522](https://www.biostars.org/p/367522)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew biostar336589
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar336589.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar336589.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar336589Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar336589Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar336589** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## input

input is a BED file (file or stdin). https://genome.ucsc.edu/FAQ/FAQformat.html#format1

  * column 1: chrom
  * column 2: start (you'll get faster results if the input is sorted on chrom/start )
  * column 3: end
  * column 4 is the name of the feature
  * column 5 is the score [0-1000] or '.'
  * column 6 strand +/-
  * column 7 ignored
  * column 8 ignored
  * column 9 is '.' or R,G,B (as in the bed specification) or it's treated as a full svg:style (e.g: `fill:red;stroke:blue;` ) 


multiple bed files are splitted into 'tracks'.

## Example


https://gist.github.com/lindenb/b6debad569dcb5112e76da893d68dd81

```
$ wget -O - -q  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz" |\
	gunzip -c | awk '{printf("%s\t%s\t%s\t%s\t%d\t+\t.\t.\t%s\n",$2,$3,$4,$8,rand()*1000,NR%20==0?"255,0,250":".");}' |\
	java -jar dist/biostar336589.jar -R src/test/resources/human_b37.dict   --url 'http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=__CHROM__%3A__START__-__END__' --title gaps -mr 300 -fh 20 > ~/jeter.svg 
```

```
$ wget -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" |\
	gunzip -c | cut -f 3,5,6 |\
	sort -t $'\t' -k1,1V -k2,2n |\
	bedtools merge |\
	java -jar dist/biostar336589.jar -md 10000 \
	 	-R src/test/resources/human_b37.dict > out.svg
```

https://gist.github.com/lindenb/5250750014441cc36586dd1f47ed0e37

## Example

```
cat src/test/resources/rotavirus_rf.fa.fai  | awk '{if(NR==1) srand(); printf("%s\t0\t%s\t%s\t%d\t%s\n",$1,$2,$1,rand()*1000,rand()<0.5?rand()<0.5?".":"+":"-");}'   > tmp1.bed
cat src/test/resources/rotavirus_rf.fa.fai  | awk '{if(NR==1) srand(); printf("%s\t0\t%s\t%s\t%d\t%s\n",$1,$2,$1,rand()*1000,rand()<0.5?rand()<0.5?".":"+":"-");}'   > tmp2.bed
cat src/test/resources/rotavirus_rf.fa.fai  | awk '{if(NR==1) srand(); printf("%s\t0\t%s\t%s\t%d\t%s\n",$1,$2,$1,rand()*1000,rand()<0.5?rand()<0.5?".":"+":"-");}'   > tmp3.bed


java -jar dist/biostar336589.jar -R src/test/resources/rotavirus_rf.fa -a 60   tmp1.bed tmp2.bed tmp3.bed > out.svg
```

https://gist.github.com/lindenb/7fd1fa1d3dedcfe38e009387d9f8579c


## Screenshot


https://twitter.com/yokofakun/status/1038060108373286912

![https://twitter.com/yokofakun/status/1038060108373286912](https://pbs.twimg.com/media/Dmft0cSXoAAp78l.jpg)

https://twitter.com/yokofakun/status/1039054467889590272

![https://pbs.twimg.com/media/Dmt2GyvWsAAHfvY.jpg](https://pbs.twimg.com/media/Dmt2GyvWsAAHfvY.jpg)


