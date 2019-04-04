# SimplePlot

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

simple figure plotter output is a R script


## Usage

```
Usage: simpleplot [options] Files
  Options:
    -chrompos, --chrom-position
      When reading a genomic file. Input is not a BED file but a CHROM(tab)POS 
      file. e.g: output of `samtools depth`
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hide-legend
      Hide Legend
      Default: false
    --min-reference-size
      When using a *.dict file, discard the contigs having a length < 'size'. 
      Useful to discard the small unused contigs like 'chrM'. -1 : ignore.
      Default: -1
    -nh, --no-header
      There is no header. Tested for Histograms
      Default: false
    -o, --out
      Output file. If defined, save the picture in this file extension:png, 
      jpg or R (experimental) and close the application.
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    -su, --sort-unique
      For PIE or SIMPLE_HISTOGRAM the input is the output of a `sort | uniq 
      -c` pipeline
      Default: false
    --title
      Chart Title
      Default: <empty string>
    --type, -t
      type
      Default: UNDEFINED
      Possible Values: [UNDEFINED, BEDGRAPH, PIE, SIMPLE_HISTOGRAM, HISTOGRAM, STACKED_HISTOGRAM, STACKED_HISTOGRAM_PIVOTED, XYV, STACKED_XYV, HEATMAP]
    --version
      print version and exit
    -xlab, -xlabel, --xlabel
      X axis label.
    -ylab, -ylabel, --ylabel
      Y axis label.

```


## Keywords

 * char
 * figure


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew simpleplot
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SimplePlot.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SimplePlot.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SimplePlotTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/SimplePlotTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **simpleplot** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Plot types:

Year  X  Y
2018  1  2
2019  3  4


## Examples

```
$ gunzip -c in.vcf.gz | grep -v "^#" | cut -f 1 | sort | uniq -c | java -jar dist/simpleplot.jar -t SIMPLE_HISTOGRAM  -su 
$ gunzip -c in.vcf.gz | grep -v "^#" | cut -f 1 | sort | uniq -c | java -jar dist/simpleplot.jar -t PIE  -su
```

### HISTOGRAM

```
Year  X  Y
2018  1  2
2019  3  4
```

```
$ echo -e "Year\tX\tY\n2018\t1\t2\n2019\t3\t4" | java -jar dist/simpleplot.jar -t HISTOGRAM
```

produces a histogram with two series in the legend (2018 and 2019). On the X-axis: 4 items X-2018, X-2019, Y-2018, Y-2019



### STACKED_HISTOGRAM

```
Year  X  Y
2018  1  2
2019  3  4
```

```
$ echo -e "Year\tX\tY\n2018\t1\t2\n2019\t3\t4" | java -jar dist/simpleplot.jar -t STACKED_HISTOGRAM
``` 

produces a histogram with two series in the legend (2018 and 2019). On the X-axis: 2 items X( 2018 under-2019), Y (2018 under 2019)


### STACKED_HISTOGRAM_PIVOTED

```
Year  X  Y
2018  1  2
2019  3  4
```

```
$  echo -e "2018\t1\t2\n2019\t3\t4\n2020\t1\t10" | java -jar dist/simpleplot.jar -t STACKED_HISTOGRAM_PIVOTED -nh
``` 

produces a histogram with two series in the legend ($1 and $2). On the X-axis: 3 items 2018, 2019, 2020. Each with the stacked $1 and $2



### Example

plot base=function(position) in a fastq:

```
gunzip -c src/test/resources/S1.R1.fq.gz | \
	awk '(NR%4==2) {L=length($0);for(i=1;i<=L;i++) printf("%d\t%s\n",i,substr($0,i,1));}' |\
	sort | uniq -c |\
	java -jar dist/simpleplot.jar -su -t STACKED_XYV --xlabel "Position"
```


### Example

HeatMap

```
echo -e "A\tA\t1\nA\tB\t2\nB\tA\t3\nB\tB\t10" | java -jar dist/simpleplot.jar  -t HEATMAP
```


## History

  * 2019: removed jfx as openjdk doesn't support it... output is now R
 
 
