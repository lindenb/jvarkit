# XContaminations

For @AdrienLeger2 : cross contamination between samples in same lane


## Usage

```
Usage: xcontaminations [options] Files
  Options:
    -filter, --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: mapqlt(1) || MapQUnavailable() || Duplicate() || FailsVendorQuality() || NotPrimaryAlignment() || SupplementaryAlignment()
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```

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
$ make xcontaminations
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/xcontamination/XContaminations.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/xcontamination/XContaminations.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **xcontaminations** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```bash
$ find . -type f -name "*.bam" > bam.list
$  head -n 10000 variant.vcf | java -jar dist/xcontaminations.jar - bam.list > out.tsv
$ verticalize out.tsv


>>> 2
$1       #Machine:FlowCell:Run:Lane-1 : HISEQ10:C3FBPACXX:0:4
$2                            sample1 : B00G5V9
$3        Machine:FlowCell:Run:Lane-2 : HISEQ10:C486PACXX:0:3
$4                            sample2 : B00G7LK
$5                          same.lane : 0
$6   reads_sample1_supporting_sample1 : 26392
$7   reads_sample1_supporting_sample2 : 70
$8    reads_sample1_supporting_others : 40
$9   reads_sample2_supporting_sample2 : 21473
$10  reads_sample2_supporting_sample1 : 39
$11    reads_sample2_supporting_other : 31
<<< 2

(...)

>>> 9
$1       #Machine:FlowCell:Run:Lane-1 : HISEQ5:C3FV0ACXX:0:7
$2                            sample1 : B00G738
$3        Machine:FlowCell:Run:Lane-2 : HISEQ5:C3FV0ACXX:0:7
$4                            sample2 : B00G754
$5                          same.lane : 1
$6   reads_sample1_supporting_sample1 : 10209
$7   reads_sample1_supporting_sample2 : 23
$8    reads_sample1_supporting_others : 15
$9   reads_sample2_supporting_sample2 : 9054
$10  reads_sample2_supporting_sample1 : 32
$11    reads_sample2_supporting_other : 9
<<< 9

```

generating in parallel:

```make
CHROMS=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22

.PHONY:all

define xcont

$$(addprefix tmp.xcont.,$$(addsuffix .tsv.gz,$(1))) :
        bcftools view -r "$(1)" unifiedgenotyper.vcf.gz -Tcapture.bed | java -Xmx1g -jar xcontaminations.jar - bal.list  | gzip --best > $$(addsuffix .tmp.gz,$$@) && mv  
$$(addsuffix .tmp.gz,$$@) $$@


endef

all: $(foreach C,${CHROMS},$(addprefix tmp.xcont.,$(addsuffix .tsv.gz,${C})))
        $(foreach I,0 1, gunzip -c $^ |  awk -F '       ' '($$5==$I)'  |awk -F '        ' 'BEGIN {T=0;N=0;} {for(i=6;i<=NF;++i) T+=int($$i); N+=int($$7); N+=int($$10);} E
ND { printf("%f\n",N/T);}'; )

$(foreach C,${CHROMS},$(eval $(call xcont,$C)))
```




