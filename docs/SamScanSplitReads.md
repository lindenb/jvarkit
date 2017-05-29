# SamScanSplitReads

scan split reads


## Usage

```
Usage: samscansplitreads [options] Files
  Options:
    --defaultSampleName
      Default Sample name if not read group
      Default: UNDEFINED
    -x, --extend
      extends interval by 'x' pb before merging.
      Default: 20
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -msr, --minSupportingReads
      Minimal number of supporting reads.
      Default: 0
    --normalize
      Optional. Normalize count to this value. e.g '1'
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * sam
 * sv
 * splitreads


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
$ make samscansplitreads
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamScanSplitReads.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamScanSplitReads.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samscansplitreads** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example:

content of circos.conf 

```
karyotype = /commun/data/packages/circos/circos-0.69-2/data/karyotype/karyotype.human.hg19.txt
chromosomes_units = 1000000


<<include ideogram.conf>>
<<include ticks.conf>>

<links>
<link>
file          = __FILENAME__
radius        = 0.99r
bezier_radius = 0r
color         = black_a4
thickness     = 2
</link>
</links>


<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>

```

content of ideogram.conf 

```
<ideogram>

<spacing>
## espace entre les contig
default = 0r
</spacing>

radius    = 0.95r
thickness = 30p
fill      = yes

## stroke du radius
stroke_color     = red
stroke_thickness = 3p

# show labels (chromosome names )
show_label       = yes
label_font       = default 
label_radius     = 1r + 75p
label_size       = 30
label_parallel   = yes


</ideogram>
```

content of ticks.conf

```
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p

# the tick label is derived by multiplying the tick position
# by 'multiplier' and casting it in 'format':
#
# sprintf(format,position*multiplier)
#

multiplier       = 1e-6

# %d   - integer
# %f   - float
# %.1f - float with one decimal
# %.2f - float with two decimals
#
# for other formats, see http://perldoc.perl.org/functions/sprintf.html

format           = %d

<tick>
spacing        = 5u
size           = 10p
</tick>

<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>

</ticks>
```

generate circos:


```
$ java -jar dist/samscansplitreads.jar -msr 2 input.bam > input.dat
$ cut -f1-7 input.dat | awk '{OFS="\t";f=$7/50.0;if(f>1.0) f=1.0;if(f<0.5) f=0.5;$1=sprintf("hs%s",$1);$4=sprintf("hs%s",$4);$7=sprintf("thickness=%s,color=%s",f,($1==$4?"blue":"red"));print;}' > tmp.txt
cat circos.conf | sed 's%__FILENAME__%tmp.txt%' >  tmp.txt.conf
circos-0.69-2/bin/circos  -outputdir ./  -outputfile output  -conf  tmp.txt.conf
```




