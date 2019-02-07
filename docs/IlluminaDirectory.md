# IlluminaDirectory

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Create a structured (**JSON** or **XML**) representation of a directory containing some Illumina FASTQs.


## Usage

```
Usage: illuminadir [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -i, --invalid
      [20180717]save invalid line/fastq names in this file.
    -J, -j, -json, --json
      Generate JSON output.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * json
 * xml
 * illumina
 * fastq
 * workflow



## See also in Biostars

 * [https://www.biostars.org/p/362767](https://www.biostars.org/p/362767)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew illuminadir
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/IlluminaDirectory.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/IlluminaDirectory.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **illuminadir** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Motivation

Illuminadir scans folders , search for FASTQs and generate a structured summary of the files (xml or json).
Currently only tested with HiSeq data.

### History

* 20180717: supports bcl2fq2
* 20171128: supports double indexing.

### Examples

```
$ find dir1 dir2 -type f -name "*.fastq.gz" |\
   java  -jar dist/illuminadir.jar | \
   xsltproc xml2script.xslt > script.bash
(...)
```

#### XML output
 
The XML ouput looks like this:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<illumina>
  <directory path="RUN62_XFC2DM8ACXX/data">
    <samples>
      <sample name="SAMPLE1">
        <pair md5="cd4b436ce7aff4cf669d282c6d9a7899" lane="8" index="ATCACG" split="2">
          <fastq md5filename="3369c3457d6603f06379b654cb78e696" side="1" path="RUN62_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R1_002.fastq.gz" file-size="359046311"/>
          <fastq md5filename="832039fa00b5f40108848e48eb437e0b" side="2" path="RUN62_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R2_002.fastq.gz" file-size="359659451"/>
        </pair>
        <pair md5="b3050fa3307e63ab9790b0e263c5d240" lane="8" index="ATCACG" split="3">
          <fastq md5filename="091727bb6b300e463c3d708e157436ab" side="1" path="RUN62_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R1_003.fastq.gz" file-size="206660736"/>
          <fastq md5filename="20235ef4ec8845515beb4e13da34b5d3" side="2" path="RUN62_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R2_003.fastq.gz" file-size="206715143"/>
        </pair>
        <pair md5="9f7ee49e87d01610372c43ab928939f6" lane="8" index="ATCACG" split="1">
          <fastq md5filename="54cb2fd33edd5c2e787287ccf1595952" side="1" path="RUN62_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R1_001.fastq.gz" file-size="354530831"/>
          <fastq md5filename="e937cbdf32020074e50d3332c67cf6b3" side="2" path="RUN62_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R2_001.fastq.gz" file-size="356908963"/>
        </pair>
        <pair md5="0697846a504158eef523c0f4ede85288" lane="7" index="ATCACG" split="2">
          <fastq md5filename="6fb35d130efae4dcfa79260281504aa3" side="1" path="RUN62_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L007_R1_002.fastq.gz" file-size="357120615"/>
(...)
      <pair md5="634cbb29ca64604174963a4fff09f37a" lane="7" split="1">
        <fastq md5filename="bc0b283a58946fd75a95b330e0aefdc8" side="1" path="RUN62_XFC2DM8ACXX/data/Undetermined_indices/Sample_lane7/lane7_Undetermined_L007_R1_001.fastq.gz" file-size="371063045"/>
        <fastq md5filename="9eab26c5b593d50d642399d172a11835" side="2" path="RUN62_XFC2DM8ACXX/data/Undetermined_indices/Sample_lane7/lane7_Undetermined_L007_R2_001.fastq.gz" file-size="372221753"/>
      </pair>
      <pair md5="bf31099075d6c3c7ea052b8038cb4a03" lane="8" split="2">
        <fastq md5filename="f229389da36a3efc20888bffdec09b80" side="1" path="RUN62_XFC2DM8ACXX/data/Undetermined_indices/Sample_lane8/lane8_Undetermined_L008_R1_002.fastq.gz" file-size="374331268"/>
        <fastq md5filename="417fd9f28d24f63ce0d0808d97543315" side="2" path="RUN62_XFC2DM8ACXX/data/Undetermined_indices/Sample_lane8/lane8_Undetermined_L008_R2_002.fastq.gz" file-size="372181102"/>
      </pair>
      <pair md5="95cab850b0608c53e8c83b25cfdb3b2b" lane="8" split="3">
        <fastq md5filename="23f5be8a962697f50e2a271394242e2f" side="1" path="RUN62_XFC2DM8ACXX/data/Undetermined_indices/Sample_lane8/lane8_Undetermined_L008_R1_003.fastq.gz" file-size="60303589"/>
        <fastq md5filename="3f39f212c36d0aa884b81649ad56630c" side="2" path="RUN62_XFC2DM8ACXX/data/Undetermined_indices/Sample_lane8/lane8_Undetermined_L008_R2_003.fastq.gz" file-size="59123627"/>
      </pair>
      <pair md5="ab108b1dda7df86f33f375367b86bfe4" lane="8" split="1">
        <fastq md5filename="14f8281cf7d1a53d29cd03cb53a45b4a" side="1" path="RUN62_XFC2DM8ACXX/data/Undetermined_indices/Sample_lane8/lane8_Undetermined_L008_R1_001.fastq.gz" file-size="371255111"/>
        <fastq md5filename="977fd388e1b3451dfcdbf9bdcbb89ed4" side="2" path="RUN62_XFC2DM8ACXX/data/Undetermined_indices/Sample_lane8/lane8_Undetermined_L008_R2_001.fastq.gz" file-size="370744530"/>
      </pair>
    </undetermined>
  </directory>
</illumina>
```

How to use that file ? here is a  example of **XSLT** stylesheet that can generate a **Makefile** to generate a **LaTex** about the number of reads per Lane/Sample/Index:


```xslt
<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	version='1.0' 
	> 
<xsl:output method="text"/>


<xsl:template match="/">
.PHONY:all clean

all: report.pdf

report.pdf: report.tex 
	pdflatex $&lt;

report.tex : all.count
	echo 'T&lt;-read.table("$&lt;",head=TRUE,sep="\t");$(foreach FTYPE,Index Sample Lane, T2&lt;-tapply(T$$count,T$$${FTYPE},sum);png("${FTYPE}.png");barplot(T2,las=3);dev.off();)' | R --no-save
	echo "\documentclass{report}" &gt; $@
	echo "\ usepackage{graphicx}" &gt;&gt; $@
	echo "\date{\today}" &gt;&gt; $@
	echo "\title{FastQ Report}" &gt;&gt; $@
	echo "\begin{document}" &gt;&gt; $@
	echo "\maketitle" &gt;&gt; $@
	$(foreach FTYPE,Index Sample Lane, echo "\section{By ${FTYPE}}#\begin{center}#\includegraphics{${FTYPE}.png}#\end{center}" | tr "#" "\n" &gt;&gt; $@ ; )
	echo "\end{document}" &gt;&gt; $@
	


all.count : $(addsuffix .count, <xsl:for-each select="//fastq" ><xsl:value-of select="@md5filename"/><xsl:text> </xsl:text></xsl:for-each>) 
	echo -e "Lane\tsplit\tside\tsize\tcount\tIndex\tSample"  &gt; $@ &amp;&amp; \
	cat $^ &gt;&gt; $@

<xsl:apply-templates select="//fastq" mode="count"/>

clean:
	rm -f all.count report.pdf report.tex $(addsuffix .count, <xsl:for-each select="//fastq" ><xsl:value-of select="@md5filename"/><xsl:text> </xsl:text></xsl:for-each>) 

</xsl:template>

<xsl:template match="fastq" mode="count">
$(addsuffix .count, <xsl:value-of select="@md5filename"/>): <xsl:value-of select="@path"/>
	gunzip -c $&lt; | awk '(NR%4==1)' | wc -l  | xargs  printf "<xsl:value-of select="../@lane"/>\t<xsl:value-of select="../@split"/>\t<xsl:value-of select="@side"/>\t<xsl:value-of select="@file-size"/>\t%s\t<xsl:choose><xsl:when test="../@index"><xsl:value-of select="../@index"/></xsl:when><xsl:otherwise>Undetermined</xsl:otherwise></xsl:choose>\t<xsl:choose><xsl:when test="../../@name"><xsl:value-of select="../../@name"/></xsl:when><xsl:otherwise>Undetermined</xsl:otherwise></xsl:choose>\n"   &gt; $@

</xsl:template>
</xsl:stylesheet>
```





```
$ xsltproc  illumina.xml illumina2makefile.xsl > Makefile
```

output:

```makefile
.PHONY:all clean

all: report.pdf

report.pdf: report.tex 
	pdflatex $<

report.tex : all.count
	echo 'T<-read.table("$<",head=TRUE,sep="\t");$(foreach FTYPE,Index Sample Lane, T2<-tapply(T$$count,T$$${FTYPE},sum);png("${FTYPE}.png");barplot(T2,las=3);dev.off();)' | R --no-save
	echo "\documentclass{report}" > $@
	echo "\ usepackage{graphicx}" >> $@
	echo "\date{\today}" >> $@
	echo "\title{FastQ Report}" >> $@
	echo "\begin{document}" >> $@
	echo "\maketitle" >> $@
	$(foreach FTYPE,Index Sample Lane, echo "\section{By ${FTYPE}}#\begin{center}#\includegraphics{${FTYPE}.png}#\end{center}" | tr "#" "\n" >> $@ ; )
	echo "\end{document}" >> $@



all.count : $(addsuffix .count, 3369c3457d6603f06379b654cb78e696 832039fa00b5f40108848e48eb437e0b 091727bb6b300e463c3d708e157436ab 20235ef4ec88....)
	echo -e "Lane\tsplit\tside\tsize\tcount\tIndex\tSample"  > $@ && \
	cat $^ >> $@


$(addsuffix .count, 3369c3457d6603f06379b654cb78e696): RUN62_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R1_002.fastq.gz
	gunzip -c $< | awk '(NR%4==1)' | wc -l  | xargs  printf "8\t2\t1\t359046311\t%s\tATCACG\tSAMPLE1\n"   > $@


$(addsuffix .count, 832039fa00b5f40108848e48eb437e0b): RUN62_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R2_002.fastq.gz
	gunzip -c $< | awk '(NR%4==1)' | wc -l  | xargs  printf "8\t2\t2\t359659451\t%s\tATCACG\tSAMPLE1\n"   > $@
(....)
```

####  JSON output

The JSON output looks like this


```json
{"directory":"RUN62_XFC2DM8ACXX/data","samples":[{"sample":"SAMPLE1","files":[{
"md5pair":"cd4b436ce7aff4cf669d282c6d9a7899","lane":8,"index":"ATCACG","split":2
,"forward":{"md5filename":"3369c3457d6603f06379b654cb78e696","path":"20131001_SN
L149_0062_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R1_002.fastq.g
z","side":1,"file-size":359046311},"reverse":{"md5filename":"832039fa00b5f401088
48e48eb437e0b","path":"20131001_SNL149_0062_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/
SAMPLE1_ATCACG_L008_R2_002.fastq.gz","side":2,"file-size":359659451}},{"md5pair"
:"b3050fa3307e63ab9790b0e263c5d240","lane":8,"index":"ATCACG","split":3,"forward
":{"md5filename":"091727bb6b300e463c3d708e157436ab","path":"20131001_SNL149_0062
_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R1_003.fastq.gz","side"
:1,"file-size":206660736},"reverse":{"md5filename":"20235ef4ec8845515beb4e13da34
b5d3","path":"20131001_SNL149_0062_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_A
TCACG_L008_R2_003.fastq.gz","side":2,"file-size":206715143}},{"md5pair":"9f7ee49
e87d01610372c43ab928939f6","lane":8,"index":"ATCACG","split":1,"forward":{"md5fi
lename":"54cb2fd33edd5c2e787287ccf1595952","path":"20131001_SNL149_0062_XFC2DM8A
CXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L008_R1_001.fastq.gz","side":1,"file-
size":354530831},"reverse":{"md5filename":"e937cbdf32020074e50d3332c67cf6b3","pa
th":"20131001_SNL149_0062_XFC2DM8ACXX/data/OUT/Sample_SAMPLE1/SAMPLE1_ATCACG_L00
8_R2_001.fastq.gz","side":2,"file-size":356908963}},{"md5pair":"0697846a504158ee
f523c0f4ede85288","lane":7,"index":"ATCACG","split":2,"forward":{"md5filename":"
```

It can be processed using a tool like [jsvelocity](https://github.com/lindenb/jsvelocity) to generate the same kind of Makefile:

The velocity template for jsvelocity (https://github.com/lindenb/jsvelocity)


```
#macro(maketarget $fastq)

$(addsuffix .count, ${fastq.md5filename}): ${fastq.path}
	gunzip -c $< | awk '(NR%4==1)' | wc -l  | xargs  printf "${fastq.parentNode.lane}\t${fastq.parentNode.split}\t${fastq.side}\t${fastq['file-size']}\t%s\t#if(${fastq.parentNode.containsKey("index")})${fastq.parentNode.index}#{else}Undetermined#{end}\t#if(${fastq.parentNode.parentNode.containsKey("name")})${fastq.parentNode.parentNode.name}#{else}Undetermined#{end}\n"   > $@

#end

.PHONY:all clean

all: report.pdf

report.pdf: report.tex 
	pdflatex $<

report.tex : all.count
	echo 'T<-read.table("$<",head=TRUE,sep="\t");$(foreach FTYPE,Index Sample Lane, T2<-tapply(T$$count,T$$${FTYPE},sum);png("${FTYPE}.png");barplot(T2,las=3);dev.off();)' | R --no-save
	echo "\documentclass{report}" > $@
	echo "\ usepackage{graphicx}" >> $@
	echo "\date{\today}" >> $@
	echo "\title{FastQ Report}" >> $@
	echo "\begin{document}" >> $@
	echo "\maketitle" >> $@
	$(foreach FTYPE,Index Sample Lane, echo "\section{By ${FTYPE}}#\begin{center}#\includegraphics{${FTYPE}.png}#\end{center}" | tr "#" "\n" >> $@ ; )
	echo "\end{document}" >> $@

all.count : $(addsuffix .count, #foreach($dir in $all) #foreach($sample in ${dir.samples})#foreach($pair in ${sample.files}) ${pair.forward.md5filename}  ${pair.reverse.md5filename} #end #end #foreach($pair in   ${dir.undetermined}) ${pair.forward.md5filename}  ${pair.reverse.md5filename} #end  #end )



#foreach($dir in $all)
#foreach($sample in ${dir.samples})
#foreach($pair in ${sample.files})
#maketarget($pair.forward)
#maketarget($pair.reverse)
#end
#end
#foreach($pair in   ${dir.undetermined})
#maketarget($pair.forward)
#maketarget($pair.reverse)
#end 
#end


clean:
	rm -f all.count  $(addsuffix .count,  #foreach($dir in $all)
#foreach($sample in ${dir.samples})
#foreach($pair in ${sample.files}) ${pair.forward.md5filename}  ${pair.reverse.md5filename}  #end #end
#foreach($pair in   ${dir.undetermined}) ${pair.forward.md5filename}  ${pair.reverse.md5filename}  #end  #end )

```

transform using jsvelocity:

```
java -jar dist/jsvelocity.jar \
     -d all illumina.json \
      illumina.vm > Makefile
```

output: same as above


