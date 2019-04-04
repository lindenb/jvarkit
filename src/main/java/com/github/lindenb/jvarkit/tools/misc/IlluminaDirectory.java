/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.tools.misc;


import java.io.BufferedReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.function.Function;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.illumina.FastQName;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.util.CloserUtil;

/**

BEGIN_DOC


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


END_DOC
*/



@Program(name="illuminadir",
	description="Create a structured (**JSON** or **XML**) representation of a directory containing some Illumina FASTQs.",
	keywords={"json","xml","illumina","fastq","workflow"},
	biostars=362767
	)
public class IlluminaDirectory
	extends Launcher
	{
	
	private static final Logger LOG = Logger.build(IlluminaDirectory.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-J","-j","-json","--json"},description="Generate JSON output.")
	private boolean JSON = false;
	@Parameter(names={"-i","--invalid"},description="[20180717]save invalid line/fastq names in this file.")
	private File invalidFile = null;

	
	private final Function<String, String> str2md5 = new StringToMd5();
	
	private int ID_GENERATOR=0;
	private PrintWriter invalidWriter = null;
    
    private class Folder
    	{
    	String projectName="Project1";
    	final SortedMap<String, Sample> sampleMap=new TreeMap<String, IlluminaDirectory.Sample>();
    	final List<Pair> undetermined=new ArrayList<Pair>();
    	void scan(final File f)
    		{
    		if(f==null) return;
    		if(!f.canRead()) return;
    		LOG.info("Scanning "+f);
    		
			FastQName fq=FastQName.parse(f);
			if(!fq.isValid())
				{
				IlluminaDirectory.this.invalidWriter.println(f.getPath());
				LOG.warn("invalid name:"+fq);
				return;
				}
			if(fq.isUndetermined())
				{
				for(int i=0;i< undetermined.size();++i)
	    			{
	    			Pair p=undetermined.get(i);
	    			if(p.complement(fq)) return;
	    			}
				undetermined.add(new Pair(fq));
				}
			else
				{
				Sample sample=this.sampleMap.get(fq.getSample());
				if(sample==null)
					{
					sample=new Sample();
					sample.name=fq.getSample();
					this.sampleMap.put(sample.name,sample);
					}
				sample.add(fq);
				
				
				final File sampleDir = f.getParentFile();
				if(sampleDir!=null && sampleDir.isDirectory() && sampleDir.getName().startsWith("Sample_"))
					{
					final File projDir = sampleDir.getParentFile();
					if(projDir!=null && projDir.isDirectory() && projDir.getName().startsWith("Project_"))
						{
						this.projectName = projDir.getName().substring(8).replace(' ', '_');
						}
					}
				}
				
    			
    		}
    	
    	void json(final JsonWriter out)  throws IOException
    		{
    		out.beginObject();
    		out.name("samples");
    		out.beginArray();

    		for(final Sample S:this.sampleMap.values())
				{
				S.json(out);
				}
    		out.endArray();
    		out.name("undetermined");
    		out.beginArray();
    		for(final Pair p:undetermined)
				{
				p.json(out);
				}
    		out.endArray();
    		out.endObject();
    		}
    	
    	void write(final XMLStreamWriter w) throws XMLStreamException
    		{
    		w.writeStartElement("project");
    		w.writeAttribute("name",this.projectName);
    		w.writeAttribute("center", "CENTER");
    		w.writeAttribute("haloplex", "false");
    		w.writeAttribute("wgs", "false");

    		for(final Sample S:this.sampleMap.values())
    			{
    			S.write(w);
    			}
    		w.writeStartElement("undetermined");
    		for(final Pair p:this.undetermined)
    			{
    			p.write(w);
    			}
    		w.writeEndElement();
    		w.writeEndElement();
    		}
    	}
    
    /** 
     * A pair of fastq , Forward, reverse
     */
    private class Pair
    	{
    	int id;
    	FastQName forward;
    	FastQName reverse;
    	
    	Pair(final FastQName fq)
    		{
    		id=++ID_GENERATOR;
    		switch(fq.getSide())
    			{
    			case Forward:forward=fq; break;
    			case Reverse:reverse=fq; break;
    			default:throw new RuntimeException("bad side "+fq);
    			}
    		}
    	
    	boolean complement(final FastQName other)
    		{
    		if(forward!=null && reverse!=null) return false;
    		if(forward!=null && forward.isComplementOf(other))
    			{
    			reverse=other;
    			return true;
    			}
    		else if(reverse!=null && reverse.isComplementOf(other))
    			{
    			forward=other;
    			return true;
    			}
    		return false;
    		}
    	
    	void json(final JsonWriter out) throws IOException
    		{
    		if(forward!=null && reverse!=null)
				{
	    		out.beginObject();
	    		out.name("id");out.value("p"+this.id);
	    		out.name("md5pair");out.value(str2md5.apply(forward.getFile().getPath()+reverse.getFile().getPath()));
	    		if(forward.hasLane()) {
	    			out.name("lane");out.value(""+forward.getLane());
	    			}
	    		out.name("index");
	    		if(forward.hasSeqIndex())
    				{
    				out.value(forward.getSeqIndex());
    				}
    			else
    				{
    				out.nullValue();
    				}
	    		out.name("split");
	    		if(forward.hasSplit()) {
	    			out.value(""+forward.getSplit());
	    			}
	    		else
					{
					out.nullValue();
					}
	    		
	    		out.name("forward");
	    		
	    		out.beginObject();
	    		out.name("md5filename");out.value(str2md5.apply(forward.getFile().getPath()));
	    		out.name("path");out.value(forward.getFile().getPath());
	    		out.name("side");out.value(forward.getSide().ordinal());
	    		out.name("file-size");out.value(forward.getFile().length());	    		
	    		out.endObject();
	    		
	    		out.name("reverse");
	    		
	    		out.beginObject();
	    		out.name("md5filename");out.value(str2md5.apply(reverse.getFile().getPath()));
	    		out.name("path");out.value(reverse.getFile().getPath());
	    		out.name("side");out.value(reverse.getSide().ordinal());
	    		out.name("file-size");out.value(reverse.getFile().length());	    		
	    		out.endObject();

	    		out.endObject();
				}
    		else
    			{
    			final FastQName F=(forward==null?reverse:forward);
    			out.beginObject();
	    		out.name("id");out.value("p"+this.id);
	    		out.name("md5filename");out.value(str2md5.apply(F.getFile().getPath()));
	    		if(F.hasLane()){
	    			out.name("lane");out.value(""+F.getLane());
	    		}
	    		out.name("index");
	    		if(F.hasSeqIndex())
    				{
    				out.value(F.getSeqIndex());
    				}
    			else
    				{
    				out.nullValue();
    				}
	    		out.name("split");
	    		if(F.hasSplit()) {
	    			out.value(""+F.getSplit());
	    			}
	    		else
	    			{
	    			out.nullValue();
	    			}
	    		out.name("path");out.value(F.getFile().getPath());
	    		out.name("side");out.value(F.getSide().ordinal());

	    		
    			out.endObject();
    			}
			}
    	
    	void write(XMLStreamWriter w,final String tagName,final FastQName fastqFile) throws XMLStreamException
    		{
			w.writeStartElement(tagName);
			w.writeAttribute("md5filename",str2md5.apply(fastqFile.getFile().getPath()));
			w.writeAttribute("file-size",String.valueOf( fastqFile.getFile().length()));
			w.writeCharacters(fastqFile.getFile().getPath());
			w.writeEndElement();
    		}
    	
    	void write(final XMLStreamWriter w) throws XMLStreamException
    		{
			w.writeStartElement("fastq");
			w.writeAttribute("id","p"+this.id);
			w.writeAttribute("md5",str2md5.apply(forward.getFile().getPath()+reverse.getFile().getPath()));
			if(forward.hasLane()) {
				w.writeAttribute("lane", String.valueOf(forward.getLane()));
			}
			if(forward.hasSeqIndex()) {
				w.writeAttribute("index", String.valueOf(forward.getSeqIndex()));
			}
			if(forward.hasSplit()) {
				w.writeAttribute("split", String.valueOf(forward.getSplit()));
				}
			w.writeAttribute("group-id", getGroupId());
			
			if(forward!=null && reverse!=null)
    			{
    			write(w,"for",forward);
    			write(w,"rev",reverse);
    			}
			else
				{
				write(w,"single",forward==null?reverse:forward);
				}
			
		
			w.writeEndElement();
    		}
    	private String getGroupId()
    		{
    		return IlluminaDirectory.this.getGroupId(this.forward);
    		}
    	}
    
    private Map<String,String> groupIdMap=new HashMap<>(); 
    private String getGroupId(final FastQName fastq)
    	{
    	final String s= str2md5.apply( fastq.getSample()+" "+
    			(fastq.hasLane()?fastq.getLane()+" ":"")+
    			(fastq.hasSeqIndex()?fastq.getSeqIndex():"")
    			);
    	String gid = this.groupIdMap.get(s);
    	if(gid==null)
    		{
    		gid = fastq.getSample()+"."+(this.groupIdMap.size()+1);
    		this.groupIdMap.put(s, gid);
    		}
    	return gid;
    	}
    
    
    private class Sample
		{
		String name;
    	final List<Pair> pairs=new ArrayList<Pair>();
    	
    	private void add(final FastQName fq)
    		{
    		for(int i=0;i< pairs.size();++i)
    			{
    			final Pair p=pairs.get(i);
    			if(p.complement(fq)) return;
    			}
    		pairs.add(new Pair(fq));
    		}
    	
    	void write(final XMLStreamWriter w) throws XMLStreamException
			{
			w.writeStartElement("sample");
			w.writeAttribute("name",this.name);
			w.writeAttribute("father","undefined");
			w.writeAttribute("mother","undefined");
			w.writeAttribute("sex","undefined");
			
			for(final Pair p:this.pairs)
				{
				p.write(w);
				}
				
			w.writeEndElement();
			}
    	
    	void json(final JsonWriter out)  throws IOException
    		{
    		out.beginObject();
    		out.name("sample");
    		out.value(this.name);
    		out.name("files");
    		out.beginArray();
    		for(final Pair p: this.pairs)
    			{
    			p.json(out);
    			}
    		out.endArray();
    		out.endObject();
    		}
    	
    	
		}
    
    @Override
    public int doWork(final List<String> args) {
    	   	BufferedReader in=null;
			try
				{
				final String inputName= oneFileOrNull(args);
				in= super.openBufferedReader(inputName);
				
				if(this.invalidFile==null)
					{
					this.invalidWriter = new PrintWriter(new NullOuputStream());
					}
				else
					{
					this.invalidWriter = new PrintWriter(this.invalidFile);
					}
				
				final Folder folder=new Folder();
				String line;
				while((line=in.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					if(!line.endsWith(".fastq.gz"))
						{
						this.invalidWriter.println(line);
						LOG.warn("ignoring "+line+" because it doesn't end with *.fastq.gz");
						continue;
						}
					final File f=new File(line);
					if(!f.exists())
						{
						this.invalidWriter.println(line);
						LOG.error("Doesn't exist:"+f);
						return -1;
						}
					if(!f.isFile())
						{
						this.invalidWriter.println(line);
						LOG.error("Not a file:"+f);
						return -1;
						}
					folder.scan(f);
					}
				in.close();
		    	
				final PrintWriter pw = this.openFileOrStdoutAsPrintWriter(outputFile);

		    	if(this.JSON)
		    		{
		    		final JsonWriter js=new JsonWriter(pw);
		    		folder.json(js);
		    		CloserUtil.close(js);
		    		}
		    	else
		    		{
		    		final XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
		    		final XMLStreamWriter w= xmlfactory.createXMLStreamWriter(pw);
	    			w.writeStartDocument("UTF-8","1.0");
	    			folder.write(w);
	    			w.writeEndDocument();
	    			w.flush();
	    			CloserUtil.close(w);
		    		}
		    	pw.flush();
		    	this.invalidWriter.flush();
		    	CloserUtil.close(pw);
		    	CloserUtil.close(this.invalidWriter);
		    	return RETURN_OK;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(in);
				CloserUtil.close(this.invalidWriter);
				}
			}
   
	public static void main(final String[] args)
		{
		new IlluminaDirectory().instanceMainWithExit(args);
		}

}
