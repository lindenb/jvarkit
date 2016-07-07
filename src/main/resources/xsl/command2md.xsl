<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns:h="http://www.w3.org/1999/xhtml"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	>
<xsl:output method='text'/>
<xsl:variable name="lowercase">abcdefghijklmnopqrstuvwxyz</xsl:variable>
<xsl:variable name="uppercase">ABCDEFGHIJKLMNOPQRSTUVWXYZ</xsl:variable>


<xsl:template match="/">
<xsl:apply-templates select="c:app"/>
</xsl:template>


<xsl:template match="c:app">

##Motivation

<xsl:apply-templates select="c:description"/>

##Compilation

### Requirements / Dependencies

Since 2016-05-30 the compilation of the "Java API for high-throughput sequencing data (HTS) formats" (htsjdk) library requires gradle http://gradle.org.

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* gradle http://gradle.org is only required to compile the "Java API for high-throughput sequencing data (HTS) formats" (htsjdk). And I think htsjdk installs it.
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make <xsl:value-of select="@jarname"/>
```

by default, the libraries are not included in the jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ). You can create a bigger but standalone executable jar by adding `standalone=yes` on the command line:


```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make <xsl:value-of select="@jarname"/> standalone=yes
```

The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```

to set the gradle user home ( https://docs.gradle.org/current/userguide/build_environment.html#sec:gradle_properties_and_system_properties )

```
gradle.user.home=/dir1/dir2/gradle_user_home
```

<xsl:if test="not(documentation/h:h3[text() = 'Synopsis'])">

##Synopsis

```bash
$ java -jar dist/<xsl:value-of select="@jarname"/>.jar  [options] (stdin|file<xsl:choose>
	<xsl:when test="c:input/@type='vcf'">.vcf|file.vcf.gz</xsl:when>
	<xsl:when test="c:input/@type='sam' or c:input/@type='bam'">.bam|file.sam</xsl:when>
	<xsl:otherwise></xsl:otherwise>
</xsl:choose>) 
```
</xsl:if>


<xsl:if test="c:snippet[@id='concatenated-vcf']">
## VCF Concatenation

This tool supports concatenated VCF on input.
For each VCf on input, a VCF will be printed in the output.
If there is more than one VCF on input, the output will be **NOT** a valid VCF but a file containing concatenated VCFs.


```bash
$ gunzip -c f1.vcf.gz f2.vcf.gz f3.vcf.gz |\
	java -jar dist/<xsl:value-of select="@jarname"/>.jar > out_3_vcfs.txt
```

If the filename is defined and ends with '.zip': each output VCF will be stored into one zip entry.

```bash
$ gunzip -c f1.vcf.gz f2.vcf.gz f3.vcf.gz |\
	java -jar dist/<xsl:value-of select="@jarname"/>.jar -o result.zip
```


</xsl:if>

<xsl:apply-templates select="c:options"/>

##Source Code

Main code is: https://github.com/lindenb/jvarkit/blob/master/src/main/java/<xsl:value-of select="translate(@package,'.','/')"/>/<xsl:value-of select="@app"/>.java

<xsl:apply-templates select="c:documentation"/>

## Contribute

- Issue Tracker: http://github.com/lindenb/jvarkit/issues
- Source Code: http://github.com/lindenb/jvarkit

## License

The project is licensed under the MIT license.

## Citing

Should you cite **<xsl:value-of select="@jarname"/>** ? https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md

The current reference is:

http://dx.doi.org/10.6084/m9.figshare.1425030

&gt; Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
&gt; http://dx.doi.org/10.6084/m9.figshare.1425030


</xsl:template>


<xsl:template match="c:documentation">
<xsl:apply-templates />
</xsl:template>


<xsl:template match="h:pre">
<xsl:text>
```</xsl:text>
<xsl:value-of select="@class"/>
<xsl:text>
</xsl:text>
<xsl:apply-templates/>
<xsl:text>
```
</xsl:text>

</xsl:template>

<xsl:template match="h:h4">
#### <xsl:apply-templates/>

</xsl:template>


<xsl:template match="h:h3">
### <xsl:apply-templates/>

</xsl:template>


<xsl:template match="h:h2">
## <xsl:apply-templates/>

</xsl:template>

<xsl:template match="h:b">
<xsl:text>**</xsl:text>
<xsl:apply-templates/>
<xsl:text>**</xsl:text>
</xsl:template>

<xsl:template match="h:i">
<xsl:text>**</xsl:text>
<xsl:apply-templates/>
<xsl:text>**</xsl:text>
</xsl:template>


<xsl:template match="h:p|h:div">
<xsl:text>
</xsl:text>
<xsl:apply-templates/>
<xsl:text>
</xsl:text>
</xsl:template>

<xsl:template match="h:ul|h:ol">
<xsl:text>
</xsl:text>
<xsl:apply-templates select="h:li"/>
<xsl:text>
</xsl:text>
</xsl:template>

<xsl:template match="h:code">
<xsl:text>`</xsl:text>
<xsl:apply-templates/>
<xsl:text>`</xsl:text>
</xsl:template>


<xsl:template match="h:span">
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="h:br">
<xsl:text>
</xsl:text>
</xsl:template>

<xsl:template match="h:li">
<xsl:text>  * </xsl:text>
<xsl:apply-templates/>
<xsl:text>
</xsl:text>
</xsl:template>


<xsl:template match="h:table">
<xsl:apply-templates select="h:tr"/>
</xsl:template>

<xsl:template match="h:tr">
<xsl:for-each select="h:th|h:td">
<xsl:if test="position()&gt;1 and position()&lt;= last()">  |  </xsl:if>
<xsl:apply-templates select="."/>
</xsl:for-each>
<xsl:text>
</xsl:text>
</xsl:template>


<xsl:template match="h:img[@src]">
<xsl:text>![</xsl:text>
<xsl:value-of select="@src"/>
<xsl:text>](</xsl:text>
<xsl:value-of select="@src"/>
<xsl:text>)</xsl:text>
</xsl:template>


<xsl:template match="h:a[not(@href)]">
<xsl:text>[</xsl:text>
<xsl:value-of select="."/>
<xsl:text>](</xsl:text>
<xsl:value-of select="."/>
<xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="h:a[@href]">
<xsl:text>[</xsl:text>
<xsl:choose>
	<xsl:when test="string-length(normalize-space(.))&gt;0">
		<xsl:apply-templates/>
	</xsl:when>
	<xsl:otherwise>
		<xsl:value-of select="@href"/>
	</xsl:otherwise>
</xsl:choose>
<xsl:value-of select="."/>
<xsl:text>](</xsl:text>
<xsl:value-of select="@href"/>
<xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="c:options">

<xsl:text>## Options
</xsl:text>
<xsl:apply-templates select="c:option"/>
<xsl:text>  * -h|--help print help
  * -version|--version show version and exit
</xsl:text>


</xsl:template>

<xsl:template match="c:option">
<xsl:text>  * -</xsl:text>
<xsl:value-of select="@opt"/>
<xsl:text>|--</xsl:text>
<xsl:value-of select="@longopt"/>
<xsl:text> </xsl:text>
<xsl:choose>
	<xsl:when test="@type='boolean'"></xsl:when>
	<xsl:when test="@argname"> (<xsl:value-of select="@argname"/>) </xsl:when>
	<xsl:when test="@arg-name"> (<xsl:value-of select="@arg-name"/>) </xsl:when>
	<xsl:otherwise> (VALUE) </xsl:otherwise>
</xsl:choose>
<xsl:apply-templates select="c:description"/>

<xsl:if test="@default">
<xsl:text> Default value : "</xsl:text>
<xsl:value-of select="@default"/>
<xsl:text>".</xsl:text>
</xsl:if>

<xsl:text>
</xsl:text>
</xsl:template>


</xsl:stylesheet>
