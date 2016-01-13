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

See also [[Compilation]].

```bash
$  make <xsl:value-of select="translate(@app,$uppercase,$lowercase)"/>
```

##Synopsis

```
$ java -jar dist/<xsl:value-of select="translate(@app,$uppercase,$lowercase)"/>.jar  [options ] (stdin|file) 
```

<xsl:apply-templates select="c:options"/>


##Source Code

Main code is: https://github.com/lindenb/jvarkit/blob/master/src/main/java/<xsl:value-of select="translate(@package,'.','/')"/>/<xsl:value-of select="@app"/>.java

<xsl:apply-templates select="c:documentation"/>

## Contribute

- Issue Tracker: http://github.com/lindenb/jvarkit/issues
- Source Code: http://github.com/lindenb/jvarkit

##History

* 2016 : Creation

## License

The project is licensed under the MIT license.

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

<xsl:template match="h:ul">
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
<xsl:if test="/c:app/c:snippet[@id='javascript']">
<xsl:text>  * -e (expression) javascript expression
  * -f (file) javascript file
</xsl:text>
</xsl:if>
<xsl:text>  * -o,--output &lt;FILENAME&gt; output file.
  * -h,--help print help
  * -version,--version show version and exit
</xsl:text>

</xsl:template>

<xsl:template match="c:option">
<xsl:text>  * -</xsl:text>
<xsl:value-of select="@opt"/>
<xsl:text> </xsl:text>
<xsl:apply-templates select="c:description"/>

</xsl:template>


</xsl:stylesheet>
