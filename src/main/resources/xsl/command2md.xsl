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
$  make <xsl:value-of select="translate(@jarname)"/>
```

by default, the libraries are not included in the jar file, so you shouldn't move them. You can create
a bigger but standalone executable jar by addinging `standalone=yes` on the command line:


```bash
$  make <xsl:value-of select="translate(@jarname)"/> standalone=yes
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
<xsl:text>  * -h,--help print help
  * -version,--version show version and exit
</xsl:text>

</xsl:template>

<xsl:template match="c:option">
<xsl:text>  * -</xsl:text>
<xsl:value-of select="@opt"/>
<xsl:text> </xsl:text>
<xsl:apply-templates select="c:description"/>
<xsl:text>
</xsl:text>
</xsl:template>


</xsl:stylesheet>
