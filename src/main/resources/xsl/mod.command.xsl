<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	>
<xsl:output method="text"/>

<xsl:template match="c:app" mode="header">/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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


<xsl:apply-templates select="c:history"/>

*/
</xsl:template>

<xsl:template match="c:app" mode="package">
<xsl:value-of select="@package"/>
</xsl:template>

<xsl:template match="c:app" mode="class-name">
<xsl:value-of select="@app"/>
</xsl:template>

<xsl:template match="c:app" mode="abstract-class-name">
<xsl:text>Abstract</xsl:text>
<xsl:apply-templates select="." mode="class-name"/>
</xsl:template>

<xsl:template match="c:app" mode="factory-class-name">
<xsl:apply-templates select="." mode="class-name"/>
<xsl:text>Factory</xsl:text>
</xsl:template>

<xsl:template match="c:option-group" mode="cli">
final OptionGroup <xsl:value-of select="generate-id()"/> = new OptionGroup();
<xsl:apply-templates select="c:option" mode="cli"/>
options.addOptionGroup(<xsl:value-of select="generate-id()"/>);
</xsl:template>

<xsl:template match="c:option" mode="cli">
options.addOption(org.apache.commons.cli.Option
	<xsl:choose>
		<xsl:when test="@opt">
			.builder("<xsl:value-of select="@opt"/>")
		</xsl:when>
		<xsl:otherwise>
			.builder("<xsl:value-of select="@name"/>")
		</xsl:otherwise>
	</xsl:choose>
	<xsl:if test="@required">
	.required(<xsl:value-of select="@required"/>)
	</xsl:if>
	<xsl:if test="@optional-arg">
	.optionalArg(<xsl:value-of select="@optional-arg"/>)
	</xsl:if>
	<xsl:if test="@longOpt">
	.longOpt("<xsl:value-of select="@longOpt"/>")
	</xsl:if>
	<xsl:if test="@value-separator">
	.valueSeparator('<xsl:value-of select="@value-separator"/>')
	</xsl:if>
	<xsl:if test="@description">
	.desc("<xsl:value-of select="@description"/>").
	</xsl:if>
	<xsl:if test="c:description">
	.desc("<xsl:value-of select="c:description"/>").
	</xsl:if>
		<xsl:choose>
		<xsl:when test="@arg-name">
			.argName("<xsl:value-of select="@arg-name"/>")
		</xsl:when>
		<xsl:otherwise>
			.argName("<xsl:value-of select="@name"/>")
		</xsl:otherwise>
	</xsl:choose>
	<xsl:choose>
		<xsl:when test="@type='int' or @type='long' or @type='short' or @type='double' or @type='float'  or @type='number'">
		.type(PatternOptionBuilder.NUMBER_VALUE)
		</xsl:when>
		<xsl:when test="@type='url'">
		.type(PatternOptionBuilder.STRING_VALUE)
		</xsl:when>
		<xsl:when test="@type='object'">
		.type(PatternOptionBuilder.OBJECT_VALUE)
		</xsl:when>
		<xsl:when test="@type='date'">
		.type(PatternOptionBuilder.DATE_VALUE)
		</xsl:when>
		<xsl:when test="@type='file'">
		.type(PatternOptionBuilder.FILE_VALUE)
		</xsl:when>
		<xsl:when test="@type='existing-file'">
		.type(PatternOptionBuilder.EXISTING_FILE_VALUE)
		</xsl:when>
		<xsl:otherwise>
			<xsl:message terminate="yes">unknown type <xsl:value-of select="@type"/></xsl:message>
		</xsl:otherwise>
		
	</xsl:choose>
	.build() );	
</xsl:template>


<xsl:template match="c:option[@type='long']" mode="validator">
LongValidator <xsl:apply-templates select="@name"/> 
</xsl:template>


<xsl:template match="c:option" mode="name">
<xsl:value-of select="@name"/>
</xsl:template>


<xsl:template match="c:option">
<xsl:variable name="cloneable">
	<xsl:apply-templates select="." mode="cloneable"/>
</xsl:variable>

/** option <xsl:apply-templates select="." mode="name"/> */
private <xsl:value-of select="@type"/><xsl:text> </xsl:text> <xsl:apply-templates select="." mode="name"/> = <xsl:choose>
		<xsl:when test="@default"><xsl:value-of select="@default"/></xsl:when>
		<xsl:otherwise>0</xsl:otherwise>
	</xsl:choose>;

/** getter for <xsl:value-of select="@name"/> */
public <xsl:value-of select="@type"/> <xsl:text> </xsl:text><xsl:apply-templates select="." mode="getter"/>()
	{
	return this.<xsl:apply-templates select="." mode="name"/>;
	}

/** setter for <xsl:value-of select="@name"/> */
public <xsl:apply-templates select="." mode="setter"/>( final <xsl:value-of select="@type"/><xsl:text> </xsl:text><xsl:apply-templates select="." mode="name"/>)
	{
	this.<xsl:apply-templates select="." mode="name"/> = <xsl:choose>
		<xsl:when test="$cloneable = 'true'">(<xsl:value-of select="@type"/>)(<xsl:apply-templates select="." mode="name"/>==null?null:<xsl:apply-templates select="." mode="name"/>.clone())</xsl:when>
		<xsl:otherwise><xsl:apply-templates select="." mode="name"/></xsl:otherwise>
		</xsl:choose>;
	}

</xsl:template>

<xsl:template match="c:option" mode="cloneable">
<xsl:variable name="nilleable">
	<xsl:apply-templates select="." mode="nilleable"/>
</xsl:variable>
<xsl:choose>
	<xsl:when test="@type='java.net.URL'">false</xsl:when>
	<xsl:when test="starts-with(@type,'java.lang')">false</xsl:when>
	<xsl:when test="$nilleable = 'true'">true</xsl:when>
</xsl:choose>
</xsl:template>

<xsl:template match="c:option" mode="nilleable">
<xsl:choose>
	<xsl:when test="@type='java.net.URL'">true</xsl:when>
    <xsl:when test="starts-with(@type,'java.lang')">true</xsl:when>
	<xsl:when test="@type='int' or @type='double'">false</xsl:when>
	<xsl:otherwise>true</xsl:otherwise>
</xsl:choose>
</xsl:template>


</xsl:stylesheet>


