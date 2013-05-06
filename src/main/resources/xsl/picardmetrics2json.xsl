<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	version='1.0' 
	xmlns:p='http://picard.sourceforge.net/'
	xml:xsi2="http://www.w3.org/2001/XMLSchema-instance"
	> 
<xsl:output method="text"/>

<xsl:template match="p:picard-metrics">
<xsl:text>{</xsl:text>
<xsl:apply-templates select="p:metrics-file"/>
<xsl:text>}</xsl:text>
</xsl:template>


<xsl:template match="p:metrics-file">
<xsl:text>"</xsl:text>
<xsl:apply-templates select="@file"/>
<xsl:text>":{</xsl:text>
<xsl:for-each select="p:headers|p:metrics|p:histogram">
<xsl:if test="position()&gt;1">,</xsl:if>
<xsl:text>"</xsl:text>
<xsl:value-of select="local-name(.)"/>
<xsl:text>":</xsl:text>
<xsl:apply-templates select="."/>
</xsl:for-each>
<xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="p:headers">
<xsl:text>[</xsl:text>
<xsl:for-each select="p:header">
<xsl:if test="position()&gt;1">,</xsl:if>
<xsl:apply-templates select="."/>
</xsl:for-each>
<xsl:text>]</xsl:text>
</xsl:template>

<xsl:template match="p:header">
<xsl:text>{"class":"</xsl:text>
<xsl:value-of select="@class"/>
<xsl:text>","value":"</xsl:text>
<xsl:value-of select="."/>
<xsl:text>"}</xsl:text>
</xsl:template>

<xsl:template match="p:metrics|p:histogram">
<xsl:variable name="ths" select="p:thead/p:th"/>
<xsl:text>[</xsl:text>
<xsl:for-each select="p:tbody/p:tr">
<xsl:variable name="tr" select="."/>
<xsl:if test="position()&gt;1">,</xsl:if>
<xsl:text>{</xsl:text>

<xsl:for-each select="$ths">
<xsl:variable name="thidx" select="position()"/>
<xsl:variable name="td" select="$tr/p:td[$thidx]"/>
<xsl:if test="$thidx&gt;1">,</xsl:if>

<xsl:text>"</xsl:text>
<xsl:value-of select="."/>
<xsl:text>":</xsl:text>
<xsl:choose>
	<xsl:when test="count($td/child::node())=0">null</xsl:when>
	<xsl:when test="number($td)=$td or $td='true' or $td='false'"><xsl:value-of select="$td"/></xsl:when>
	<xsl:otherwise>
		<xsl:text>"</xsl:text>
		<xsl:value-of select="$td"/>
		<xsl:text>"</xsl:text>
	</xsl:otherwise>
</xsl:choose>


</xsl:for-each>
<xsl:text>}</xsl:text>
</xsl:for-each>
<xsl:text>]</xsl:text>
</xsl:template>

</xsl:stylesheet>

