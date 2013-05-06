<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	version='1.0' 
	xmlns:p='http://picard.sourceforge.net/'
	xml:xsi2="http://www.w3.org/2001/XMLSchema-instance"
	> 
<xsl:output method="html"/>

<xsl:template match="p:picard-metrics">
<html>
<head>
 <script type="text/javascript" src="https://www.google.com/jsapi"></script>
 <script type="text/javascript">
 <xsl:apply-templates select="//p:histogram" mode="google"/>
 </script>
</head>
<body>
<xsl:apply-templates select="p:metrics-file"/>
</body>
</html>
</xsl:template>


<xsl:template match="p:metrics-file">
<div>
<h1><xsl:apply-templates select="@file"/></h1>
<xsl:apply-templates select="p:headers|p:metrics|p:histogram"/>
</div>
</xsl:template>

<xsl:template match="p:headers">
<h2>Headers</h2><dl>
<xsl:apply-templates select="p:header"/>
</dl>
</xsl:template>

<xsl:template match="p:header">
<dt><xsl:value-of select="@class"/></dt>
<dd><xsl:value-of select="."/></dd>
</xsl:template>



<xsl:template match="p:metrics">
<h2><xsl:value-of select="p:thead/@class"/></h2>
<xsl:choose>
<xsl:when test="count(p:tbody/p:tr)=1">
<table border="1">
<xsl:for-each select="p:thead/p:th">
<xsl:variable name="idx" select="position()"/>
<tr>
 <th><xsl:value-of select="."/></th>
 <td><xsl:value-of select="../../p:tbody/p:tr[1]/p:td[$idx]"/></td>
</tr>
</xsl:for-each>
</table>
</xsl:when>
<xsl:otherwise>
<table border="1">
<thead>
	<xsl:for-each select="p:thead/p:th">
		<th><xsl:value-of select="."/></th>
	</xsl:for-each>
</thead>

<tbody>
<xsl:for-each select="p:tbody/p:tr">
<tr>
<xsl:for-each select="p:td">
<td>
<xsl:value-of select="."/>
</td>
</xsl:for-each>
</tr>

</xsl:for-each>
</tbody>
</table>
</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="p:histogram">
<h2>Histogram</h2>
<div style="width: 900px; height: 500px;">
<xsl:attribute name="id"><xsl:value-of select="generate-id(.)"/></xsl:attribute>
</div>
</xsl:template>

<xsl:template match="p:histogram" mode="google">
 google.load("visualization", "1", {packages:["corechart"]});
      google.setOnLoadCallback(<xsl:value-of select="concat('drawChart',generate-id(.))"/>);
      function <xsl:value-of select="concat('drawChart',generate-id(.))"/>() {
        var <xsl:value-of select="concat('data',generate-id(.))"/> = google.visualization.arrayToDataTable([
          <xsl:text>[</xsl:text><xsl:for-each select="p:thead/p:th"><xsl:if test="position()&gt;1">,</xsl:if>'<xsl:value-of select="."/>'</xsl:for-each><xsl:text>]</xsl:text>
          <xsl:for-each select="p:tbody/p:tr">
          <xsl:text>,</xsl:text>
          <xsl:text>[</xsl:text><xsl:for-each select="p:td"><xsl:if test="position()&gt;1">,</xsl:if><xsl:value-of select="."/></xsl:for-each><xsl:text>]</xsl:text>
          </xsl:for-each>
        ]);

        var <xsl:value-of select="concat('opt',generate-id(.))"/> = {
          title: '<xsl:value-of select="../@file"/>'
        };

        var <xsl:value-of select="concat('chart',generate-id(.))"/> = new google.visualization.LineChart(document.getElementById('<xsl:value-of select="generate-id(.)"/>'));
        <xsl:value-of select="concat('chart',generate-id(.))"/>.draw(<xsl:value-of select="concat('data',generate-id(.))"/>, <xsl:value-of select="concat('opt',generate-id(.))"/>);
        }
</xsl:template>

</xsl:stylesheet>

