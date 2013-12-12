<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	version='1.0' > 
<xsl:output method="html"/>

<xsl:template match="/">
<html >
<xsl:apply-templates/>
</html>
</xsl:template>

<xsl:template match="vcf-statistics">

<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8"/>
<title>VCF statistics: <xsl:value-of select="@input"/></title>
 <script type="text/javascript" src="http://www.google.com/jsapi"></script>
 <script type="text/javascript">
  google.load('visualization', '1',{packages:["corechart", "table"]});
  
 function drawVisualization()
 	{ 
	<xsl:apply-templates select="//data-table" mode="g1"/>
	<xsl:apply-templates select="//chart" mode="g1"/>
	}
 google.setOnLoadCallback(drawVisualization);
 </script>


</head>
<body style="font-family: Arial;border: 0 none;">
<div>
Chromosomes : <xsl:for-each select="//statistics[@chromosome]"> [<a><xsl:attribute name="href"><xsl:value-of select="concat('#',generate-id(.))"/></xsl:attribute><xsl:value-of select="@chromosome"/></a>]</xsl:for-each>
</div>
<div>
Samples : <xsl:for-each select="//statistics[@sample]"> [<a><xsl:attribute name="href"><xsl:value-of select="concat('#',generate-id(.))"/></xsl:attribute><xsl:value-of select="@sample"/></a>]</xsl:for-each>
</div>
<xsl:apply-templates/>
</body>
</xsl:template>

<xsl:template match="statistics">
<div>
<a><xsl:attribute name="name"><xsl:value-of select="generate-id(.)"/></xsl:attribute></a>
<h3><xsl:value-of select="@name"/></h3>
<p><xsl:value-of select="@description"/></p>
<xsl:apply-templates/>
</div>
</xsl:template>

<xsl:template match="section">
<div>
<a><xsl:attribute name="name"><xsl:value-of select="generate-id(.)"/></xsl:attribute></a>
<h2><xsl:value-of select="@name"/></h2>
<p><xsl:value-of select="@description"/></p>
<xsl:apply-templates/>
</div>
</xsl:template>


<xsl:template match="data-table|chart">
<h4><xsl:value-of select="@name"/></h4>
<p><xsl:value-of select="@description"/></p>
<xsl:apply-templates select="." mode="g2"/>
</xsl:template>

<xsl:template match="counts">
<h3><xsl:value-of select="@name"/></h3>
<p><xsl:value-of select="@description"/></p>
<table>
	<thead>
		<th>Key</th>
		<th>Count</th>
	</thead>
	<tbody>
		<xsl:for-each select="property">
			<tr>
				<th><xsl:value-of select="@key"/></th>
				<td><xsl:value-of select="."/></td>
			</tr>
		</xsl:for-each>
	</tbody>
</table>
</xsl:template>




<xsl:template match="chart|data-table" mode="g1">
<xsl:variable name="hid" select="generate-id(.)"/>


  var wrapper<xsl:value-of select="$hid"/> = new google.visualization.ChartWrapper({
    chartType: 'ColumnChart',
    dataTable: [[''<xsl:for-each select="property">,'<xsl:choose>
    	<xsl:when test="@label">
    		<xsl:value-of select="@label"/>
    	</xsl:when>
    	<xsl:otherwise>
    		<xsl:value-of select="@key"/>
    	</xsl:otherwise>
    	</xsl:choose>'</xsl:for-each>],
                [''<xsl:for-each select="property">,<xsl:value-of select="."/></xsl:for-each>]],
    options: {'title': '<xsl:value-of select="@description"/>' },
    containerId: 'visualization<xsl:value-of select="$hid"/>'
  });
  wrapper<xsl:value-of select="$hid"/>.draw();
</xsl:template>


<xsl:template match="chart|data-table" mode="g2">
<xsl:variable name="hid" select="generate-id(.)"/>
<div   style="width: 800px; height: 500px;">
	<xsl:attribute name="id">visualization<xsl:value-of select="$hid"/></xsl:attribute>
	<xsl:text> </xsl:text>
</div>
</xsl:template>

</xsl:stylesheet>

