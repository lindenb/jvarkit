<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:fx="http://javafx.com/fxml"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	>
<xsl:output method="html" indent="no" encoding="UTF-8" omit-xml-declaration="yes"/>
<xsl:param name="name"></xsl:param>

<xsl:template match="/">
<xsl:apply-templates select="command"/>
</xsl:template>

<xsl:template match="command">
<tr><th>
<a>
 <xsl:attribute name="href">
	<xsl:value-of select="$name"/>
	<xsl:text>.jnlp</xsl:text>
 </xsl:attribute>
  <xsl:value-of select="$name"/>
 </a>

 </th>
 <td> 
 <xsl:apply-templates select="description"/>
 </td>
</tr>
</xsl:template>

</xsl:stylesheet>

