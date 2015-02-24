<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	version='1.0' 
	> 
<xsl:output method="xml"/>

<xsl:param name="class"/>
<xsl:param name="version">undefined</xsl:param>
<xsl:param name="name">undefined</xsl:param>

<xsl:template match="/">
	<xsl:choose>
		<xsl:when test="tools/tool[@id= $class]">
			<xsl:apply-templates select="tools/tool[@id= $class]"/>
		</xsl:when>
		<xsl:otherwise>
			<xsl:message terminate="yes"></xsl:message>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>


<xsl:template match="tool">
<tool hidden="false">
<xsl:attribute name="id"><xsl:value-of select="@id"/></xsl:attribute>
<xsl:attribute name="version"><xsl:value-of select="$version"/></xsl:attribute>
<xsl:attribute name="name"><xsl:value-of select="$name"/></xsl:attribute>
</tool>
</xsl:template>

</xsl:stylesheet>

