<?xml version='1.0' ?>
<xsl:stylesheet
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns="http://www.w3.org/1999/xhtml"
	version='1.0'
	>
<xsl:output method="xml" indent="no" />


<xsl:template match="/">
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="@*">
<xsl:copy select="."/>
</xsl:template>

<xsl:template match="c:options">
<c:options>
	<xsl:if test="../c:snippet[@id='javascript']">
	<c:option name="javascriptFile" type="input-file" argname="SCRIPT.JS" opt="f" longopt="jsfile">
		<c:description>Javascript file</c:description>
	</c:option>
	<c:option name="javascriptExpr" type="string" argname="SCRIPT" opt="e" longopt="jsexpr"
		multiline="true">
		<c:description>Javascript expression</c:description>
	</c:option>
	</xsl:if>
	<xsl:apply-templates select="@*"/>
	<xsl:apply-templates select="*|text()"/>
</c:options>
</xsl:template>

<xsl:template match="*">
<xsl:copy select=".">
<xsl:apply-templates select="@*"/>
<xsl:apply-templates select="*|text()"/>
</xsl:copy>
</xsl:template>

</xsl:stylesheet>

