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

<xsl:template match="c:app">
<c:app >
	<xsl:apply-templates select="@*"/>
	
	<xsl:if test="not(c:options)">
		<c:options>
			<xsl:call-template name="make-output-option"/>
			<xsl:call-template name="make-sam-output-option"/>
		</c:options>
	</xsl:if>
	
	<xsl:apply-templates select="*|text()"/>
</c:app>
</xsl:template>


<xsl:template match="c:input[not(@type)]">
<c:input type="stdin-or-many">
	<xsl:apply-templates select="@*"/>
	<xsl:apply-templates select="*|text()"/>
</c:input>
</xsl:template>


<xsl:template match="c:options">
<c:options>
	<xsl:apply-templates select="@*"/>
	
	<xsl:call-template name="make-output-option"/>
	<xsl:call-template name="make-sam-output-option"/>
	
	<xsl:if test="../c:snippet[@id='http.proxy']">
		<c:option name="http_proxy_str" type="string" argname="HOST:PORT" opt="http_proxy" longopt="http_proxy">
			<c:description>set the http and the https proxy ( HOST:PORT ) </c:description>
		</c:option>
	</xsl:if>

	<xsl:if test="../c:snippet[@id='javascript']">
	<c:option name="javascriptFile" type="input-file" argname="SCRIPT.JS" opt="f" longopt="jsfile">
		<c:description>Javascript file</c:description>
	</c:option>
	<c:option name="javascriptExpr" type="string" argname="SCRIPT" opt="e" longopt="jsexpr"  multiline="true">
		<c:description>Javascript expression</c:description>
	</c:option>
	</xsl:if>
	
	
	<xsl:if test="../c:snippet[@id='ref.faidx']">
	<xsl:variable name="snippet" select="../c:snippet[@id='ref.faidx']"/>
	<c:option type="input-file" argname="FASTA">
		<xsl:attribute name="name">
			<xsl:choose>
				<xsl:when test="$snippet/@name"><xsl:value-of select="$snippet/@name"/></xsl:when>
				<xsl:otherwise><xsl:message terminate="yes">No Name for snpippet faidx</xsl:message></xsl:otherwise>
			</xsl:choose>
		</xsl:attribute>
		<xsl:attribute name="opt">
			<xsl:choose>
				<xsl:when test="$snippet/@opt"><xsl:value-of select="$snippet/@opt"/></xsl:when>
				<xsl:otherwise>R</xsl:otherwise>
			</xsl:choose>
		</xsl:attribute>
		<xsl:attribute name="longopt">
			<xsl:choose>
				<xsl:when test="$snippet/@longopt"><xsl:value-of select="$snippet/@longopt"/></xsl:when>
				<xsl:otherwise>REF</xsl:otherwise>
			</xsl:choose>
		</xsl:attribute>
		<c:description>indexed Fasta sequence</c:description>
	</c:option>
	</xsl:if>
	
	
	<xsl:apply-templates select="*|text()"/>
</c:options>
</xsl:template>

<xsl:template match="*">
<xsl:copy select=".">
<xsl:apply-templates select="@*"/>
<xsl:apply-templates select="*|text()"/>
</xsl:copy>
</xsl:template>


<xsl:template name="make-output-option">
	<xsl:if test="/c:app[not(@generate-output-option='false')]">
		<c:option name="outputFile" type="output-file" arg-name="OUTPUT-FILE" opt="o" longopt="output">
			<c:description>Output file. Default:stdout</c:description>
		</c:option>
	</xsl:if>
</xsl:template>


<xsl:template name="make-sam-output-option">
	<xsl:if test="/c:app/c:output[@type='sam'] or /c:app/c:output[@type='bam']">
		<c:option name="formatout" type="string" arg-name="FORMAT" opt="formatout" longopt="formatout" default="sam">
			<c:regex>(sam|bam)</c:regex>
			<c:description>output format : sam or bam.  if stdout is used</c:description>
		</c:option>
	</xsl:if>
</xsl:template>


</xsl:stylesheet>

