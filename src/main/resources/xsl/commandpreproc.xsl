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
			<xsl:apply-templates select="." mode="otheroptions"/>
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
	<xsl:comment>Generate options</xsl:comment>
	<xsl:apply-templates select="/c:app" mode="otheroptions"/>
	<xsl:comment>user options</xsl:comment>
	
	
	<xsl:apply-templates select="*"/>
</c:options>
</xsl:template>

<xsl:template match="*">
<xsl:copy select=".">
<xsl:apply-templates select="@*"/>
<xsl:apply-templates select="*|text()"/>
</xsl:copy>
</xsl:template>


<xsl:template match="c:app" mode="otheroptions">

	<xsl:if test="/c:app[not(@generate-output-option='false')]">
		<c:option name="outputFile" type="output-file" arg-name="OUTPUT-FILE" opt="o" longopt="output" label="Output">
			<c:description>Output file. Default:stdout</c:description>
		</c:option>
	</xsl:if>

	<xsl:if test="/c:app/c:output[@type='sam'] or /c:app/c:output[@type='bam']">
		<c:option name="formatout" type="string" arg-name="FORMAT" opt="formatout" longopt="formatout" default="sam" label="Sam Format">
			<c:regex>(sam|bam)</c:regex>
			<c:description>output format : sam or bam.  if stdout is used</c:description>
		</c:option>
	</xsl:if>

	<xsl:if test="/c:app/c:snippet[@id='sorting-collection']">
		<c:option name="maxRecordsInRam" type="int" arg-name="NUMBER" label="Max Records in RAM" opt="maxRecordsInRam" longopt="maxRecordsInRam" default="500000">
			<c:description>When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.</c:description>
		</c:option>
	</xsl:if>
	
	<xsl:if test="/c:app/c:snippet[@id='sorting-collection'] or /c:app/c:snippet[@id='tmp-dir']" >
		<c:option name="tmpdir" type="input-directory" arg-name="TMPDIR" label="Temporary directory" opt="tmpdir" longopt="tmpdir" >
			<c:description>Set tmp directory</c:description>
		</c:option>
	</xsl:if>
		<xsl:if test="/c:app/c:snippet[@id='http.proxy']">
		<c:option name="http_proxy_str" type="string" argname="HOST:PORT" opt="http_proxy" longopt="http_proxy">
			<c:description>set the http and the https proxy ( HOST:PORT ) </c:description>
		</c:option>
	</xsl:if>

	<xsl:if test="/c:app/c:snippet[@id='javascript']">
	<c:option name="javascriptFile" type="input-file" argname="SCRIPT.JS" opt="f" longopt="jsfile">
		<c:description>Javascript file</c:description>
	</c:option>
	<c:option name="javascriptExpr" type="string" argname="SCRIPT" opt="e" longopt="jsexpr"  multiline="true">
		<c:description>Javascript expression</c:description>
	</c:option>
	</xsl:if>
	
	<xsl:if test="/c:app/c:snippet[@id='berkeleydb']">
	<c:option name="berkeleyDbHome" type="input-directory" argname="BDB.HOME" opt="bdbHome" longopt="bdbHome" >
		<c:description>BerkeleyDB home directory used to store data</c:description>
	</c:option>
	</xsl:if>
	
	
	<xsl:if test="/c:app/c:snippet[@id='ref.faidx']">
	<xsl:variable name="snippet" select="/c:app/c:snippet[@id='ref.faidx']"/>
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

	
</xsl:template>


</xsl:stylesheet>

