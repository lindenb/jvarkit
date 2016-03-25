<?xml version='1.0' ?>
<xsl:stylesheet
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns="http://www.w3.org/1999/xhtml"
	xmlns:galaxy="https://usegalaxy.org/"
	version='1.0'
	>
<xsl:output method="xml" indent="no" />
<xsl:param name="githash">undefined</xsl:param>
<xsl:param name="jarname">jarname</xsl:param>
<xsl:template match="/">
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="@*">
<xsl:copy select="."/>
</xsl:template>

<xsl:template match="c:app">
<c:app >
	<xsl:attribute name="githash"><xsl:value-of select="$githash"/></xsl:attribute>
	<xsl:attribute name="jarname"><xsl:value-of select="$jarname"/></xsl:attribute>
	<xsl:apply-templates select="@*"/>
	
	<xsl:if test="not(c:options)">
		<c:options>
			<xsl:apply-templates select="." mode="otheroptions"/>
		</c:options>
	</xsl:if>
	
	<xsl:if test="not(c:elixir)">
		<c:elixir>
			<xsl:apply-templates select="." mode="elixir"/>
		</c:elixir>
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

<xsl:template match="c:elixir">
<c:elixir>
	<xsl:apply-templates select="@*"/>
	<xsl:apply-templates select="*"/>
	<xsl:apply-templates select="/c:app" mode="elixir"/>
</c:elixir>
</xsl:template>

<xsl:template match="*">
<xsl:copy select=".">
<xsl:apply-templates select="@*"/>
<xsl:apply-templates select="*|text()"/>
</xsl:copy>
</xsl:template>

<xsl:template match="c:app" mode="elixir">
</xsl:template>

<xsl:template match="c:description">
<c:description>
<xsl:apply-templates select="@*"/>
<xsl:apply-templates select="*|text()"/>

<xsl:if test="../@type='input-file-set' or ../@type='string-list' or ../@type='string-set' or ../@type='uri-set'">
 <xsl:text>Multiple calls to this option should end with double hyphen : --</xsl:text>
</xsl:if>
</c:description>
</xsl:template>


<xsl:template match="c:app" mode="otheroptions">

	<xsl:if test="/c:app[not(@generate-output-option='false')]">
		<c:option name="outputFile" type="output-file" arg-name="OUTPUT-FILE" opt="o" longopt="output" label="Output">
			<c:description>Output file. Default:stdout</c:description>
		</c:option>
	</xsl:if>

	<xsl:if test="/c:app/c:output[@type='sam'] or /c:app/c:output[@type='bam'] or  /c:app/c:snippet[@id='write-sam']">
		<c:option name="formatout" type="string" arg-name="FORMAT" opt="formatout" longopt="formatout" default="sam" label="Sam Format">
			<c:regex>(sam|bam)</c:regex>
			<c:description>output format : sam or bam.  if stdout is used</c:description>
		</c:option>
		<c:option name="bam_compression_level" type="int" arg-name="LEVEL" opt="bam_compression_level" longopt="bam_compression_level" default="9">
			<c:description>BAM Compression level (0-9)</c:description>
		</c:option>
	</xsl:if>

	<xsl:if test="/c:app/c:snippet[@id='sorting-collection']">
		<c:option name="maxRecordsInRam" type="int" arg-name="NUMBER" label="Max Records in RAM" opt="maxRecordsInRam" longopt="maxRecordsInRam" default="500000" galaxy:ignore="true">
			<c:description>When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM/VCF/... file, and increases the amount of RAM needed.</c:description>
		</c:option>
	</xsl:if>
	
	<xsl:if test="/c:app/c:snippet[@id='sorting-collection'] or /c:app/c:snippet[@id='tmp-dir']" >
		<c:option name="tmpdir" type="input-directory" arg-name="TMPDIR" label="Temporary directory" opt="tmpdir" longopt="tmpdir" galaxy:ignore="true" >
			<c:description>Set temporary directory</c:description>
		</c:option>
	</xsl:if>
		<xsl:if test="/c:app/c:snippet[@id='http.proxy']">
		<c:option name="http_proxy_str" type="string" argname="HOST:PORT" opt="http_proxy" longopt="http_proxy" galaxy:ignore="true">
			<c:description>set the http and the https proxy ( HOST:PORT ) </c:description>
		</c:option>
	</xsl:if>

	<xsl:if test="/c:app/c:snippet[@id='javascript']">
	<c:option  name="javascriptFile" type="input-file" argname="SCRIPT.JS" opt="f" longopt="jsfile" galaxy:optional="true">
		<c:description>Javascript file. Use either javascript file or javascript expression.</c:description>
	</c:option>
	<c:option name="javascriptExpr"  galaxy:optional="true" galaxy:area="true" galaxy:size="20x80" type="string" argname="SCRIPT" opt="e" longopt="jsexpr"  multiline="true">
		<c:description>Javascript expression. Use either javascript file or javascript expression.</c:description>
	</c:option>
	</xsl:if>
	
	<xsl:if test="/c:app/c:snippet[@id='berkeleydb']">
	<c:option name="berkeleyDbHome" type="input-directory" argname="BDB.HOME" opt="bdbHome" longopt="bdbHome" galaxy:ignore="true">
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
		
	    <galaxy:conditional>
	     <xsl:attribute name="name"><xsl:value-of select="$snippet/@name"/>_source</xsl:attribute>
	      <galaxy:param name="reference_source_selector" type="select" label="Load reference genome from">
	        <galaxy:option value="cached">Local cache</galaxy:option>
	        <galaxy:option value="history">History</galaxy:option>
	      </galaxy:param>
	      <galaxy:when value="cached">
	        <galaxy:param name="ref_file" type="select" label="Using reference genome" help="REFERENCE_SEQUENCE">
	          <galaxy:options from_data_table="all_fasta">
	          </galaxy:options>
	          <galaxy:validator type="no_options" message="A built-in reference genome is not available for the build associated with the selected input file"/>
	        </galaxy:param>
	      </galaxy:when>
	      <galaxy:when value="history">
	        <galaxy:param name="ref_file" type="data" format="fasta" label="Use the folloing dataset as the reference sequence" help="REFERENCE_SEQUENCE; You can upload a FASTA sequence to the history and use it as reference" />
	      </galaxy:when>
	    </galaxy:conditional>
	    <galaxy:command>
 #if $<xsl:value-of select="$snippet/@name"/>_source.reference_source_selector != "history":
     -<xsl:value-of select="$snippet/@opt"/> "${<xsl:value-of select="$snippet/@name"/>_source.ref_file.fields.path}"
 #else:
     -<xsl:value-of select="$snippet/@opt"/> "${<xsl:value-of select="$snippet/@name"/>_source.ref_file}"
 #end if
</galaxy:command>
	    
	    
	</c:option>
	</xsl:if>

	<xsl:if test="/c:app/c:snippet[@id='jetty-server']">
		<c:option name="serverPort" type="int" argname="PORT" opt="port" longopt="port" default="8080">
			<c:description>server port</c:description>
		</c:option>
	</xsl:if>
	
	
</xsl:template>

<xsl:template match="c:doc[@id='pedigree-file']">
<xsl:text>A pedigree file is a text file containing the following columns: FAMILY-ID INDIVIDUAL-ID FATHER-ID(or 0) MOTHER-ID(or 0) SEX(1 male, 2 female, 0 unknown) PHENOTYPE (1 affected, 0 unaffected, 9 not available)</xsl:text>
</xsl:template>

<xsl:template match="c:doc[@id='file.list']">
<xsl:text>If the filename ends with .list , it will be interpreted as a text file containing the path(s) of the file(s) to be read.</xsl:text>
</xsl:template>




</xsl:stylesheet>

