<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:date="http://exslt.org/dates-and-times" 
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns:h="http://www.w3.org/1999/xhtml"
	xmlns:galaxy="https://usegalaxy.org/"
	version='1.0'
	exclude-result-prefixes="date c xsl h galaxy"
	>
<xsl:import href="mod.command.xsl"/>
<!-- Motivation: generate tools.xml for galaxy -->
<xsl:output method="xml" indent="yes"/>


<xsl:template match="/">
 <xsl:apply-templates select="c:app"/>
</xsl:template>


<xsl:template match="c:app">
<tool hidden="false">
<xsl:attribute name="id"><xsl:value-of select="@app"/></xsl:attribute>
<xsl:attribute name="version"><xsl:value-of select="@githash"/></xsl:attribute>
<xsl:attribute name="name"><xsl:value-of select="@app"/></xsl:attribute>

	<description><xsl:apply-templates select="c:description"/></description>
	<requirements>
		 <requirement type="binary">java</requirement>
		 <requirement type="set_environment">JVARKIT_DIST</requirement> 
	</requirements>
	<command>
		<xsl:choose>
			<xsl:when test="galaxy:galaxy/galaxy:command">
				<xsl:value-of select="galaxy:galaxy/galaxy:command"/>
			</xsl:when>
			<xsl:otherwise>
			 	<xsl:apply-templates select="." mode="command"/>
			</xsl:otherwise>
		</xsl:choose>
	</command>
	<inputs>
		 <xsl:choose>
			<xsl:when test="c:input/@type='vcf'">
				<param format="vcf,vcf_bgzip" name="input" type="data">
					 <xsl:attribute name="label">VCF input" for '<xsl:value-of select="@app"/>'</xsl:attribute>
				</param>
			</xsl:when>
			<xsl:otherwise>
				<xsl:message terminate="yes">unknown input for '<xsl:value-of select="@app"/>'</xsl:message>
			</xsl:otherwise>
		</xsl:choose>
		<xsl:apply-templates select="c:options/c:option" mode="input"/>
		<xsl:apply-templates select="galaxy:galaxy/galaxy:inputs/galaxy:param"/>
		<xsl:apply-templates select="galaxy:galaxy/galaxy:inputs/galaxy:conditional"/>
		
	</inputs>
	<outputs>
		<xsl:apply-templates select="c:options/c:option" mode="output"/>
		<xsl:apply-templates select="galaxy:galaxy/galaxy:outputs/galaxy:data"/>
	</outputs>
	<stdio>
		<exit_code range="1:" />
		<exit_code range=":-1" />
	</stdio>
	<help>

<xsl:apply-templates select="c:documentation" mode="rst"/>

**Author**

Pierre Lindenbaum PhD. @yokofakun. Institut du Thorax, U1087, 44000 Nantes.

**Source Code**

<xsl:value-of select="concat('https://github.com/lindenb/jvarkit/blob/master/src/main/java/',translate(@package,'.','/'),'/',@app,'.java')"/>


**Report bugs / Contribute**

https://github.com/lindenb/jvarkit/issues

**License**

The project is licensed under the MIT license.



Date: <xsl:value-of select="date:date-time()"/>
  </help>
  <citations>
    <citation type="doi">10.6084/m9.figshare.1425030.v1</citation>
  </citations>
</tool>
</xsl:template>

<xsl:template match="c:option" mode="input">
<xsl:choose>
	<xsl:when test="@galaxy:ignore='true'"></xsl:when>
	<xsl:when test="@type='string-list'">
		<repeat>
			<xsl:attribute name="name"><xsl:value-of select="@name"/>_list</xsl:attribute>
			<xsl:attribute name="title"><xsl:value-of select="@name"/></xsl:attribute>
			<xsl:attribute name="help"><xsl:value-of select="c:description"/></xsl:attribute>
			<param type="text">
				<xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
				<xsl:attribute name="label"><xsl:value-of select="c:description"/></xsl:attribute>
				<xsl:attribute name="value"></xsl:attribute>
			</param>
		</repeat>
	</xsl:when>
	<xsl:when test="@type='string'">
		<param type="text">
			<xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
			<xsl:attribute name="label"><xsl:value-of select="c:description"/></xsl:attribute>
			<xsl:attribute name="value"><xsl:value-of select="@default"/></xsl:attribute>
		</param>
	</xsl:when>
	<xsl:when test="@type='boolean'">
		<param type="boolean">
			<xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
			<xsl:attribute name="label"><xsl:value-of select="c:description"/></xsl:attribute>
			<xsl:attribute name="checked"><xsl:value-of select="@default"/></xsl:attribute>
			<xsl:attribute name="truevalue"><xsl:value-of select="concat('-',@opt)"/></xsl:attribute>
			<xsl:attribute name="falsevalue"></xsl:attribute>
		</param>
	</xsl:when>
	<xsl:when test="@type='int'">
		<param type="integer">
			<xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
			<xsl:attribute name="label"><xsl:value-of select="c:description"/></xsl:attribute>
			<xsl:attribute name="value"><xsl:value-of select="@default"/></xsl:attribute>
		</param>
	</xsl:when>
	<xsl:when test="@type='double'">
		<param type="float">
			<xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
			<xsl:attribute name="label"><xsl:value-of select="c:description"/></xsl:attribute>
			<xsl:attribute name="value"><xsl:value-of select="@default"/></xsl:attribute>
		</param>
	</xsl:when>
	<xsl:when test="@type='input-file' and galaxy:conditional">
		<xsl:apply-templates select="galaxy:conditional"/>
	</xsl:when>
	<xsl:when test="@type='input-file'">
		<param type="data">
			<xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
			<xsl:attribute name="label"><xsl:value-of select="c:description"/></xsl:attribute>
			<xsl:attribute name="format">
				<xsl:choose>
					<xsl:when test="@galaxy:format">
						<xsl:value-of select="@galaxy:format"/>
					</xsl:when>
					<xsl:otherwise></xsl:otherwise>
				</xsl:choose>
			</xsl:attribute>
		</param>
	</xsl:when>
	<xsl:when test="@opt='o' and /c:app/c:output/@type='vcf'"></xsl:when>
	<xsl:otherwise>
		<xsl:message terminate="yes">c:option/mode=input : unknown input for '<xsl:value-of select="@name"/>/<xsl:value-of select="@type"/>'</xsl:message>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>


<xsl:template match="c:option" mode="output">
<xsl:choose>
	<xsl:when test="@galaxy:ignore='true'"></xsl:when>
	<xsl:when test="@type='string-list'"></xsl:when>
	<xsl:when test="@type='int'"></xsl:when>
	<xsl:when test="@type='input-file'"></xsl:when>
	<xsl:when test="@type='double'"></xsl:when>
	<xsl:when test="@type='string'"></xsl:when>
	<xsl:when test="@type='boolean'"></xsl:when>
	<xsl:when test="@opt='o' and /c:app/c:output/@type='vcf'">
		<data format="vcf_bgzip">
			<xsl:attribute name="name"><xsl:value-of select="@name"/></xsl:attribute>
		</data>
	</xsl:when>
	<xsl:otherwise>
		<xsl:message terminate="yes">c:option/mode=output : unknown output for '<xsl:value-of select="@name"/>/<xsl:value-of select="@type"/>'</xsl:message>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:app" mode="command">
<xsl:choose>
	<xsl:when test="c:input[@type='vcf']">
		<xsl:text>(gunzip -c ${input} 2&gt; /dev/null || cat ${input})</xsl:text>
	</xsl:when>
	<xsl:otherwise>
		<xsl:message terminate="yes">unknown c:input/type</xsl:message>
	</xsl:otherwise>
</xsl:choose>

<!-- call java -->
<xsl:text> |  java -jar \$JVARKIT_DIST/</xsl:text>
<xsl:value-of select="@jarname"/>
<xsl:text>.jar </xsl:text>
<xsl:for-each select="c:options/c:option">
	<xsl:text> </xsl:text>
	<xsl:apply-templates select="." mode="command"/>
	<xsl:text> </xsl:text>
</xsl:for-each>

<!-- move output -->
<xsl:choose>
 <xsl:when test='c:output/@type="vcf"'>
 	<xsl:text> &amp;&amp; cp '${outputFile}.vcf.gz' '${outputFile}' &amp;&amp; rm '${outputFile}.vcf.gz' </xsl:text>
  </xsl:when>
  <xsl:otherwise>
	<xsl:message terminate="yes">unknown c:output/type '<xsl:value-of select="c:output/@type"/>'</xsl:message>
  </xsl:otherwise>
</xsl:choose>


</xsl:template>


<xsl:template match="c:option" mode="command">
<xsl:choose>
	<xsl:when test="@galaxy:ignore='true'"></xsl:when>
	<xsl:when test='@type="boolean"'>
		<xsl:text> ${</xsl:text>
		<xsl:value-of select="@name"/>
		<xsl:text>} </xsl:text>
	</xsl:when>
	<xsl:when test='@type="string-list"'> -<xsl:value-of select="@opt"/>
#for $<xsl:value-of select="generate-id(.)"/>, in enumerate($<xsl:value-of select="@name"/>_list):
            '${<xsl:value-of select="generate-id(.)"/>.<xsl:value-of select="@name"/>}'
#end for
 --  </xsl:when>
	<xsl:when test='@type="output-file" and @opt="o" and /c:app/c:output/@type="vcf"'>
		<xsl:text>-o '${outputFile}.vcf.gz'</xsl:text>
	</xsl:when>
	<xsl:when test='@type="input-file" and galaxy:conditional[@id="vcf"]'>
<xsl:text>
#if $</xsl:text><xsl:value-of select="@name"/><xsl:text>.index_type == "fromHistory":
  ${</xsl:text><xsl:value-of select="@name"/>.<xsl:value-of select="@name"/><xsl:text>_vcf_file}
#else
  ${</xsl:text><xsl:value-of select="@name"/>.<xsl:value-of select="@name"/><xsl:text>_index.fields.path}
#end if
</xsl:text>
	</xsl:when>
	
	<xsl:when test='@type="int" or @type="input-file" or @type="double" or @type="string"'>
		<xsl:text>-</xsl:text>
		<xsl:value-of select="@opt"/>
		<xsl:text> '${</xsl:text>
		<xsl:value-of select="@name"/>
		<xsl:text>}'</xsl:text>
	</xsl:when>
	<xsl:otherwise>
		<xsl:message terminate="yes">unknown option type '<xsl:value-of select="@type"/>'</xsl:message>
  </xsl:otherwise>
</xsl:choose>
</xsl:template>


<!-- copy galaxy -->
<xsl:template match="galaxy:*">
<xsl:variable name="N" select="local-name(.)"/>
<xsl:element name="{$N}">
	<xsl:copy-of select="@*"/>
	<xsl:apply-templates select="galaxy:*|text()"/>
</xsl:element>
</xsl:template>

<xsl:template match="galaxy:conditional[@id='vcf']">
<xsl:variable name="varname">
<xsl:choose>
  	<xsl:when test="@name"><xsl:value-of select="@name"/></xsl:when>
  	<xsl:when test="local-name(..) ='option'"><xsl:value-of select="../@name"/></xsl:when>
	<xsl:otherwise test="not(@name)">
		<xsl:message terminate="yes">no @name in galaxy:conditional.</xsl:message>
	</xsl:otherwise>
</xsl:choose>
</xsl:variable>
<xsl:if test="not(@data-table)">
	<xsl:message terminate="yes">no @data-table in galaxy:conditional</xsl:message>
</xsl:if>
<conditional>
  <xsl:attribute name="name"><xsl:value-of select="$varname"/></xsl:attribute>
  <param name="index_type" type="select" label="Will you select a VCF from your history or use a built-in index?">
    <option value="builtin">Use a built-in VCF</option>
    <option value="history">Use one from your history</option>
  </param>
  <when value="builtin">
    <param type="select" label="Select a VCF">
      <xsl:attribute name="name"><xsl:value-of select="$varname"/>_index</xsl:attribute>
      <options>
      	<xsl:attribute name="from_data_table"><xsl:value-of select="@data-table"/></xsl:attribute>
        <validator type="no_options" message="No file are available" />
      </options>
    </param>
  </when>
  <when value="history">
    <param type="data" format="vcf,vcf_bgzip" metadata_name="dbkey" label="Select a VCF from your history" >
   	 <xsl:attribute name="name"><xsl:value-of select="$varname"/>_vcf_file</xsl:attribute>
    </param>
  </when>
</conditional>
</xsl:template>


<!-- restructuredtext  http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html -->

<xsl:template match="description" mode="rst">
 <xsl:apply-templates mode="rst"/>
</xsl:template>

<xsl:template match="c:inputs"  mode="rst">
<xsl:text>

**Inputs**

</xsl:text>
<xsl:apply-templates select="c:param" mode="rst"/>
</xsl:template>


<xsl:template match="c:outputs"  mode="rst">
<xsl:text>

**Outputs**

</xsl:text>
<xsl:apply-templates select="c:data" mode="rst"/>
</xsl:template>





<xsl:template match="h:a"  mode="rst">
<xsl:apply-templates mode="rst"/>
</xsl:template>

<xsl:template match="h:i"  mode="rst">
<xsl:text>*</xsl:text>
<xsl:apply-templates mode="rst"/>
<xsl:text>*</xsl:text>
</xsl:template>

<xsl:template match="h:b"  mode="rst">
<xsl:text>**</xsl:text>
<xsl:apply-templates mode="rst"/>
<xsl:text>**</xsl:text>
</xsl:template>

<xsl:template match="h:div"  mode="rst">
<xsl:text>
</xsl:text>
<xsl:apply-templates mode="rst"/>
<xsl:text>
</xsl:text>
</xsl:template>


<xsl:template match="h:h3"  mode="rst">
<xsl:text>

**</xsl:text>
<xsl:apply-templates mode="rst"/>
<xsl:text>**
</xsl:text>
</xsl:template>

<xsl:template match="h:li"  mode="rst">
<xsl:text>
- </xsl:text>
<xsl:apply-templates mode="rst"/>
</xsl:template>

<xsl:template match="h:ul"  mode="rst">
<xsl:text>
</xsl:text>
<xsl:apply-templates select="*" mode="rst"/>
<xsl:text>
</xsl:text>
</xsl:template>


<xsl:template match="h:pre"  mode="rst">
<xsl:text>

.. code:: bash
    
    </xsl:text>
<xsl:variable name="newline" select="'&#10;'" />
<xsl:call-template name="string-replace">
	<xsl:with-param name="string">
		<xsl:value-of select="."/>
	</xsl:with-param>
	<xsl:with-param name="from" select="$newline"/>
	<xsl:with-param name="to" select="concat($newline,'    ')"/>
</xsl:call-template>
<xsl:text>

</xsl:text>
</xsl:template>





</xsl:stylesheet>

