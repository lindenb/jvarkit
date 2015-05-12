<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:date="http://exslt.org/dates-and-times" 
	version='1.0' 
	>
<!-- Motivation: generate tools.xml for galaxy -->
<xsl:output method="xml" indent="yes"/>

<xsl:param name="class"/>
<xsl:param name="version">undefined</xsl:param>
<xsl:param name="name">undefined</xsl:param>
<xsl:param name="classpath">undefined</xsl:param>

<xsl:template match="/">
	<xsl:choose>
		<xsl:when test="tools/tool[@main-class= $class]">
			<xsl:apply-templates select="tools/tool[@main-class= $class]"/>
		</xsl:when>
		<xsl:otherwise>
			<xsl:message terminate="yes"></xsl:message>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>


<xsl:template match="tool">
<tool hidden="false">
<xsl:attribute name="id"><xsl:value-of select="@main-class"/></xsl:attribute>
<xsl:attribute name="version"><xsl:value-of select="$version"/></xsl:attribute>
<xsl:attribute name="name"><xsl:value-of select="$name"/></xsl:attribute>

<xsl:comment>Date: <xsl:value-of select="date:date-time()"/></xsl:comment>	
	<description><xsl:apply-templates select="." mode="description"/></description>
	<requirements>
		<requirement type="binary">java</requirement>
	</requirements>
	 <xsl:apply-templates select="command"/>
	 <xsl:apply-templates select="inputs"/>
	 <xsl:apply-templates select="outputs"/>
	 <stdio>
		<exit_code range="1:" />
		<exit_code range=":-1" />
	 </stdio>
	 <help>
<xsl:apply-templates select="description" mode="rst"/>
<xsl:apply-templates select="inputs" mode="rst"/>
<xsl:apply-templates select="outputs" mode="rst"/>

<xsl:if test="wiki">
<xsl:text>

**Wiki**

</xsl:text>
<xsl:apply-templates select="wiki"/>
</xsl:if>

<xsl:text>

-----


**Author**

Pierre Lindenbaum PhD @yokofakun

**Citation**

If you use this Galaxy tool in work leading to a scientific publication please
cite: Pierre Lindenbaum PhD  https://github.com/lindenb/jvarkit

**Contribute**

- Issue Tracker: http://github.com/lindenb/jvarkit/issues`
- Source Code: http://github.com/lindenb/jvarkit

**License**

The project is licensed under the MIT license.

**Compilation**

Date: </xsl:text>
<xsl:value-of select="date:date-time()"/>
<xsl:text>
Version: </xsl:text>
<xsl:value-of select="$version"/>
<xsl:text>
</xsl:text>
	 
	 </help>
</tool>
</xsl:template>

<xsl:template match="classpath">
<xsl:call-template name="print.classpath">
  <xsl:with-param name="cp" select="normalize-space($classpath)"/>
  <xsl:with-param name="index" select="number(0)"/>
</xsl:call-template>
</xsl:template>

<xsl:template match="main-class">
<xsl:value-of select="../@main-class"/>
</xsl:template>

<xsl:template match="command">
<command>
<xsl:apply-templates/>
</command>
</xsl:template>



<xsl:template match="inputs">
<inputs>
 <xsl:apply-templates/>
</inputs>
</xsl:template>

<xsl:template match="outputs">
<outputs>
<xsl:apply-templates select="data"/>
</outputs>
</xsl:template>



<xsl:template match="param|repeat|data">
<xsl:copy-of select="."/>
</xsl:template>




<!-- restructuredtext  http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html -->

<xsl:template match="description" mode="rst">
 <xsl:apply-templates mode="rst"/>
</xsl:template>

<xsl:template match="inputs"  mode="rst">
<xsl:text>

**Inputs**

</xsl:text>
<xsl:apply-templates select="param" mode="rst"/>
</xsl:template>


<xsl:template match="outputs"  mode="rst">
<xsl:text>

**Outputs**

</xsl:text>
<xsl:apply-templates select="data" mode="rst"/>
</xsl:template>


<xsl:template match="param"  mode="rst">
<xsl:text>
- **</xsl:text>
<xsl:apply-templates select="." mode="label.rst" />
<xsl:text>** :  </xsl:text>
<xsl:apply-templates select="." mode="description.rst"/>
<xsl:text>
</xsl:text>
</xsl:template>

<xsl:template match="data"  mode="rst">
<xsl:text>
- **</xsl:text>
<xsl:apply-templates select="." mode="label.rst" />
<xsl:text>** :  </xsl:text>
<xsl:apply-templates select="." mode="description.rst"/>
<xsl:text>
</xsl:text>
</xsl:template>


<xsl:template match="param|data"  mode="label.rst">
<xsl:choose>
  <xsl:when test="label"><xsl:apply-templates select="label"/></xsl:when>
  <xsl:when test="@label"><xsl:value-of select="@label"/></xsl:when>
  <xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="param|data"  mode="description.rst">
<xsl:choose>
  <xsl:when test="description"><xsl:apply-templates select="description"/></xsl:when>
  <xsl:when test="@description"><xsl:value-of select="@description"/></xsl:when>
  <xsl:otherwise><xsl:apply-templates select="." mode="label.rst" /></xsl:otherwise>
</xsl:choose>
</xsl:template>


<xsl:template match="tool"  mode="description">
<xsl:choose>
	 <xsl:when test="label"><xsl:apply-templates select="label"/></xsl:when>
  <xsl:when test="@label"><value-of select="@label"/></xsl:when>
  <xsl:otherwise><value-of select="@name"/></xsl:otherwise>
</xsl:choose>
</xsl:template>


<xsl:template match="a"  mode="rst">
<xsl:apply-templates mode="rst"/>
</xsl:template>

<xsl:template match="i"  mode="rst">
<xsl:text>*</xsl:text>
<xsl:apply-templates mode="rst"/>
<xsl:text>*</xsl:text>
</xsl:template>

<xsl:template match="b"  mode="rst">
<xsl:text>**</xsl:text>
<xsl:apply-templates mode="rst"/>
<xsl:text>**</xsl:text>
</xsl:template>

<xsl:template match="div"  mode="rst">
<xsl:text>
</xsl:text>
<xsl:apply-templates mode="rst"/>
<xsl:text>
</xsl:text>
</xsl:template>


<xsl:template match="h3"  mode="rst">
<xsl:text>

**</xsl:text>
<xsl:apply-templates mode="rst"/>
<xsl:text>**
</xsl:text>
</xsl:template>

<xsl:template name="print.classpath">
<xsl:param name="cp"/>
<xsl:param name="index"/>
<xsl:choose>
   <xsl:when test="string-length(normalize-space($cp))=0">
  </xsl:when>
  <xsl:when test="contains(normalize-space($cp),' ')">
  		<xsl:if test="$index &gt;0"><xsl:text>:</xsl:text></xsl:if>
  		<xsl:text>$__tool_directory__/</xsl:text>
  		<xsl:value-of select="substring-before($cp,' ')"/>
		<xsl:call-template name="print.classpath">
		  <xsl:with-param name="cp" select="substring-after($cp,' ')"/>
		  <xsl:with-param name="index" select="number(1)"/>
		</xsl:call-template>
  </xsl:when>
  <xsl:otherwise>
  		<xsl:if test="$index &gt;0"><xsl:text>:</xsl:text></xsl:if>
  		<xsl:text>$__tool_directory__/</xsl:text>
  		<xsl:value-of select="$cp"/>
  </xsl:otherwise>
</xsl:choose>
</xsl:template>


</xsl:stylesheet>

