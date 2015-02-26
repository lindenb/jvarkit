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

<xsl:comment>Date: <xsl:value-of select="date:date-time()"/></xsl:comment>	
	
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
<xsl:value-of select="$classpath"/>
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
<xsl:apply-templates />
</outputs>
</xsl:template>



<xsl:template match="param">
<xsl:copy-of select="."/>
</xsl:template>

<xsl:template match="main-class">
<xsl:value-of select="$class"/>
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
  <xsl:when test="@label"><xsl:apply-templates select="@label"/></xsl:when>
  <xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="param|data"  mode="description.rst">
<xsl:choose>
  <xsl:when test="description"><xsl:apply-templates select="description"/></xsl:when>
  <xsl:when test="@description"><xsl:apply-templates select="@description"/></xsl:when>
  <xsl:otherwise><xsl:apply-templates select="." mode="label.rst" /></xsl:otherwise>
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

</xsl:stylesheet>

