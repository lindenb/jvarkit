<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:fx="http://javafx.com/fxml"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns="http://javafx.com/fxml"
	>
<xsl:output method="xml" indent="yes" encoding="UTF-8"/>
<xsl:param name="name"></xsl:param>
<xsl:param name="codebase"></xsl:param>
<xsl:param name="mainclass"></xsl:param>

<xsl:template match="/">
<xsl:apply-templates select="command"/>
</xsl:template>

<xsl:template match="command">
<jnlp spec="6.0+">
	<xsl:attribute name="codebase">
		<xsl:value-of select="$codebase"/>
	</xsl:attribute>
	<xsl:attribute name="href">
		<xsl:value-of select="$name"/>
		<xsl:text>.jnlp</xsl:text>
	</xsl:attribute>
   <information>
      <title><xsl:value-of select="$name"/></title>
      <vendor>Pierre Lindenbaum PhD</vendor>
      <homepage href="http://www.umr1087.univ-nantes.fr"/>
      <description><xsl:value-of select="description"/></description>
      <description kind="short"><xsl:value-of select="description"/></description>
      <icon>
      	 <xsl:attribute name="href">
			<xsl:value-of select="$name"/>
			<xsl:text>.png</xsl:text>
		</xsl:attribute>
      </icon>
      <icon kind="splash">
      	 <xsl:attribute name="href">
			<xsl:value-of select="$name"/>
			<xsl:text>.png</xsl:text>
		</xsl:attribute>
      </icon>
      <offline-allowed/>
	<shortcut online="true">
	  <desktop/>
	  <menu>
	  	<xsl:attribute name="submenu">
			<xsl:value-of select="$name"/>
		</xsl:attribute>
	  </menu>
	</shortcut>
   </information>
   <security>
      <all-permissions/>
   </security>
   <update check="always" policy="prompt-run" />
   <resources>
      <j2se version="1.8+" >
      	<xsl:if test="contains($mainclass,'gatk')">
      		<xsl:attribute name="max-heap-size">3g</xsl:attribute>
      	</xsl:if>
      </j2se>
      <xsl:for-each select="libraries/library">
      	<jar download="null">
      	<xsl:attribute name="href">
			<xsl:value-of select="@href"/>
		</xsl:attribute>
		<xsl:if test="@main = 'true'">
			<xsl:attribute name="main">true</xsl:attribute>
		</xsl:if>
      	</jar>
      </xsl:for-each>
   </resources>
   <application-desc width="1000" height="800">
   		<xsl:attribute name="name">
			<xsl:value-of select="$name"/>
		</xsl:attribute>
   		<xsl:attribute name="main-class">
			<xsl:value-of select="$mainclass"/>
		</xsl:attribute>
   </application-desc>
   <update check="background"/>
</jnlp>
</xsl:template>

</xsl:stylesheet>

