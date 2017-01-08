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

<xsl:template match="/">
<xsl:apply-templates select="command"/>
</xsl:template>

<xsl:template match="command">
<jnlp spec="6.0+" codebase="__WEBSTART__/ped" href="ped.jnlp" >
	<xsl:attribute name="codebase">
		<xsl:value-of select="$codebase"/>
	</xsl:attribute>
   <information>
      <title><xsl:value-of select="$name"/></title>
      <vendor>http://www.inserm.f</vendor>
      <description><xsl:value-of select="descritpion"/></description>
      <description kind="short">Pedigree Drawer</description>
   </information>
   <security>
      <all-permissions/>
   </security>
   <update check="always" policy="prompt-run" />
   <resources>
      <j2se version="1.6+" />
      <jar href="ped.jar" />
      <jar href="codec.jar" />
      <jar href="httpclient.jar" />
      <jar href="logging.jar" />
   </resources>
   <application-desc main-class="fr.inserm.umr915.bomcat.ped.PedigreeDrawer"/>
</jnlp>
</xsl:template>

</xsl:stylesheet>

