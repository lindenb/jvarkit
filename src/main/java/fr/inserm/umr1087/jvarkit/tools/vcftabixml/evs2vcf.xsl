<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform' version='1.0' > 
<xsl:output method="text"/>

<xsl:param name="vcfchrom"/>
<xsl:param name="vcfpos"/>
<xsl:param name="vcfref"/>
<xsl:param name="vcfalt"/>




<xsl:template match="/">
<xsl:apply-templates select="evsData|snpList"/>
</xsl:template>



<xsl:template match="evsData">
<xsl:apply-templates select="snpList"/>
</xsl:template>


<xsl:template match="snpList">
<xsl:choose>
 <xsl:when test="chromosome=$vcfchrom and chrPosition=$vcfpos and refAllele=$vcfref">
<xsl:apply-templates select="clinicalLink|rsIds|uaMAF|aaMAF|totalMAF|avgSampleReadDepth|geneList|conservationScore|conservationScoreGERP|gwasPubmedIds|onExomeChip|gwasPubmedIds"/>
  <xsl:if test="altAlleles!=$vcfalt"><xsl:value-of select="concat('EVS_CONFLICTALT=',altAlleles,';')"/></xsl:if>
  </xsl:when>
  <xsl:otherwise/>
</xsl:choose>
</xsl:template>


<xsl:template match="clinicalLink|rsIds|uaMAF|aaMAF|totalMAF|avgSampleReadDepth|geneList|conservationScore|conservationScoreGERP|onExomeChip|gwasPubmedIds">
<xsl:if test="string-length(normalize-space(text()))&gt;0">
<xsl:value-of select="concat('EVS_',translate(name(.),'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'))"/>
<xsl:text>=</xsl:text>
<xsl:value-of select="text()"/>
<xsl:text>;</xsl:text>
</xsl:if>
</xsl:template>


</xsl:stylesheet>
