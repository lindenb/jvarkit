<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0"
      xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
      xmlns:fo="http://www.w3.org/1999/XSL/Format">
  <xsl:output method="xml" indent="yes"/>
  <xsl:template match="/">
    <fo:root>
      <fo:layout-master-set>
        <fo:simple-page-master master-name="A4-portrait"
              page-height="21.0cm" page-width="40.7cm" margin="1cm">
          <fo:region-body/>
        </fo:simple-page-master>
      </fo:layout-master-set>
	<xsl:apply-templates select="html"/>
    </fo:root>
  </xsl:template>


<xsl:template match="html">
	<xsl:apply-templates select="body"/>
</xsl:template>

<xsl:template match="body">
	<xsl:apply-templates select="div/div/div[h3]" mode="variant"/>
</xsl:template>

<xsl:template match="div" mode="variant">
      <fo:page-sequence master-reference="A4-portrait">
        <fo:flow flow-name="xsl-region-body">
	   <xsl:apply-templates/>
        </fo:flow>
      </fo:page-sequence>
</xsl:template>

  <xsl:template match="table">
       <xsl:apply-templates select="thead/caption"/>
      <fo:table border="solid 0.1mm black">
	<xsl:choose>
	   <xsl:when test="thead/caption/text()='Genotype Types'"><xsl:attribute name="width">15cm</xsl:attribute></xsl:when>
	   <xsl:when test="thead/caption/text()='FILTERS'"><xsl:attribute name="width">20cm</xsl:attribute></xsl:when>
	   <xsl:when test="thead/caption/text()='Alleles'"><xsl:attribute name="width">20cm</xsl:attribute></xsl:when>
	   <xsl:when test="thead/caption/text()='SpliceAI'"><xsl:attribute name="width">20cm</xsl:attribute></xsl:when>
	   <xsl:when test="thead/caption/text()='ANN'"><xsl:attribute name="width">50cm</xsl:attribute></xsl:when>
	</xsl:choose>
    	<xsl:apply-templates select="thead"/>
    	<xsl:apply-templates select="tbody"/>
      </fo:table>
	<fo:block>
    		<fo:leader />
	</fo:block>
  </xsl:template>


  <xsl:template match="caption">
   <fo:block>
	 <xsl:apply-templates/>
   </fo:block>
  </xsl:template>

  <xsl:template match="thead">
    <fo:table-header>
	<xsl:apply-templates select="tr[1]"/>
    </fo:table-header>
  </xsl:template>


  <xsl:template match="tbody">
	<fo:table-body  border="inherit">
	<xsl:apply-templates select="tr"/>
	</fo:table-body>
  </xsl:template>



  <xsl:template match="tr">
	<fo:table-row  border="inherit">
		<xsl:choose>
		<xsl:when test="local-name(..)='thead'">
			<xsl:attribute name="background-color">#F2F2F2</xsl:attribute>
		</xsl:when>
		<xsl:when test="count(preceding-sibling::tr) mod 2 = 0 ">
			<xsl:attribute name="background-color">#F3F2E3</xsl:attribute>
		</xsl:when>
		</xsl:choose>
		<xsl:apply-templates select="th|td"/>
	</fo:table-row>
  </xsl:template>

  <xsl:template match="th">
	<fo:table-cell border="inherit"  wrap-option="wrap">
		<fo:block>
		<xsl:apply-templates/>
		</fo:block>
	</fo:table-cell>
  </xsl:template>

  <xsl:template match="td">
	<fo:table-cell>
		<fo:block>
		<xsl:apply-templates/>
		</fo:block>
	</fo:table-cell>
  </xsl:template>


  <xsl:template match="h1">
    <fo:block>
	<xsl:apply-templates/>
    </fo:block>
  </xsl:template>


  <xsl:template match="h2">
    <fo:block>
	<xsl:apply-templates/>
    </fo:block>
  </xsl:template>


  <xsl:template match="code">
    <fo:block>
	<xsl:apply-templates/>
    </fo:block>
  </xsl:template>

  <xsl:template match="span">
    <fo:block>
	<xsl:choose>
	  <xsl:when test="@style='color:green'"><xsl:attribute name="color">green</xsl:attribute></xsl:when>
	  <xsl:when test="@style='color:red'"><xsl:attribute name="color">red</xsl:attribute></xsl:when>
	  <xsl:when test="@style='color:yellow'"><xsl:attribute name="color">yellow</xsl:attribute></xsl:when>
	</xsl:choose>
	<xsl:apply-templates/>
    </fo:block>
  </xsl:template>



  <xsl:template match="h3">
    <fo:block font-weight="bold" font-size="1.17em">
	<xsl:apply-templates/>
    </fo:block>
  </xsl:template>

  <xsl:template match="text()">
    <fo:block><xsl:value-of select="."/></fo:block>
  </xsl:template>

  <xsl:template match="a[@name]">
  </xsl:template>

  <xsl:template match="a[@href]">
    <fo:block><fo:basic-link color="blue"  text-decoration="underline">
    <xsl:choose>
      <xsl:when test="starts-with(@href,'#')">
        <xsl:attribute name="internal-destination">
          <xsl:value-of select="substring-after(@href,'#')"/>
        </xsl:attribute>
      </xsl:when>
      <xsl:otherwise>
        <xsl:attribute name="external-destination">
          <xsl:text>url('</xsl:text>
          <xsl:value-of select="@href"/>
          <xsl:text>')</xsl:text>
        </xsl:attribute>
      </xsl:otherwise>
    </xsl:choose>
    <xsl:if test="@title">
      <xsl:attribute name="role">
        <xsl:value-of select="@title"/>
      </xsl:attribute>
    </xsl:if>
    <xsl:apply-templates/>
    </fo:basic-link></fo:block>
  </xsl:template>

</xsl:stylesheet>
