<?xml version='1.0' ?>
<xsl:stylesheet
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns="http://www.w3.org/1999/xhtml"
	version='1.0'
	xmlns:x="http://www.ibm.com/xmlns/prod/2009/jsonx"
	>
<xsl:import href="mod.command.xsl"/>
<xsl:output method="xml" indent="yes" />


<xsl:template match="/">
<xsl:apply-templates/>
</xsl:template>



<xsl:template match="c:app">
<jsonx:object xmlns:jsonx="http://www.ibm.com/xmlns/prod/2009/jsonx">
  <jsonx:string name="accessibility">Public</jsonx:string>
  <jsonx:string name="affiliation">univ-nantes.fr</jsonx:string>
  <jsonx:string name="cost">Free</jsonx:string>
  <jsonx:array name="platform">
    <jsonx:string>Linux</jsonx:string>
    <jsonx:string>Mac</jsonx:string>
  </jsonx:array>
  <jsonx:string name="version">1.0</jsonx:string>
  <jsonx:string name="homepage"><xsl:apply-templates select="." mode="wikiurl"/></jsonx:string>
  <jsonx:array name="function">
    <jsonx:object>
      <jsonx:array name="input">
      	<xsl:if test="c:snippet[@id='ref.faidx']">
      		<jsonx:object>
		      <jsonx:object name="dataType">
		        <jsonx:string name="term">File name</jsonx:string>
		        <jsonx:string name="uri">http://edamontology.org/data_1050</jsonx:string>
		      </jsonx:object>
		      <jsonx:array name="dataFormat">
		      	 <jsonx:object>
		          <jsonx:string name="term">FASTA</jsonx:string>
		          <jsonx:string name="uri">http://edamontology.org/format_1929</jsonx:string>
		         </jsonx:object>
		      </jsonx:array>
		    </jsonx:object>
      	</xsl:if>
      	<xsl:if test="c:input[@type='sam'] or c:input[@type='bam']">
       		<jsonx:object>
		      <jsonx:object name="dataType">
		        <jsonx:string name="term">File name</jsonx:string>
		        <jsonx:string name="uri">http://edamontology.org/data_1050</jsonx:string>
		      </jsonx:object>
		      <jsonx:array name="dataFormat">
		        <jsonx:object>
		          <jsonx:string name="term">BAM</jsonx:string>
		          <jsonx:string name="uri">http://edamontology.org/format_2572</jsonx:string>
		        </jsonx:object>
		      </jsonx:array>
		    </jsonx:object>
		</xsl:if>
       	<xsl:if test="c:input[@type='vcf']">
       		<jsonx:object>
		      <jsonx:object name="dataType">
		        <jsonx:string name="term">File name</jsonx:string>
		        <jsonx:string name="uri">http://edamontology.org/data_1050</jsonx:string>
		      </jsonx:object>
		      <jsonx:array name="dataFormat">
		        <jsonx:object>
		          <jsonx:string name="term">VCF</jsonx:string>
		          <jsonx:string name="uri">http://edamontology.org/format_3016</jsonx:string>
		        </jsonx:object>
		      </jsonx:array>
		    </jsonx:object>
		</xsl:if>
      </jsonx:array>
      <jsonx:array name="output">
        <jsonx:object>
          <jsonx:object name="dataType">
            <jsonx:string name="term">File name</jsonx:string>
            <jsonx:string name="uri">http://edamontology.org/data_1050</jsonx:string>
          </jsonx:object>
          <jsonx:string name="dataDescription">any format</jsonx:string>
          <jsonx:array name="dataFormat">
          	<xsl:if test="c:output/@type='sam' or c:output/@type='bam'">
            <jsonx:object>
              <jsonx:string name="term">BAM</jsonx:string>
              <jsonx:string name="uri">http://edamontology.org/format_2572</jsonx:string>
            </jsonx:object>
            </xsl:if>
          	<xsl:if test="c:output/@type='vcf'">
            <jsonx:object>
              <jsonx:string name="term">VCF</jsonx:string>
              <jsonx:string name="uri">http://edamontology.org/format_3016</jsonx:string>
            </jsonx:object>
            </xsl:if>
            <xsl:if test="not(c:output/@type='sam' or c:output/@type='bam' or c:output/@type='vcf')">
            <xsl:comment>no elixir output</xsl:comment>
            <jsonx:object>
              <jsonx:string name="term">Textual format</jsonx:string>
              <jsonx:string name="uri">http://edamontology.org/format_2330</jsonx:string>
            </jsonx:object>
            </xsl:if>
          </jsonx:array>
        </jsonx:object>
      </jsonx:array>
      <jsonx:array name="functionName">
        <jsonx:object>
          <jsonx:string name="term">Formatting</jsonx:string>
          <jsonx:string name="uri">http://edamontology.org/operation_0335</jsonx:string>
        </jsonx:object>
      </jsonx:array>
      <jsonx:string name="functionDescription"><xsl:apply-templates select="c:description"/></jsonx:string>
    </jsonx:object>
  </jsonx:array>
  <jsonx:string name="description"><xsl:apply-templates select="c:description"/></jsonx:string>
  <jsonx:object name="docs">
    <jsonx:string name="docsTermsOfUse">https://opensource.org/licenses/MIT</jsonx:string>
    <jsonx:string name="docsGithub"><xsl:apply-templates select="." mode="wikiurl"/></jsonx:string>
    <jsonx:string name="docsHome"><xsl:apply-templates select="." mode="wikiurl"/></jsonx:string>
    <jsonx:string name="docsCitationInstructions">https://github.com/lindenb/jvarkit/wiki/Citing</jsonx:string>
    <jsonx:string name="docsDownloadSource">https://github.com/lindenb/jvarkit/archive/master.zip</jsonx:string>
    <jsonx:string name="docsDownload">https://github.com/lindenb/jvarkit/archive/master.zip</jsonx:string>
  </jsonx:object>
  <jsonx:array name="collection">
    <jsonx:string>jvarkit</jsonx:string>
  </jsonx:array>
  <jsonx:object name="credits">
    <jsonx:array name="creditsInstitution">
      <jsonx:string>Institut du Thorax, Nantes, France</jsonx:string>
    </jsonx:array>
    <jsonx:array name="creditsDeveloper">
      <jsonx:string>Pierre Lindenbaum</jsonx:string>
    </jsonx:array>
  </jsonx:object>
  <jsonx:array name="interface">
    <jsonx:object>
      <jsonx:string name="interfaceType">Command line</jsonx:string>
    </jsonx:object>
  </jsonx:array>
  <jsonx:string name="name"><xsl:value-of select="@app"/></jsonx:string>
  <jsonx:array name="topic">
    <jsonx:object>
      <jsonx:string name="term">Omics</jsonx:string>
      <jsonx:string name="uri">http://edamontology.org/topic_3391</jsonx:string>
    </jsonx:object>
  </jsonx:array>
  <jsonx:string name="license">MIT License</jsonx:string>
  <jsonx:array name="language">
    <jsonx:string>Java</jsonx:string>
    <xsl:if test="c:snippet[@id='javascript']">
     <jsonx:string>Javascript</jsonx:string>
    </xsl:if>
  </jsonx:array>
  <jsonx:array name="resourceType">
    <jsonx:string>Tool</jsonx:string>
  </jsonx:array>
  <jsonx:string name="maturity">Stable</jsonx:string>
  <jsonx:array name="contact">
    <jsonx:object>
      <jsonx:string name="contactURL">https://github.com/lindenb</jsonx:string>
      <jsonx:string name="contactName">Pierre Lindenbaum</jsonx:string>
      <jsonx:array name="contactRole">
        <jsonx:string>Developer</jsonx:string>
        <jsonx:string>Maintainer</jsonx:string>
        <jsonx:string>Helpdesk</jsonx:string>
      </jsonx:array>
    </jsonx:object>
  </jsonx:array>
  <jsonx:object name="publications">
    <jsonx:string name="publicationsPrimaryID">doi:10.6084/m9.figshare.1425030.v1</jsonx:string>
  </jsonx:object>
</jsonx:object>
</xsl:template>




</xsl:stylesheet>

