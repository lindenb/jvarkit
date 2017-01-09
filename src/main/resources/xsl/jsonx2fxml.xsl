<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns:j="http://www.ibm.com/xmlns/prod/2009/jsonx"
	xmlns:fx="http://javafx.com/fxml"
	version="1.0"
	>
<xsl:output method="xml" indent="yes" encoding="UTF-8"/>

<xsl:template match="/">

<xsl:apply-templates select="/j:object/j:array[@name='arguments']/j:object" mode="java"/>
<xsl:apply-templates select="/j:object/j:array[@name='arguments']/j:object" mode="java2"/>
<command>
<xsl:apply-templates select="/j:object/j:array[@name='arguments']/j:object" mode="xml"/>
</command>
</xsl:template>

<xsl:template match="j:object[j:string[@name='type'] = 'Boolean' or j:string[@name='type'] = 'boolean' ]" mode="java2">
 			new OptionBuilder(<xsl:value-of select="substring(j:string[@name='name']/text(),3)"/>,"<xsl:value-of select="j:string[@name='name']/text()"/>").fill(args);		
</xsl:template>

<xsl:template match="j:object[j:string[@name='type'] = 'Boolean' or j:string[@name='type'] = 'boolean' ]" mode="java">
 	@FXML
 	private CheckBox <xsl:value-of select="substring(j:string[@name='name']/text(),3)"/>;
</xsl:template>

<xsl:template match="j:object[j:string[@name='type'] = 'Boolean' or j:string[@name='type'] = 'boolean' ]" mode="xml">	
 	<CheckBox>
 			<xsl:attribute name="fx:id">
 				<xsl:value-of select="substring(j:string[@name='name']/text(),3)"/>
 			</xsl:attribute>
 			<xsl:attribute name="selected">
 				<xsl:value-of select="j:string[@name='defaultValue']/text()"/>
 			</xsl:attribute>
 		<label><xsl:value-of select="j:string[@name='summary']/text()"/></label>
 		<description><xsl:value-of select="j:string[@name='fulltext']/text()"/></description>
 	</CheckBox>
</xsl:template>

<xsl:template match="j:object" mode="xml">

</xsl:template>

<xsl:template match="j:object" mode="java">

</xsl:template>

<xsl:template match="j:object" mode="java2">

</xsl:template>
</xsl:stylesheet>
