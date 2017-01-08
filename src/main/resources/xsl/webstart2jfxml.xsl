<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:fx="http://javafx.com/fxml"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns="http://javafx.com/fxml"
	>
<xsl:output method="xml" indent="yes"/>

<xsl:template match="/">
<xsl:apply-templates select="command"/>
</xsl:template>

<xsl:template match="command">
<xsl:processing-instruction name="import">javafx.scene.layout.*</xsl:processing-instruction>
<xsl:processing-instruction name="import">javafx.scene.paint.*</xsl:processing-instruction>
<xsl:processing-instruction name="import">javafx.scene.control.*</xsl:processing-instruction>
<xsl:processing-instruction name="import">javafx.geometry.*</xsl:processing-instruction>

<VBox  xmlns="http://javafx.com/fxml">
  <MenuBar>
    <Menu text="File">
      <MenuItem text="Quit"/>
    </Menu>
  </MenuBar>
  <BorderPane>
  	<top>
  		<Label  style="-fx-font-weight:bold;-fx-font-size:150%"  maxWidth="1000" wrapText="true">
  			<xsl:attribute name="text">
  				<xsl:value-of select="description"/>
			</xsl:attribute>
  		</Label>
  	</top>
  	<center>
	  <ScrollPane>
	    <VBox alignment="CENTER_RIGHT" style="-fx-padding: 15px;">
			<xsl:for-each select="options/*">
				<BorderPane style="-fx-padding: 5 5 5 5;" xmlns:fx="http://javafx.com/fxml">
					<top>
						<xsl:apply-templates select="." mode="labelt"/>
					</top>
					<center>
						 <FlowPane style="-fx-padding: 5 0 5 20;">
							<xsl:apply-templates select="."/>
						</FlowPane>
					</center>
					 <right>
					 	<xsl:apply-templates select="." mode="labelr"/>
					 </right>
 				</BorderPane>
				<Separator/>
			</xsl:for-each>
	   </VBox>
	 </ScrollPane>
	</center>
	<bottom>
		<BorderPane>
			<top>
				<HBox>
					<Button text="Run" onAction="#doCommandStart" id="runbutton"/>
				</HBox>
			</top>
			<center>
				<ScrollPane>
					<TextArea fx:id="console" text="Hello"/>
				</ScrollPane>
			</center>
		</BorderPane>
	</bottom>
	</BorderPane>
</VBox>
</xsl:template>


<xsl:template match="*" mode="labelr">
<Label maxWidth="500" wrapText="true" textAlignment="JUSTIFY">
	<xsl:attribute name="text">
		<xsl:choose>
			<xsl:when test="description"><xsl:value-of select="description"/></xsl:when>
			<xsl:when test="@description"><xsl:value-of select="@description"/></xsl:when>
			<xsl:when test="label"><xsl:value-of select="label"/></xsl:when>
			<xsl:otherwise><xsl:value-of select="@label"/></xsl:otherwise>
		</xsl:choose>
	</xsl:attribute>
</Label>
</xsl:template>

<xsl:template match="*" mode="labelt">
<Label style="-fx-font-weight:bold;">
	<xsl:attribute name="text">
			<xsl:choose>
				<xsl:when test="label"><xsl:value-of select="label"/></xsl:when>
				<xsl:otherwise><xsl:value-of select="@label"/></xsl:otherwise>
			</xsl:choose>
	</xsl:attribute>
	<tooltip>
	  <Tooltip>
		<xsl:attribute name="text">
			<xsl:choose>
				<xsl:when test="description"><xsl:value-of select="description"/></xsl:when>
				<xsl:when test="@description"><xsl:value-of select="@description"/></xsl:when>
				<xsl:when test="label"><xsl:value-of select="label"/></xsl:when>
				<xsl:otherwise><xsl:value-of select="@label"/></xsl:otherwise>
			</xsl:choose>
		</xsl:attribute>
	  </Tooltip>
	</tooltip>
</Label>
</xsl:template>

<xsl:template match="FileChooserPane|com.github.lindenb.jvarkit.jfx.components.FileChooserPane">
 <com.github.lindenb.jvarkit.jfx.components.FileChooserPane>
	 <xsl:copy-of select="@id|@fx:id|@filter|@saveKey|@required|@open" />
 </com.github.lindenb.jvarkit.jfx.components.FileChooserPane>
</xsl:template>

<xsl:template match="Spinner">
 <Spinner>
	 <xsl:copy-of select="@fx:id" />
	 <xsl:copy-of select="valueFactory"/>
 </Spinner>
</xsl:template>
			

			
			




</xsl:stylesheet>

