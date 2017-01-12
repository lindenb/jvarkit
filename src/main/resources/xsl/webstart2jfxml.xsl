<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:fx="http://javafx.com/fxml"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns="http://javafx.com/fxml"
	>
<xsl:output method="xml" indent="yes"/>
<xsl:param name="mainclass"></xsl:param>

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
      <MenuItem text="About..." onAction="#doMenuAbout"/>
      <MenuItem text="Quit" onAction="#doMenuQuit"/>
    </Menu>
  </MenuBar>
  <BorderPane>
	<padding>
			<Insets bottom="5" left="5" right="5" top="5" />
	</padding>
  	<top>
  		<Label  style="-fx-font-weight:bold;-fx-font-size:150%; -fx-text-fill: white; -fx-background-color:darkgray;"  wrapText="true">
  			<xsl:attribute name="text">
  				<xsl:value-of select="description"/>
  				
  				<xsl:choose>
		  			<xsl:when test="contains($mainclass,'gatk')">
			  			<xsl:text> This is a wrapper for a GATK tool licensed under the Broad Institute software license agreement.</xsl:text>
			  		</xsl:when>
			  		<xsl:when test="contains($mainclass,'picard')">
			  			<xsl:text> This is a wrapper for a picard tool licensed under the MIT license .</xsl:text>
			  		</xsl:when>
      			</xsl:choose>
  				
			</xsl:attribute>
			<padding>
				<Insets bottom="5" left="5" right="5" top="5" />
			</padding>
  		</Label>
  	</top>
  	<center>
	  <ScrollPane>
	    <VBox alignment="CENTER_RIGHT" style="-fx-padding: 15px;">
			<xsl:for-each select="options/*">
				<xsl:apply-templates select="." mode="frame"/>
				<Separator/>
			</xsl:for-each>
	   </VBox>
	 </ScrollPane>
	</center>
	<bottom>
		<BorderPane>
			<top>
				<FlowPane  alignment="CENTER_RIGHT">
					<Button style="-fx-font-weight:bold;-fx-font-size:150%; -fx-background-color: red; -fx-text-fill: white;" text="STOP" fx:id="cancelCommandButton"/>
					<Button style="-fx-font-weight:bold;-fx-font-size:150%; -fx-background-color: green; -fx-text-fill: white;" text="GO  !"  fx:id="runCommandButton"/>
					<padding>
						<Insets bottom="5" left="5" right="5" top="5" />
					</padding>
				</FlowPane>
			</top>
			<center>
				<BorderPane>
					<padding>
						<Insets bottom="5" left="5" right="5" top="5" />
					</padding>
					<top>
						<Label text="Console:"/>
					</top>
					<center>
						<ScrollPane>
							<TextArea fx:id="console" text="" editable="false"  prefColumnCount="1000" 
									prefRowCount="15" 
									promptText="console"
									style="-fx-text-fill: white; -fx-background-color:lightgray; -fx-text-fill: black; -fx-font-size:10px; -fx-font-family:monospaced;"/>
						</ScrollPane>
					</center>
				</BorderPane>
			</center>
		</BorderPane>
	</bottom>
	</BorderPane>
</VBox>
</xsl:template>

<xsl:template match="*" mode="frame">
<BorderPane style="-fx-padding: 5 5 5 5;" xmlns:fx="http://javafx.com/fxml">
<xsl:if test="local-name(.)!='CheckBox'">
<top>
	<xsl:apply-templates select="." mode="labelt"/>
</top>
</xsl:if>
<center>
	 <FlowPane style="-fx-padding: 5 0 5 20;">
		<xsl:apply-templates select="."/>
	</FlowPane>
</center>
 <right>
 	<xsl:apply-templates select="." mode="labelr"/>
 </right>
</BorderPane>
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
		<xsl:if test="local-name(.)='ComboBox'">
			<xsl:text> Available options: </xsl:text>
			<xsl:for-each select="options/option">
				<xsl:text>**</xsl:text>
	 	 		<xsl:value-of select="@value"/>
	 	 		<xsl:text>** : </xsl:text>
	 	 		<xsl:value-of select="text()"/>
	 	 		<xsl:text>. </xsl:text>
	 	 	</xsl:for-each>
		</xsl:if>
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
	 <xsl:copy-of select="@id|@fx:id|@filter|@saveKey|@required|@open|@directory|@remember" />
 </com.github.lindenb.jvarkit.jfx.components.FileChooserPane>
</xsl:template>

<xsl:template match="FilesChooserPane|com.github.lindenb.jvarkit.jfx.components.FilesChooserPane">
 <com.github.lindenb.jvarkit.jfx.components.FilesChooserPane>
	 <xsl:copy-of select="@id|@fx:id|@filter|@saveKey|@minCardinality|@maxCardinality" />
 </com.github.lindenb.jvarkit.jfx.components.FilesChooserPane>
</xsl:template>

<xsl:template match="verbatim">
	  <xsl:copy-of select="*[not(name()='label' or name()='description')]"/>
</xsl:template>


<xsl:template match="GatkResource">
	<HBox spacing="5">
		<Label>
			<xsl:attribute name="text">
				<xsl:text>Resource </xsl:text>
				<xsl:value-of select="@suffix"/>
				<xsl:text> name:</xsl:text>
			</xsl:attribute>
		</Label>
		<TextField >
				<xsl:attribute name="fx:id">
				<xsl:value-of select="@fx:id"/>
				<xsl:text>_name</xsl:text>
				<xsl:value-of select="@suffix"/>
			</xsl:attribute>		
		</TextField>
		<com.github.lindenb.jvarkit.jfx.components.FileChooserPane >
				<xsl:attribute name="fx:id">
				<xsl:value-of select="@fx:id"/>
				<xsl:text>_file</xsl:text>
				<xsl:value-of select="@suffix"/>
			</xsl:attribute>		
			<xsl:copy-of select="@filter|@saveKey|@required|@remember" />
		</com.github.lindenb.jvarkit.jfx.components.FileChooserPane>
	</HBox>
</xsl:template>




<xsl:template match="Spinner">
 <Spinner>
	 <xsl:copy-of select="@fx:id" />
	 <xsl:copy-of select="valueFactory"/>
 </Spinner>
</xsl:template>
			
<xsl:template match="TextArea">
<ScrollPane>
	 <TextArea>
		 <xsl:copy-of select="@id|@fx:id|@text|@promptText|@prefColumnCount|@prefRowCount|@wrapText" />
	 </TextArea>
</ScrollPane>
</xsl:template>

<xsl:template match="TextField">
 <TextField>
	 <xsl:copy-of select="@id|@fx:id|@text|@promptText|@prefColumnCount" />
 </TextField>
</xsl:template>
			
<xsl:template match="CheckBox|Checkbox">
 <CheckBox>
	 <xsl:copy-of select="@id|@fx:id|@selected|@allowIndeterminate" />
	 <xsl:attribute name="text">
	 		<xsl:choose>
				<xsl:when test="label"><xsl:value-of select="label"/></xsl:when>
				<xsl:otherwise><xsl:value-of select="@label"/></xsl:otherwise>
			</xsl:choose>
	 </xsl:attribute>
 </CheckBox>
</xsl:template>

<xsl:template match="Combobox|ComboBox">
 <ComboBox>
	 <xsl:copy-of select="@id|@fx:id" />
	 <items>
	 	 <javafx.collections.FXCollections fx:factory="observableArrayList">
	 	 	<xsl:for-each select="options/option">
	 	 		<java.lang.String>
	 	 			<xsl:attribute name="fx:value">
	 	 				<xsl:value-of select="@value"/>
	 	 			</xsl:attribute>
	 	 		</java.lang.String>
	 	 	</xsl:for-each>
	 	 </javafx.collections.FXCollections>
	 </items>
	 <value>
	 	<java.lang.String>
	 		<xsl:attribute name="fx:value">
			 <xsl:choose>
			 	<xsl:when test="options/option[@selected='true']">
			 		<xsl:value-of select="options/option[@selected='true']/@value"/>
			 	</xsl:when>
			 	<xsl:otherwise>
			 		<xsl:value-of select="options/option[1]/@value"/>
			 	</xsl:otherwise>
			 </xsl:choose>
	 	</xsl:attribute>
	 	 </java.lang.String>
	</value>
 </ComboBox>
</xsl:template>



</xsl:stylesheet>

