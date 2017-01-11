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

<!-- ============================================================================== -->


<xsl:template match="j:object" mode="java2">
 			new OptionBuilder(this.<xsl:value-of select="substring(j:string[@name='name']/text(),3)"/>,"<xsl:value-of select="j:string[@name='name']/text()"/>").fill(args);		
</xsl:template>

<!-- ============================================================================== -->

<xsl:template match="j:object[j:string[@name='type'] = 'Boolean' or j:string[@name='type'] = 'boolean' ]" mode="java">
 	@FXML private CheckBox <xsl:value-of select="substring(j:string[@name='name']/text(),3)"/>;
</xsl:template>

<xsl:template match="j:object[j:string[@name='type'] = 'Boolean' or j:string[@name='type'] = 'boolean' ]" mode="xml">	
 	<CheckBox>
 			<xsl:attribute name="selected">
 				<xsl:value-of select="j:string[@name='defaultValue']/text()"/>
 			</xsl:attribute>
 			<xsl:apply-templates select="." mode="inner"/>
 	</CheckBox>
</xsl:template>

<!-- ============================================================================== -->

<xsl:template match="j:object[j:string[@name='type'] = 'String' or j:string[@name='type'] = 'string' ]" mode="java">
 	@FXML private TextField <xsl:value-of select="substring(j:string[@name='name']/text(),3)"/>;
</xsl:template>

<xsl:template match="j:object[j:string[@name='type'] = 'String' or j:string[@name='type'] = 'string' ]" mode="xml">	
 	<TextField>
 			<xsl:attribute name="text">
 				<xsl:value-of select="j:string[@name='defaultValue']/text()"/>
 			</xsl:attribute>
 			<xsl:apply-templates select="." mode="inner"/>
 	</TextField>
</xsl:template>

<!-- ============================================================================== -->

<!-- ============================================================================== -->

<xsl:template match="j:object[j:string[@name='type'] = 'Integer' or j:string[@name='type'] = 'int' ]" mode="java">
 	@FXML private Spinner&lt;Integer&gt; <xsl:value-of select="substring(j:string[@name='name']/text(),3)"/>;
</xsl:template>

<xsl:template match="j:object[j:string[@name='type'] = 'Integer' or j:string[@name='type'] = 'int' ]" mode="xml">	
 	<Spinner>
 			<xsl:apply-templates select="." mode="inner"/>
 			<valueFactory>
		        <SpinnerValueFactory.IntegerSpinnerValueFactory  >
		        		<xsl:attribute name="initialValue">
		        			<xsl:value-of select="j:string[@name='defaultValue']/text()"/>
		        		</xsl:attribute>
		        		<xsl:if test="j:string[@name='minValue']/text() != '-Infinity' ">
		        			<xsl:attribute name="min">
		        				<xsl:value-of select="j:string[@name='minValue']/text()"/>
		        			</xsl:attribute>
		        		</xsl:if>
		        		<xsl:if test="j:string[@name='maxValue']/text() != 'Infinity' ">
		        			<xsl:attribute name="max">
		        				<xsl:value-of select="j:string[@name='maxValue']/text()"/>
		        			</xsl:attribute>
		        		</xsl:if>
		        		<xsl:if test="j:string[@name='minValue']/text() = '-Infinity' ">
		        			<min>
		        					<java.lang.Integer fx:constant="MIN_VALUE"/>
		        			</min>
		        		</xsl:if>
		        		<xsl:if test="j:string[@name='maxValue']/text() = 'Infinity' ">
		        			<max>
		        				<java.lang.Integer fx:constant="MAX_VALUE"/>
		        			</max>
		        		</xsl:if>
		        </SpinnerValueFactory.IntegerSpinnerValueFactory>
		    </valueFactory>
 	</Spinner>
</xsl:template>

<!-- ============================================================================== -->
<xsl:template match="j:object[j:string[@name='type'] = 'Double' or j:string[@name='type'] = 'double' ]" mode="java">
 	@FXML private Spinner&lt;Double&gt; <xsl:value-of select="substring(j:string[@name='name']/text(),3)"/>;
</xsl:template>

<xsl:template match="j:object[j:string[@name='type'] = 'Double' or j:string[@name='type'] = 'double' ]" mode="xml">	
 	<Spinner>
 			<xsl:apply-templates select="." mode="inner"/>
 			<valueFactory>
		        <SpinnerValueFactory.DoubleSpinnerValueFactory  >
		        		<xsl:attribute name="initialValue">
		        			<xsl:value-of select="j:string[@name='defaultValue']/text()"/>
		        		</xsl:attribute>
		        		<xsl:if test="j:string[@name='minValue']/text() != '-Infinity' ">
		        			<xsl:attribute name="min">
		        				<xsl:value-of select="j:string[@name='minValue']/text()"/>
		        			</xsl:attribute>
		        		</xsl:if>
		        		<xsl:if test="j:string[@name='maxValue']/text() != 'Infinity' ">
		        			<xsl:attribute name="max">
		        				<xsl:value-of select="j:string[@name='maxValue']/text()"/>
		        			</xsl:attribute>
		        		</xsl:if>
		        		<xsl:if test="j:string[@name='minValue']/text() = '-Infinity' ">
		        			<min>
		        					<java.lang.Double fx:constant="MIN_VALUE"/>
		        			</min>
		        		</xsl:if>
		        		<xsl:if test="j:string[@name='maxValue']/text() = 'Infinity' ">
		        			<max>
		        				<java.lang.Double fx:constant="MAX_VALUE"/>
		        			</max>
		        		</xsl:if>
		        </SpinnerValueFactory.DoubleSpinnerValueFactory>
		    </valueFactory>
 	</Spinner>
</xsl:template>

<!-- ============================================================================== -->


<xsl:template match="j:object" mode="xml">
	
	<TODOTODOTODOTODO>
		<xsl:apply-templates select="." mode="inner"/>
	</TODOTODOTODOTODO>
	
</xsl:template>

<xsl:template match="j:object" mode="java">

</xsl:template>



<xsl:template match="j:object" mode="inner">	
 	<xsl:attribute name="fx:id">
 				<xsl:value-of select="substring(j:string[@name='name']/text(),3)"/>
		</xsl:attribute>
	<label><xsl:value-of select="j:string[@name='summary']/text()"/></label>
	<description><xsl:value-of select="normalize-space(j:string[@name='fulltext']/text())"/></description>
</xsl:template>


</xsl:stylesheet>
