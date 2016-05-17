<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	>
<xsl:output method="text"/>

<xsl:template match="c:app" mode="header">/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


<xsl:apply-templates select="c:history"/>

*/
</xsl:template>

<xsl:template match="c:app" mode="package">
<xsl:value-of select="@package"/>
</xsl:template>

<xsl:template match="c:app" mode="class-name">
<xsl:value-of select="@app"/>
</xsl:template>

<xsl:template match="c:app" mode="abstract-class-name">
<xsl:text>Abstract</xsl:text>
<xsl:apply-templates select="." mode="class-name"/>
</xsl:template>

<xsl:template match="c:app" mode="swing-name">
<xsl:apply-templates select="." mode="class-name"/>
<xsl:text>SwingUI</xsl:text>
</xsl:template>

<xsl:template match="c:app" mode="factory-class-name">
<xsl:message terminate="yes">Factory is deprecated</xsl:message>
</xsl:template>

<xsl:template match="c:app" mode="abstract-command-name">
<xsl:text>Abstract</xsl:text>
<xsl:apply-templates select="." mode="class-name"/>
<xsl:text>Command</xsl:text>
</xsl:template>

<xsl:template match="c:app" mode="label">
	<xsl:choose>
		<xsl:when test="c:label">
		@Override
		public String getLabel()
			{
			return "<xsl:value-of select="c:label"/>";
			}
		</xsl:when>
		<xsl:when test="@label">
		@Override
		public String getLabel()
			{
			return "<xsl:value-of select="@label"/>";
			}
		</xsl:when>
		<xsl:otherwise>
		@Override
		public String getLabel()
			{
			return "<xsl:apply-templates select="." mode="class-name"/>";
			}
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>
	
<xsl:template match="c:app" mode="description">
	<xsl:choose>
		<xsl:when test="c:description">
		@Override
		public String getDescription()
			{
			return "<xsl:value-of select="c:description"/>";
			}
		</xsl:when>
		<xsl:when test="@description">
		@Override
		public String getDescription()
			{
			return "<xsl:value-of select="@description"/>";
			}
		</xsl:when>
	</xsl:choose>
</xsl:template>

<xsl:template match="c:app" mode="jfx">
public static class <xsl:apply-templates select="." mode="class-name"/>Application extends javafx.application.Application
	{
    public void start(javafx.stage.Stage stage)
    	{
    	stage.setTitle("<xsl:apply-templates select="." mode="class-name"/>");
    	javafx.scene.layout.VBox parent = new javafx.scene.layout.VBox();
    	<xsl:apply-templates select=".//c:option" mode="jfx-ctrl"/>
    	javafx.scene.Scene scene = new javafx.scene.Scene(parent);
    	stage.setScene(scene);
        stage.show();
    	}
    public static void main(String args[])
    	{
    	javafx.application.Application.launch(<xsl:apply-templates select="." mode="class-name"/>Application.class,args);
    	}
    }
</xsl:template>


<xsl:template match="c:app" mode="html">
<xsl:if test="c:documentation">
	@Override
	protected void writeHtmlDoc(final javax.xml.stream.XMLStreamWriter w)
		throws javax.xml.stream.XMLStreamException
	{
	w.writeStartElement("div");
	
	w.writeStartElement("h3");
	w.writeCharacters("Description");
	w.writeEndElement();
	w.writeStartElement("div");
	w.writeEndElement();
	
	
	w.writeEndElement();
	}
</xsl:if>
</xsl:template>

<xsl:template match="c:option-group" mode="cli">
final OptionGroup <xsl:value-of select="generate-id()"/> = new OptionGroup();
<xsl:apply-templates select="c:option" mode="cli"/>
options.addOptionGroup(<xsl:value-of select="generate-id()"/>);
</xsl:template>

<xsl:template match="c:option" mode="label">
<xsl:choose>
	<xsl:when test="c:label"><xsl:value-of select="c:label"/></xsl:when>
	<xsl:when test="@label"><xsl:value-of select="@label"/></xsl:when>
	<xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:option" mode="description">
<xsl:choose>
	<xsl:when test="c:description"><xsl:apply-templates select="c:description"/></xsl:when>
	<xsl:when test="@description"><xsl:value-of select="@description"/></xsl:when>
	<xsl:otherwise><xsl:apply-templates select="." mode="label"/></xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:option" mode="OPTION_OPT">
<xsl:value-of select="concat('OPTION_',translate(@name,'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'))"/>
</xsl:template>

<xsl:template match="c:option" mode="cli">
options.addOption(org.apache.commons.cli.Option
	<xsl:choose>
		<xsl:when test="@opt">
			.builder(<xsl:apply-templates select="." mode="OPTION_OPT"/>)
		</xsl:when>
		<xsl:otherwise>
			.builder("<xsl:value-of select="@name"/>")
		</xsl:otherwise>
	</xsl:choose>
	<xsl:if test="@required">
	.required(<xsl:value-of select="@required"/>)
	</xsl:if>
	<xsl:if test="@optional-arg">
	.optionalArg(<xsl:value-of select="@optional-arg"/>)
	</xsl:if>
	<xsl:if test="@longOpt">
	.longOpt("<xsl:value-of select="@longOpt"/>")
	</xsl:if>
	<xsl:if test="@longopt">
	.longOpt("<xsl:value-of select="@longopt"/>")
	</xsl:if>
	<xsl:if test="@value-separator">
	.valueSeparator('<xsl:value-of select="@value-separator"/>')
	</xsl:if>
	<xsl:if test="@description">
	.desc("<xsl:value-of select="@description"/>"
		<xsl:if test="@default">
		+ ". default: <xsl:value-of select="@default"/>"
		</xsl:if>
		<xsl:if test="@type='input-file-set' or @type='string-list' or @type='string-set' or @type='uri-set'">
		+ ". Multiple calls to this option should end with double hyphen : --."
		</xsl:if>
		
		)
	</xsl:if>
	<xsl:if test="c:description">
	.desc("<xsl:apply-templates select="c:description"/>"
		<xsl:if test="@default">
		+ ". default: <xsl:value-of select="@default"/>"
		</xsl:if>
		)
	</xsl:if>
		<xsl:choose>
		<xsl:when test="@arg-name">
			.argName("<xsl:value-of select="@arg-name"/>")
		</xsl:when>
		<xsl:when test="@argname">
			.argName("<xsl:value-of select="@argname"/>")
		</xsl:when>
		<xsl:otherwise>
			.argName("<xsl:value-of select="@name"/>")
		</xsl:otherwise>
	</xsl:choose>
	<xsl:choose>
		<xsl:when test="@type='bool' or @type='boolean' or @type='Boolean' or @type='java.lang.Boolean'">
		.hasArg(false)
		</xsl:when>
		<xsl:when test="@type='int' or @type='java.lang.Integer' or @type='long' or @type='short' or @type='double' or @type='float'  or @type='number'">
		.hasArg(true)
		.type(org.apache.commons.cli.PatternOptionBuilder.NUMBER_VALUE)
		</xsl:when>
		<xsl:when test="@type='url' or @type='string' or @type='String'">
		.hasArg(true)
		.type(org.apache.commons.cli.PatternOptionBuilder.STRING_VALUE)
		</xsl:when>
		<xsl:when test="@type='object'">
		.hasArg(true)
		.type(org.apache.commons.cli.PatternOptionBuilder.OBJECT_VALUE)
		</xsl:when>
		<xsl:when test="@type='date'">
		.hasArg(true)
		.type(org.apache.commons.cli.PatternOptionBuilder.DATE_VALUE)
		</xsl:when>
		<xsl:when test="@type='file' or @type='output-file'">
		.hasArg(true)
		.type(org.apache.commons.cli.PatternOptionBuilder.FILE_VALUE)
		</xsl:when>
		<xsl:when test="@type='existing-file' or @type='input-file' or @type='input-directory'">
		.hasArg(true)
		.type(org.apache.commons.cli.PatternOptionBuilder.EXISTING_FILE_VALUE)
		</xsl:when>
		<xsl:when test="@type='input-file-set'">
		.hasArgs() //unlimited
		.type(org.apache.commons.cli.PatternOptionBuilder.EXISTING_FILE_VALUE)
		</xsl:when>
		<xsl:when test="@type='string-set' or @type='uri-set'">
		.hasArgs() //unlimited
		.type(org.apache.commons.cli.PatternOptionBuilder.STRING_VALUE)
		</xsl:when>
		<xsl:when test="@type='string-list'">
		.hasArgs() //unlimited
		.type(org.apache.commons.cli.PatternOptionBuilder.STRING_VALUE)
		</xsl:when>
		<xsl:otherwise>
			<xsl:message terminate="yes">option:cli unknown type <xsl:value-of select="@type"/></xsl:message>
		</xsl:otherwise>

	</xsl:choose>
	.build() );	
</xsl:template>


<xsl:template match="c:option[@type='long']" mode="validator">
LongValidator <xsl:apply-templates select="@name"/> 
</xsl:template>


<xsl:template match="c:option" mode="name">
<xsl:value-of select="@name"/>
</xsl:template>


<xsl:template match="c:option">
<xsl:variable name="cloneable">
	<xsl:apply-templates select="." mode="cloneable"/>
</xsl:variable>
<xsl:variable name="nilleable">
	<xsl:apply-templates select="." mode="nilleable"/>
</xsl:variable>

/** option <xsl:apply-templates select="." mode="name"/> */

public static final String <xsl:apply-templates select="." mode="OPTION_OPT"/> = "<xsl:value-of select="@opt"/>";

protected <xsl:apply-templates select="." mode="java-type"/><xsl:text> </xsl:text> <xsl:apply-templates select="." mode="name"/> = <xsl:choose>
		<xsl:when test="@type='input-file-set'"> new java.util.HashSet&lt;java.io.File&gt;()</xsl:when>
		<xsl:when test="@type='string-set' or @type='uri-set'"> new java.util.HashSet&lt;java.lang.String&gt;()</xsl:when>
		<xsl:when test="@type='string-list'"> new java.util.ArrayList&lt;java.lang.String&gt;()</xsl:when>
		<xsl:when test="@default and (not(@type) or @type='string' or @type='String' or @type='java.lang.String')">"<xsl:value-of select="@default"/>"</xsl:when>
		<xsl:when test="@default"><xsl:value-of select="@default"/></xsl:when>
		<xsl:when test="$nilleable = 'true'">null</xsl:when>
		<xsl:when test="@type='boolean' or @type='bool'"> false</xsl:when>
		<xsl:otherwise>0</xsl:otherwise>
	</xsl:choose>;

/** getter for <xsl:value-of select="@name"/> */
public <xsl:apply-templates select="." mode="java-type"/>
		<xsl:text> </xsl:text>
		<xsl:apply-templates select="." mode="getter"/>()
	{
	return this.<xsl:apply-templates select="." mode="name"/>;
	}

/** setter for <xsl:value-of select="@name"/> */
public  void  <xsl:apply-templates select="." mode="setter"/>( final <xsl:apply-templates select="." mode="java-type"/><xsl:text> </xsl:text><xsl:apply-templates select="." mode="name"/>)
	{
	this.<xsl:apply-templates select="." mode="name"/> = <xsl:choose>
		<xsl:when test="@type='input-file-set'">(<xsl:apply-templates select="." mode="java-type"/>)(<xsl:apply-templates select="." mode="name"/>==null?null: new java.util.HashSet&lt;java.io.File&gt;(<xsl:apply-templates select="." mode="name"/>))</xsl:when>
		<xsl:when test="@type='string-set' or @type='uri-set'">(<xsl:apply-templates select="." mode="java-type"/>)(<xsl:apply-templates select="." mode="name"/>==null?null: new java.util.HashSet&lt;java.lang.String&gt;(<xsl:apply-templates select="." mode="name"/>))</xsl:when>
		<xsl:when test="@type='string-list'">(<xsl:apply-templates select="." mode="java-type"/>)(<xsl:apply-templates select="." mode="name"/>==null?null: new java.util.ArrayList&lt;java.lang.String&gt;(<xsl:apply-templates select="." mode="name"/>))</xsl:when>
		<xsl:when test="$cloneable = 'true'">(<xsl:apply-templates select="." mode="java-type"/>)(<xsl:apply-templates select="." mode="name"/>==null?null:<xsl:apply-templates select="." mode="name"/>.clone())</xsl:when>
		<xsl:otherwise><xsl:apply-templates select="." mode="name"/></xsl:otherwise>
		</xsl:choose>;
	}

</xsl:template>

<xsl:template match="c:option" mode="copy">
this.<xsl:apply-templates select="." mode="setter"/>(factory.<xsl:apply-templates select="." mode="getter"/>());
</xsl:template>

<xsl:template match="c:option" mode="cloneable">
<xsl:variable name="nilleable">
	<xsl:apply-templates select="." mode="nilleable"/>
</xsl:variable>
<xsl:choose>
	<xsl:when test="@type='java.net.URL'">false</xsl:when>
	<xsl:when test="@type='string' or @type='String' or @type='java.lang.String'">false</xsl:when>
	<xsl:when test="@type='output-file'">false</xsl:when>
	<xsl:when test="@type='input-file'">false</xsl:when>
	<xsl:when test="@type='input-file-set' or @type='string-set'or @type='uri-set'">true</xsl:when>
	<xsl:when test="starts-with(@type,'java.lang')">false</xsl:when>
	<xsl:when test="@type='bool' or @type='boolean'">false</xsl:when>
	<xsl:when test="@type='int'">false</xsl:when>
	<xsl:when test="@type='input-directory'">false</xsl:when>
	<xsl:when test="@type='input-file-set' or @type='string-set' or @type='string-list'">true</xsl:when>
	<xsl:when test="starts-with(@type,'java.lang')">false</xsl:when>
	<xsl:when test="@type='bool' or @type='boolean'">false</xsl:when>
	<xsl:when test="@type='int'">false</xsl:when>
	<xsl:when test="@type='long'">false</xsl:when>
	<xsl:when test="@type='double'">false</xsl:when>
	<xsl:when test="@type='float'">false</xsl:when>
	<xsl:when test="$nilleable = 'true'">true</xsl:when>
	<xsl:otherwise>
		<xsl:message terminate='yes'>cloneable: unknown type <xsl:value-of select="@type"/>.</xsl:message>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:option" mode="nilleable">
<xsl:choose>
	<xsl:when test="@type='java.net.URL'">true</xsl:when>
	<xsl:when test="@type='output-file'">true</xsl:when>
	<xsl:when test="@type='input-file'">true</xsl:when>
	<xsl:when test="@type='input-file-set' or @type='string-set' or @type='uri-set'">true</xsl:when>
	<xsl:when test="@type='input-directory'">true</xsl:when>
	<xsl:when test="@type='input-file-set' or @type='string-set' or @type='string-list'">true</xsl:when>
	<xsl:when test="@type='java.net.URL'">true</xsl:when>
	<xsl:when test="@type='bool' or @type='boolean'">false</xsl:when>
	<xsl:when test="@type='string' or @type='String' or @type='java.lang.String'">true</xsl:when>
    <xsl:when test="starts-with(@type,'java.lang')">true</xsl:when>
	<xsl:when test="@type='int' or @type='double'or @type='long'">false</xsl:when>
	<xsl:when test="@type='float'">false</xsl:when>
	<xsl:otherwise>
		<xsl:message terminate='yes'>nilleable: unknown type <xsl:value-of select="@type"/> name=<xsl:value-of select="@name"/> opt=<xsl:value-of select="@opt"/>.</xsl:message>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>



<xsl:template match="c:option" mode="setter">
<xsl:variable name="s">
<xsl:call-template name="titleize">
	<xsl:with-param name="s" select="@name"/>
</xsl:call-template>
</xsl:variable>
<xsl:value-of select="concat('set',$s)"/>
</xsl:template>

<xsl:template match="c:option" mode="getter">
<xsl:variable name="s">
<xsl:call-template name="titleize">
	<xsl:with-param name="s" select="@name"/>
</xsl:call-template>
</xsl:variable>
<xsl:choose>
		<xsl:when test="@type='bool' or @type='boolean'  or @type='Boolean'  or @type='java.lang.Boolean'"><xsl:value-of select="concat('is',$s)"/></xsl:when>
		<xsl:otherwise><xsl:value-of select="concat('get',$s)"/></xsl:otherwise>
</xsl:choose>

</xsl:template>


<xsl:template match="c:option" mode="java-type">
<xsl:choose>
	<xsl:when test="@type='output-file'">java.io.File</xsl:when>
	<xsl:when test="@type='input-file'">java.io.File</xsl:when>
	<xsl:when test="@type='input-directory'">java.io.File</xsl:when>
	<xsl:when test="@type='input-file-set'">java.util.Set&lt;java.io.File&gt;</xsl:when>
	<xsl:when test="@type='string-set' or @type='uri-set'">java.util.Set&lt;java.lang.String&gt;</xsl:when>
	<xsl:when test="@type='string-list'">java.util.List&lt;java.lang.String&gt;</xsl:when>
	<xsl:when test="@type='int'">int</xsl:when>
	<xsl:when test="@type='long'">long</xsl:when>
	<xsl:when test="@type='double'">double</xsl:when>
	<xsl:when test="@type='float'">float</xsl:when>
	<xsl:when test="@type='bool' or @type='boolean'">boolean</xsl:when>
	<xsl:when test="@type='string' or @type='String' or @type='java.lang.String'">java.lang.String</xsl:when>
	<xsl:when test="@type='Integer' or @type='java.lang.Integer'">java.lang.Integer</xsl:when>
	<xsl:otherwise>
		<xsl:message terminate='yes'>java-type: unknown type "<xsl:value-of select="@type"/>".</xsl:message>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:option" mode="jfx-ctrl">
/* begin : javafx control for <xsl:value-of select="@name"/> */
final javafx.scene.control.Tooltip <xsl:value-of select="concat('tooltip',generate-id())"/> = new javafx.scene.control.Tooltip("<xsl:apply-templates select="." mode="description"/>");
final javafx.scene.layout.HBox  <xsl:value-of select="concat('hbox',generate-id())"/> = new javafx.scene.layout.HBox();
parent.getChildren().add( <xsl:value-of select="concat('hbox',generate-id())"/>);
final javafx.scene.control.Label <xsl:value-of select="concat('lbl',generate-id())"/> = new javafx.scene.control.Label("<xsl:apply-templates select="." mode="label"/>");
<xsl:value-of select="concat('hbox',generate-id())"/>.getChildren().add(<xsl:value-of select="concat('lbl',generate-id())"/>);
<xsl:choose>
	<xsl:when test="@type='input-file'">
	final javafx.scene.control.Button <xsl:value-of select="generate-id()"/> = new javafx.scene.control.Button();
	<xsl:value-of select="concat('hbox',generate-id())"/>.getChildren().add(<xsl:value-of select="generate-id()"/>);
	<xsl:value-of select="concat('lbl',generate-id())"/>.setLabelFor(<xsl:value-of select="generate-id()"/>);
	<xsl:value-of select="generate-id()"/>.setTooltip(<xsl:value-of select="concat('tooltip',generate-id())"/>);
	<xsl:value-of select="generate-id()"/>.setOnAction(new javafx.event.EventHandler&lt;javafx.event.ActionEvent&gt;() {
	    @Override public void handle(final javafx.event.ActionEvent evt)
	    	{
	        javafx.stage.FileChooser <xsl:value-of select="concat('fc',generate-id())"/> = new javafx.stage.FileChooser();
			<xsl:value-of select="concat('fc',generate-id())"/>.setTitle("<xsl:apply-templates select="." mode="label"/>");
			java.io.File selectedFile = <xsl:value-of select="concat('fc',generate-id())"/>.showOpenDialog(<xsl:value-of select="generate-id()"/>.getScene().getWindow());
			if( selectedFile == null ) return;
			<xsl:value-of select="generate-id()"/>.getProperties().put("file",selectedFile);
			<xsl:value-of select="generate-id()"/>.setText(selectedFile.getName());
	    	}
		});
	
	</xsl:when>
	<xsl:when test="@type='bool' or @type='boolean' or @type='java.lang.Boolean' or @type='Boolean'">
		final javafx.scene.control.CheckBox <xsl:value-of select="generate-id()"/> = new javafx.scene.control.CheckBox("<xsl:apply-templates select="." mode="label"/>");
		<xsl:value-of select="generate-id()"/>.setTooltip(<xsl:value-of select="concat('tooltip',generate-id())"/>);
		<xsl:value-of select="concat('hbox',generate-id())"/>.getChildren().remove(<xsl:value-of select="concat('lbl',generate-id())"/>);
		<xsl:value-of select="concat('hbox',generate-id())"/>.getChildren().add(<xsl:value-of select="generate-id()"/>);
	</xsl:when>
	
	<xsl:when test="@type='string' or @type='String' or @type='java.lang.String'">
		final javafx.scene.control.TextField <xsl:value-of select="generate-id()"/> = new javafx.scene.control.TextField();
		<xsl:value-of select="generate-id()"/>.setTooltip(<xsl:value-of select="concat('tooltip',generate-id())"/>);
		<xsl:value-of select="concat('hbox',generate-id())"/>.getChildren().add(<xsl:value-of select="generate-id()"/>);
		<xsl:if test="@default">
		<xsl:value-of select="generate-id()"/>.setText("<xsl:value-of select="@default"/>");
		</xsl:if>
	</xsl:when>
	<xsl:message terminate='yes'>jfx-ctrl:unknown type <xsl:value-of select="@type"/>.</xsl:message>
</xsl:choose>
/* end : javafx control for <xsl:value-of select="@name"/> */
</xsl:template>


<xsl:template match="c:option" mode="visit">
<xsl:variable name="process_option">
<xsl:choose>
		<xsl:when test="@name='formatout'  and ((//c:output/@type='sam' or //c:output/@type='bam' or //c:snippet[@id='write-sam']))">
			<xsl:text>false</xsl:text>
		</xsl:when>
	<xsl:otherwise>true</xsl:otherwise>
</xsl:choose>
</xsl:variable>
<xsl:choose>
<xsl:when test="$process_option = 'false'">
/** option  <xsl:value-of select="@name"/> will be processed elsewhere */
</xsl:when>
<xsl:otherwise>if(opt.getOpt().equals(<xsl:apply-templates select="." mode="OPTION_OPT"/>))
	{
	/* <xsl:value-of select="@name"/> : <xsl:value-of select="@type"/> */
	<xsl:choose>
		<xsl:when test="@name='http_proxy_str' and ../../c:snippet[@id='http.proxy']">
		
			final String <xsl:value-of select="generate-id()"/> = opt.getValue();
			final int colon = <xsl:value-of select="generate-id()"/>.lastIndexOf(':');
			if(colon==-1)
				{
				LOG.error("bad proxy \""+<xsl:value-of select="generate-id()"/>+"\"");
				return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;
				}
			LOG.debug("setting proxy "+<xsl:value-of select="generate-id()"/>);
			System.setProperty("http.proxyHost", <xsl:value-of select="generate-id()"/>.substring(0,colon));
			System.setProperty("http.proxyPort", <xsl:value-of select="generate-id()"/>.substring(colon+1));
			System.setProperty("https.proxyHost", <xsl:value-of select="generate-id()"/>.substring(0,colon));
			System.setProperty("https.proxyPort", <xsl:value-of select="generate-id()"/>.substring(colon+1));
		</xsl:when>

	
		<xsl:when test="@type='bool' or @type='boolean' or @type='Boolean' or @type='java.lang.Boolean' ">
			boolean <xsl:value-of select="generate-id()"/>;
			<xsl:choose>
				<xsl:when test="@default='true'">
				 <xsl:value-of select="generate-id()"/> = false;
				</xsl:when>
				<xsl:when test="@default='false'">
				 <xsl:value-of select="generate-id()"/> = true;
				</xsl:when>
				<xsl:otherwise>
					<xsl:message terminate='true'>c:option visit no @default value set for  <xsl:value-of select="@name"/></xsl:message>
				</xsl:otherwise>
			</xsl:choose>
		</xsl:when>
	
		<xsl:when test="not(@type) or @type='string' or @type='String' or @type='java.lang.String'">
		java.lang.String <xsl:value-of select="generate-id()"/> = opt.getValue();
		</xsl:when>
	
		<xsl:when test="@type='int'">
		int <xsl:value-of select="generate-id()"/> = 0;
		try { <xsl:value-of select="generate-id()"/> = Integer.parseInt(opt.getValue());}
		catch(Exception err) { LOG.error("Cannot cast "+opt.getValue()+" to integer",err); return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;}
		</xsl:when>
		
		
		<xsl:when test="@type='long'">
		long <xsl:value-of select="generate-id()"/> = 0L;
		try { <xsl:value-of select="generate-id()"/> = Long.parseLong(opt.getValue());}
		catch(Exception err) { LOG.error("Cannot cast "+opt.getValue()+" to long",err); return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;}
		</xsl:when>
		
		<xsl:when test="@type='float'">
		float <xsl:value-of select="generate-id()"/> = 0f;
		try { <xsl:value-of select="generate-id()"/> = Float.parseFloat(opt.getValue());}
		catch(Exception err) { LOG.error("Cannot cast "+opt.getValue()+" to float",err); return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;}
		</xsl:when>
		
		<xsl:when test="@type='double'">
		double <xsl:value-of select="generate-id()"/> = 0.0;
		try { <xsl:value-of select="generate-id()"/> = Double.parseDouble(opt.getValue());}
		catch(Exception err) { LOG.error("Cannot cast "+opt.getValue()+" to double",err); return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;}
		</xsl:when>
		
		<xsl:when test="@type='java.lang.Integer'">
		java.lang.Integer <xsl:value-of select="generate-id()"/> = null;
		try { <xsl:value-of select="generate-id()"/> =new Integer(opt.getValue());}
		catch(Exception err) { LOG.error("Cannot cast "+opt.getValue()+" to integer",err); return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;}
		</xsl:when>
		
		<xsl:when test="@type='output-file'">
		java.io.File <xsl:value-of select="generate-id()"/> =  null;
		try { <xsl:value-of select="generate-id()"/> = new java.io.File(opt.getValue());}
		catch(Exception err) { LOG.error("Cannot cast "+opt.getValue()+" to File",err); return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;}
		</xsl:when>
		<xsl:when test="@type='input-file' or @type='input-directory'">
		java.io.File <xsl:value-of select="generate-id()"/> =  null;
		try { <xsl:value-of select="generate-id()"/> = new java.io.File(opt.getValue());}
		catch(Exception err) { LOG.error("Cannot cast "+opt.getValue()+" to File",err); return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;}
		if(!<xsl:value-of select="generate-id()"/>.exists())
			{
			LOG.error("option -"+opt.getOpt()+": file "+<xsl:value-of select="generate-id()"/>+" doesn't exists");
			return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;
			}
		<xsl:choose>
			<xsl:when test="@type='input-directory'">
			if(!<xsl:value-of select="generate-id()"/>.isDirectory())
				{
				LOG.error("option -"+opt.getOpt()+": path "+<xsl:value-of select="generate-id()"/>+" is not a directory.");
				return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;
				}
			</xsl:when>
			<xsl:otherwise>
			if(!<xsl:value-of select="generate-id()"/>.isFile())
				{
				LOG.error("option -"+opt.getOpt()+": file "+<xsl:value-of select="generate-id()"/>+" is not a file.");
				return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;
				}
			if(!<xsl:value-of select="generate-id()"/>.canRead())
				{
				LOG.error("option -"+opt.getOpt()+": file "+<xsl:value-of select="generate-id()"/>+" is not readeable.");
				return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;
				}
			</xsl:otherwise>
		</xsl:choose>
		</xsl:when>
		

		<xsl:when test="@type='input-file-set'">
		final <xsl:apply-templates select="." mode="java-type"/> <xsl:text> </xsl:text> <xsl:value-of select="generate-id()"/> = new java.util.HashSet&lt;java.io.File&gt;();
		for(final String <xsl:value-of select="concat('s_',generate-id())"/>: opt.getValues())
			{
			java.io.File <xsl:value-of select="concat('f_',generate-id())"/> =  null;
			try { <xsl:value-of select="concat('f_',generate-id())"/> = new java.io.File( <xsl:value-of select="concat('s_',generate-id())"/>);}
			catch(Exception err) { LOG.error("Cannot cast "+ <xsl:value-of select="concat('s_',generate-id())"/>+" to File",err); return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;}
			if(!<xsl:value-of select="concat('f_',generate-id())"/>.exists())
				{
				LOG.error("option -"+opt.getOpt()+": file "+<xsl:value-of select="concat('f_',generate-id())"/>+" doesn't exists");
				return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;
				}
			if(!<xsl:value-of select="concat('f_',generate-id())"/>.isFile())
				{
				LOG.error("option -"+opt.getOpt()+": file "+<xsl:value-of select="generate-id()"/>+" is not a file.");
				return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;
				}
			if(!<xsl:value-of select="concat('f_',generate-id())"/>.canRead())
				{
				LOG.error("option -"+opt.getOpt()+": file "+<xsl:value-of select="concat('f_',generate-id())"/>+" is not readeable.");
				return com.github.lindenb.jvarkit.util.command.Command.Status.EXIT_FAILURE;
				}
			<xsl:value-of select="generate-id()"/>.add(<xsl:value-of select="concat('f_',generate-id())"/>);
			}
		</xsl:when>

		<xsl:when test="@type='string-list'">
		final <xsl:apply-templates select="." mode="java-type"/> <xsl:text> </xsl:text> <xsl:value-of select="generate-id()"/> = opt.getValuesList();
		</xsl:when>
		<xsl:when test="@type='string-set' or @type='uri-set'">
		final java.util.Set&lt;String&gt; <xsl:value-of select="generate-id()"/> = new java.util.HashSet&lt;String&gt;(opt.getValuesList());
		</xsl:when>



		<xsl:otherwise>
			<xsl:message terminate='yes'>visit: unknown type <xsl:value-of select="@type"/>.</xsl:message>
		</xsl:otherwise>
	</xsl:choose>
	this.<xsl:apply-templates select="." mode="setter"/>(<xsl:value-of select="generate-id()"/>);
	return com.github.lindenb.jvarkit.util.command.Command.Status.OK;
	}
</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:option" mode="height">
<xsl:choose>
	<xsl:otherwise>
		<xsl:text>DEFAULT_HEIGHT</xsl:text>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>


<xsl:template match="c:option" mode="swing-declare">
/** <xsl:value-of select="@name"/> */
<xsl:choose>
	<xsl:when test="@type='boolean'">
		private javax.swing.JCheckBox  <xsl:value-of select="generate-id(.)"/> = null;
	</xsl:when>
	<xsl:when test="(@type='string' and @multiline='true') or @type='string-set' or @type='string-list'">
		private javax.swing.JTextArea <xsl:value-of select="generate-id(.)"/> = null;
	</xsl:when>
	<xsl:when test="@type='string' or @type='int' or @type='long' or @type='double' or @type='float'">
		private javax.swing.JTextField <xsl:value-of select="generate-id(.)"/> = null;
	</xsl:when>
	<xsl:when test="@type='input-file' or @type='input-directory'">
		private com.github.lindenb.jvarkit.util.swing.InputChooser <xsl:value-of select="generate-id(.)"/> = null;
	</xsl:when>
	<xsl:when test="@type='output-file'">
		private com.github.lindenb.jvarkit.util.swing.OutputChooser <xsl:value-of select="generate-id(.)"/> = null;
	</xsl:when>
	<xsl:when test="@type='uri-set' or @type='input-file-set'">
		private  com.github.lindenb.jvarkit.util.swing.MultipleInputChooser <xsl:value-of select="generate-id(.)"/> = null;
	</xsl:when>
	<xsl:otherwise>
		<xsl:message terminate='yes'>swing-declare: unknown type <xsl:value-of select="@type"/>.</xsl:message>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:option" mode="swing-build">
/* BEGIN option <xsl:value-of select="@name"/> ****/
<xsl:choose>
	<xsl:when test="@type='boolean'">
	pane.add(new javax.swing.JLabel(""));
	</xsl:when>
	<xsl:otherwise>
		final javax.swing.JLabel <xsl:value-of select="generate-id(.)"/>lbl = new javax.swing.JLabel(makeLabel("<xsl:apply-templates select="." mode="label"/>")+":",javax.swing.JLabel.RIGHT);
		pane.add(<xsl:value-of select="generate-id(.)"/>lbl);
		<xsl:value-of select="generate-id(.)"/>lbl.setToolTipText("<xsl:apply-templates select="." mode="description"/>");
	</xsl:otherwise>
</xsl:choose>

<xsl:choose>
	<xsl:when test="(@type='string' and @multiline='true') or @type='string-list' or @type='string-set'">
		this.<xsl:value-of select="generate-id(.)"/> = new javax.swing.JTextArea(5,20);
		<xsl:if test="@default">
		this.<xsl:value-of select="generate-id(.)"/>.setText("<xsl:value-of select="@default"/>");
		</xsl:if>
		final javax.swing.JScrollPane <xsl:value-of select="generate-id(.)"/>scroll = new javax.swing.JScrollPane(this.<xsl:value-of select="generate-id(.)"/>);
		pane.add(<xsl:value-of select="generate-id(.)"/>scroll);
	</xsl:when>
	<xsl:when test="@type='string'">
		this.<xsl:value-of select="generate-id(.)"/> = new javax.swing.JTextField();
		<xsl:if test="@default">
		this.<xsl:value-of select="generate-id(.)"/>.setText("<xsl:value-of select="@default"/>");
		</xsl:if>
		pane.add(<xsl:value-of select="generate-id(.)"/>);

	</xsl:when>
	<xsl:when test="@type='input-file' or @type='input-directory'">
		this.<xsl:value-of select="generate-id(.)"/> = new com.github.lindenb.jvarkit.util.swing.InputChooser();
		<xsl:if test="@type='input-directory'">
		this.<xsl:value-of select="generate-id(.)"/>.setSelectType(com.github.lindenb.jvarkit.util.swing.InputChooser.SelectType.SELECT_DIRECTORY);
		</xsl:if>
		pane.add(<xsl:value-of select="generate-id(.)"/>);
	</xsl:when>
	<xsl:when test="@type='output-file' ">
		this.<xsl:value-of select="generate-id(.)"/> = new com.github.lindenb.jvarkit.util.swing.OutputChooser();
		pane.add(<xsl:value-of select="generate-id(.)"/>);
	</xsl:when>
	<xsl:when test="@type='int' or @type='long'  or @type='double'  or @type='float'">
		this.<xsl:value-of select="generate-id(.)"/> = new javax.swing.JTextField(20);
		this.<xsl:value-of select="generate-id(.)"/>.setText("<xsl:value-of select="@default"/>");
		pane.add(<xsl:value-of select="generate-id(.)"/>);
	</xsl:when>
	
	<xsl:when test="@type='boolean'">
		this.<xsl:value-of select="generate-id(.)"/> = new javax.swing.JCheckBox("<xsl:apply-templates select="." mode="label"/>");
		<xsl:choose>
			<xsl:when test="@default='true'">this.<xsl:value-of select="generate-id(.)"/>.setSelected(true);</xsl:when>
			<xsl:when test="@default='false'">this.<xsl:value-of select="generate-id(.)"/>.setSelected(false);</xsl:when>
			<xsl:otherwise>
				<xsl:message terminate='yes'>no default for '<xsl:value-of select="@name"/>'.</xsl:message>
			</xsl:otherwise>
		</xsl:choose>
		pane.add(<xsl:value-of select="generate-id(.)"/>);
	</xsl:when>
	<xsl:when test="@type='uri-set' or @type='input-file-set'">
		this.<xsl:value-of select="generate-id(.)"/> = new  com.github.lindenb.jvarkit.util.swing.MultipleInputChooser();
		pane.add(<xsl:value-of select="generate-id(.)"/>);
	 </xsl:when>
	<xsl:otherwise>
		<xsl:message terminate='yes'>swing-build: unknown type '<xsl:value-of select="@type"/>'.</xsl:message>
	</xsl:otherwise>
</xsl:choose>

<xsl:choose>
	<xsl:when test="@type='boolean'">
	</xsl:when>
	<xsl:otherwise>
		<xsl:value-of select="generate-id(.)"/>lbl.setLabelFor(this.<xsl:value-of select="generate-id(.)"/>);
	</xsl:otherwise>
</xsl:choose>

/* END option <xsl:value-of select="@name"/> ****/

</xsl:template>


<xsl:template match="c:option" mode="swing-validation">
<xsl:choose>
	<xsl:when test="@type='output-file' and @opt='o'">
	if( <xsl:value-of select="generate-id(.)"/>.isEmpty())
		{
		return "undefined output file";
		}
	</xsl:when>
	<xsl:when test="@type='string'">
	
	<xsl:if test="c:regex">
	final java.util.regex.Pattern  <xsl:value-of select="generate-id(.)"/>regex = java.util.regex.Pattern.compile("<xsl:value-of select="c:regex/text()"/>");
	if(! <xsl:value-of select="generate-id(.)"/>regex.matcher(<xsl:value-of select="generate-id(.)"/>.getText().trim()).matches())
		{
		this.<xsl:value-of select="generate-id(.)"/>.requestFocus();
		return " <xsl:value-of select="@name"/> doesn't match "+  <xsl:value-of select="generate-id(.)"/>regex.pattern();
		}
	</xsl:if>
	
	
	
	
	</xsl:when>
	<xsl:when test='@type="int"'>
		{
		try
			{
			Integer.parseInt(<xsl:value-of select="generate-id(.)"/>.getText());
			}
		catch(Exception err)
			{
			this.<xsl:value-of select="generate-id(.)"/>.requestFocus();
			return "Not an integer number: <xsl:value-of select="@name"/>";
			}
		}
	</xsl:when>
	
	<xsl:when test='@type="long"'>
		{
		try
			{
			Long.parseLong(<xsl:value-of select="generate-id(.)"/>.getText());
			}
		catch(Exception err)
			{
			this.<xsl:value-of select="generate-id(.)"/>.requestFocus();
			return "Not a Long number: <xsl:value-of select="@name"/>";
			}
		}
	</xsl:when>
	
	<xsl:when test='@type="double"'>
		{
		try
			{
			Double.parseDouble(<xsl:value-of select="generate-id(.)"/>.getText());
			}
		catch(Exception err)
			{
			this.<xsl:value-of select="generate-id(.)"/>.requestFocus();
			return "Not an double number: <xsl:value-of select="@name"/>";
			}
		}
	</xsl:when>
	
	<xsl:when test='@type="float"'>
		{
		try
			{
			Float.parseFloat(<xsl:value-of select="generate-id(.)"/>.getText());
			}
		catch(Exception err)
			{
			this.<xsl:value-of select="generate-id(.)"/>.requestFocus();
			return "Not an float number: <xsl:value-of select="@name"/>";
			}
		}
	</xsl:when>
	
	
	<xsl:when test="@type='input-file' or @type='input-directory'">
		<xsl:if test="@required='true'">
			if(<xsl:value-of select="generate-id(.)"/>.isEmpty())
				{
				return "Undefined input file for <xsl:value-of select="@name"/>";
				}
		</xsl:if>
	</xsl:when>
	<xsl:when test="@type='output-file'">
		<xsl:if test="@required='true'">
			if(<xsl:value-of select="generate-id(.)"/>.isEmpty())
				{
				return "Undefined output file for <xsl:value-of select="@name"/>";
				}
		</xsl:if>
	</xsl:when>
	<xsl:when test="@type='uri-set'">

	</xsl:when>
	<xsl:when test="@type='boolean'">

	</xsl:when>
	
	<xsl:when test="@type='string-set'">

	</xsl:when>
	
	<xsl:when test="@type='string-list'">

	</xsl:when>
	
	<xsl:when test="@type='input-file-set'">

	</xsl:when>
	
	
	<xsl:otherwise>
		<xsl:message terminate='yes'>swing-validationvalidation: unknown type <xsl:value-of select="@type"/>.</xsl:message>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>


<xsl:template match="c:option" mode="swing-fill-command">
<xsl:choose>
	<xsl:when test="@type='string' or @type='int' or @type='long' or @type='double' or @type='float'">
	if( !this.<xsl:value-of select="generate-id(.)"/>.getText().trim().isEmpty())
		{
		command.add("-<xsl:value-of select="@opt"/>");
		command.add(this.<xsl:value-of select="generate-id(.)"/>.getText().trim());
		}
	</xsl:when>
	
	<xsl:when test="@type='input-file' or  @type='output-file' or  @type='input-directory'">
	if( !this.<xsl:value-of select="generate-id(.)"/>.isEmpty())
		{
		command.add("-<xsl:value-of select="@opt"/>");
		command.add(this.<xsl:value-of select="generate-id(.)"/>.getText().trim());
		}
	</xsl:when>
	
	<xsl:when test="@type='uri-set' or @type='input-file-set'">
	if(!this.<xsl:value-of select="generate-id(.)"/>.getAsList().isEmpty())
		{
		command.add("-<xsl:value-of select="@opt"/>");
		for(final String s: this.<xsl:value-of select="generate-id(.)"/>.getAsList())
			{
			command.add(s);
			}
		}
	</xsl:when>
	
	<xsl:when test="@type='boolean'">
	if( this.<xsl:value-of select="generate-id(.)"/>.isSelected())
		{
		command.add("-<xsl:value-of select="@opt"/>");
		}
	</xsl:when>
	
	<xsl:when test="@type='string-set'">
	
		{
		final java.util.Set&lt;String&gt; set = new java.util.LinkedHashSet&lt;String&gt;();
		for( String s: this.<xsl:value-of select="generate-id(.)"/>.getText().trim().split("[\n\r]"))
			{
			s=s.trim();
			if(s.isEmpty()) continue;
			set.add(s);
			}
		if( !set.isEmpty())
			{
			command.add("-<xsl:value-of select="@opt"/>");
			for(String s:set)
				{
				command.add(s);
				}
			}
		}
	</xsl:when>
	
	
	<xsl:when test="@type='string-list'">
	
		{
		final java.util.List&lt;String&gt; list = new java.util.ArrayList&lt;String&gt;();
		for( String s: this.<xsl:value-of select="generate-id(.)"/>.getText().trim().split("[\n\r]"))
			{
			s=s.trim();
			if(s.isEmpty()) continue;
			list.add(s);
			}
		if( !list.isEmpty())
			{
			command.add("-<xsl:value-of select="@opt"/>");
			for(String s:list)
				{
				command.add(s);
				}
			}
		}
	</xsl:when>
	
	<xsl:otherwise>
		<xsl:message terminate='yes'>swing-fill-command: unknown type <xsl:value-of select="@type"/>.</xsl:message>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:description[@id='faidx']">
<xsl:text>Indexed Reference genome</xsl:text>
</xsl:template>


<xsl:template match="c:description">
<xsl:value-of select="text()"/>
</xsl:template>


<xsl:template match="c:app" mode="online-urls">

	@Override
	public String getOnlineSrcUrl()
		{
		return "https://github.com/lindenb/jvarkit/blob/master/src/main/java/<xsl:value-of select="translate(@package,'.','/')"/>/<xsl:value-of select="@app"/>.java";
		}
	
	@Override
	public String getOnlineDocUrl()
		{
		return "<xsl:apply-templates select="." mode="wikiurl"/>";
		}
	
</xsl:template>

<xsl:template match="c:app" mode="wikiurl">
<xsl:text>https://github.com/lindenb/jvarkit/wiki/</xsl:text>
<xsl:value-of select="@app"/>
</xsl:template>

<xsl:template name="titleize">
<xsl:param name="s"/>
<xsl:value-of select="concat(translate(substring($s,1,1),'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'),substring($s,2))"/>
</xsl:template>



<xsl:template match="c:history">
* xsl TODO
</xsl:template>



<xsl:template name="string-replace">
    <xsl:param name="string"/>
    <xsl:param name="from"/>
    <xsl:param name="to"/>
    
    <xsl:if test="string-length($from)=0"><xsl:message terminate="yes">BOUM:'<xsl:value-of select="$from"/>' vs '<xsl:value-of select="$to"/>'</xsl:message></xsl:if>
    
    <xsl:choose>
      <xsl:when test="contains($string,$from)">
        <xsl:value-of select="substring-before($string,$from)"/>
        <xsl:value-of select="$to"/>
        <xsl:call-template name="string-replace">
          <xsl:with-param name="string" select="substring-after($string,$from)"/>
          <xsl:with-param name="from" select="$from"/>
          <xsl:with-param name="to" select="$to"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$string"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

</xsl:stylesheet>



