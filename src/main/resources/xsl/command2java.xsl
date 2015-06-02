<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:c="https://github.com/lindenb/jvarkit/"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	>
<xsl:output method="text"/>

<xsl:template match="/">
 <xsl:apply-templates select="c:command"/>
</xsl:template>

<xsl:template match="c:command">/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package <xsl:apply-templates select="." mode="package"/>;

@Generated("xslt")
public abstract class <xsl:apply-templates select="." mode="abstract-class-name"/>
	extends <xsl:choose>
		<xsl:when test="@extends"><xsl:value-of select="@extends"/></xsl:when>
		<xsl:otherwise>com.github.lindenb.jvarkit.util.AbstractCommandLineProgram</xsl:otherwise>
	</xsl:choose>
	{
	protected <xsl:apply-templates select="." mode="abstract-class-name"/>()
		{
		}
	public String getName()
		{
		return "<xsl:apply-templates select="." mode="class-name"/>";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"<xsl:apply-templates select="." mode="class-name"/>";
		}
	
	<xsl:apply-templates select="c:option" />
	
	@Override
	public void printOptions(PrintStream out)
		{
		<xsl:apply-templates select="c:option" mode="print.options"/>
		}
	
	protected int parseOptions(String args[])
		{
		int optind=0;
		while(optind &lt; args.length)
			{
			<xsl:for-each select="c:option" >
			<xsl:if test="position()&gt;1"> else </xsl:if> if(args[optind].equals("-c"))
				{
				
				}
			</xsl:for-each>
			}
		return optind;
		}
	
	}
</xsl:template>

<xsl:template match="c:command" mode="package">
<xsl:value-of select="@package"/>
</xsl:template>

<xsl:template match="c:command" mode="class-name">
<xsl:value-of select="@class"/>
</xsl:template>

<xsl:template match="c:command" mode="abstract-class-name">
<xsl:text>Abstract</xsl:text>
<xsl:apply-templates select="." mode="class-name"/>
</xsl:template>

<xsl:template match="c:option">
/** */
private <xsl:apply-templates select="." mode="type"/><xsl:text> </xsl:text><xsl:apply-templates select="." mode="name"/> = <xsl:apply-templates select="." mode="default-value"/>;



public <xsl:apply-templates select="." mode="type"/><xsl:text> </xsl:text><xsl:apply-templates select="." mode="getter"/>()
	{
	return this.<xsl:apply-templates select="." mode="name"/>;
	}

public <xsl:apply-templates select="." mode="setter"/>( <xsl:apply-templates select="." mode="type"/><xsl:text> </xsl:text><xsl:apply-templates select="." mode="name"/>)
	{
	this.<xsl:apply-templates select="." mode="name"/> = <xsl:apply-templates select="." mode="name"/>;
	}
</xsl:template>


<xsl:template match="c:option" mode="default-value">
<xsl:choose>
  <xsl:when test="@type='int'">-1</xsl:when>
  <xsl:otherwise>null</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:option" mode="type">
<xsl:choose>
  <xsl:when test="@type='int'">int</xsl:when>
  <xsl:otherwise></xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="c:option" mode="name">
<xsl:value-of select="@name"/>
</xsl:template>

<xsl:template match="c:option" mode="getter">
<xsl:text>get</xsl:text>
<xsl:value-of select="@name"/>
</xsl:template>

<xsl:template match="c:option" mode="setter">
<xsl:text>set</xsl:text>
<xsl:value-of select="@name"/>
</xsl:template>

<xsl:template match="c:option" mode="print.options">
out.println("-<xsl:value-of select="@optopt">");
</xsl:template>

</xsl:stylesheet>


