<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	>
<xsl:output method="text"/>
<xsl:param name="githash">undefined</xsl:param>

<xsl:template match="/">

 <xsl:apply-templates select="c:app"/>
</xsl:template>

<xsl:template match="c:app">/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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


<xsl:apply-templates select="history"/>

*/
package <xsl:apply-templates select="." mode="package"/>;



@Generated("xslt")
public abstract class <xsl:apply-templates select="." mode="abstract-class-name"/>
	extends <xsl:choose>
		<xsl:when test="@extends"><xsl:value-of select="@extends"/></xsl:when>
		<xsl:otherwise>com.github.lindenb.jvarkit.util.AbstractCommandLineProgram</xsl:otherwise>
	</xsl:choose>
	{
	/** error stream */
	protected java.io.PrintStream _errStream = System.err;
	/** stdout stream */
	protected java.io.PrintStream _outStream = System.out;	
	/** stdin stream */
	protected java.io.InputStream _inStream  = System.in;	
	/** options */
	protected org.apache.commons.cli.Options options = new org.apache.commons.cli.Options();
	<xsl:for-each select="c:options/c:option" >
	/** option <xsl:apply-templates select="." mode="name"/> */
	protected org.apache.commons.cli.Option Option = org.apache.commons.cli.Option.Builder("<xsl:apply-templates select="." mode="name"/>")
		.desc("<xsl:apply-templates select="." mode="description"/>")
		.build();
	</xsl:for-each>
	
	<xsl:if test="not(@generate-constructor='false')">
	/** Constructor */
	protected <xsl:apply-templates select="." mode="abstract-class-name"/>()
		{
		}
	</xsl:if>
	
	
	
	public java.io.PrintStream getStderr()
		{
		return this._errStream;
		}
	
	public void setStderr(final java.io.PrintStream stream)
		{
		this._errStream = stream;
		}
	
	public java.io.PrintStream getStdout()
		{
		return this._outStream;
		}
	
	public void setStdout(final java.io.PrintStream stream)
		{
		this._outStream = stream;
		}
		

	public java.io.InputStream getStdin()
		{
		return this._inStream;
		}
	
	public void setStdin(final java.io.InputStream stream)
		{
		this._inStream = stream;
		}	
	
	
	/** return application Name */
	public String getName()
		{
		return "<xsl:apply-templates select="." mode="class-name"/>";
		}
	
	/** return application Description */
	public String getDescription()
		{
		return "<xsl:apply-templates select="c:description"/>";
		}
	
	
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"<xsl:apply-templates select="." mode="class-name"/>";
		}
	
	<xsl:apply-templates select="c:options" />
	
	
	
	protected int parseOptions(final java.util.List&lt;String&gt; args)
		{
		int optind=0;
		while(optind &lt; args.size() )
			{
			<xsl:for-each select="c:options/c:option" >
			<xsl:if test="position()&gt;1"> else </xsl:if>
				<xsl:apply-templates select="." mode="parse"/>
			</xsl:for-each>
			else if(args.get(optind).equals("--"))
				{
				args.remove(optind);
				break;
				}
			else if(args.get(optind).startsWith("-"))
				{
				error("Unknown Option "+args.get(optind));
				return -1;
				}
			else
				{
				++optind;
				continue;	
				}
			}
		return 0;
		}
	
	}
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

<xsl:template match="c:options">
	public void printOptions(PrintStream out)
		{
		<xsl:apply-templates select="c:options/c:optio" mode="print.options"/>
		}
</xsl:template>

<xsl:template match="c:option[@type='int' and not(@multiple='true')]">
/** option <xsl:apply-templates select="." mode="name"/> */
protected int <xsl:apply-templates select="." mode="name"/> = <xsl:choose>
		<xsl:when test="@default"><xsl:value-of select="@default"/></xsl:when>
		<xsl:otherwise>0</xsl:otherwise>
	</xsl:choose>;

public int <xsl:text> </xsl:text><xsl:apply-templates select="." mode="getter"/>()
	{
	return this.<xsl:apply-templates select="." mode="name"/>;
	}

public <xsl:apply-templates select="." mode="setter"/>( final int <xsl:apply-templates select="." mode="name"/>)
	{
	this.<xsl:apply-templates select="." mode="name"/> = <xsl:apply-templates select="." mode="name"/>;
	}
</xsl:template>

<xsl:template match="c:option[@type='outFile' and not(@multiple='true')]">
/** option <xsl:apply-templates select="." mode="name"/> */
protected java.io.File <xsl:apply-templates select="." mode="name"/> = <xsl:choose>
		<xsl:when test="@default">new java.io.File("<xsl:value-of select="@default"/>");</xsl:when>
		<xsl:otherwise>null</xsl:otherwise>
	</xsl:choose>;

public java.io.File <xsl:text> </xsl:text><xsl:apply-templates select="." mode="getter"/>()
	{
	return this.<xsl:apply-templates select="." mode="name"/>;
	}

public <xsl:apply-templates select="." mode="setter"/>( final java.io.File <xsl:apply-templates select="." mode="name"/>)
	{
	this.<xsl:apply-templates select="." mode="name"/> = <xsl:apply-templates select="." mode="name"/>;
	}
</xsl:template>


<xsl:template match="c:option[@type='version']" mode="parse">
<xsl:apply-templates select="." mode="parsearg"/>
	{
	
	return 0;
	}
</xsl:template>

<xsl:template match="c:option[@type='help']" mode="parse">
<xsl:apply-templates select="." mode="parsearg"/>
	{
	
	return 0;
	}
</xsl:template>

<xsl:template match="c:option[@type='int' and not(@multiple='true')]" mode="parse">
<xsl:apply-templates select="." mode="parsearg"/>
	{
	if( optind+1 &gt;= args.size())
		{
		error("Option \"<xsl:value-of select="@name"/>\" : argument missing.");
		return -1;
		}
	final String _v = args.get(optind + 1);
	try
		{
		final int _i = Integer.parseInt(_v);
		<xsl:if test="@min-inclusive">
		if( _i &lt; <xsl:value-of select="@min-inclusive"/>)
			{
			error("<xsl:value-of select="@name"/> should be greater of equal to <xsl:value-of select="@min-inclusive"/>");
			return -1; 
			}
		</xsl:if>
		<xsl:if test="@max-inclusive">
		if( _i &gt; <xsl:value-of select="@max-inclusive"/>)
			{
			error("<xsl:value-of select="@name"/> should be lower of equal to <xsl:value-of select="@max-inclusive"/>");
			return -1; 
			}
		</xsl:if>
		<xsl:if test="@min-exclusive">
		if( _i &lt;= <xsl:value-of select="@min-exclusive"/>)
			{
			error("<xsl:value-of select="@name"/> should be greater than <xsl:value-of select="@min-exclusive"/>");
			return -1; 
			}
		</xsl:if>
		<xsl:if test="@max-exclusive">
		if( _i &gt;= <xsl:value-of select="@max-exclusive"/>)
			{
			error("<xsl:value-of select="@name"/> should be lower than <xsl:value-of select="@max-exclusive"/>");
			return -1; 
			}
		</xsl:if>
		
		args.remove(optind+1);
		args.remove(optind);
		this.<xsl:apply-templates select="." mode="setter"/>(_i);
		}
	catch(NumberFormatException err)
		{
		error("Option \"\" : Cannot convert "+_v+" to an integer");
		return -1;
		}
	continue;
	}
</xsl:template>


<xsl:template match="c:option[@type='string' and not(@multiple='true')]" mode="parse">
<xsl:apply-templates select="." mode="parsearg"/>
	{
	if( optind+1 &gt;= args.size())
		{
		error("Option \"<xsl:value-of select="@name"/>\" : argument missing.");
		return -1;
		}
	final String s = args.get(optind + 1);
	try
		{
		<xsl:if test="@length">
		
		if( s.length() != <xsl:value-of select="@length"/>)
			{
			error("<xsl:value-of select="@name"/> should contain <xsl:value-of select="@length"/> characters");
			return -1; 
			}
		</xsl:if>
		
		<xsl:if test="@min-length">
		if( s.length() &lt; <xsl:value-of select="@min-length"/>)
			{
			error("<xsl:value-of select="@name"/> should longer than <xsl:value-of select="@min-length"/> characters");
			return -1; 
			}
		</xsl:if>
		
		<xsl:if test="@max-length">
		if( s.length() &gt; <xsl:value-of select="@max-length"/>)
			{
			error("<xsl:value-of select="@name"/> should longer than <xsl:value-of select="@max-length"/> characters");
			return -1; 
			}
		</xsl:if>
		
		args.remove(optind+1);
		args.remove(optind);
		this.<xsl:apply-templates select="." mode="setter"/>(_i);
		}
	catch(Exception err)
		{
		error("Option \"<xsl:value-of select="@name"/>\" : Error");
		return -1;
		}
	continue;
	}
</xsl:template>


<xsl:template match="c:option[@type='outFile' and not(@multiple='true')]" mode="parse">
<xsl:apply-templates select="." mode="parsearg"/>
	{
	if( optind+1 &gt;= args.size())
		{
		error("Option \"<xsl:value-of select="@name"/>\" : argument missing.");
		return -1;
		}
	final String _v = args.get(optind + 1);
	try
		{
		final java.io.File _f = new java.io.File(_v);
		args.remove(optind+1);
		args.remove(optind);
		this.<xsl:apply-templates select="." mode="setter"/>(_f);
		}
	catch(Exception err)
		{
		error("Option \"\" : Cannot convert "+_v+" to an file");
		return -1;
		}
	continue;
	}
</xsl:template>



<xsl:template match="c:option" mode="parsearg">
<xsl:text>if ( </xsl:text>
<xsl:if test="@opt">args.get(optind).equals("-<xsl:value-of select="@opt"/>")</xsl:if>
<xsl:if test="@longopt"><xsl:if test="@opt"> || </xsl:if> args.get(optind).equals("--<xsl:value-of select="@longopt"/>") </xsl:if>
<xsl:text>)
	</xsl:text>
</xsl:template>

<!-- option type not handled -->
<xsl:template match="c:option"  mode="parse">
<xsl:message terminate="yes"><xsl:apply-templates select="." mode="name"/> : Unknown Option type '<xsl:value-of select="@type"/>'</xsl:message>
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
out.println("-<xsl:value-of select="@optopt"/>");
</xsl:template>

<xsl:template match="c:history" >
<xsl:text>
History:
</xsl:text>
<xsl:apply-templates/>
</xsl:template>



</xsl:stylesheet>


