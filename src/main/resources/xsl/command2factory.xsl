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

<xsl:template match="c:app">
<xsl:apply-templates select="." mode="header"/>

package <xsl:apply-templates select="." mode="package"/>;
import com.github.lindenb.jvarkit.util.command.Command;


@Generated("xslt")	
public class <xsl:apply-templates select="." mode="factory-class-name"/> extends
	<xsl:choose>
	<xsl:when test="@extends-factory"><xsl:value-of select="@extends-factory"/></xsl:when>
	<xsl:otherwise> com.github.lindenb.jvarkit.util.command.AbstractCommandFactory </xsl:otherwise>
	</xsl:choose>
	 {
	 private static final Log LOG=LogFactory.getLog(<xsl:apply-templates select="." mode="abstract-class-name"/>.class);
	 public <xsl:apply-templates select="." mode="factory-class-name"/>()
	 	{
	 	/* constructor */
	 	public <xsl:apply-templates select="." mode="factory-class-name"/>()
	 		{
	 		}
	 	
	 	
	 	/** return application Name */
		@Override
		public String getName()
			{
			return "<xsl:apply-templates select="." mode="class-name"/>";
			}
		
		public <xsl:apply-templates select="." mode="class-name"/> createCommand()
			{
			<xsl:apply-templates select="." mode="class-name"/> cmd = new <xsl:apply-templates select="." mode="class-name"/>();
			cmd.setFactory(this);
			return cmd;
			}

		@Override
		protected void fillOptions(final java.util.List&lt;org.apache.commons.cli.Option&gt; opts)
			{
			<xsl:apply-templates select="//c:option" mode="fill-options"/>
			super.fillOptions(opts);
			}
		
		<xsl:choose>
			<xsl:when test="c:label">
			/** get application label */
			@Override
			public String getLabel()
				{
				return "<xsl:apply-templates select="c:label"/>";
				}
			</xsl:when>
			<xsl:when test="@label">
			/** get application label */
			@Override
			public String getLabel()
				{
				return "<xsl:value-of select="@label"/>";
				}
			</xsl:when>
		</xsl:choose>
		
		<xsl:choose>
			<xsl:when test="c:description">
			/** get application description */
			@Override
			public String getDescription()
				{
				return "<xsl:apply-templates select="c:description"/>";
				}
			</xsl:when>
			<xsl:when test="@description">
			/** get application description */
			@Override
			public String getDescription()
				{
				return "<xsl:value-of select="@description"/>";
				}
			</xsl:when>
		</xsl:choose>
	 }
</xsl:template>



<xsl:template match="c:option" mode="fill-options">
<xsl:variable name="nargs">
</xsl:variable>

options.addOption(org.apache.commons.cli.Option
	<xsl:choose>
		<xsl:when test="@opt">
			.builder("<xsl:value-of select="@opt"/>")
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
	<xsl:if test="@value-separator">
	.valueSeparator('<xsl:value-of select="@value-separator"/>')
	</xsl:if>
	<xsl:if test="@description">
	.desc("<xsl:value-of select="@description"/>").
	</xsl:if>
	<xsl:if test="c:description">
	.desc("<xsl:value-of select="c:description"/>").
	</xsl:if>
		<xsl:choose>
		<xsl:when test="@arg-name">
			.argName("<xsl:value-of select="@arg-name"/>")
		</xsl:when>
		<xsl:otherwise>
			.argName("<xsl:value-of select="@name"/>")
		</xsl:otherwise>
	</xsl:choose>
	<xsl:choose>
		<xsl:when test="@type='bool' or @type='boolean'">
		.hasArg(false)
		</xsl:when>
		<xsl:when test="not(@type) or @type='string' or @type='char'">
		.type(PatternOptionBuilder.STRING_VALUE)
		</xsl:when>
		<xsl:when test="@type='int' or @type='long' or @type='short' or @type='double' or @type='float'  or @type='number'">
		.type(PatternOptionBuilder.NUMBER_VALUE)
		</xsl:when>
		<xsl:when test="@type='url'">
		.type(PatternOptionBuilder.STRING_VALUE)
		</xsl:when>
		<xsl:when test="@type='object'">
		.type(PatternOptionBuilder.OBJECT_VALUE)
		</xsl:when>
		<xsl:when test="@type='date'">
		.type(PatternOptionBuilder.DATE_VALUE)
		</xsl:when>
		<xsl:when test="@type='file'">
		.type(PatternOptionBuilder.FILE_VALUE)
		</xsl:when>
		<xsl:when test="@type='existing-file'">
		.type(PatternOptionBuilder.EXISTING_FILE_VALUE)
		</xsl:when>
		<xsl:otherwise>
			<xsl:message terminate="yes">unknown type <xsl:value-of select="@type"/></xsl:message>
		</xsl:otherwise>
	</xsl:choose>
	<xsl:choose>
		<xsl:when test="not(@count) and not(@type='bool' or @type='boolean')">
			.numberOfArgs(1)
		</xsl:when>
		<xsl:when test="@count='unlimited'">
			.numberOfArgs(Option.UNLIMITED_VALUES)
		</xsl:when>
		<xsl:otherwise>
			.numberOfArgs(<xsl:value-of select="@count"/>)
		</xsl:otherwise>
	</xsl:choose>
	build()
	);	

</xsl:template>



<xsl:template match="c:option">
<xsl:variable name="basetype">
	<xsl:choose>
		<xsl:when test="@type='bool' or @type='boolean'">bool</xsl:when>
		<xsl:when test="not(@type) or @type='string' or @type='char'">java.lang.String</xsl:when>
		<xsl:when test="@type='short'">java.lang.Short</xsl:when>
		<xsl:when test="@type='int'">java.lang.Integer</xsl:when>
		<xsl:when test="@type='long'">java.lang.Long</xsl:when>
		<xsl:when test="@type='float'">java.lang.Float</xsl:when>
		<xsl:when test="@type='double'">java.lang.Float</xsl:when>
		<xsl:when test="@type='int' or @type='long' or @type='short' or @type='double' or @type='float'  or @type='number'">
		.type(PatternOptionBuilder.NUMBER_VALUE)
		</xsl:when>
		<xsl:when test="@type='url'">
		.type(PatternOptionBuilder.STRING_VALUE)
		</xsl:when>
		<xsl:when test="@type='object'">
		.type(PatternOptionBuilder.OBJECT_VALUE)
		</xsl:when>
		<xsl:when test="@type='date'">
		.type(PatternOptionBuilder.DATE_VALUE)
		</xsl:when>
		<xsl:when test="@type='file'">
		.type(PatternOptionBuilder.FILE_VALUE)
		</xsl:when>
		<xsl:when test="@type='existing-file'">
		.type(PatternOptionBuilder.EXISTING_FILE_VALUE)
		</xsl:when>
		<xsl:otherwise>
			<xsl:message terminate="yes">unknown type <xsl:value-of select="@type"/></xsl:message>
		</xsl:otherwise>
	</xsl:choose>
</xsl:variable>


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


