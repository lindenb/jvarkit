<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	version='1.0'
	xmlns:c="http://github.com/lindenb/jvarkit/"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns:j="http://github.com/lindenb/jvarkit/"
	xmlns:h="http://www.w3.org/1999/xhtml"
	>

<xsl:output method="text"/>


<xsl:template match="/">

<xsl:apply-templates select="//j:documentation"/>

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="XXXXX",description="<xsl:value-of select="/j:app/j:description"/>")

private static final Logger LOG = Logger.build(XXXXXX.class).make();


@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
private File outputFile = null;

<xsl:apply-templates select="//j:option"/>

<xsl:if test="//j:snippet[@id='sorting-collection']">

	@Parameter(names={"-T","--tmpDir"},description="mp directory")
	private File tmpDir = new File(System.getProperty("java.io.tmpdir"));

	@Parameter(names={"-maxRecordsInRam","--maxRecordsInRam"},description="Max records in RAM")
	private int maxRecordsInRam =50000;
</xsl:if>


<xsl:for-each select="//j:snippet[@id='ref.faidx']">
@Parameter(names={"-<xsl:value-of select="@opt"/>","--reference"},description="Indexed fasta Reference")
private File <xsl:value-of select="@name"/> = null;
</xsl:for-each>

</xsl:template>


<xsl:template match="j:option[@type='boolean']">
@Parameter(names={"-<xsl:value-of select="@opt"/>","--<xsl:value-of select="@longopt"/>"},description="<xsl:value-of select="normalize-space(j:description)"/>")
private boolean <xsl:value-of select="@name"/> = false;
</xsl:template>


<xsl:template match="j:option[@type='int']">
@Parameter(names={"-<xsl:value-of select="@opt"/>","--<xsl:value-of select="@longopt"/>"},description="<xsl:value-of select="normalize-space(j:description)"/>")
private int <xsl:value-of select="@name"/> = <xsl:value-of select="@default"/> ;
</xsl:template>

<xsl:template match="j:option[@type='long']">
@Parameter(names={"-<xsl:value-of select="@opt"/>","--<xsl:value-of select="@longopt"/>"},description="<xsl:value-of select="normalize-space(j:description)"/>")
private long <xsl:value-of select="@name"/> = <xsl:value-of select="@default"/> ;
</xsl:template>


<xsl:template match="j:option[@type='double']">
@Parameter(names={"-<xsl:value-of select="@opt"/>","--<xsl:value-of select="@longopt"/>"},description="<xsl:value-of select="normalize-space(j:description)"/>")
private double <xsl:value-of select="@name"/> = <xsl:value-of select="@default"/> ;
</xsl:template>


<xsl:template match="j:option[@type='string']">
@Parameter(names={"-<xsl:value-of select="@opt"/>","--<xsl:value-of select="@longopt"/>"},description="<xsl:value-of select="normalize-space(j:description)"/>")
private String <xsl:value-of select="@name"/> = <xsl:choose>
		<xsl:when test="@default">"<xsl:value-of select="@default"/>"</xsl:when>
		<xsl:otherwise>null</xsl:otherwise>
	</xsl:choose>;
</xsl:template>

<xsl:template match="j:option[@type='string-list']">
@Parameter(names={"-<xsl:value-of select="@opt"/>","--<xsl:value-of select="@longopt"/>"},description="<xsl:value-of select="normalize-space(j:description)"/>")
private List&lt;String&gt; <xsl:value-of select="@name"/> = new ArrayList&lt;&gt;();
</xsl:template>


<xsl:template match="j:option[@type='string-set']">
@Parameter(names={"-<xsl:value-of select="@opt"/>","--<xsl:value-of select="@longopt"/>"},description="<xsl:value-of select="normalize-space(j:description)"/>")
private Set&lt;String&gt; <xsl:value-of select="@name"/> = new HashSet&lt;&gt;();

</xsl:template>


<xsl:template match="j:option[@type='output-file' or @type='input-file']">
@Parameter(names={"-<xsl:value-of select="@opt"/>","--<xsl:value-of select="@longopt"/>"},description="<xsl:value-of select="normalize-space(j:description)"/>")
private File <xsl:value-of select="@name"/> = <xsl:choose>
		<xsl:when test="@default">"<xsl:value-of select="@default"/></xsl:when>
		<xsl:otherwise>null</xsl:otherwise>
	</xsl:choose>;
</xsl:template>



<xsl:template match="j:option">
<xsl:message >###############BOUM <xsl:value-of select="@name"/>:<xsl:value-of select="@type"/></xsl:message>
</xsl:template>

<xsl:template match="j:documentation">
/**

BEGIN_DOC

<xsl:apply-templates/>

END_DOC
*/
</xsl:template>


<xsl:template match="h:pre">

```
<xsl:apply-templates/>
```

</xsl:template>

<xsl:template match="h:ul">
<xsl:text>
</xsl:text>
<xsl:apply-templates select="h:li"/>
<xsl:text>
</xsl:text>
</xsl:template>


<xsl:template match="h:li">
<xsl:text> *  </xsl:text>
<xsl:apply-templates/>
<xsl:text>
</xsl:text>
</xsl:template>


<xsl:template match="h:img">![img](<xsl:value-of select="@src"/>)</xsl:template>



<xsl:template match="h:h3">

### <xsl:apply-templates/><xsl:text>
</xsl:text>

</xsl:template>

<xsl:template match="h:h4">

#### <xsl:apply-templates/><xsl:text>
</xsl:text>

</xsl:template>


</xsl:stylesheet>
