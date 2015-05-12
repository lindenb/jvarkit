<?xml version='1.0'  encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:date="http://exslt.org/dates-and-times" 
	version='1.0' 
	>
<xsl:param name="param">Main</xsl:param>

<xsl:template match="/">

import htsjdk.variant.vcf.*;
import htsjdk.variant.vcf.writer.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.tribble.readers.*;
import java.util.*;
import java.util.logging.*;


public class <xsl:value-of select="$mainClass"/>
	{
    private static final Logger LOG= Logger.getLogger("<xsl:value-of select="$mainClass"/>");
	
	private interface VariantFilter
		{
		public void initialize(VCFHeader header);
		public void dispose();
		}
	
	private class VariantFilter implements VariantFilter
		{
		
		public void initialize(VCFHeader header)
			{
			
			}
		public void dispose()
			{
			
			}
		}
	
	<xsl:apply-templates select="filter"/>
	
	
	private VariantFilter[] filters=new VariantFilter[]{
		
		};
	
	private void doWork(String args[])
		{
		VCFCodec codec= new VCFCodec();
		LineReader r= LineReaderUtil.fromBufferedStream(System.in);
		LineIteratorImpl t= new LineIteratorImpl(r);
		VCFHeader header = (VCFHeader)codec.readActualHeader(t);
		/* copy header */
		final VCFHeader header2 = new VCFHeader(header);
		
		VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
		vcwb.setOutputVCFStream(System.out);
		VariantContextWriter out =  VariantContextWriter build(); 
		out.writeHeader(header2);
		
		while(t.hasNext())
			{
			VariantContext ctx = codec.decode(t.next());
			
			 public void add(VariantContext vc);
			}
		r.close();
		}
	public static void main(String args[]) throws Exception
		{
		<xsl:value-of select="$mainClass"/> app = new <xsl:value-of select="$mainClass"/>();
		app.run();
		}
	}

</xsl:template>

<xsl:template match="x">
<xsl:for-each select="filter">
header.
</xsl:for-each>
</xsl:template>


<xsl:template match="filter">


header2.add(new VCFFilterHeaderLine <xsl:value-of select="generate-id()"/> = new
	VCFFilterHeaderLine(
		&quot;<xsl:value-of select="name"/>&quot;,
		&quot;<xsl:value-of select="description"/>&quot;
		));
</xsl:template>

<xsl:template match="info">
header2.add(new VCFFilterHeaderLine <xsl:value-of select="generate-id()"/> = new
	VCFFilterHeaderLine(
		&quot;<xsl:value-of select="name"/>&quot;,
		&quot;<xsl:value-of select="description"/>&quot;
		));
</xsl:template>



</xsl:stylesheet>

