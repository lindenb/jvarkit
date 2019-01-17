/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

BEGIN_DOC


4th column of the BED indexed with TABIX is a XML string.
It will be processed with the xslt-stylesheet and should procuce a xml result <properties><entry key='key1'>value1</property><property key='key2'>values1</property></properies>
INFO fields. Carriage returns will be removed." 
Parameters to be passed to the stylesheet: vcfchrom (string) vcfpos(int) vcfref(string) vcfalt(string). 

END_DOC
*/
package com.github.lindenb.jvarkit.tools.vcftabixml;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlValue;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Templates;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
/**
BEGIN_DOC

## Example

The header of the tabix-index BED+XML file:
```
1	69427	69428	<snpList><positionString>1:69428</positionString><chrPosition>69428</chrPosition><alleles>G/T</alleles><ua..
1	69475	69476	<snpList><positionString>1:69476</positionString><chrPosition>69476</chrPosition><alleles>C/T</alleles><ua..
1	69495	69496	<snpList><positionString>1:69496</positionString><chrPosition>69496</chrPosition><alleles>A/G</alleles><ua..
1	69510	69511	<snpList><positionString>1:69511</positionString><chrPosition>69511</chrPosition><alleles>G/A</alleles><ua...
1	69589	69590	<snpList><positionString>1:69590</positionString><chrPosition>69590</chrPosition><alleles>A/T</alleles><ua...
1	69593	69594	<snpList><positionString>1:69594</positionString><chrPosition>69594</chrPosition><alleles>C/T</alleles><ua..
1	69619	69620	<snpList><positionString>1:69620</positionString><chrPosition>69620</chrPosition><alleles>T/TA</alleles><u...
1	69744	69745	<snpList><positionString>1:69745</positionString><chrPosition>69745</chrPosition><alleles>CA/C</alleles><u..
1	69760	69761	<snpList><positionString>1:69761</positionString><chrPosition>69761</chrPosition><alleles>T/A</alleles><ua...
1	801942	801943	<snpList><positionString>1:801943</positionString><chrPosition>801943</chrPosition><alleles>T/C</alleles...
(...)
```
The content of the tags.txt file:
```
##INFO=<ID=EVS_AAMAF,Number=1,Type=String,Description="EVS_AAMAF from EVS">
##INFO=<ID=EVS_AVGSAMPLEREADDEPTH,Number=1,Type=String,Description="EVS_AVGSAMPLEREADDEPTH from EVS">
##INFO=<ID=EVS_CLINICALLINK,Number=1,Type=String,Description="EVS_CLINICALLINK from EVS">
##INFO=<ID=EVS_CONFLICTALT,Number=1,Type=String,Description="EVS_CONFLICTALT from EVS">
##INFO=<ID=EVS_CONSERVATIONSCOREGERP,Number=1,Type=String,Description="EVS_CONSERVATIONSCOREGERP from EVS">
##INFO=<ID=EVS_CONSERVATIONSCORE,Number=1,Type=String,Description="EVS_CONSERVATIONSCORE from EVS">
##INFO=<ID=EVS_GENELIST,Number=1,Type=String,Description="EVS_GENELIST from EVS">
##INFO=<ID=EVS_GWASPUBMEDIDS,Number=1,Type=String,Description="EVS_GWASPUBMEDIDS from EVS">
##INFO=<ID=EVS_ONEXOMECHIP,Number=1,Type=String,Description="EVS_ONEXOMECHIP from EVS">
##INFO=<ID=EVS_RSIDS,Number=1,Type=String,Description="EVS_RSIDS from EVS">
##INFO=<ID=EVS_TOTALMAF,Number=1,Type=String,Description="EVS_TOTALMAF from EVS">
##INFO=<ID=EVS_UAMAF,Number=1,Type=String,Description="EVS_UAMAF from EVS">
```
The XSLT file:

```xml
<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"> 
<xsl:output method="xml"/>
<xsl:param name="vcfchrom"/>
<xsl:param name="vcfpos"/>
<xsl:param name="vcfref"/>
<xsl:param name="vcfalt"/>
<xsl:template match="/">
<properties>
<xsl:apply-templates select="evsData|snpList"/>
</properties>
</xsl:template>
<xsl:template match="evsData">
<xsl:a* 
 * @author lindenb
 *pply-templates select="snpList"/>
</xsl:template>
<xsl:template match="snpList">
<xsl:choose>
 <xsl:when test="chromosome=$vcfchrom and chrPosition=$vcfpos and refAllele=$vcfref">
<xsl:apply-templates select="clinicalLink|rsIds|uaMAF|aaMAF|totalMAF|avgSampleReadDepth|geneList|conservationScore|conservationScoreGERP|gwasPubmedIds|onExomeChip|gwasPubmedIds"/>
  <xsl:if test="altAlleles!=$vcfalt">
  	<entry key="EVS_CONFLICTALT">
  		<xsl:value-of select="altAlleles"/>
  	</entry>		
  	</xsl:if>
  </xsl:when>
  <xsl:otherwise/>
</xsl:choose>
</xsl:template>
<xsl:template match="clinicalLink|rsIds|uaMAF|aaMAF|totalMAF|avgSampleReadDepth|geneList|conservationScore|conservationScoreGERP|onExomeChip|gwasPubmedIds">
<xsl:if test="string-length(normalize-space(text()))&gt;0">
<entry>
<xsl:attribute name="key"><xsl:value-of select="concat('EVS_',translate(name(.),'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'))"/></xsl:attribute>
<xsl:value-of select="text()"/>
</entry>
</xsl:if>
</xsl:template>
</xsl:stylesheet>
```

Running:

```bash
java -jar dist/vcftabixml.jar \
	IN=input.vcf.gz \
	XSL=evs2vcf.xsl \
	BEDFILE=evs.data.gz \
	TAGS=tags.txt | grep -v "##" | head

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	U0925
A	613326	.	G	A	182	.	AC1=1;AF1=0.5;DP=327;DP4=0,147,0,178;FQ=163;MQ=47;PV4=1,0.21,1.2e-35,1;RPB=6.999638e+00;VDB=8.105711e-05	GT:PL:DP:GQ	0/1:212,0,190:325:99
A	614565	.	A	G	222	.	AC1=2;AF1=1;DP=910;DP4=0,1,71,836;EVS_AAMAF=23.2659;EVS_AVGSAMPLEREADDEPTH=70;EVS_CLINICALLINK=unknown;EVS_CONSERVATIONSCORE=0.0;EVS_CONSERVATIONSCOREGERP=-3.5;EVS_GENELIST=KCNAB2;EVS_GWASPUBMEDIDS=unknown;EVS_ONEXOMECHIP=false;EVS_RSIDS=rs26;EVS_TOTALMAF=30.5519;EVS_UAMAF=33.7209;FQ=-282;MQ=41;PV4=1,0.26,0.032,1;RPB=1.705346e+00;VDB=1.656917e-35	GT:PL:DP:GQ	1/1:255,255,0:908:99
A	614379	.	C	T	225	.	AC1=1;AF1=0.5;DP=979;DP4=33,440,37,456;EVS_AAMAF=2.8179;EVS_AVGSAMPLEREADDEPTH=59;EVS_CLINICALLINK=unknown;EVS_CONSERVATIONSCORE=0.0;EVS_CONSERVATIONSCOREGERP=-4.7;EVS_GENELIST=KCNAB2;EVS_GWASPUBMEDIDS=unknown;EVS_ONEXOMECHIP=false;EVS_RSIDS=rs249;EVS_TOTALMAF=8.8261;EVS_UAMAF=11.4393;FQ=225;MQ=42;PV4=0.8,1,3.8e-152,1;RPB=1.317662e+01;VDB=3.857882e-21	GT:PL:DP:GQ	0/1:255,0,255:966:99
A	614565	.	T	G	222	.	AC1=2;AF1=1;DP=209;DP4=0,0,0,187;FQ=-282;MQ=46;VDB=2.410569e-12	GT:PL:DP:GQ	1/1:255,255,0:187:99
A	614953	.	C	T	225	.	AC1=1;AF1=0.5;DP=810;DP4=197,214,194,195;EVS_AAMAF=2.4603;EVS_AVGSAMPLEREADDEPTH=45;EVS_CLINICALLINK=unknown;EVS_CONSERVATIONSCORE=0.0;EVS_CONSERVATIONSCOREGERP=-3.1;EVS_GENELIST=KCNAB2;EVS_GWASPUBMEDIDS=unknown;EVS_ONEXOMECHIP=false;EVS_RSIDS=rs733;EVS_TOTALMAF=7.0506;EVS_UAMAF=9.4205;FQ=225;MQ=47;PV4=0.62,1,9.2e-58,0.015;RPB=4.560133e+00;VDB=6.880316e-03	GT:PL:DP:GQ	0/1:255,0,255:800:99
A	614922	.	G	A	130	.	AC1=1;AF1=0.5;DP=183;DP4=0,94,0,86;EVS_AAMAF=2.3053;EVS_AVGSAMPLEREADDEPTH=43;EVS_CLINICALLINK=unknown;EVS_CONSERVATIONSCORE=0.0;EVS_CONSERVATIONSCOREGERP=-3.0;EVS_GENELIST=KCNAB2;EVS_GWASPUBMEDIDS=unknown;EVS_ONEXOMECHIP=false;EVS_RSIDS=rs202;EVS_TOTALMAF=7.1215;EVS_UAMAF=9.6386;FQ=133;MQ=44;PV4=1,0.0065,3.7e-51,1;RPB=3.500959e+00;VDB=9.915205e-29	GT:PL:DP:GQ	0/1:160,0,184:180:99
A	614986	.	G	C	188	.	AC1=2;AF1=1;DP=176;DP4=0,0,0,175;FQ=-282;MQ=46;VDB=0.000000e+00	GT:PL:DP:GQ	1/1:221,255,0:175:99
A	615009	.	T	A	125	.	AC1=1;AF1=0.5;DP=103;DP4=45,0,56,0;FQ=120;MQ=45;PV4=1,0.14,2.8e-19,1;RPB=1.520268e+00;VDB=1.539079e-06	GT:PL:DP:GQ	0/1:155,0,148:101:99
A	615037	.	C	T	161	.	AC1=1;AF1=0.5;DP=353;DP4=0,164,0,165;FQ=110;MQ=48;PV4=1,1,1.1e-23,1;RPB=5.549816e+00;VDB=1.486773e-11	GT:PL:DP:GQ	0/1:191,0,137:329:99
```

END_DOC

 */
@Program(name="vcftabixml",
	description=" annotate a value from a vcf+xml file",
	keywords={"vcf","xml"})
public class VCFTabixml extends Launcher
	{
	
	private static final Logger LOG=Logger.build(VCFTabixml.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File fileout = null;

	
	@Parameter(names="-B",description=" BED file indexed with tabix. The 4th column *is* a XML string.)",required=true)
	private String BEDFILE=null;
	
	@Parameter(names="-xsl",description="x xslt-stylesheet. REQUIRED. Should produce a valid set of INFO field.",required=true)
	private File STYLESHEET=null;
	
	@Parameter(names="-F",description="file containing extra INFO headers line to add version: 4.1",required=true)
	public File TAGS=null;
	
	
	private Templates stylesheet=null;
	
	@XmlRootElement(name="entry")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class Property
		{
		@XmlAttribute(name="key")
		public String key;
		@XmlValue
		public String value;
		
		
		
		@Override
		public String toString() {
			return ""+key+"="+value+";";
			}
		}
	


	@XmlRootElement(name="properties")
	@XmlAccessorType(XmlAccessType.PROPERTY)
	public static class Properties
		{
		private List<Property> property=new ArrayList<Property>();
		
		@XmlElement(name="entry")
		public List<Property> getProperty() {
			return property;
			}
		public void setProperty(List<Property> property)
			{
			this.property = property;
			}
		@Override
		public String toString() {
			return property.toString();
			}
		}

	@Override
	protected int doVcfToVcf(String inputName, VCFIterator r, VariantContextWriter w) {
		TabixReader tabixReader =null;

		try {
			LOG.info("opening BED"+BEDFILE);

			tabixReader=new TabixReader(this.BEDFILE);
			
			Pattern tab=Pattern.compile("[\t]");
			LOG.info("loading xslt "+STYLESHEET);
			this.stylesheet=TransformerFactory.newInstance().newTemplates(new StreamSource(STYLESHEET));
			Transformer transformer=this.stylesheet.newTransformer();
			transformer.setOutputProperty(OutputKeys.METHOD,"xml");

			final VCFHeader header=r.getHeader();
			final VCFHeader h2=new VCFHeader(header);			
			
			LOG.info("reading Tags "+TAGS);
			BufferedReader rT=IOUtils.openFileForBufferedReading(TAGS);
			String line;
			while((line=rT.readLine())!=null)
				{
				if(!line.startsWith(VCFHeader.METADATA_INDICATOR))
					{
					throw new RuntimeException("should start with "+ VCFHeader.METADATA_INDICATOR +":"+line);
					}
				 if (!line.startsWith(VCFConstants.INFO_HEADER_START) )
				 	{
					throw new RuntimeException("should start with "+ VCFConstants.INFO_HEADER_START +":"+line);
				 	}
				VCFInfoHeaderLine hi=new VCFInfoHeaderLine(line.substring(7), VCFHeaderVersion.VCF4_1);
				if(hi.getCount()!=1)
					{
					throw new IllegalArgumentException("VCFHeaderLineCount not supported : "+hi);
					}
				switch(hi.getType())
					{
					case String:break;
					default: throw new IllegalArgumentException("VCFHeaderLineTyoe not supported : "+hi);
					}
				LOG.info(hi.toString());
				h2.addMetaDataLine(hi);
				}		
			rT.close();
			LOG.info("writing header");
			w.writeHeader(h2);
			
			
			
			JAXBContext jaxbCtx=JAXBContext.newInstance(Properties.class,Property.class);
			Unmarshaller unmarshaller=jaxbCtx.createUnmarshaller();
			
			while(r.hasNext())
				{
				VariantContext ctx=r.next();
				
				HashMap<String, Set<String>> insert=new LinkedHashMap<String,Set<String>>();
				int[] array = tabixReader.parseReg(ctx.getContig()+":"+(ctx.getStart())+"-"+(ctx.getEnd()+1));
				TabixReader.Iterator iter=null;
				
				if(array!=null && array.length==3 && array[0]!=-1 && array[1]>=0 && array[2]>=0)
					{
					iter=tabixReader.query(array[0],array[1],array[2]);
					}
				else
					{
					LOG.info("Cannot get "+ctx.getContig()+":"+(ctx.getStart())+"-"+(ctx.getEnd()+1));
					}
				String line2=null;
				
				while(iter!=null && (line2=iter.next())!=null)
					{
					
					String tokens2[]=tab.split(line2,5);
					
					if(tokens2.length<4)
						{
						LOG.error("[VCFTabixml] VCF. Error not enough columns in tabix.line "+line2);
						return -1;
						}
					
					int chromStart=Integer.parseInt(tokens2[1]);
					int chromEnd=Integer.parseInt(tokens2[2]);
					if(chromStart+1!=chromEnd)
						{
						LOG.error("Error in "+this.BEDFILE+" extected start+1=end int "+tokens2[0]+":"+tokens2[1]+"-"+tokens2[2]);
						continue;
						}
					
					

					if(ctx.getStart()-1!=chromStart) continue;

					
					transformer.setParameter("vcfchrom",ctx.getContig());
					transformer.setParameter("vcfpos",ctx.getStart());
					transformer.setParameter("vcfref",ctx.getReference().getBaseString());
					transformer.setParameter("vcfalt",ctx.getAltAlleleWithHighestAlleleCount().getBaseString());
					
					
					try {
						StringWriter sw=new StringWriter();
						StreamSource src=new StreamSource(new StringReader(tokens2[3]));
						StreamResult rez=new StreamResult(sw);
						transformer.transform(src, rez);
						Properties props=unmarshaller.unmarshal(new StreamSource(new StringReader(sw.toString())),Properties.class).getValue();
												
						for(Property p:props.getProperty())
							{
							
							if(p.key.isEmpty()) continue;
							if(h2.getInfoHeaderLine(p.key)==null) 
								{
								LOG.info("ignoring key "+p.key+" you could set it to:\n"+
										"##INFO=<ID="+p.key+",Number=1,Type=String,Description=\""+p.key+" from "+BEDFILE+"\">"
										);
								continue;
								}
							Set<String> x=insert.get(p.key);
							if(x==null)
								{
								x=new LinkedHashSet<String>();
								insert.put(p.key,x);
								}
							x.add(p.value);
							}
						}
					catch (Exception e)
						{
						e.printStackTrace();
						throw new RuntimeException("error",e);
						}
					
					}
				
				if(insert.isEmpty())
					{
					w.add(ctx);
					continue;
					}
				VariantContextBuilder b=new VariantContextBuilder(ctx);
				for(String key:insert.keySet())
					{
					for(String s2:insert.get(key))
						{
						b.attribute(key,s2);
						break;//limit 1
						}
					}
				w.add(b.make());
				}
			return 0;
			}
		catch (IOException err)
			{
			err.printStackTrace();
			return -1;
			}
		catch (Throwable err)
			{
			err.printStackTrace();
			return -1;
			}
		}
	

	@Override
	public int doWork(List<String> args) {
		return doVcfToVcf(args, this.fileout);
		}
	
	public static void main(String[] args) throws Exception
		{
		new VCFTabixml().instanceMainWithExit(args);
		}
	}
