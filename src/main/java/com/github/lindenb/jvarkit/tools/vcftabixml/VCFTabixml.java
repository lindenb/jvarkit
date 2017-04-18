/*
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
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
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

@Program(name="vcftabixml",description=" annotate a value from a vcf+xml file")
public class VCFTabixml extends Launcher
	{
	
	private static final Logger LOG=Logger.build(VCFTabixml.class).make();
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
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
	protected int doVcfToVcf(String inputName, VcfIterator r, VariantContextWriter w) {
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
