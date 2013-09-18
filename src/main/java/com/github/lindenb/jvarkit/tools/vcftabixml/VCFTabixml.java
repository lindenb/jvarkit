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
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlValue;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Templates;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.picard.vcf.VcfIterator;

import org.broad.tribble.readers.TabixReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderVersion;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.picard.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;


public class VCFTabixml extends AbstractVCFFilter
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" annotate a value from a vcf+xml file."+
		("4th column of the BED indexed with TABIX is a XML string." +
		"It will be processed with the xslt-stylesheet and should procuce a xml result <properties><property key='key1'>value1</property><property key='key2'>values1</property></properies>" +
		" INFO fields. Carriage returns will be removed." +
		"Parameters to be passed to the stylesheet: vcfchrom (string) vcfpos(int) vcfref(string) vcfalt(string). "
		);
	
	private static final Log LOG=Log.getInstance(VCFTabixml.class);

	
	@Option(shortName="BED",doc=" BED file indexed with tabix. The 4th column *is* a XML string.)",optional=false)
	public String BEDFILE=null;
	
	@Option(shortName="XSL",doc="x xslt-stylesheet. REQUIRED. Should produce a valid set of INFO field.",optional=false)
	public File STYLESHEET=null;
	
	@Option(shortName="F",doc="file containing extra INFO headers line to add version: 4.1",optional=false)
	public File TAGS=null;
	
	
	private Templates stylesheet=null;
	
	@XmlRootElement(name="property")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class Property
		{
		@XmlAttribute(name="key")
		public String key;
		@XmlValue
		public String value;
		}
	


	@XmlRootElement(name="properties")
	@XmlAccessorType(XmlAccessType.PROPERTY)
	public static class Properties
		{
		private List<Property> property=new ArrayList<Property>();
		public List<Property> getProperty() {
			return property;
			}
		public void setProperty(List<Property> property)
			{
			this.property = property;
			}
		}

	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		TabixReader tabixReader =null;

		try {
			LOG.info("opening BED"+BEDFILE);

			tabixReader=new TabixReader(this.BEDFILE);
			
			Pattern tab=Pattern.compile("[\t]");
			LOG.info("loading xslt "+STYLESHEET);
			this.stylesheet=TransformerFactory.newInstance().newTemplates(new StreamSource(STYLESHEET));
			Transformer transformer=this.stylesheet.newTransformer();
			transformer.setOutputProperty(OutputKeys.METHOD,"xml");

			VCFHeader header=r.getHeader();
			VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());			
			
			LOG.info("reading Tags "+TAGS);
			BufferedReader rT=IOUtils.openFileForBufferedReading(TAGS);
			String line;
			while((line=rT.readLine())!=null)
				{
				if(!line.startsWith(VCFHeader.METADATA_INDICATOR))
					{
					throw new PicardException("should start with "+ VCFHeader.METADATA_INDICATOR +":"+line);
					}
				 if (!line.startsWith(VCFConstants.INFO_HEADER_START) )
				 	{
					throw new PicardException("should start with "+ VCFConstants.INFO_HEADER_START +":"+line);
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
				
				TabixReader.Iterator iter= tabixReader.query(
						ctx.getChr()+":"+(ctx.getStart())+"-"+(ctx.getEnd()+1));
				String line2=null;
				
				while(iter!=null && (line2=iter.next())!=null)
					{

					String tokens2[]=tab.split(line2,5);
					
					if(tokens2.length<4)
						{
						LOG.error("[VCFTabixml] VCF. Error not enough columns in tabix.line "+line2);
						return;
						}
					
					int chromStart=Integer.parseInt(tokens2[1]);
					int chromEnd=Integer.parseInt(tokens2[2]);
					if(chromStart+1!=chromEnd)
						{
						LOG.error("Error in "+this.BEDFILE+" extected start+1=end int "+tokens2[0]+":"+tokens2[1]+"-"+tokens2[2]);
						continue;
						}
					
					

					if(ctx.getStart()-1!=chromStart) continue;

					
					transformer.setParameter("vcfchrom",ctx.getChr());
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
						throw new PicardException("error",e);
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
			}
		catch (IOException err)
			{
			throw err;
			}
		catch (Exception err)
			{
			throw new IOException(err);
			}
		}
	

	public static void main(String[] args) throws Exception
		{
		new VCFTabixml().instanceMainWithExit(args);
		}
	}
