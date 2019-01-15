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


*/
package com.github.lindenb.jvarkit.util.vcf;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;


/**
 * XML VCF writer
 * @author lindenb
 *
 */
public class XMLVcfWriterFactory 
	{
	
	private static class XMLVcfWriter implements VariantContextWriter
		{
		private XMLStreamWriter writer=null;
		private OutputStream delegateOut=null;
		private VCFHeader  header=null;
		private Map<String,XMLInfoHandler> info2handler= new HashMap<String,XMLInfoHandler>();
		private Map<String,XMLFormatHandler> format2handler= new HashMap<String,XMLFormatHandler>();


		@Override
		public boolean checkError()
			{
			if( delegateOut==null) return false;
			if( !(delegateOut instanceof java.io.PrintStream) ) return false;
			return  java.io.PrintStream.class.cast(delegateOut).checkError();
			}

		
		private XMLVcfWriter()
			{
			info2handler.put("DP4",new DP4Handler());	
			info2handler.put("PV4",new PV4Handler());	
			}
		
				
		@SuppressWarnings("unused")
		public void putInfoHandler(XMLInfoHandler handler)
			{
			this.info2handler.put(handler.getKey(), handler);
			}
		
		
		protected void start(String tag) throws XMLStreamException
			{
			this.writer.writeStartElement(tag);
			}
		
		protected void attribute(String tag,Object o) throws XMLStreamException
			{
			this.writer.writeAttribute(tag, String.valueOf(o));
			}
		
		protected void end() throws XMLStreamException
			{
			this.writer.writeEndElement();
			}
		
		protected void element(String tag,Object content) throws XMLStreamException
			{
			if(content==null)
				{
				this.writer.writeEmptyElement(tag);
				return;
				}
			start(tag);
			characters(content);
			end();
			}
	
		protected void characters(Object content)
				throws XMLStreamException
			{
			if(content==null) return;
			this.writer.writeCharacters(String.valueOf(content));
			}
		
		@Override
		public void setHeader(final VCFHeader header) {
			throw new UnsupportedOperationException("setHeader shouldn't be called"); 
			}
		
		@Override
		public void writeHeader(VCFHeader header)
			{
			if(this.header!=null) throw new RuntimeException("Header was already written");
			this.header=header;
	
			try
				{
				start("vcf");
				start("header");
				if(header.getInfoHeaderLines()!=null)
					{
					start("infos");
					for(final VCFInfoHeaderLine h:header.getInfoHeaderLines())
						{
						start("info");
						attribute("key",h.getID());
						attribute("countType",h.getCountType());
						if(h.getCountType()==VCFHeaderLineCount.INTEGER)
							{
							attribute("count",h.getCount());
							}
						if(h.getValue()!=null && !h.getValue().isEmpty())
							attribute("value",h.getValue());
						characters(h.getDescription());
						end();
						if(!info2handler.containsKey(h.getID()))
							{
							XMLInfoHandler handler= new DefaultXMLInfoHandler(h);
							info2handler.put(h.getID(),handler);
							}
						}
					end();
					}
				
				if(header.getFormatHeaderLines()!=null)
					{
					start("formats");
					for(final VCFFormatHeaderLine h: header.getFormatHeaderLines())
						{
						start("format");
						attribute("key",h.getID());
						attribute("countType",h.getCountType());
						if(h.getCountType()==VCFHeaderLineCount.INTEGER)
							{
							attribute("count",h.getCount());
							}
						if(h.getValue()!=null && !h.getValue().isEmpty())
							attribute("value",h.getValue());
						characters(h.getDescription());
						end();
						
						if(!format2handler.containsKey(h.getID()))
							{
							XMLFormatHandler handler= new DefaultXMLFormatHandler(h);
							format2handler.put(h.getID(),handler);
							}
						
						}
					end();
					}
				
				if(header.getFilterLines()!=null)
					{
					start("filters");
					for(VCFFilterHeaderLine h: header.getFilterLines())
						{
						start("filter");
						element("key",h.getID());
						characters(h.getValue());
						end();
						}
					end();
					}
				
				if(header.getContigLines()!=null)
					{
					start("contigs");
					for(VCFContigHeaderLine h:header.getContigLines())
						{
						start("contig");
						attribute("tid", h.getContigIndex());
						element("key",h.getID());
						characters(h.getValue());
						end();
						}
					end();
					}
				
				if(header.getSampleNamesInOrder()!=null)
					{
					start("samples");
					for(String name:header.getSampleNamesInOrder())
						{
						start("sample");
						attribute("index", header.getSampleNameToOffset().get(name));
						characters(name);
						end();
						}
					end();
					}
				if(header.getMetaDataInInputOrder()!=null)
					{
					start("metas");
					for(VCFHeaderLine meta:header.getMetaDataInInputOrder())
						{
						if(meta.getKey().equals("INFO"))continue;
						if(meta.getKey().equals("FORMAT"))continue;
						if(meta.getKey().equals("contig"))continue;
						if(meta.getKey().equals("FILTER"))continue;
						if(meta.getKey().equals("fileformat"))continue;
						start("meta");
						attribute("key", meta.getKey());
						characters(meta.getValue());
						end();
						}
					end();
					}
				end();//header
				start("variations");
				characters("\n");
				}
			catch (XMLStreamException e)
				{
				e.printStackTrace();
				throw new RuntimeException(String.valueOf(e.getMessage()),e);
				}
			}
		
	    // should we write genotypes or just sites?
		private boolean doNotWriteGenotypes=false;
		
		@Override
		public void add(VariantContext variant)
			{
			if(this.header==null) throw new RuntimeException("No header was written.");
	
			try
				{
				 if ( doNotWriteGenotypes )
					 variant = new VariantContextBuilder(variant).noGenotypes().make();
							 
				start("variation");
				element("chrom",variant.getContig());
				element("start",variant.getStart());
				element("end",variant.getEnd());
				if(variant.hasID())
					{
					element("id",variant.getID());	
					}
				element("ref",variant.getReference().getDisplayString());
				
				if ( variant.isVariant() )
					{
					for(Allele a:variant.getAlternateAlleles())
						{
						element("alt",a.getDisplayString());
						}
					}
				
				if(variant.hasLog10PError())
					{
					element("qual", variant.getPhredScaledQual());
					}
				
				if(variant.isFiltered() || variant.filtersWereApplied())
					{
					start("filters");
					if(variant.isFiltered())
						{
						for(String s: variant.getFilters())
							{
							element("filter",s);
							}
						}
					else if(variant.filtersWereApplied())
						{
						element("filter",VCFConstants.PASSES_FILTERS_v4);
						}
					end();
					}
					
				if(variant.getAttributes()!=null)
					{
					start("infos");
					Map<String,Object> infos=variant.getAttributes();
					for(String key:infos.keySet())
						{
						XMLInfoHandler infoHandler=this.info2handler.get(key);
						if(infoHandler==null) continue;
						infoHandler.handle(this.header,this.writer, variant);
						}
					end();
					}
				
				
				
				if(variant.hasGenotypes())
					{
					
					start("genotypes");
					for(String sample:variant.getSampleNames())
						{
						Genotype g=variant.getGenotype(sample);
						if(g==null) continue;
						start("genotype");
						attribute("available",g.isAvailable());
						attribute("called",g.isCalled());
						attribute("het",g.isHet());
						attribute("hom",g.isHom());
						attribute("homRef",g.isHomRef());
						attribute("homVar",g.isHomVar());
						attribute("mixed",g.isMixed());
						attribute("noCall",g.isNoCall());
						attribute("nonInformative",g.isNonInformative());
						attribute("filtered",g.isFiltered());
						attribute("phased",g.isPhased());
						attribute("sample",g.getSampleName());
						if(g.hasAD())
							{
							start("AD");
							for(int ad:g.getAD())
								{
								element("value", ad);
								}
							end();
							}
						if(g.hasDP())
							{
							element("DP", g.getDP());
							}
						if(g.hasGQ())
							{
							element("GQ", g.getGQ());
							}
						if(g.hasPL())
							{
							start("PL");
							int index=0;
							for(int v:g.getPL())
								{
								start("value");
								attribute("index", ++index);
								characters(v);
								end();
								}
							end();
							}
						
						
						
						start("alleles");
						for(Allele a:g.getAlleles())
							{
							if(a.isNoCall()) continue;
							if(a.getBaseString().isEmpty()) continue;
							if(a.getBaseString().equals(".")) continue;
							start("allele");
							if(a.isReference()) attribute("ref", a.isReference());
							if(a.isSymbolic()) attribute("symbolic",true);
							characters(a.getBaseString());
							end();
							}
						end();
						
						
						Map<String,Object> xatt=g.getExtendedAttributes();
						if(xatt!=null)
							{
							for(String key:xatt.keySet())
								{
								XMLFormatHandler fmtHandler=this.format2handler.get(key);
								if(fmtHandler==null) continue;
								fmtHandler.handle(this.writer, variant,g);
								}
							}
						
						end();
						}
					end();
					}
				end();//variation
				characters("\n");
				}
			catch (XMLStreamException e)
				{
				throw new RuntimeException(e);
				}
			}
		
		@Override
		public void close()
			{
			if(this.writer==null) return;
			if(this.header==null) throw new RuntimeException("No header was written.");
			try {
				end();//variations
				end();//vcf
				writer.close();
				if(this.delegateOut!=null) {delegateOut.flush(); delegateOut.close();}
				this.writer=null;
				} 
			catch (Exception e)
				{
				e.printStackTrace();
				throw new RuntimeException("close failed",e);
				}
			}
	
		
	    public static interface XMLFormatHandler
			{
			public String getKey();
			public void handle(
					XMLStreamWriter w,
					final  VariantContext ctx,
					final Genotype g
					) throws XMLStreamException;
			}
	    
		public static interface XMLInfoHandler
			{
			public String getKey();
			public void handle(
					VCFHeader header,
					XMLStreamWriter w,
					final  VariantContext ctx
					) throws XMLStreamException;
			}
	
	public static abstract class AbstractXMLInfoHandler
		implements XMLInfoHandler
		{
		
		protected void handleObject(VCFHeader header,XMLStreamWriter w,int index,Object o)
				throws XMLStreamException
			{
			w.writeStartElement(getKey());
			if(index>=0) w.writeAttribute("index", String.valueOf(index));
			w.writeCharacters(String.valueOf(o));
			w.writeEndElement();
			}
		
		@SuppressWarnings("rawtypes")
		protected void handleArray(VCFHeader header,XMLStreamWriter w,Collection array)
				throws XMLStreamException
			{
			int index=0;
			for(Object o2:array) handleObject(header,w,++index,o2);
			}
	
		
		@SuppressWarnings("rawtypes")
		@Override
		public void handle(
				VCFHeader header,
				XMLStreamWriter w,
				final  VariantContext ctx
				)
				throws XMLStreamException
			{
			Object o=ctx.getAttribute(this.getKey());
			if(o==null) return;
			if(o.getClass().isArray())
				{
				Object array[]=(Object[])o;
				if(array.length==0) return;
				handleArray(header,w,Arrays.asList(array));
				}
			else if(o instanceof Collection)
				{
				Collection array=(Collection)o;
				if(array.isEmpty()) return;
				handleArray(header,w,array);
				}
			else
				{
				handleObject(header,w,-1,o);
				}
			}
		}
	
	
	public static class DefaultXMLInfoHandler extends AbstractXMLInfoHandler
		{
		private VCFInfoHeaderLine vihl;
		public DefaultXMLInfoHandler(VCFInfoHeaderLine vihl)
			{
			this.vihl=vihl;
			}
		
		@Override
		public String getKey()
			{
			return vihl.getID();
			}
		}
	
	    
	public static abstract class AbstractXMLFormatHandler
	implements XMLFormatHandler
		{
		
		protected void handleObject(XMLStreamWriter w,Object o)
				throws XMLStreamException
			{
			w.writeStartElement(getKey());
			w.writeCharacters(String.valueOf(o));
			w.writeEndElement();
			}
		
		
		@SuppressWarnings("rawtypes")
		@Override
		public void handle(
				XMLStreamWriter w,
				final  VariantContext ctx,
				final Genotype g
				)
				throws XMLStreamException
			{
			Object o=g.getExtendedAttribute(this.getKey());
			if(o==null) return;
			if(o.getClass().isArray())
				{
				Object array[]=(Object[])o;
				if(array.length==0) return;
				for(Object o2:array) handleObject(w,o2);
				}
			else if(o instanceof Collection)
				{
				Collection array=(Collection)o;
				if(array.isEmpty()) return;
				for(Object o2:array) handleObject(w,o2);
				}
			else
				{
				handleObject(w,o);
				}
			}
	}
	
	
	public static class DefaultXMLFormatHandler extends AbstractXMLFormatHandler
		{
		private VCFFormatHeaderLine vfhl;
		public DefaultXMLFormatHandler(VCFFormatHeaderLine vfhl)
			{
			this.vfhl=vfhl;
			}
		
		@Override
		public String getKey()
			{
			return vfhl.getID();
			}
		}
	
		
		private static class DP4Handler extends AbstractXMLInfoHandler
			{
			@Override
			protected void handleObject(VCFHeader header,
						XMLStreamWriter w, int index, Object o)
						throws XMLStreamException {
				}
			@SuppressWarnings("rawtypes")
			@Override
			protected void handleArray(VCFHeader header, XMLStreamWriter w,
						Collection array) throws XMLStreamException {
				if(array.size()!=4) return;
				w.writeStartElement(getKey());
				int i=0;
				for(Object o:array)
					{
					switch(i)
						{
						case 0: w.writeStartElement("ref-forward"); break;
						case 1: w.writeStartElement("ref-reverse"); break;
						case 2: w.writeStartElement("alt-forward"); break;
						case 3: w.writeStartElement("alt-reverse"); break;
						default: throw new XMLStreamException("bad index");
						}
					
					w.writeCharacters(String.valueOf(o));
					w.writeEndElement();
					i++;
					}
				w.writeEndElement();
				}
			@Override
			public String getKey()
				{
				return "DP4";
				}
			}
	
		private static class PV4Handler extends AbstractXMLInfoHandler
			{
			@Override
			protected void handleObject(VCFHeader header,
						XMLStreamWriter w, int index, Object o)
						throws XMLStreamException {
				}
			@Override
			protected void handleArray(VCFHeader header, XMLStreamWriter w,
						@SuppressWarnings("rawtypes") Collection array) throws XMLStreamException {
				if(array.size()!=4) return;
				w.writeStartElement(getKey());
				int i=0;
				for(Object o:array)
					{
					switch(i)
						{
						case 0: w.writeStartElement("strand-bias"); break;
						case 1: w.writeStartElement("baseQ-bias"); break;
						case 2: w.writeStartElement("mapQ-bias"); break;
						case 3: w.writeStartElement("tail-distance-bias"); break;
						default: throw new XMLStreamException("bad index");
						}
					
					w.writeCharacters(String.valueOf(o));
					w.writeEndElement();
					i++;
					}
				w.writeEndElement();
				}
			@Override
			public String getKey()
				{
				return "PV4";
				}
			}
		}
	
	private File outputFile=null;
	private XMLVcfWriterFactory()
		{
		
		}
	
	public void setOutputFile(File out)
		{
		this.outputFile=out;
		}
	
	
	public static XMLVcfWriterFactory newInstance()
		{
		return new XMLVcfWriterFactory();
		}	
	
	
	
	public VariantContextWriter createVariantContextWriter() throws IOException
		{
		XMLVcfWriter w=new XMLVcfWriter();
		try {
			XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
			if(this.outputFile!=null)
				{
				w.delegateOut=IOUtils.openFileForWriting(this.outputFile);
				}
			else
				{
				w.delegateOut=System.out;
				}
			w.writer = xmlfactory.createXMLStreamWriter(w.delegateOut);
			} 
		catch (XMLStreamException e)
			{
			CloserUtil.close(w.writer);;
			CloserUtil.close(w.delegateOut);;
			throw new IOException(e);
			}
		
		
		return w;
		}
	}
