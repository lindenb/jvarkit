package com.github.lindenb.jvarkit.util.vcf;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.CommonInfo;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFContigHeaderLine;
import org.broadinstitute.variant.vcf.VCFFilterHeaderLine;
import org.broadinstitute.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

public class XMLVcfWriter implements VariantContextWriter
	{
	private XMLStreamWriter writer;
	public XMLVcfWriter(XMLStreamWriter writer)
		{
		this.writer=writer;
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
	public void writeHeader(VCFHeader header)
		{
		try
			{
			start("vcf");
			start("header");
			start("infos");
			for(VCFInfoHeaderLine h:header.getInfoHeaderLines())
				{
				start("info");
				attribute("ID", h.getID());
				element("key",h.getKey());
				element("count",h.getCount());
				element("description",h.getDescription());
				element("countType",h.getCountType());
				element("value",h.getValue());
				h.getCountType();
				end();
				}
			end();
			
			start("formats");
			for(VCFFormatHeaderLine h: header.getFormatHeaderLines())
				{
				start("format");
				attribute("ID", h.getID());
				element("key",h.getKey());
				element("count",h.getCount());
				element("description",h.getDescription());
				element("countType",h.getCountType());
				element("value",h.getValue());
				end();
				}
			end();
			
			start("filters");
			for(VCFFilterHeaderLine h: header.getFilterLines())
				{
				start("filter");
				attribute("ID", h.getID());
				element("key",h.getKey());
				element("value",h.getValue());
				end();
				}
			end();
			
			
			start("contigs");
			for(VCFContigHeaderLine h:header.getContigLines())
				{
				start("contig");
				attribute("ID", h.getID());
				attribute("tid", h.getContigIndex());
				element("key",h.getKey());
				element("value",h.getValue());
				end();
				}
			end();
			
			start("samples");
			for(String name:header.getSampleNamesInOrder())
				{
				start("sample");
				attribute("index", header.getSampleNameToOffset().get(name));
				characters(name);
				end();
				}
			end();
			
			for(VCFHeaderLine meta:header.getMetaDataInInputOrder())
				{
				
				}
			
			end();//header
			start("variations");
			}
		catch (XMLStreamException e)
			{
			// TODO: handle exception
			}
		}
	@Override
	public void add(VariantContext variant)
		{
		try
			{
			start("variation");
			element("chrom",variant.getChr());
			element("start",variant.getStart());
			element("end",variant.getEnd());
			element("id",variant.getID());
			element("ref",variant.getReference().getDisplayString());
			
			for(Allele a:variant.getAlternateAlleles())
				{
				element("alt",a.getDisplayString());
				}
			element("qual",variant.getLog10PError());
			
			
			start("filters");
			for(String s: variant.getFilters())
				{
				element("filter",s);
				}
			end();
			
			start("infos");
			CommonInfo info=variant.getCommonInfo();
			for(String key:info.getAttributes().keySet())
				{
				start("info");
				attribute("key", key);
				Object v=info.getAttribute(key);
				
				end();
				}
			end();
			
			if(variant.hasGenotypes())
				{
				start("genotype");
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
					element("sample",g.getSampleName());
					if(g.hasAD())
						{
						start("ad");
						for(int ad:g.getAD())
							{
							element("value", ad);
							}
						g.getAD();
						end();
						}
					if(g.hasDP())
						{
						element("dp", g.getDP());
						}
					if(g.hasGQ())
						{
						element("gq", g.getGQ());
						}
					if(g.hasPL())
						{
						start("pl");
						for(int v:g.getPL())
							{
							element("value", v);
							}
						g.getAD();
						end();
						}
					if(g.hasLikelihoods())
						{
						
						}
					
					start("alleles");
					for(Allele a:g.getAlleles())
						{
						element("allele",a.getBaseString());
						}
					end();
					
					
					end();
					}
				end();
				}
			end();
			}
		catch (XMLStreamException e)
			{
			e.printStackTrace();
			}
		}

	@Override
	public void close()
		{
		try
			{
			end();//variations
			end();//vcf
			writer.close();
			}
		catch (XMLStreamException e)
			{
			e.printStackTrace();
			}
		
		}
	}
