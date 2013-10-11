package com.github.lindenb.jvarkit.util.vcf;

import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.vcf.VCFHeader;

public class XMLVcfReader 
	{
	private XMLEventReader in;
	
	
	
	private Allele parseAllele(StartElement root) throws XMLStreamException
		{
		boolean isRef=false;
		Attribute att=root.getAttributeByName(new QName("isRef"));
		if(att!=null && att.getValue().equals("true"))
			{
			isRef=true;
			}
		return Allele.create(in.getElementText(), isRef);
		}
	
	private VariantContext variant(StartElement root) throws XMLStreamException
		{
		String chrom=null;
		List<Allele> alleles=new ArrayList<Allele>();
		while(in.hasNext())
			{
			XMLEvent evt=in.nextEvent();
			if(evt.isStartElement())
				{
				StartElement SE=evt.asStartElement();
				if(SE.getName().getLocalPart().equals("chrom"))
					{
					chrom=in.getElementText().trim();
					}
				else if(SE.getName().getLocalPart().equals("allele"))
					{
					alleles.add(parseAllele(SE));
					}
				}
			}
	
		VariantContextBuilder b=new VariantContextBuilder("", chrom, 0, 0, alleles);
		return b.make();
		}
	
	
	
	private List<Object> _list(StartElement root) throws XMLStreamException
		{
		List<Object> L=new ArrayList<Object>();
		while(in.hasNext())
			{
			XMLEvent evt=in.nextEvent();
			if(evt.isStartElement())
				{
				StartElement E=evt.asStartElement();
				L.add(_any(E));
				}
			else if(evt.isEndElement())
				{
				EndElement E=evt.asEndElement();
				if(E.getName().equals(root.getName()))
					{
					return L;
					}
				}
			}
		throw new IllegalStateException("unclosed list");
		}
	
	private Map<String,Object> _object(StartElement root) throws XMLStreamException
		{
		Map<String,Object> hash=new LinkedHashMap<String,Object>();
		while(in.hasNext())
			{
			XMLEvent evt=in.nextEvent();
			if(evt.isStartElement())
				{
				StartElement E=evt.asStartElement();
				Attribute att=E.getAttributeByName(new QName("key"));
				if(att==null) throw new XMLStreamException("key missing");
				if(hash.containsKey(att.getValue()))
					{
					throw new XMLStreamException("duplicate key "+att);
					}
				hash.put(att.getValue(), _any(E));
				}
			else if(evt.isEndElement())
				{
				EndElement E=evt.asEndElement();
				if(E.getName().equals(root.getName()))
					{
					return hash;
					}
				}
			}
		throw new IllegalStateException("unclosed object");
		}
	
	private Object _any(StartElement root) throws XMLStreamException
		{
		String name=root.getName().getLocalPart();
		if(name.equals("list"))
			{
			return _list(root);
			}
		else if(name.equals("object"))
			{
			return _object(root);
			}
		throw new IllegalStateException("undefined "+name);
		}
	
	
	public XMLVcfReader(InputStream in) throws IOException
		{
		try
			{
			XMLInputFactory xif=XMLInputFactory.newFactory();
			this.in=xif.createXMLEventReader(in);
			}
		catch(XMLStreamException err)
			{
			throw new IOException(err);
			}
		}
	
	
	public VCFHeader readHeader() throws IOException {
		return null;
		}

	
	public void close() throws IOException {
		if(in!=null)
			{
			try {in.close();}
			catch(XMLStreamException err) {throw new IOException(err);}
			in=null;
			}
		}
	
	public VariantContext next()throws IOException
		{
		if(in==null) return null;
		try
			{
			while(in.hasNext())
				{
				XMLEvent evt=in.nextEvent();
				if(evt.isStartElement())
					{
					StartElement E=evt.asStartElement();
					if(E.getName().getLocalPart().equals("variant"))
						{
						return variant(E);
						}
					}
				}
			close();
			return null;
			}
		catch(XMLStreamException err)
			{
			throw new IOException(err);
			}
		}
	}
