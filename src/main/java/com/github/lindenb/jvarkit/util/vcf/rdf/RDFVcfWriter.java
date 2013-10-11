package com.github.lindenb.jvarkit.util.vcf.rdf;

import java.io.IOException;
import java.io.OutputStream;

import javax.xml.XMLConstants;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import net.sf.picard.PicardException;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;

public class RDFVcfWriter
	implements VariantContextWriter
	{
	private static final String XSD="http://www.w3.org/2001/XMLSchema#";
	private static final String RDF="http://www.w3.org/1999/02/22-rdf-syntax-ns#";
	private static final String DC="http://purl.org/dc/elements/1.1/";
	private static final String NS="http://github.com/lindenb/jvarkit/";
	private static final String PFX="vcf";
	private XMLStreamWriter w;
	private VCFHeader  header;
	private long id_generator=0L;
	private OutputStream delegateOut;
	
	public RDFVcfWriter(XMLStreamWriter writer)
		{
		this.w=writer;
		}
	
	public RDFVcfWriter(OutputStream delegateOut) throws IOException
		{
		try {
			this.delegateOut=delegateOut;
			XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
			this.w= xmlfactory.createXMLStreamWriter(delegateOut,"UTF-8");
			} 
		catch (XMLStreamException e)
			{
			throw new IOException(e);
			}
		}
	
	private void datatype(String t) throws XMLStreamException
		{
		this.w.writeAttribute("rdf", RDF, "datatype", "xsd:"+t);
		}
	
	@Override
	public void writeHeader(VCFHeader  header)
		{
		if(this.header!=null) throw new PicardException("Header was already written");
		this.header=header;
		try {
			w.writeStartDocument("UTF-8","1.0");
			this.w.writeStartElement("rdf", "RDF", RDF);
			w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "rdf",RDF);
			w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "dc", DC);
			w.writeAttribute("xmlns", XMLConstants.XML_NS_URI,  PFX, NS);
			w.writeAttribute("xmlns", XMLConstants.XML_NS_URI,  "xsd", XSD);

			
			SAMSequenceDictionary dict=header.getSequenceDictionary();
			for(SAMSequenceRecord ssr:dict.getSequences())
				{
				this.w.writeStartElement(PFX, "Chromosome", NS);
				this.w.writeAttribute("rdf",RDF,"about","urn:chromosome/"+ssr.getSequenceName());
				
				this.w.writeStartElement("dc","title",DC);
				this.w.writeCharacters(ssr.getSequenceName());
				this.w.writeEndElement();//dc:title
				
				this.w.writeStartElement(PFX,"length",NS);
				datatype("int");
				this.w.writeCharacters(String.valueOf(ssr.getSequenceLength()));
				this.w.writeEndElement();//length

				this.w.writeStartElement(PFX,"index",NS);
				datatype("int");
				this.w.writeCharacters(String.valueOf(ssr.getSequenceIndex()));
				this.w.writeEndElement();//length

				
				this.w.writeEndElement();//rdf:RDF
				}
			
			//Sample
			for(String sample:header.getSampleNamesInOrder())
				{
			
				this.w.writeStartElement(PFX, "Sample", NS);
				this.w.writeAttribute("rdf",RDF,"about","urn:sample/"+sample);
				
				this.w.writeStartElement("dc","title",DC);
				this.w.writeCharacters(sample);
				this.w.writeEndElement();//dc:title
				
				this.w.writeEndElement();//rdf:RDF
				}
			
			}
		catch(Exception e) {
			throw new PicardException("close failed",e);
			}
		}
	
	@Override
	public void add(VariantContext ctx)
		{
		if(this.header==null) throw new PicardException("No header was written.");
		try {
			++id_generator;
			this.w.writeStartElement(PFX, "Variant", NS);
			this.w.writeAttribute("rdf",RDF,"about","urn:variant/"+id_generator);

			
			this.w.writeEmptyElement(PFX, "chromosome",NS);
			this.w.writeAttribute("rdf",RDF,"resource","urn:chromosome/"+ctx.getChr());
			
			this.w.writeStartElement(PFX,"start",NS);
			datatype("int");
			this.w.writeCharacters(String.valueOf(ctx.getStart()));
			this.w.writeEndElement();
			
			this.w.writeStartElement(PFX,"end",NS);
			datatype("int");
			this.w.writeCharacters(String.valueOf(ctx.getEnd()));
			this.w.writeEndElement();

			if(ctx.hasID())
				{
				this.w.writeStartElement(PFX,"ID",NS);
				this.w.writeCharacters(ctx.getID());
				this.w.writeEndElement();
				}
			
			this.w.writeStartElement(PFX,"ref",NS);
			this.w.writeCharacters(ctx.getReference().getBaseString());
			this.w.writeEndElement();
			
			for(Allele a:ctx.getAlleles())
				{
				this.w.writeStartElement(PFX,"alt",NS);
				this.w.writeCharacters(a.getBaseString());
				this.w.writeEndElement();
				}
			
			if(ctx.hasLog10PError())
				{
				this.w.writeStartElement(PFX,"qual",NS);
				datatype("double");
				this.w.writeCharacters(String.valueOf(ctx.getPhredScaledQual()));
				this.w.writeEndElement();
				}
			
			for(String filt:ctx.getFilters())
				{
				this.w.writeStartElement(PFX,"filter",NS);
				this.w.writeCharacters(filt);
				this.w.writeEndElement();
				}
			//INFO
			
			this.w.writeEndElement();//Variant
			
			for(Genotype g:ctx.getGenotypes())
				{
				this.w.writeStartElement(PFX, "Call", NS);
				
				this.w.writeEmptyElement(PFX, "sample",NS);
				this.w.writeAttribute("rdf",RDF,"resource","urn:sample/"+g.getSampleName());

				
				this.w.writeEmptyElement(PFX, "variant",NS);
				this.w.writeAttribute("rdf",RDF,"resource","urn:variant/"+id_generator);
				
				this.w.writeEmptyElement("rdf", "type",RDF);
				this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/"+(g.isAvailable()?"available":"unavaliable"));
				
				if(g.isCalled())
					{
					this.w.writeEmptyElement("rdf", "type",RDF);
					this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/called");
					}
				
				if(g.isFiltered())
					{
					this.w.writeEmptyElement("rdf", "type",RDF);
					this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/filtered");
					}
				
				if(g.isHom())
					{
					this.w.writeEmptyElement("rdf", "type",RDF);
					this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/hom");
					}
				
				if(g.isHet())
					{
					this.w.writeEmptyElement("rdf", "type",RDF);
					this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/het");
					}

				
				for(Allele a:g.getAlleles())
					{
					this.w.writeStartElement(PFX,"allele",NS);
					this.w.writeCharacters(a.getBaseString());
					this.w.writeEndElement();
					}
				
				
				
				
				this.w.writeEndElement();
				}
			
			
			
			}
		catch(XMLStreamException e) {
			throw new PicardException("add failed",e);
			}
		}
	@Override
	public void close()
		{
		if(this.w==null) return;
		if(this.header==null) throw new PicardException("No header was written.");
		try {
			this.w.writeEndElement();//rdf:RDF
			this.w.writeEndDocument();
			this.w.flush();
			this.w.close();
			if(this.delegateOut!=null) {delegateOut.flush(); delegateOut.close();}
			this.w=null;
			} 
		catch (Exception e)
			{
			throw new PicardException("close failed",e);
			}
		}
	}
