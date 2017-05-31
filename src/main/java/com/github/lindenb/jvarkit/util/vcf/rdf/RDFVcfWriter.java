package com.github.lindenb.jvarkit.util.vcf.rdf;

import java.io.IOException;
import java.io.OutputStream;
import java.net.URI;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.xml.XMLConstants;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.tools.vcfannot.VCFPredictions;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree.Term;
import com.github.lindenb.jvarkit.util.vcf.predictions.MyPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.Prediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.PredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

/** 
 * RDFVcfWriter
 * @author lindenb
 *
 */
public class RDFVcfWriter
	implements VariantContextWriter
	{
	private static final Logger LOG=Logger.build(RDFVcfWriter.class).make();

	private static final String XSD="http://www.w3.org/2001/XMLSchema#";
	private static final String RDF=com.github.lindenb.jvarkit.util.ns.RDF.NS;
	private static final String DC="http://purl.org/dc/elements/1.1/";
	private static final String NS="http://github.com/lindenb/jvarkit/";
	private static final String PFX="vcf";
	private XMLStreamWriter w;
	private VCFHeader  header;
	private long id_generator=0L;
	private OutputStream delegateOut;
	private Map<String,RDFVcfInfoHandler> key2infoHandler=new HashMap<String,RDFVcfInfoHandler>();
	private URI source=null;
	private boolean _xmlHeaderPrinted=false;

	
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
	
	
	public void addInfoHandler(RDFVcfInfoHandler handler)
		{
		this.key2infoHandler.put(handler.getKey(), handler);
		}
	
	
	private void datatype(String t) throws XMLStreamException
		{
		this.w.writeAttribute("rdf", RDF, "datatype", "xsd:"+t);
		}
	
	protected RDFVcfInfoHandler createDefaultRdfVcfInfoHandlerFor(VCFInfoHeaderLine h)
		{
		return new DefaultInfoHandler(h);
		}
	
	public void writeStartDocument() throws XMLStreamException
		{
		if(_xmlHeaderPrinted) return;
		w.writeStartDocument("UTF-8","1.0");
		this.w.writeStartElement("rdf", "RDF", RDF);
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "rdf",RDF);
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "dc", DC);
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI,  PFX, NS);
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI,  "xsd", XSD);
		_xmlHeaderPrinted=true;
		}
	
	
	@Override
	public void writeHeader(VCFHeader  header)
		{
		this.writeHeader(header, null);
		}
	public void writeHeader(VCFHeader  header,URI source)
		{
		if(this.header!=null) throw new RuntimeException("Header was already written");
		this.header=header;
		this.source=source;
		if(this.source==null) this.source=URI.create("urn:source/id"+(++id_generator));

		
		try {
			writeStartDocument();
			
			this.w.writeStartElement(PFX, "Source", NS);
			this.w.writeAttribute("rdf",RDF,"about",this.source.toString());
				this.w.writeStartElement("dc","title",DC);
				this.w.writeCharacters(this.source.toString());
				this.w.writeEndElement();//dc:title
			this.w.writeEndElement();
			
			
			SAMSequenceDictionary dict=header.getSequenceDictionary();
			if(dict!=null)
				{
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
	
					
					this.w.writeEndElement();//Chromosome
					}
				}
			key2infoHandler.put(SnpEffPredictionParser.getDefaultTag(), new SnpEffHandler());
			key2infoHandler.put(VepPredictionParser.getDefaultTag(), new VepHandler());
			key2infoHandler.put(VCFPredictions.TAG, new MyPredictionHandler());
			for(VCFInfoHeaderLine h:header.getInfoHeaderLines())
				{
				RDFVcfInfoHandler handler= key2infoHandler.get(h.getID());
				if(handler==null)
					{
					LOG.info("creating default handler for INFO:"+h.getID());
					handler=createDefaultRdfVcfInfoHandlerFor(h);
					key2infoHandler.put(handler.getKey(), handler);
					}
				handler.init( h);
				
				}
			
			
			for(VCFFilterHeaderLine h:header.getFilterLines())
				{
				this.w.writeStartElement(PFX, "Filter", NS);
				this.w.writeAttribute("rdf",RDF,"about","urn:filter/"+h.getKey());
				
				this.w.writeStartElement("dc","title",DC);
				this.w.writeCharacters(h.getKey());
				this.w.writeEndElement();//dc:title


				this.w.writeStartElement("dc","description",DC);
				this.w.writeCharacters(h.getValue());
				this.w.writeEndElement();//dc:title
				
				this.w.writeEndElement();//Filter
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
			throw new RuntimeException("close failed",e);
			}
		}
	
	@Override
	public void add(VariantContext ctx)
		{
		if(this.header==null) throw new RuntimeException("No header was written.");
		try {
			long variant_id=++id_generator;
			this.w.writeStartElement(PFX, "Variant", NS);
			this.w.writeAttribute("rdf",RDF,"about","urn:variant/"+variant_id);
			
			this.w.writeEmptyElement(PFX, "source",NS);
			this.w.writeAttribute("rdf",RDF,"resource",this.source.toString());
			
			this.w.writeEmptyElement(PFX, "chromosome",NS);
			this.w.writeAttribute("rdf",RDF,"resource","urn:chromosome/"+ctx.getContig());
			
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
				if(ctx.getID().matches("rs[0-9]+"))
					{
					this.w.writeEmptyElement(PFX,"ID",NS);
					this.w.writeAttribute("rdf",RDF,"resource",
							"http://www.ncbi.nlm.nih.gov/snp/"+
							ctx.getID().substring(2)
							);
					}
				else
					{
					this.w.writeStartElement(PFX,"ID",NS);
					this.w.writeCharacters(ctx.getID());
					this.w.writeEndElement();
					}
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
				this.w.writeEmptyElement(PFX,"filter",NS);
				this.w.writeAttribute("rdf",RDF,"resource",
						"urn:filter/"+ filt
						);
				}
			
			//INFO
			for(String key:ctx.getAttributes().keySet())
				{
				RDFVcfInfoHandler handler=this.key2infoHandler.get(key);
				if(handler==null) continue;
				handler.handle(ctx);
				
				}
			
			this.w.writeEndElement();//Variant
			
			for(Genotype g:ctx.getGenotypes())
				{
				if(!g.isAvailable()) continue;
				if(!g.isCalled()) continue;
				long genotype_id=++id_generator;
				this.w.writeStartElement(PFX, "Genotype", NS);
				this.w.writeAttribute("rdf",RDF,"about","urn:genotype/"+genotype_id);

				
				
				this.w.writeEmptyElement(PFX, "sample",NS);
				this.w.writeAttribute("rdf",RDF,"resource","urn:sample/"+g.getSampleName());

				
				this.w.writeEmptyElement(PFX, "variant",NS);
				this.w.writeAttribute("rdf",RDF,"resource","urn:variant/"+variant_id);
				
				/*
				this.w.writeEmptyElement("rdf", "type",RDF);
				this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/"+(g.isAvailable()?"available":"unavaliable"));
				*/
				
				/*
				if(g.isCalled() )
					{
					this.w.writeEmptyElement("rdf", "type",RDF);
					this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/called");
					}*/
				
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
				
				if(g.isHomRef())
					{
					this.w.writeEmptyElement("rdf", "type",RDF);
					this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/homRef");
					}
				if(g.isHomVar())
					{
					this.w.writeEmptyElement("rdf", "type",RDF);
					this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/homVar");
					}
				
				if(g.isMixed())
					{
					this.w.writeEmptyElement("rdf", "type",RDF);
					this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/mixed");
					}
				
				if(g.isPhased())
					{
					this.w.writeEmptyElement("rdf", "type",RDF);
					this.w.writeAttribute("rdf",RDF,"resource","urn:genotype/phased");
					}
				
				Set<String> seen=new HashSet<String>();
				for(Allele a:g.getAlleles())
					{
					if(a.isNoCall() || seen.contains(a.getBaseString())) continue;
					this.w.writeStartElement(PFX,"allele",NS);
					this.w.writeCharacters(a.getBaseString());
					this.w.writeEndElement();
					seen.add(a.getBaseString());
					}
				
				if(g.hasDP())
					{
					this.w.writeStartElement(PFX,"dp",NS);
					datatype("int");
					this.w.writeCharacters(String.valueOf(g.getDP()));
					this.w.writeEndElement();
					}
				
				if(g.hasGQ())
					{
					this.w.writeStartElement(PFX,"gq",NS);
					datatype("int");
					this.w.writeCharacters(String.valueOf(g.getGQ()));
					this.w.writeEndElement();
					}
				if(g.hasPL())
					{
					int pl[]=g.getPL();
					this.w.writeStartElement(PFX,"pl",NS);
					for(int i=0;i<pl.length;++i)
						{
						if(i>0) this.w.writeCharacters(",");
						this.w.writeCharacters(String.valueOf(pl[i]));
						}
					this.w.writeEndElement();
					}
				
				this.w.writeEndElement();
				}
			
			
			
			}
		catch(XMLStreamException e) {
			throw new RuntimeException("add failed",e);
			}
		}
	@Override
	public void close()
		{
		if(this.w==null) return;
		if(this.header==null) throw new RuntimeException("No header was written.");
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
			e.printStackTrace();
			throw new RuntimeException("close failed",e);
			}
		}
	
	

		public interface RDFVcfInfoHandler
			{
			public String getKey();
			public void init(
					VCFInfoHeaderLine line
					)  throws XMLStreamException;
			public void handle(
					VariantContext ctx
					) throws XMLStreamException;
			}

	  private abstract class AbstractInfoHandler
	  		implements RDFVcfInfoHandler
			{
			protected VCFInfoHeaderLine info;
			public AbstractInfoHandler(VCFInfoHeaderLine info)
				{
				this.info=info;
				}
			
			@Override
			public void init(VCFInfoHeaderLine line)
					throws XMLStreamException
				{
				w.writeStartElement(PFX,"InfoHeader",NS);
				w.writeAttribute("rdf", RDF, "about", "urn:info/"+line.getID());
				
				w.writeStartElement(PFX,"id",NS);
				w.writeCharacters(line.getID());
				w.writeEndElement();
				
				w.writeStartElement(PFX,"description",NS);
				w.writeCharacters(line.getDescription());
				w.writeEndElement();
				
				w.writeStartElement(PFX,"countType",NS);
				w.writeCharacters(line.getCountType().name());
				w.writeEndElement();
				
				w.writeStartElement(PFX,"type",NS);
				w.writeCharacters(line.getType().name());
				w.writeEndElement();
				
				
				w.writeEndElement();
				}
			
			protected abstract void handleObject(Object o)
					throws XMLStreamException;
			
			@Override
			public void handle(VariantContext ctx)
					throws XMLStreamException
				{
				Object o=ctx.getAttribute(this.getKey());
				if(o==null) return;
				if(o.getClass().isArray())
					{
					Object array[]=(Object[])o;
					for(Object o2:array) handleObject(o2);
					}
				else if(o instanceof Collection)
					{
					@SuppressWarnings("rawtypes")
					Collection array=(Collection)o;
					for(Object o2:array) handleObject(o2);
					}
				else
					{
					handleObject(o);
					}
				}
			
			@Override
			public String getKey()
				{
				return info.getID();
				}
			};
	  		
	private abstract class AbstractPredHandler
		implements RDFVcfInfoHandler
	{	
	public AbstractPredHandler()
		{
		}
	
	abstract PredictionParser getPredictionParser();
	abstract String getLocalName();
	
	@Override
	public void init(VCFInfoHeaderLine line) throws XMLStreamException
		{
		// TODO Auto-generated method stub
		
		}
	
	@Override
	public void handle(VariantContext ctx) throws XMLStreamException
		{
		w.writeComment("X"+ctx.getAttribute(getPredictionParser().getTag()));
		
		for(final Prediction pred:this.getPredictionParser().getPredictions(ctx))
			{
			/*
			w.writeStartElement(PFX,"prediction",NS);
			w.writeStartElement(PFX,getLocalName(),NS);
			Integer i=pred.getAminoAcidPosition();
			if(i!=null)
				{
				w.writeStartElement(PFX,"aminoAcidPosition",NS);
				datatype("int");
				w.writeCharacters(String.valueOf(i));
				w.writeEndElement();
				}
			
			String s=pred.getEnsemblTranscript();
			if(s!=null)
				{
				w.writeEmptyElement(PFX,"enst",NS);
				w.writeAttribute("rdf",RDF,"resource","http://www.ensembl.org/"+s);
				}
			s=pred.getEnsemblProtein();
			if(s!=null)
				{
				w.writeEmptyElement(PFX,"ensp",NS);
				w.writeAttribute("rdf",RDF,"resource","http://www.ensembl.org/"+s);
				}
			s=pred.getGeneName();
			
			if(s!=null)
				{
				w.writeStartElement(PFX,"geneName",NS);
				w.writeCharacters(s);
				w.writeEndElement();
				}
			
			s=pred.getReferenceAminoAcid();
			if(s!=null)
				{
				w.writeStartElement(PFX,"refAA",NS);
				w.writeCharacters(s);
				w.writeEndElement();
				}
			s=pred.getAltAminoAcid();
			if(s!=null)
				{
				w.writeStartElement(PFX,"altAA",NS);
				w.writeCharacters(s);
				w.writeEndElement();
				}
			for(Term term:pred.getSOTerms())
				{
				w.writeStartElement(PFX,"so",NS);
				w.writeCharacters(term.getAcn());
				w.writeEndElement();
				}
			
			w.writeEndElement();
			w.writeEndElement();
			*/
			}
		
		}
	}
			
			
			
	private class SnpEffHandler
		  		extends AbstractPredHandler
		{	
		private SnpEffPredictionParser predFactory;
		public SnpEffHandler()
			{
			this.predFactory=new SnpEffPredictionParserFactory(RDFVcfWriter.this.header).get();
			}
		
		@Override
		PredictionParser getPredictionParser()
			{
			return predFactory;
			}
		
		@Override
		public String getKey()
			{
			return predFactory.getTag();
			}
		@Override
		String getLocalName()
			{
			return "snpEff";
			}
		}
			
	private class VepHandler
				extends AbstractPredHandler
		{	
		private VepPredictionParser predFactory;
		public VepHandler()
			{
			this.predFactory=new VepPredictionParserFactory(RDFVcfWriter.this.header).get();
			}
		
		@Override
		PredictionParser getPredictionParser()
			{
			return predFactory;
			}
		
		@Override
		public String getKey()
			{
			return predFactory.getTag();
			}
		@Override
		String getLocalName()
			{
			return "vep";
			}
		}
	
	private class MyPredictionHandler
	extends AbstractPredHandler
		{	
		private MyPredictionParser predFactory;
		public MyPredictionHandler()
			{
			this.predFactory=new MyPredictionParser(RDFVcfWriter.this.header);
			}
		
		@Override
		PredictionParser getPredictionParser()
			{
			return predFactory;
			}
		
		@Override
		public String getKey()
			{
			return predFactory.getTag();
			}
		@Override
		String getLocalName()
			{
			return "myprediction";
			}
		}
	
  private class DefaultInfoHandler
  		extends AbstractInfoHandler
		{
		public DefaultInfoHandler(VCFInfoHeaderLine info)
			{
			super(info);
			}

		@Override
		protected void handleObject(Object o)
				throws XMLStreamException
			{
			w.writeStartElement(PFX, getKey(), NS);
			if(o.getClass()==Double.class)
				{
				datatype("double");
				}
			else if(o.getClass()==Float.class)
				{
				datatype("float");
				}
			else if(o.getClass()==Integer.class)
				{
				datatype("int");
				}
			w.writeCharacters(String.valueOf(o));
			w.writeEndElement();
			}
		};
	
	@Override
	public boolean checkError() {
		return false;
		}
	}
