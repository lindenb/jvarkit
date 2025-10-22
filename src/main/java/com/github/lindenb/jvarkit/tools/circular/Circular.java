package com.github.lindenb.jvarkit.tools.circular;

import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.transform.dom.DOMSource;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.yaml.snakeyaml.Yaml;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.svg.SVGDocument;
import com.github.lindenb.jvarkit.xml.ExtendedXmlStreamWriter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class Circular extends Launcher {
	private static final Logger LOG = Logger.of(Circular.class);
	@Parameter(names="-o",description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output = null;

	
	private XPath xPathInstance=null;
	
	private static String format(Object s) {
		return  String.valueOf(s);
		}
	
	private static class ContigInfo {
		final SAMSequenceRecord ssr;
		int tid;
		double angle_start=0;
		double angle_end=0;
		ContigInfo(final SAMSequenceRecord ssr) {
			this.ssr = ssr;
			}
		}
	private final List<ContigInfo> contigs= new ArrayList<>(); 
	private final Map<String,ContigInfo> name2contig = new HashMap<>(); 
	
	
	private double pos2angle(String contig,int pos) {
		final ContigInfo ci= this.name2contig.get(contig);
		pos = Math.max(1, pos);
		pos = Math.min(pos, ci.ssr.getLengthOnReference());
		return ci.angle_start+ (pos/(double)ci.ssr.getLengthOnReference())*(ci.angle_end-ci.angle_end);
		}
	
	private String pointToStr(final Point2D.Double p)
		{
		return format(p.x)+" "+format(p.y);
		}

	
		
	private Point2D.Double polarToCartestian(double radius,double angle)
		{
		return new Point2D.Double(
				Math.cos(angle)*radius,
				Math.sin(angle)*radius
				);
		}

	private Map<String,Object> loadConfig(BufferedReader br) throws IOException {
		final List<Object> configs = new ArrayList<>();
		final Yaml ymal=new Yaml();
		ymal.loadAll(br).forEach(O->configs.add(O));
		if(configs.size()!=1) throw new IOException("expectedn one object");
		if(!(configs.get(0) instanceof Map)) throw new IOException("object is not a map");
		return ( Map<String,Object>)configs.get(0);
		}
	
	
	
	private Optional<String> get(Node n,String path) throws Exception {
		String s=(String)this.xPathInstance.evaluate(path, new DOMSource(n), XPathConstants.STRING);
		if(StringUtils.isBlank(s)) return Optional.empty();
		return Optional.of(s);
		}
	
	
	@Override
	public int doWork(List<String> args) {
		try {
			this.xPathInstance = XPathFactory.newInstance().newXPath();
			final String input = super.oneFileOrNull(args);
			final Document config;
			
			final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			dbf.setNamespaceAware(false);
			dbf.setXIncludeAware(true);
			dbf.setCoalescing(true);
			final DocumentBuilder db = dbf.newDocumentBuilder();
			
			if(input==null || input.equals("-") ) {
				config= db.parse(stdin());
				}
			else
				{
				config = db.parse(IOUtils.openURIForReading(input));
				}
			
			if(this.xPathInstance.evaluate("/config/reference", config,XPathConstants.NODE)==null) {
				LOG.error("reference missing");
				return -1;
				}
			if(this.xPathInstance.evaluate("/config/image", config,XPathConstants.NODE)==null) {
				LOG.error("image missing");
				return -1;
				}
			
			final Element refNode = (Element)this.xPathInstance.evaluate("/config/reference", config,XPathConstants.NODE);
			final Element imageNode = (Element)this.xPathInstance.evaluate("/config/image",config, XPathConstants.NODE);
				
			
			String dictPath =(String)this.xPathInstance.evaluate("fasta/text()",refNode, XPathConstants.STRING);
			if(dictPath==null) dictPath = (String)this.xPathInstance.evaluate("dict/text()",refNode, XPathConstants.STRING);
			if(dictPath==null) {
				LOG.error("no fasta or dict reference");
				}
			SAMSequenceDictionary dict0 = new SequenceDictionaryExtractor().extractRequiredDictionary(dictPath.trim());
			String s= (String)this.xPathInstance.evaluate("regex",refNode,XPathConstants.STRING);
			if(!StringUtils.isBlank(s)) {
				final Pattern regex = Pattern.compile(s.trim());
				dict0 = new SAMSequenceDictionary(dict0
						.getSequences()
						.stream()
						.filter(C->regex.matcher(C.getContig()).matches())
						.collect(Collectors.toList()));
				}
			s= (String)this.xPathInstance.evaluate("min-length",refNode,XPathConstants.STRING);
			if(!StringUtils.isBlank(s)) {
				final int minLen = Integer.parseInt(s);
				dict0 = new SAMSequenceDictionary(dict0
						.getSequences()
						.stream()
						.filter(C->C.getSequenceLength()>=minLen)
						.collect(Collectors.toList()));
				}
			final SAMSequenceDictionary dict = dict0;
			final long genomeLen = dict.getReferenceLength();
			if(genomeLen < 1) {
				LOG.error("genome is empty");
				return -1;
				}
			//draw chromosomes
			
			long x=0;
			for(SAMSequenceRecord ssr:dict.getSequences()) {
				final ContigInfo ci=new ContigInfo(ssr);
				ci.tid = this.contigs.size();
				ci.angle_start = (x/(double)genomeLen)*Math.PI*2;
				x += ssr.getSequenceLength();
				ci.angle_end = (x/(double)genomeLen)*Math.PI*2;
				this.contigs.add(ci);
				this.name2contig.put(ssr.getSequenceName(), ci);
				}
			
			final double imageRadius = 1000;
			
			try(OutputStream os= super.openPathOrStdoutAsStream(this.output)) {
				final XMLOutputFactory xof = XMLOutputFactory.newFactory();
				final XMLStreamWriter w=  xof.createXMLStreamWriter(os, "UTF-8");
				w.writeStartDocument("UTF-8", "1.0");
				w.writeStartElement("svg");
				w.writeDefaultNamespace(SVG.NS);
				w.writeAttribute("width", format(imageRadius));
				w.writeAttribute("height", format(imageRadius));
				w.writeStartElement("style");
				w.writeEndElement();//style
				
				w.writeStartElement("g");

				w.writeStartElement("g");
				w.writeAttribute("transform", "translate("+format(imageRadius)+","+format(imageRadius)+")");
				w.writeStartElement("g");

				
				
				w.writeEndElement();//g
				w.writeEndElement();//g
				w.writeEndElement();//g
				w.writeEndElement();//svg
				w.writeEndDocument();
				w.close();
				w.flush();
				os.flush();
				}	
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new Circular().instanceMainWithExit(args);
		}
}
