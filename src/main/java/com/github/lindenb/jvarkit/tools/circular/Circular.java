/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.circular;

import java.awt.Dimension;
import java.awt.Insets;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.transform.dom.DOMSource;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathException;
import javax.xml.xpath.XPathFactory;

import org.apache.commons.math3.geometry.spherical.oned.Arc;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.ExtendedLocatable;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.tools.genome2svg.beans.BamCoverageTrack;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
@Program(name="circular",
description="circular genome browser renderer",
keywords= {"genome","browser","circular","vcf","svg"},
modificationDate="20251022",
creationDate="20220522",
jvarkit_amalgamion = true
)
public class Circular extends Launcher {
	private static final Logger LOG = Logger.of(Circular.class);
	@Parameter(names="-o",description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output = null;

	
	private XPath xPathInstance=null;
	
	private static String format(Object s) {
		return  String.valueOf(s);
		}
	
	
	
	private abstract class Track {
		double radius_max;
		double radius_min;
		abstract void plot(XMLStreamWriter w) throws IOException,XMLStreamException;
		
		void writeTitle(XMLStreamWriter w,Object t)throws IOException,XMLStreamException {
			if(t==null) return;
			w.writeStartElement("title");
			w.writeCharacters(String.valueOf(t));
			w.writeEndElement();
			}
		double getWeigth() {
			return 1.0;
			}
		}
	
	private  class PlotContig extends Track{
		@Override
		void plot(XMLStreamWriter w) throws IOException, XMLStreamException {
			w.writeStartElement("g");
			for(ContigInfo ci:Circular.this.contigs) {
				w.writeStartElement("g");
				w.writeStartElement("path");
				w.writeAttribute("p",ci.arc(this.radius_min,this.radius_max,1,ci.getLengthOnReference()));
				w.writeEndElement();
				
				w.writeStartElement("text");
				w.writeStartElement("textPath");
				w.writeAttribute("path", "");
				w.writeCharacters(ci.getContig());
				w.writeEndElement();

				w.writeEndElement();
				
				writeTitle(w, ci.getContig());
				w.writeEndElement();
				}
			w.writeEndElement();
			}
		}
	private abstract class ConfigurableTrack extends Track{
		protected Element root;
		ConfigurableTrack(final Element root) {
			this.root=root;
			}
		double getWeigth() {
			return 1.0;
			}
		}
	
	private  class GroupTrack extends ConfigurableTrack{
		List<Track> subTracks=null;
		GroupTrack(final Element root) {
			super(root);
			}
		
		List<Track> getTracks() {
			if(this.subTracks==null) {
				this.subTracks=new ArrayList<>();
				for(Node c=root.getFirstChild();c!=null;c=c.getNextSibling()) {
					if(c.getNodeType()!=Node.ELEMENT_NODE) continue;
					}
				}
			return this.subTracks;
		}
		
		double getWeigth() {
			return 1.0;
			}
		@Override
		void plot(XMLStreamWriter w) throws IOException, XMLStreamException {
			
			
			
			 
			}
		}

	
	private abstract class AbstractBamTrack extends ConfigurableTrack{
		AbstractBamTrack(final Element root) {
			super(root);
			}
		 Predicate<SAMRecord> createReadFilter() {
			Predicate<SAMRecord> p =R->!R.getReadUnmappedFlag();
			Optional<String> opt=Circular.this.evalString(this.root, "mapq");
			if(opt.isPresent()) {
				int mapq = Integer.parseInt(opt.get());
				p =p.and(R->R.getMappingQuality()>=mapq);
				}
			return p;
		 	}
		 SamReader open() throws IOException {
			SamReaderFactory srf = SamReaderFactory.make();
			Optional<String> opt=Circular.this.evalString(this.root, "fasta");
			if(!opt.isPresent()) opt=Circular.this.evalString(this.root, "/config/reference/fasta");
			Path fasta=null;
			if(opt.isPresent()) fasta= Paths.get(opt.get().trim());
			String bamPath=Circular.this.evalString(this.root, "bam").orElseThrow(()->new IllegalArgumentException("bam undefined"));
			final SamInputResource si = SamInputResource.of(bamPath);
			srf.validationStringency(ValidationStringency.LENIENT);
			if(fasta!=null) srf.referenceSequence(fasta);
			return srf.open(si);
			}
		}
	
	private class CoverageBamTrack extends AbstractBamTrack{
		CoverageBamTrack(Element root) {
			super(root);
			}
		double extract(int[] cov,int start,int end) {
			return Arrays.stream(cov, start, end).average().orElseThrow();
			}
		
		void plot(XMLStreamWriter w) throws IOException,XMLStreamException {
			try(SamReader sr=open()) {
				for(ContigInfo ci:Circular.this.contigs) {
					if(sr.getFileHeader().getSequenceDictionary().getSequence(ci.getContig())==null) continue;
					final int[] cov= new int[ci.getLengthOnReference()];
					final Predicate<SAMRecord> filter = createReadFilter();
					try(CloseableIterator<SAMRecord> iter = sr.query(ci.getContig(), 1, ci.getLengthOnReference(), true)) {
						while(iter.hasNext()) {
							final SAMRecord rec = iter.next();
							if(!filter.test(rec)) continue;
							for(AlignmentBlock ab: rec.getAlignmentBlocks()) {
								for(int i=0;i< ab.getLength();i++) {
									int pos1 = ab.getReferenceStart()+i;
									if(pos1<1 || pos1> cov.length) continue;
									cov[pos1-1]++;
									}
								}
							final StringBuilder sb=new StringBuilder();
							sb.append("M "+pointToStr( ci.toPoint(1, radius_min)));
							
							double max_value = 40;
							double arc_len = radius_max*(ci.angle_end-ci.angle_start);
							int dx= (int)Math.max(1,(ci.getLengthOnReference()/arc_len));
							int x=0;
							while(x<cov.length) {
								int x2= Math.min(ci.getLengthOnReference(),x+dx);
								double value =  Math.min(max_value,extract(cov,x,x2));
								double r = this.radius_min+ (value/max_value)*(radius_max-radius_min);
								x=x2;
								sb.append("L "+pointToStr( ci.toPoint(x+dx/2, r)));
								}
							sb.append("L "+pointToStr( ci.toPoint(ci.getLengthOnReference(), radius_min)));
							sb.append(ci.arcTo(radius_min,1,0));
							sb.append("Z");
							
							w.writeStartElement("path");
							w.writeAttribute("p",sb.toString());
							w.writeEndElement();
							}
						
						}
					}
				}
			}
		}
	
	private class ContigInfo implements ExtendedLocatable{
		final SAMSequenceRecord ssr;
		int tid;
		double angle_start=0;
		double angle_end=0;
		ContigInfo(final SAMSequenceRecord ssr) {
			this.ssr = ssr;
			}
		
		double pos2angle(final int genomic_pos1) {
			return angle_start+ (genomic_pos1/(double)getLengthOnReference())*(angle_end-angle_start);
			}
		
		Point2D toPoint(int genomic_pos1,double radius) {
			final double a = pos2angle(genomic_pos1);
			return new Point2D.Double(
					radius * Math.cos(a),
					radius * Math.sin(a)
					);
			}
		
		String arc(double r_min,double r_max,int start,int end) {
			return String.join(" ",
					arcTo(r_min, start, end),
					arcTo(r_max, start, end)
					);
			}
		@Override
		public String getContig() {
			return ssr.getContig();
			}
		@Override
		public int getStart() {
			return 1;
			}
		@Override
		public int getEnd() {
			return ssr.getSequenceLength();
			}
		private String arcTo(double radius,int start,int end) {
			return String.join(" ",
					"A",
					format(radius), // radius X
					format(radius), // radius Y
					format(0), //rotate
					format(0), //small
					(start<end?"1":"0"),//direction
					pointToStr(toPoint(end,radius))
					);
			}
		}
	private final List<ContigInfo> contigs= new ArrayList<>(); 
	private final Map<String,ContigInfo> name2contig = new HashMap<>(); 
	
	
	private String pointToStr(final Point2D p)
		{
		return format(p.getX())+" "+format(p.getY());
		}

	
	
	private Optional<String> evalString(final Node node,final String path) {
		try 
			{
			final String s=(String)this.xPathInstance.evaluate(path, node, XPathConstants.STRING);
			if(StringUtils.isBlank(s)) return Optional.empty();
			return Optional.of(s);
			} 
		catch(XPathException err) {
			throw new RuntimeException(err);
			}
		}
	
	private Track createTrack(Element root) {
		if(root.getNodeName().equals("g") || root.getNodeName().equals("group")) {
			return new GroupTrack(root);
			}
		if(root.getNodeName().equals("coverage")) {
			return new CoverageBamTrack(root);
			}
		throw new IllegalArgumentException("cannot create track for "+root.getNodeName());
		}
	
	@Override
	public int doWork(final List<String> args) {
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
			
			
			
			final Element refNode = (Element)this.xPathInstance.evaluate("/config/reference", config,XPathConstants.NODE);
			if(refNode==null) {
				LOG.error("reference missing");
				return -1;
				}
			final Element imageNode = (Element)this.xPathInstance.evaluate("/config/image",config, XPathConstants.NODE);
			if(imageNode==null) {
				LOG.error("image missing");
				return -1;
				}
			final Element tracks = (Element)this.xPathInstance.evaluate("/config/tracks",config, XPathConstants.NODE);
			if(tracks==null) {
				LOG.warning("/config/tracks missing");
				}
			
			Optional<String> opt  =  this.evalString(refNode,"fasta/text()");
			if(!opt.isPresent()) opt = this.evalString(refNode,"dict/text()");
			if(!opt.isPresent()) {
				LOG.error("no fasta or dict reference");
				}
			SAMSequenceDictionary dict0 = new SequenceDictionaryExtractor().extractRequiredDictionary(opt.get());
			opt = this.evalString(refNode,"regex");
			if(opt.isPresent()) {
				final Pattern regex = Pattern.compile(opt.get().trim());
				dict0 = new SAMSequenceDictionary(dict0
						.getSequences()
						.stream()
						.filter(C->regex.matcher(C.getContig()).matches())
						.collect(Collectors.toList()));
				}
			opt = this.evalString(refNode,"min-length");
			if(opt.isPresent()) {
				final int minLen = Integer.parseInt(opt.get());
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
			
			final int default_margin = Integer.parseInt(this.evalString(imageNode,"margin").orElse("10"));
			final Insets margin=new Insets(default_margin, default_margin, default_margin, default_margin);
			this.evalString(imageNode,"margin-top").ifPresent(S->margin.top=Integer.parseInt(S));
			this.evalString(imageNode,"margin-bottom").ifPresent(S->margin.bottom=Integer.parseInt(S));
			this.evalString(imageNode,"margin-left").ifPresent(S->margin.left=Integer.parseInt(S));
			this.evalString(imageNode,"margin-right").ifPresent(S->margin.right=Integer.parseInt(S));
			final int imageRadius=  Integer.parseInt(this.evalString(imageNode,"radius").orElse("1000"));
			final Dimension imageSize= new Dimension(
					2*imageRadius+margin.left+margin.right,
					2*imageRadius+margin.top+margin.bottom
					);
			
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
			
			double radius = imageRadius;
			final double min_radius = imageRadius*0.5;
			try(OutputStream os= super.openPathOrStdoutAsStream(this.output)) {
				final XMLOutputFactory xof = XMLOutputFactory.newFactory();
				final XMLStreamWriter w=  xof.createXMLStreamWriter(os, "UTF-8");
				w.writeStartDocument("UTF-8", "1.0");
				w.writeStartElement("svg");
				w.writeDefaultNamespace(SVG.NS);
				w.writeAttribute("width", format(imageSize.width+1));
				w.writeAttribute("height", format(imageSize.height+1));
				w.writeStartElement("style");
				w.writeEndElement();//style
				
				w.writeStartElement("g");

				// write background panel
				opt=this.evalString(imageNode,"background-color");
				if(opt.isPresent()) {
					w.writeEmptyElement("rect");
					w.writeAttribute("x", format(0));
					w.writeAttribute("y", format(0));
					w.writeAttribute("width", format(imageSize.width));
					w.writeAttribute("height", format(imageSize.height));
					w.writeAttribute("style", "fill:"+opt.get());
					}

				
				
				w.writeStartElement("g");
				w.writeAttribute("transform", "translate("+format(imageRadius+margin.top)+","+format(imageRadius+margin.left)+")");
				w.writeStartElement("g");
				
				final double contig_height=10;
				PlotContig plotContigs = new PlotContig();
				plotContigs.radius_max = radius;
				plotContigs.radius_min = radius - contig_height;
				plotContigs.plot(w);
				radius = plotContigs.radius_min;
				
				
				if(tracks!=null ) {
					GroupTrack group= new GroupTrack(tracks);
					group.radius_max= radius;
					group.radius_min= min_radius;
					group.plot(w);
					}
				
				
				
				w.writeEndElement();//g
				w.writeEndElement();//g
				
				//at the end, write frame
				// write background panel
				w.writeEmptyElement("rect");
				w.writeAttribute("x", format(0));
				w.writeAttribute("y", format(0));
				w.writeAttribute("width", format(imageSize.width));
				w.writeAttribute("height", format(imageSize.height));
				w.writeAttribute("style", "fill:none;stroke:black;");
					
				
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
