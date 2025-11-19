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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.namespace.QName;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.ExtendedLocatable;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.tools.structvar.indexcov.IndexCovUtils;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;
import com.github.lindenb.jvarkit.xml.DOMUtils;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;


/**
BEGIN_DOC

## input



END_DOC
 */
@Program(name="circular",
	description="circular genome browser renderer",
	keywords= {"genome","browser","circular","vcf","svg"},
	modificationDate="20251109",
	creationDate="20220522",
	generate_doc = false,
	jvarkit_hidden = true,
	jvarkit_amalgamion = true
	)
public class Circular extends Launcher {
	private static final Logger LOG = Logger.of(Circular.class);
	@Parameter(names="-o",description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output = null;
	@Parameter(names="-R",description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path override_faidx = null;
	
	private XPath xPathInstance=null;
	
	private static String format(final Object s) {
		return  String.valueOf(s);
		}
	
	
	private class RadialAndLocatable implements ExtendedLocatable {
		final ContigInfo contigInfo;
		final int start;
		int end; //updatable
		double outer_radius =0;
		double inner_radius =0;
		
		RadialAndLocatable(ContigInfo contigInfo,int start,int end) {
			this.contigInfo = contigInfo;
			this.start = start;
			this.end = end;
			}
		@Override
		public String getContig() {
			return contigInfo.getContig();
			}
		@Override
		public int getStart() {
			return start;
			}
		@Override
		public int getEnd() {
			return start;
			}
		}
	
	/**
	 * 
	 * A Track in the browser
	 *
	 */
	private abstract class Track {
		double radius_max;
		double radius_min;
		abstract void plot(XMLStreamWriter w) throws IOException,XMLStreamException,XPathException;
		
		void writeTitle(XMLStreamWriter w,Object t)throws IOException,XMLStreamException {
			if(t==null) return;
			w.writeStartElement("title");
			w.writeCharacters(String.valueOf(t));
			w.writeEndElement();
			}
		 void pack() {
			 
		 }
		double getWeigth() {
			return 1.0;
			}
		}
	
	private  class PlotContigTrack extends Track{
		@Override
		void plot(XMLStreamWriter w) throws IOException, XMLStreamException,XPathException {
			w.writeStartElement("g");
			for(ContigInfo ci:Circular.this.contigs) {
				w.writeStartElement("g");
				
				w.writeStartElement("path");
				w.writeAttribute("style","fill:"+(ci.tid%2==0?"rgb(250,250,250)":"rgb(240,240,240)")+";stroke:rgb(230,230,230)");
				w.writeAttribute("d",ci.arc(0,this.radius_max,1,ci.getLengthOnReference()));
				w.writeEndElement();
				
				
				w.writeStartElement("path");
				w.writeAttribute("style","fill:"+(ci.tid%2==0?"red":"blue")+";stroke:black");
				w.writeAttribute("d",ci.arc(this.radius_min,this.radius_max,1,ci.getLengthOnReference()));
				w.writeEndElement();
				
				w.writeStartElement("text");
				w.writeStartElement("textPath");
				w.writeAttribute("path", ci.arc(this.radius_min,this.radius_max,ci.getLengthOnReference()/2,ci.getLengthOnReference()/2));
				w.writeCharacters(ci.getContig());
				w.writeEndElement();

				w.writeEndElement();
				
				writeTitle(w, ci.getContig());
				w.writeEndElement();//g
				}
			w.writeEndElement();
			}
		}
	private abstract class ConfigurableTrack extends Track{
		protected final Element root;
		ConfigurableTrack(final Element root) {
			this.root=root;
			}
		@Override
		double getWeigth() {
			return 1.0;
			}
		protected Element resolvePath(String...path) {
			return Circular.this.resolvePath(this.root,path);
			}
		protected Element resolveNode(final List<String> path) {
			return Circular.this.resolveNode(this.root,path);
			}
		}
	
	
	private  class AbstractPileupTrack extends ConfigurableTrack {
		AbstractPileupTrack(final Element root) {
			super(root);
			}
		@Override
		void plot(XMLStreamWriter w) throws IOException, XMLStreamException, XPathException {
			w.writeStartElement("g");
			for(ContigInfo ci:Circular.this.contigs) {
				
				}
			w.writeEndElement();
			}
		}
	
	
	private  class GroupTrack extends ConfigurableTrack{
		final List<Track> subTracks=new ArrayList<>();
		GroupTrack(final Element root) throws XPathException {
			super(root);
			for(Node c=root.getFirstChild();c!=null;c=c.getNextSibling()) {
				if(c.getNodeType()!=Node.ELEMENT_NODE) continue;
				Element e= Element.class.cast(c);
				Track track=createTrack(e);
				
				if(track==null) {
					LOG.warn("track is null");
					continue;
					}
				this.subTracks.add(track);
				}
			}
		
		@Override
		void pack() {
			final double total_w = subTracks.stream().mapToDouble(T->T.getWeigth()).sum();
			double r = super.radius_max;
			LOG.info("radius p:"+this.radius_max+"-"+this.radius_min);

			for(Track sub: subTracks) {
				sub.radius_max = r;
				sub.radius_min = r - (sub.getWeigth()/total_w)*( super.radius_max- super.radius_min);
				r= sub.radius_min;
				
				double dr= (sub.radius_max-sub.radius_min)*0.05;
				sub.radius_max-=dr;
				sub.radius_min+=dr;
				
				LOG.info("radius "+sub.radius_max+"-"+sub.radius_min);
				}
			
			super.pack();
			}
		
		
		double getWeigth() {
			return 1.0;
			}
		@Override
		void plot(final XMLStreamWriter w) throws IOException, XMLStreamException,XPathException {
			LOG.info("plotxx"+this.subTracks);
			for(Track t:this.subTracks) {
				LOG.info("plot");
				t.plot(w);
				}
			
			
			 
			}
		}

	
	private abstract class AbstractIndexCov extends ConfigurableTrack {
		class Item extends RadialAndLocatable {
			IndexCovUtils.SvType type;
			Item(ContigInfo ci,int start,int end) {
				super(ci,start,end);
				}
			}
		AbstractIndexCov(final Element root) {
			super(root);
			}
		Path getPath() {
			String p  = valueOf(super.resolvePath("bed")).orElse(null);
			if(StringUtils.isBlank(p)) throw new IllegalArgumentException("vcf missing in "+DOMUtils.getNodePath(root));
			if(!(p.endsWith(".bed") || p.endsWith(".bed.gz"))) throw new IllegalArgumentException("doesnt look like a bed "+DOMUtils.getNodePath(root));
			return Paths.get(p);
			}
		Predicate<BedLine> createFilter() {
			return V->true;
			}
		@Override
		void plot(XMLStreamWriter w) throws IOException, XMLStreamException, XPathException {
			Path bedPath = getPath();
			FileHeader header = null;
			try(BufferedReader br= IOUtils.openPathForBufferedReading(bedPath) ) {
				String line= br.readLine();
				if(StringUtils.isBlank(line)) throw new IllegalArgumentException("Cannot read first line of "+bedPath);
				header  = new FileHeader(line, CharSplitter.TAB);
				}
			double treshold = Double.parseDouble(valueOf(this.resolvePath("treshold")).orElse(String.valueOf(IndexCovUtils.DEFAULT_TRESHOLD)));
			int  min_size = Integer.parseInt(valueOf(this.resolvePath("min-size")).orElse(String.valueOf("1000")));
			int n_samples = header.size()-3;
			double d_radius = (this.radius_max - this.radius_min)/n_samples;
			double outer_radius = this.radius_max;
			IndexCovUtils utils = new IndexCovUtils(treshold);
			for(int col=3;col< header.size();col++) {
				String sn = header.get(col);
				w.writeStartElement("g");
				final List<Item> items= new ArrayList<>();
				try(BufferedReader br= IOUtils.openPathForBufferedReading(bedPath) ) {
					for(;;) {
						String line= br.readLine();
						if(line==null) break;
						if(line.startsWith("#")) continue;
						List<String> tokens = header.split(line);
						final ContigInfo ci =  name2contig.get(tokens.get(0));
						if(ci==null) continue;
						float value  = Float.parseFloat(tokens.get(col));
						IndexCovUtils.SvType type= utils.getType(value);
						if(type.isAmbigous() || type.isReference()) continue;
						if(!utils.isVariant(value)) continue;
						// looks like centromer
						final double hom_del  = 0.01;
						if(value< hom_del && n_samples>10 &&  tokens.subList(3, tokens.size()).stream().mapToDouble(Float::parseFloat).allMatch(V->V< hom_del)) {
							continue;
							}
						
						Item item = new Item(ci,Integer.parseInt(tokens.get(1)),Integer.parseInt(tokens.get(2)));
						
						item.type = type;
						if(!items.isEmpty()) {
							Item last = items.get(items.size()-1);
							if(last.contigInfo==item.contigInfo && last.type==item.type && last.end+IndexCovUtils.BAI_BLOCK_SIZE>=item.start) {
								last.end = item.end;
								continue;
								}
							}
						items.add(item);
						}
					}
		 		double inner_radius = outer_radius - d_radius * 0.95;

				items.removeIf(V->V.getLengthOnReference() < min_size);
				for(Item item:items) {
					String fill;
					switch(item.type) {
						case  HET_DEL : fill="red";break;
						case HOM_DEL : fill="black";break;
						case HET_DUP : fill="green";break;
						case HOM_DUP : fill="blue"; break;
						default : fill="gray"; break;
						}
					w.writeStartElement("path");
					w.writeAttribute("style", "fill:"+fill+";stroke:"+(java.lang.Math.random()<0.5?"blue":"red")+";");
					w.writeAttribute("d", item.contigInfo.arc(inner_radius,outer_radius,item.start,item.start));
					w.writeEndElement();
					}
				w.writeEndElement();
				outer_radius -= d_radius;
				}
			
			}

		}
	
	private abstract class AbstractSVGTrack extends ConfigurableTrack{
		AbstractSVGTrack(final Element root) {
			super(root);
			}
		Path getPath() {
			String p = (String)evalXPath("vcf", super.root, XPathConstants.STRING);
			if(StringUtils.isBlank(p)) throw new IllegalArgumentException("vcf missing in "+DOMUtils.getNodePath(root));
			return Paths.get(p);
			}
		Predicate<VariantContext> createFilter() {
			return V->true;
			}
		@Override
		void plot(XMLStreamWriter w) throws IOException,XMLStreamException,XPathException {
			Path p = getPath();
			Predicate<VariantContext> filter=createFilter();
			try(VCFReader r = VCFReaderFactory.makeDefault().open(p,true)) {
				for(int ntry=0;ntry<2;++ntry) {
					
					for(ContigInfo ci:Circular.this.contigs) {
						try(CloseableIterator<VariantContext> iter=r.query(ci)) {
							while(iter.hasNext()) {
								VariantContext ctx = iter.next();
								if(!filter.test(ctx)) continue;
								}
							}
						}
					if(ntry==1) {
						w.writeStartElement("g");
						w.writeEndElement();
						}
					}
				}
			}
		}
		
	
	private abstract class AbstractBamTrack extends ConfigurableTrack{
		AbstractBamTrack(final Element root) {
			super(root);
			}
		 Predicate<SAMRecord> createReadFilter() {
			Predicate<SAMRecord> p =R->!R.getReadUnmappedFlag();
			Optional<String> opt= valueOf(this.resolvePath("mapq"));
			if(opt.isPresent()) {
				int mapq = Integer.parseInt(opt.get());
				p =p.and(R->R.getMappingQuality()>=mapq);
				}
			return p;
		 	}
		 
		 Optional<Path> getReferencePath()  {
			String s= (String)Circular.this.evalXPath("fasta",this.root,XPathConstants.STRING);
			if(!StringUtils.isBlank(s)) return Optional.of(Paths.get(s.trim()));
			if(Circular.this.override_faidx!=null)  return Optional.of(Circular.this.override_faidx);
			s= (String)Circular.this.evalXPath("/confif/reference/fasta",this.root,XPathConstants.STRING);
			if(!StringUtils.isBlank(s)) return Optional.of(Paths.get(s));
			return Optional.empty();
		 	}
		 
		 SamInputResource getSamInputResource() {
			String s= (String)Circular.this.evalXPath("bam",this.root,XPathConstants.STRING);
			if(!StringUtils.isBlank(s)) return SamInputResource.of(s.trim());
			throw new IllegalArgumentException("bam undefined under "+ DOMUtils.getNodePath(this.root));
		 	}
		 
		 SamReader open() throws IOException {
			final SamReaderFactory srf = SamReaderFactory.make();
			final SamInputResource si = getSamInputResource();
			final Optional<Path> ref = getReferencePath();
			srf.validationStringency(ValidationStringency.LENIENT);
			if(!ref.isEmpty()) srf.referenceSequence(ref.get());
			return srf.open(si);
			}
		}
	
	private class CoverageBamTrack extends AbstractBamTrack {
		private class Count {
			long total=0L;
			}
		CoverageBamTrack(Element root) {
			super(root);
			}

		
		@Override
		void plot(XMLStreamWriter w) throws IOException,XMLStreamException,XPathException {
			LOG.info("write here");
			final long genome_length = Circular.this.contigs.stream().mapToLong(CI->CI.getLengthOnReference()).sum();
			try(SamReader sr=open()) {
				final SAMFileHeader hdr= sr.getFileHeader();
				final String sn  = hdr.getReadGroups().
						stream().
						filter(RG->!StringUtils.isBlank(RG.getSample())).
						map(RG->RG.getSample()).
						findAny().
						orElse("BAM")
						;
				
				w.writeStartElement("g");
				for(ContigInfo ci:Circular.this.contigs) {
					LOG.info("arc_length="+ci.getLengthOnReference());
					if(hdr.getSequenceDictionary().getSequence(ci.getContig())==null) continue;
					final double contig_arc_length = ((Math.PI*2*this.radius_max)/genome_length)*ci.getLengthOnReference();
					final double unit_arc_length = 10;
					final int bin_size_bp = (int) Math.max(1,(unit_arc_length/contig_arc_length)*ci.getLengthOnReference());
					final List<Count> counts = new ArrayList<>(1+ci.getLengthOnReference()/bin_size_bp);
					
					final Predicate<SAMRecord> filter = createReadFilter();
					try(CloseableIterator<SAMRecord> iter = sr.query(ci.getContig(), 1, ci.getLengthOnReference(), true)) {
						while(iter.hasNext()) {
							final SAMRecord rec = iter.next();
							if(!filter.test(rec)) continue;
							for(AlignmentBlock ab: rec.getAlignmentBlocks()) {
								for(int i=0;i< ab.getLength();i++) {
									final int pos1 = ab.getReferenceStart()+i;
									if(pos1>= ci.getLengthOnReference()) break;
									final int bin_idx = pos1/bin_size_bp;
									while(bin_idx>= counts.size()) {
										counts.add(new Count());
										}
									final Count c = counts.get(bin_idx);
									c.total++;
									}
								}
							final double max_cov= 10;
							for(int i=0;i< counts.size();++i) {
								final Count count = counts.get(i);
								double cov = Math.min(max_cov, count.total/(double)bin_size_bp);
								w.writeStartElement("path");
								double h = (cov/max_cov)*(this.radius_max-this.radius_min);
								w.writeAttribute("style", "fill:"+(i%2==0?"pink":"green")+";stroke:orange;");
								w.writeAttribute("d", ci.arc(this.radius_min, this.radius_min+h, i*bin_size_bp, (i+1)*bin_size_bp));
								
								w.writeStartElement("title");
								w.writeCharacters(String.valueOf(count.total/(double)bin_size_bp));
								w.writeEndElement();
								
								w.writeEndElement();
								}
								
							
							}
						
						}
					}
				w.writeStartElement("title");
				w.writeCharacters(sn);
				w.writeEndElement();
				
				w.writeEndElement();
				}//try samReader
			}//plot
		}
	
	private class PileupBamTrack extends AbstractBamTrack {
		
		PileupBamTrack(Element root) {
			super(root);
			}
		
		
		@Override
		void plot(XMLStreamWriter w) throws IOException,XMLStreamException,XPathException {
			LOG.info("write here");
			final long genome_length = Circular.this.contigs.stream().mapToLong(CI->CI.getLengthOnReference()).sum();
			try(SamReader sr=open()) {
				final SAMFileHeader hdr= sr.getFileHeader();
				final String sn  = hdr.getReadGroups().
						stream().
						filter(RG->!StringUtils.isBlank(RG.getSample())).
						map(RG->RG.getSample()).
						findAny().
						orElse("BAM")
						;
				
				w.writeStartElement("g");
				final Predicate<SAMRecord> filter = createReadFilter();
				
				w.writeStartElement("text");
				w.writeAttribute("style", "alignment-baseline:middle");
				w.writeStartElement("textPath");
				ContigInfo k1=Circular.this.contigs.get(0);
				w.writeAttribute("path", k1.arc(this.radius_max,this.radius_max,1,k1.getLengthOnReference()));
				w.writeCharacters("LEGEND!!");
				w.writeEndElement();
				w.writeEndElement();
				
				
				for(ContigInfo ci:Circular.this.contigs) {
					 final List<RadialAndLocatable> items= new ArrayList<>();
					 try(CloseableIterator<SAMRecord> iter = sr.query(ci.getContig(), 1, ci.getLengthOnReference(), true)) {
						 while(iter.hasNext()) {
							 final SAMRecord rec = iter.next();
								if(!filter.test(rec)) continue;
								final RadialAndLocatable item = new RadialAndLocatable(ci,rec.getStart(),rec.getEnd());
								items.add(item);
						 		}
							}
					 	
					 if(!items.isEmpty()) {
						w.writeStartElement("g");
					 	Collections.sort(items,(A,B)->Integer.compare(A.getStart(), B.getStart()));
					 	Pileup<RadialAndLocatable> pileup = new Pileup<>((left,right)->{
					 		//use the smallest radius for this track
					 		double left_a = ci.pos2angle(left.getEnd()+1);
					 		double right_a = ci.pos2angle(right.getStart());
					 		double radial_dist = PileupBamTrack.this.radius_min*(right_a-left_a);
					 		return radial_dist > 1 ;
					 		});
					 	final List<List<RadialAndLocatable>> rows= pileup.addAll(items).getRows();
					 	
					 	double dr = (this.radius_max-this.radius_min)/rows.size();
					 	double r = this.radius_max;
					 	for(List<RadialAndLocatable> row: rows) {
					 		double outer_radius = r;
					 		double inner_radius = outer_radius - dr * 0.95;
					 		for(RadialAndLocatable item:row) {
					 			item.outer_radius = outer_radius;
					 			item.inner_radius = inner_radius;
					 			w.writeStartElement("path");
								w.writeAttribute("style", "fill:"+(java.lang.Math.random()<0.5?"pink":"green")+";stroke:"+(java.lang.Math.random()<0.5?"blue":"red")+";");
								w.writeAttribute("d", ci.arc(inner_radius,outer_radius,item.getStart(),item.getEnd()));
								
								
								
								w.writeEndElement();
					 			}
					 		
					 		r  -= dr;
					 		}
					 	w.writeEndElement();
					 	}
						
					}
				w.writeStartElement("title");
				w.writeCharacters(sn);
				w.writeEndElement();
				
				w.writeEndElement();
				}//try samReader
			}//plot
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
			return (angle_start+ (genomic_pos1/(double)getLengthOnReference())*(angle_end-angle_start)) - Math.PI/2.0;
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
					"M",
					pointToStr(toPoint(start,r_min)),
					"L",
					pointToStr(toPoint(start,r_max)),
					arcTo(r_max, start, end),
					"L",
					pointToStr(toPoint(end,r_min)),
					arcTo(r_min, end,start)
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

	private Object evalXPath(final String path,Node node, QName what) {
		try 
			{
			
			return  this.xPathInstance.evaluate(path, node, what);
			} 
		catch(XPathException err) {
			throw new RuntimeException(err);
			}
		}
	
	private Optional<String> valueOf(Element root) {
		if(root==null) return Optional.empty();
		String s= root.getTextContent();
		if(StringUtils.isBlank(s)) return Optional.empty();
		return Optional.of(s);
		}
	
	private Element resolvePath(Element root,String...paths) {
		List<String> array = Arrays.asList(paths);
		if(array.isEmpty()) throw new IllegalArgumentException();
		String path= array.get(0);
		if(path.startsWith("/")) throw new IllegalArgumentException();
		if(array.size()==1 && path.contains("/")) {
			array = CharSplitter.of("/").splitAsStringList(path);
			}
		return resolveNode(root, array);
		}
	
	private Element resolveNode(Element root,final List<String> path) {
		if(path.isEmpty()) return root;
		Attr att = root.getAttributeNode("ref");
		if(att==null)  att = root.getAttributeNode("href");
		if(att!=null) {
			String node_id= att.getValue();
			if(node_id.startsWith("#")) node_id=node_id.substring(1).trim();
			if(root.hasChildNodes()) {
				LOG.warn("node has @ref but has children "+DOMUtils.getNodePath(root));
				}
			final Element defs = (Element)this.evalXPath("/config/defs", root.getOwnerDocument(), XPathConstants.NODE);
			if(defs==null) {
				LOG.warn("Cannot find /config/defs for "+ DOMUtils.getNodePath(att));
				return null;
				}
			for(Node c=defs.getFirstChild();c!=null;c=c.getNextSibling()) {
				if(c.getNodeType()!=Node.ELEMENT_NODE) continue;
				final Element E = Element.class.cast(c);
				if(!E.getNodeName().equals(path.get(0))) continue;
				if(!E.hasAttribute("id")) continue;
				if(!E.getAttribute("id").equals(node_id)) continue;
				if(path.size()==1) return E;
				return resolveNode(E, path.subList(1, path.size()));
				}
			LOG.warn("Cannot find <"+path.get(0)+" @id="+node_id+"/> under /config/defs for "+DOMUtils.getNodePath(att));
			return null;
			}
		for(Node c=root.getFirstChild();c!=null;c=c.getNextSibling()) {
			if(c.getNodeType()!=Node.ELEMENT_NODE) continue;
			final Element E = Element.class.cast(c);
			if(!E.getNodeName().equals(path.get(0))) continue;
			return resolveNode(E, path.subList(1, path.size()));
			}
		return null;
		}


	
	private Track createTrack(final Element root) throws XPathException {
		if(!root.getNodeName().equals("track")) {
			throw new IllegalArgumentException("exepected <track> but got "+DOMUtils.getNodePath(root));
			}
		final String track_type =(String)this.evalXPath( "type", root,XPathConstants.STRING);
		if(StringUtils.isBlank(track_type)) {
			throw new IllegalArgumentException("track is missing a type in "+ DOMUtils.getNodePath(root));
			}		
		
		
		if(track_type.equals("group")) {
			return new GroupTrack(root);
			}
		else if(track_type.equals("coverage")) {
			return new PileupBamTrack(root);
			//return new CoverageBamTrack(root);
			}
		throw new IllegalArgumentException("cannot create track for track_type="+track_type+" "+DOMUtils.getNodePath(root));
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
			
			
			
		
			
			final Element tracks = (Element)this.xPathInstance.evaluate("/config/tracks",config, XPathConstants.NODE);
			if(tracks==null) {
				LOG.warning("/config/tracks missing");
				}
			
			final Element defs = (Element)this.xPathInstance.evaluate("/config/defs",config, XPathConstants.NODE);
			if(defs==null) {
				LOG.warning("/config/defs missing");
				}
			
			
			SAMSequenceDictionary dict0 = null;
			if(this.override_faidx==null) {
				String fasta_path = (String)this.evalXPath("/config/reference/fasta/text()",config,XPathConstants.STRING);
				if(StringUtils.isBlank(fasta_path) ) fasta_path = (String)this.evalXPath("/config/reference/dict/text()",config,XPathConstants.STRING);
				if(StringUtils.isBlank(fasta_path)) {
					LOG.error("no /config/reference/fasta or /contig/reference/dict");
					return -1;
					}
				dict0 = new SequenceDictionaryExtractor().extractRequiredDictionary(fasta_path);
				}
			else
				{
				dict0 = new SequenceDictionaryExtractor().extractRequiredDictionary(this.override_faidx);
				}
			String opt;
			
			opt  =  (String)this.evalXPath("/config/reference/regex",config,XPathConstants.STRING);
			if(!StringUtils.isBlank(opt) ) {
				LOG.info("applying regex "+opt+" to chromosome names");
				final Pattern regex = Pattern.compile(opt.trim());
				dict0 = new SAMSequenceDictionary(dict0
						.getSequences()
						.stream()
						.filter(C->regex.matcher(C.getContig()).matches())
						.collect(Collectors.toList()));
				}
			opt = (String)this.evalXPath("/config/reference/min-length",config,XPathConstants.STRING);
			if(!StringUtils.isBlank(opt) ) {
				LOG.info("applying minLen "+opt+" to chromosome length");
				final int minLen = Integer.parseInt(opt);
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
			
			opt = (String)this.evalXPath("/contig/reference/margin",config,XPathConstants.STRING);
			
			final int default_margin = Integer.parseInt( StringUtils.ifBlank(opt, "10"));
			final Insets margin=new Insets(default_margin, default_margin, default_margin, default_margin);
			
			opt = (String)this.evalXPath("/contig/reference/margin-top",config,XPathConstants.STRING);
			if(!StringUtils.isBlank(opt)) margin.top=Integer.parseInt(opt);
			
			opt = (String)this.evalXPath("/contig/reference/margin-bottom",config,XPathConstants.STRING);
			if(!StringUtils.isBlank(opt)) margin.bottom=Integer.parseInt(opt);

			opt = (String)this.evalXPath("/contig/reference/margin-left",config,XPathConstants.STRING);
			if(!StringUtils.isBlank(opt)) margin.left=Integer.parseInt(opt);
			
			opt = (String)this.evalXPath("/contig/reference/margin-right",config,XPathConstants.STRING);
			if(!StringUtils.isBlank(opt)) margin.right=Integer.parseInt(opt);
			
			opt = (String)this.evalXPath("/contig/image/radius",config,XPathConstants.STRING);
			final int imageRadius=  Integer.parseInt( StringUtils.ifBlank(opt, "1000"));
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
				
				w.writeStartElement("defs");
				
				w.writeEndElement();//def
				
				
				w.writeStartElement("g");

				// write background panel
				opt = (String)this.evalXPath("/contig/image/background-color",config,XPathConstants.STRING);
				if(!StringUtils.isBlank(opt)) {
					w.writeEmptyElement("rect");
					w.writeAttribute("x", format(0));
					w.writeAttribute("y", format(0));
					w.writeAttribute("width", format(imageSize.width));
					w.writeAttribute("height", format(imageSize.height));
					w.writeAttribute("style", "fill:"+opt);
					}

				
				
				w.writeStartElement("g");
				w.writeAttribute("transform", "translate("+format(imageRadius+margin.top)+","+format(imageRadius+margin.left)+")");
				w.writeStartElement("g");
				
				final double contig_height=10;
				final PlotContigTrack plotContigs = new PlotContigTrack();
				plotContigs.radius_max = radius;
				plotContigs.radius_min = radius - contig_height;
				plotContigs.plot(w);
				radius = plotContigs.radius_min;
				
				
				if(tracks!=null ) {
					final GroupTrack group= new GroupTrack(tracks);
					group.radius_max= plotContigs.radius_min;
					group.radius_min= (group.radius_max)/3;
					group.pack();
					group.plot(w);
					LOG.info("done");
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
				w.flush();
				w.close();
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
