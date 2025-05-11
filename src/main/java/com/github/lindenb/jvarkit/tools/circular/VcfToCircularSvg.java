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

import java.awt.Color;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.rdf.ns.XLINK;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

## input



END_DOC
 */
@Program(name="vcf2circularsvg",
description="displays a VCF as a circular SVG",
keywords= {"genome","browser","circular","vcf","svg"},
modificationDate="20220522",
creationDate="20220522"
)

public class VcfToCircularSvg extends Launcher{
	private static final Logger LOG = Logger.of(VcfToCircularSvg.class);

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names="--include",description="include chromosomes matching the following expression")
	private String includeContigsRegex = "(chr)?[0-9XY]+";

	@Parameter(names="--handler",description="handler's name")
	private String handlerName = "default";
	@DynamicParameter(names="-D",description="handler's name")
	private Map<String,String> __properties = new HashMap<>();

	private AttributeMap properties = null;
	private final ColorUtils colorUtils = new ColorUtils();
	private final DecimalFormat niceIntFormat = new DecimalFormat("###,###");
	private Color scoreColorStart = Color.WHITE;
	private Color scoreColorEnd = Color.BLACK;

	
	
	private static final Pattern RGB_PATTERN = Pattern.compile("\\d+,\\d+,\\d+");

	private class Arc extends SimpleInterval{
		Arc(VariantContext ctx) {
			super(ctx);
			}
		}
	private class ArcScore extends Arc{
		double score;
		ArcScore(VariantContext ctx,double score) {
			super(ctx);
			this.score = score;
			}
		}
	
	
	private byte toStrand(final String strandStr) {
		if(strandStr.equals("+")) {
			return (byte)1;
			}
		else if(strandStr.equals("-")) {
			return  (byte)-1;
		}else  {
			return (byte)0;
			} 
	}
	
	private String toCss(final String s) {
		if(StringUtil.isBlank(s)) return null;
		if(s.contains(":") || s.contains(";"))
			{
			return s.trim();
			}
		if(!RGB_PATTERN.matcher(s).matches()) return null;
		final String tokens[] = CharSplitter.COMMA.split(s);
		if(tokens.length!=3) return null;
		
		try 
			{
			final int r = Integer.parseInt(tokens[0]);
				final int g = Integer.parseInt(tokens[1]);
				final int b = Integer.parseInt(tokens[2]);
				if(r>=0 && r<256 && g>=0 && g<256 && b>=0 && b<256)
					{
					return "fill:rgb("+r+","+g+","+b+")";
					}
			}
		catch(final NumberFormatException err)
			{
			return null;
			}
		return null;
		}
	
	private Color between(final Color i,final Color j,float f)
		{
		if(f<0) f=0f;
		if(f>1) f=1f;
		final int r= i.getRed()+(int)((j.getRed()-i.getRed())*f);
		final int g= i.getGreen()+(int)((j.getGreen()-i.getGreen())*f);
		final int b= i.getBlue()+(int)((j.getBlue()-i.getBlue())*f);
		return new Color(r,g,b);
		}
	
	private static class Point2D {
		double x;
		double y;
		Point2D(double x,double y) {
			this.x = x;
			this.y = y;
		}
	}
	
	private class Sample {
		final String sn;
		final int index;
		Color color = Color.BLACK;
		final List<Arc> arcs = new ArrayList<>();
		final Map<String,List<List<Arc>>> contig2rows = new HashMap<>();
		Sample(final String sn,int index) {
			this.sn = sn;
			this.index = index;
			}
		
		int getMaxRows() {
			return this.contig2rows.values().stream().mapToInt(L->L.size()).max().orElse(0);
			}
		
		void pileup() {
			int within_distance = 10;
			final Set<String> contigs = this.arcs.stream().map(S->S.getContig()).collect(Collectors.toSet());
			for(final String ctg: contigs) {
				final List<Arc> L = this.arcs.stream().filter(S->S.getContig().equals(ctg)).collect(Collectors.toCollection(ArrayList::new));
				Collections.sort(L,(A,B)->Integer.compare(A.getStart(), B.getStart()));
				final List<List<Arc>> rows = new ArrayList<>();
				while(!L.isEmpty()) {
					final Arc curr = L.remove(0);
					int y=0;
					for(y=0;y < rows.size();++y) {
						final List<Arc> row = rows.get(y);
						if(row.stream().anyMatch(A->A.withinDistanceOf(curr,within_distance))) continue;
						row.add(curr);
						break;
						}
					if(y==rows.size()) {
						final List<Arc> row = new ArrayList<>();
						row.add(curr);
						rows.add(row);
						}
					}
				this.contig2rows.put(ctg,rows);
				this.arcs.removeIf(S->S.getContig().equals(ctg));
				}
			}
		}
	
	private abstract class Handler {
		private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
		final SAMSequenceDictionary dict;
		final long reference_length;
		private long tid_to_start[];
		final List<Sample> samples = new ArrayList<>();
		final Map<String,Sample> name2sample = new HashMap<>();
		final Pattern contigRegex;
		XMLStreamWriter w;
		
		protected class ArcDrawer {
			String title= "";
			String url = "";
			final Locatable loc;
			int strand = 0;
			double radius1;
			double radius2;
			ArcDrawer(final Locatable loc) {
				this.loc = loc;
				}
			
			private Point2D[] points() {
				final double r_start = loc2rad(loc.getContig(),loc.getStart());
				final double r_end = loc2rad(loc.getContig(),loc.getEnd());
				final double arc_length = (r_end-r_start)*Math.PI;
				final double angle_arrow= 1;//this.arrow_size / mid_radius;

				final Point2D p1 = polarToCartestian(radius1, r_start);
				final Point2D p2 = polarToCartestian(radius1, r_end);
				final Point2D p3 = polarToCartestian(radius2, r_end);
				final Point2D p4 = polarToCartestian(radius2, r_start);
				return new Point2D[] {p1,p2,p3,p4};
				}
			
			String makePath() {
				switch(this.strand) {
					case -1: return makeNegativePath();
					case  1: return makePositivePath();
					default: return makeUnstrandedPath();
					}
				}
			
			private String makeNegativePath() {
				final double mid_radius= (radius1+radius2)/2.0;
				final Point2D[] p = points();
				final StringBuilder sb = new StringBuilder();
				final double r_start = loc2rad(loc.getContig(),loc.getStart());
				final double angle_arrow= 1;

				sb.append("M ").append(pointToStr(polarToCartestian(mid_radius, r_start)));
				sb.append("L ").append(pointToStr(polarToCartestian(radius1, r_start+angle_arrow)));
				
				sb.append(" A ").
					append(format(radius1)).
					append(" ").
					append(format(radius1)).
					append(" 0"). //X axis rotation
					append(" 0").// large arc
					append(" 1 ").// sweep flag (positive angle direction)
					append(pointToStr(p[1]));

				sb.append("L ").append(pointToStr(p[2]));	
				
				sb.append(" A ").
					append(format(radius2)).
					append(" ").
					append(format(radius2)).
					append(" 0"). //X axis rotation
					append(" 0").// large arc
					append(" 0 ").// sweep flag (positive angle direction)
					append(pointToStr(polarToCartestian(radius2, r_start+angle_arrow)));
				
				return sb.append(" Z").toString();
				}
			private String makePositivePath() {
				final double mid_radius= (radius1+radius2)/2.0;
				final Point2D[] p = points();
				final StringBuilder sb = new StringBuilder();
				final double r_end = loc2rad(loc.getContig(),loc.getEnd());
				final double angle_arrow= 1;

				sb.append("M ").append(pointToStr(p[0]));
				
				sb.append(" A ").
					append(format(radius1)).
					append(" ").
					append(format(radius1)).
					append(" 0"). //X axis rotation
					append(" 0").// large arc
					append(" 1 ").// sweep flag (positive angle direction)
					append(pointToStr(polarToCartestian(radius1, r_end - angle_arrow)));
				
				sb.append("L ").append(pointToStr(polarToCartestian(mid_radius, r_end)));
				sb.append("L ").append(pointToStr(polarToCartestian(radius2, r_end-angle_arrow)));	
				
				sb.append(" A ").
					append(format(radius2)).
					append(" ").
					append(format(radius2)).
					append(" 0"). //X axis rotation
					append(" 0").// large arc
					append(" 0 ").// sweep flag (positive angle direction)
					append(pointToStr(p[3]));
				
				return sb.append(" Z").toString();
				}
			private String makeUnstrandedPath() {
				final Point2D[] p = points();
				final StringBuilder sb = new StringBuilder();
				sb.append("M ").append(pointToStr(p[0]));
				
				sb.append(" A ").
					append(format(radius1)).
					append(" ").
					append(format(radius1)).
					append(" 0"). //X axis rotation
					append(" 0").// large arc
					append(" 1 ").// sweep flag (positive angle direction)
					append(pointToStr(p[1]));
				
				sb.append(" L").append(pointToStr(p[2]));
				
				sb.append(" A ").
					append(format(radius2)).
					append(" ").
					append(format(radius2)).
					append(" 0"). //X axis rotation
					append(" 0").// large arc
					append(" 0 ").// sweep flag (positive angle direction)
					append(pointToStr(p[3]));
				
				sb.append(" L ").append(pointToStr(p[0]));
				return sb.append(" Z").toString();
				}
			
			void write() throws XMLStreamException{
				final String d = makePath();
				if(StringUtils.isBlank(d)) return;
				if(!StringUtils.isBlank(this.url)) {
					w.writeStartElement("a");
					w.writeAttribute("href", this.url);
					}
				
				w.writeStartElement("path");
				w.writeAttribute("d",d);
				w.writeEndElement();//path;
				if(!StringUtils.isBlank(this.url)) {
					w.writeEndElement();
					}
				}
			} /** end ArcDrawer */
		
		Handler(final VCFHeader header) {
			final SAMSequenceDictionary dict0  = SequenceDictionaryUtils.extractRequired(header);
			this.contigRegex = Pattern.compile(VcfToCircularSvg.this.includeContigsRegex);
			
			this.dict = new SAMSequenceDictionary(dict0.getSequences().stream().
				filter(SR->this.contigRegex.matcher(SR.getSequenceName()).matches()).
				collect(Collectors.toList()));
			this.reference_length = this.dict.getReferenceLength();
			
			this.tid_to_start = new long[this.dict.size()];
			Arrays.fill(this.tid_to_start,0L);
			
			
			long n = 0;
			for(int i=0;i< this.dict.size();i++) {
				this.tid_to_start[i] = n;
				n += this.dict.getSequence(i).getSequenceLength();
				}
			
			final List<String> snL = header.getGenotypeSamples();
			for(int i=0;i< snL.size();i++) {
				final Sample sn = new Sample(snL.get(i),i);
				sn.color = i%2==0?Color.BLUE:Color.CYAN;
				this.samples.add(sn);
				name2sample.put(sn.sn,sn);
				}
			
			
			}
		
		void beginDocument() throws XMLStreamException {
			this.w.writeStartDocument("UTF-8", "1.0");
			this.w.writeStartElement("svg");
			this.w.writeDefaultNamespace(SVG.NS);
			String r= String.valueOf(Math.ceil(getImageRadius()));
			w.writeAttribute("width",r);
			w.writeAttribute("height",r);
			
			writeStyle();
			writeScript();
			writeTitle();
			}
		
		void writeStyle() throws XMLStreamException {
			}
		void writeScript() throws XMLStreamException {
			}
		void writeTitle() throws XMLStreamException {
			}
		
		void visit(final VariantContext ctx) {
			
			}
		
		double getMinRadius() {
			return 100;
			}
		double getDataRadius() {
			return 1000;
			}
		
		double getImageRadius() {
			return getMinRadius()+getDataRadius()+100;
			}
		
		void title(String title) throws XMLStreamException {
			if(StringUtils.isBlank(title)) return;
			w.writeStartElement("title");
			w.writeCharacters(title);
			w.writeEndElement();
			}
		
		void finish() throws XMLStreamException {
			this.w.writeEndElement();//svg
			this.w.writeEndDocument();
			this.w.flush();
			this.w.close();
			}
		

		private String pointToStr(final Point2D p)
			{
			return format(p.x)+" "+format(p.y);
			}
		
		private double loc2rad(final String contig,final int pos) {
			final SAMSequenceRecord ssr = this.dict.getSequence(contig);
			if(ssr==null) throw new IllegalStateException();
			final int tid  = ssr.getSequenceIndex();
			final long index_at_start = this.tid_to_start[tid];
			return ((index_at_start+pos)/(double)this.reference_length)*(2.0*Math.PI);
			}
		
		
		
		private Point2D polarToCartestian(double radius,double angle)
			{
			return new Point2D(
					Math.cos(angle)*radius,
					Math.sin(angle)*radius
					);
			}

		private String format(double v) {
			return this.decimalFormater.format(v);
			}

		private void writeContigs() throws XMLStreamException {
			w.writeStartElement("g");
			radius+=this.feature_height;
			for(int tid=0;tid< dict.size();++tid)
				{
				final SAMSequenceRecord ssr = this.dict.getSequence(tid);
				w.writeStartElement("path");
				w.writeAttribute("class", "contig"+(tid%2));
				w.writeAttribute("d", arc(tid,radius,radius+this.contig_height,0,ssr.getSequenceLength(),(byte)0));

				w.writeStartElement("title");
				w.writeCharacters(ssr.getSequenceName());
				w.writeEndElement();
				
				w.writeEndElement();//path
				
				w.writeStartElement("text");
				w.writeAttribute("class", "contig");
				w.writeAttribute("x", format(radius+this.distance_between_arc+this.contig_height));
				w.writeAttribute("y", "0");
				w.writeAttribute("transform","rotate("+format(Math.toDegrees(
						loc2rad(ssr.getContig(),ssr.getSequenceLength()/2)
						))+")");
				w.writeCharacters(ssr.getSequenceName());
				w.writeEndElement();
				}
			w.writeEndElement();//g
			}
		abstract void finish(XMLStreamWriter w) throws XMLStreamException;
		}
	
	private class DellyCNVHandler extends Handler {
		DellyCNVHandler(final VCFHeader header){
				super(header);
				

				
				}
		@Override
		void visit(final VariantContext ctx) {
			for(Sample sn:this.samples) {
				final Genotype gt = ctx.getGenotype(sn.sn);
				if(gt==null) continue;
				if(!gt.hasAnyAttribute("CN")) continue;
				int cn = gt.getAttributeAsInt("CN", 2);
				if(cn==2) continue;
				sn.arcs.add(new ArcScore(ctx,cn));
				}
			}
		void finish(XMLStreamWriter w) throws XMLStreamException {
			this.w = w;
			for(final Sample sn:this.samples) {
				sn.pileup();
				}
			int numRows = this.samples.stream().mapToInt(SN->SN.getMaxRows()).sum();
			double rowHeight = getDataRadius()/numRows;
			double x0 = getMinRadius();
			for(int i=0;i< this.samples.size();i++)  {
				final Sample sample = this.samples.get(i);
				int nRows = sample.getMaxRows();
				if(nRows==0) continue;
				}
			beginDocument();
			}
		}
	
	private Handler createHandler(final VCFHeader header) {
		if(this.handlerName.equalsIgnoreCase("dellycnv1")) {
			return new DellyCNVHandler(header);
			}
		throw new IllegalArgumentException("cannot create handler: "+this.handlerName);
		}
	
	@Override
	public int doWork(final List<String> args) {
		this.properties = AttributeMap.wrap(this.__properties);
		try{
			final String input = oneFileOrNull(args);
			try(VCFIterator iter = super.openVCFIterator(input)) {
				final VCFHeader header = iter.getHeader();
				final Handler handler = createHandler(header);
				while(iter.hasNext()) {
					final VariantContext ctx = iter.next();
					if(!ctx.getContig().matches(this.includeContigsRegex)) continue;
					handler.visit(ctx);
					}
				try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					final XMLOutputFactory xof = XMLOutputFactory.newInstance();
					final XMLStreamWriter w = xof.createXMLStreamWriter(out);
					handler.finish(w);
					out.flush();
					}
				}
			return 0;
			} 
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
	}
	
	

	
	
	
	public static void main(final String[] args) {
		new VcfToCircularSvg().instanceMainWithExit(args);
	}
}
