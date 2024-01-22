package com.github.lindenb.jvarkit.tools.genome2svg;

import java.awt.geom.Point2D;
import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.BiFunction;
import java.util.function.BiPredicate;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.springframework.context.ApplicationContext;
import org.springframework.context.support.FileSystemXmlApplicationContext;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.rdf.ns.XLINK;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptCodec;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.structure.Exon;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

@Program(name="genome2svg",
description="static genome browser as SVG",
keywords={"vcf","svg","xml","visualization"},
modificationDate="20240120",
creationDate="20240120",
generate_doc = false
)
public class GenomeToSvg extends Launcher {
	private static final Logger LOG=Logger.build(GenomeToSvg.class).make();
	private static int ID_GENERATOR = 0;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputPath=null;
	@Parameter(names={"-r","--interval","--region"},description=IntervalParserFactory.OPT_DESC,required = true)
	private String intervalStr = null;
	@Parameter(names={"-w","--width"},description="Image width")
	private double image_width = 1000;
	
	private class SVGContext {
		double y=0;
		int left_margin = 100;
		Document svgDom;
		Element defsNode;
		Element styleNode;
		Element tracksNode;
		ReferenceBean referenceBean;
		Locatable loc;
		double image_width;
		
		Element element(final String name) {
			return svgDom.createElementNS(SVG.NS,name);
		}
		Element element(final String name,final String content) {
			Element e = element(name);
			e.appendChild(text(content));
			return e;
			}
		
		double pos2pixel(double pos1) {
			return ((pos1 - loc.getStart())/(double)loc.getLengthOnReference()) * image_width;
		}
		
		double pixel2genomic(double pixX) {
			 return loc.getStart()+(pixX/image_width)*loc.getLengthOnReference();
		}
		
		double trimpos(double pos1) {
			return Math.max(loc.getStart(), Math.min(loc.getEnd(), pos1));
		}
		Text text(final String s) {
			return svgDom.createTextNode(s);
		}
		Element title(final String s) {
			return element("title",s);
		}
	}
	
	public static abstract class AbstractPathBean implements Closeable {
		private String path;
		protected AbstractPathBean() {}
		public void setPath(String path) {
			this.path = path;
			}
		public String getPath() {
			return path;
			}
		String requiredPath() {
			return path;
			}
		@Override
		public void close() {
			
		}
	}
	

	
	public static class ReferenceBean extends AbstractPathBean {
		private ReferenceSequenceFile refseqfile;
		private SAMSequenceDictionary dict;
		public ReferenceBean() {
			}
		private void open() {
			if(refseqfile!=null) return ;
			refseqfile = ReferenceSequenceFileFactory.getReferenceSequenceFile(Paths.get(requiredPath()));
			dict = SequenceDictionaryUtils.extractRequired(this.refseqfile);
			return;
			}
		public SAMSequenceDictionary getSequenceDictionary() {
			open(); return dict;
			}
		@Override
		public void close() {
			if(refseqfile!=null) try{refseqfile.close();}catch(Throwable err) {}
			refseqfile=null;
		}
	}
	
	public static class TabixBean extends AbstractPathBean {
		private TabixFileReader tabix = null;
		public TabixBean() {
			
			}
		private void open() {
			if(tabix!=null) return ;
			try {
				tabix = new TabixFileReader(requiredPath());
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() {
			if(tabix!=null) try{tabix.close();}catch(Throwable err) {}
			tabix=null;
		}
	}
	
	public static abstract class Track {
		private final String id = "id"+(++ID_GENERATOR);
		private String shortDesc = null;
		private String longDesc = null;
		private String display = "";
		private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
		protected Track() {
			}
		
		public String getId() {
			return id;
			}
		
		public void setShortDesc(String shortDesc) {
			this.shortDesc = shortDesc;
		}
		
		public String getShortDesc() {
			return StringUtils.isBlank(this.shortDesc)?getClass().getName():this.shortDesc;
			}
		
		public void setLongDesc(String longDesc) {
			this.longDesc = longDesc;
			}
		public String getLongDesc() {
			return StringUtils.isBlank(this.longDesc)?getShortDesc():this.longDesc;
			}
		
		public void setDisplay(String display) {
			this.display = display;
			}
		public String getDisplay() {
			return display;
		}
		public boolean isVisible() {
			return !getDisplay().equals("hidden");
			}
		abstract void paint(SVGContext ctx);
		
		protected void insertTitle(Element root,SVGContext ctx) {
			String t = getShortDesc();
			if(StringUtils.isBlank(t)) return;
			int fontSize = 10;
			ctx.y += 2;
			Element txt = ctx.element("text");
			txt.setAttribute("x", format(ctx.image_width/2.0));
			txt.setAttribute("y", format(ctx.y+fontSize));
			txt.setAttribute("style", "stroke:none;fill:darkgray;stroke-width:1px;text-anchor:middle;font-size:"+fontSize+";");
			txt.appendChild(ctx.text(t));
			t = getLongDesc();
			if(!StringUtils.isBlank(t)) {
				txt.appendChild(ctx.title(t));
				}
			root.appendChild(txt);
			ctx.y += fontSize;
			ctx.y += 5;
			}
		
		/** convert double to string */
		protected String format(final double v)
			{
			return this.decimalFormater.format(v);
			}
	}
	
	public static class BasesTrack extends Track {
		@Override
		void paint(SVGContext ctx) {
		
			double featureWidth=  ctx.image_width/(double)(ctx.loc.getLengthOnReference()); 
			double featureHeight= Math.min(Math.max(5.0,featureWidth),30); 			//base
			if(!isVisible() || featureWidth<5) return;
			final Hershey hershey = new Hershey();
			
			for(char base : new char[] {'A','C','G','T','N'}) {
				ctx.styleNode.appendChild(ctx.text(".b"+base+" {fill:none;stroke:"+AcidNucleics.cssColor(base)+";stroke-width:"+(featureWidth<5?0.5:1.0)+"px}\n"));
				}
			
			for(String base:new String[]{"a","A","t","T","g","G","c","C","n","N"})
				{
				final Element path  = ctx.element("path");
				ctx.defsNode.appendChild(path);
				
				path.setAttribute("id","b"+base);
				
				path.setAttribute("class","b"+base.toUpperCase());
				path.setAttribute("d",hershey.svgPath(
						base,
						0,
						0,
						featureWidth*0.95,
						featureHeight*0.95
						));
				path.appendChild(ctx.title(base));
				}
			
			final GenomicSequence genomicSequence = new GenomicSequence(ctx.referenceBean.refseqfile, ctx.loc.getContig());
			Element g = ctx.element("g");
			for(int pos= ctx.loc.getStart();
					featureWidth>5 && //ignore if too small
					pos<= ctx.loc.getEnd() && pos<=genomicSequence.length();++pos)
				{
				final char c= genomicSequence.charAt(pos-1);
				Element use = ctx.element("use");
				g.appendChild(use);
				use.setAttribute("x",format( ctx.pos2pixel(pos)));
				use.setAttribute("y",format(ctx.y));
				use.setAttribute("href", "#b"+c);
				use.appendChild(ctx.title(String.valueOf(c)+" "+StringUtils.niceInt(pos)));
				}
			ctx.tracksNode.appendChild(g);
			ctx.y+=featureHeight+1;
		}
	}
	
	public static class KnownGeneTrack extends Track {
		private TabixBean tabix;
		public KnownGeneTrack() {
			
			}
		public void setTabix(TabixBean tabix) {
			this.tabix = tabix;
			}
		public TabixBean getTabix() {
			return tabix;
			}
		@Override void paint(SVGContext ctx) {
			
			
			ctx.styleNode.appendChild(ctx.text("."+getId()+"strand {fill:none;stroke:gray;stroke-width:0.5px;}\n"));
			ctx.styleNode.appendChild(ctx.text("."+getId()+"kgtr {fill:none;stroke:gray;stroke-width:1px;}\n"));
			ctx.styleNode.appendChild(ctx.text("."+getId()+"kgexon {fill:lightgray;stroke:gray;stroke-width:1px;}\n"));

			
			final Element g0 = ctx.element("g");
			ctx.tracksNode.appendChild(g0);
			insertTitle(g0,ctx);
			tabix.open();
			Iterator<String> iter=tabix.tabix.iterator(
					ctx.loc.getContig(),
					ctx.loc.getStart(),
					ctx.loc.getEnd()
					);
			final UcscTranscriptCodec codec = new UcscTranscriptCodec();
			Pileup<UcscTranscript> pileup = new Pileup<>();
			
			while(iter.hasNext()) {
				UcscTranscript tr=codec.decode(iter.next());
				LOG.debug(tr);
				pileup.add(tr);
				}
			if(!pileup.isEmpty()) {
				Element polyline = ctx.element("polyline");
				ctx.defsNode.appendChild(polyline);
				polyline.setAttribute("id", "strandF");
				polyline.setAttribute("class", getId()+"strand");
				polyline.setAttribute("points", "-5,-5 0,0 -5,5");
				
				polyline = ctx.element("polyline");
				ctx.defsNode.appendChild(polyline);
				polyline.setAttribute("id", "strandR");
				polyline.setAttribute("class", getId()+"strand");
				polyline.setAttribute("points", "5,-5 0,0 5,5");
				}
			final double exonHeight = 10;
			
			for(List<UcscTranscript> row:pileup.getRows()) {
				double midY= ctx.y + exonHeight/2.0;
				for(UcscTranscript tr:row) {
					final Element g = ctx.element("g");
					g0.appendChild(g);
					
					/* transcript line */
					Element line = ctx.element("line");
					g.appendChild(line);
					line.setAttribute("class",getId()+"kgtr");
					line.setAttribute("x1",String.valueOf(ctx.pos2pixel(ctx.trimpos(tr.getStart()))));
					line.setAttribute("y1",String.valueOf(midY));
					line.setAttribute("x2",String.valueOf(ctx.pos2pixel(ctx.trimpos(tr.getEnd()))));
					line.setAttribute("y2",String.valueOf(midY));
					

					
					/* strand symbols */
					for(double pixX=0;
						pixX< ctx.image_width;
						pixX+=30)
						{
						double pos1= ctx.pixel2genomic(pixX);
						if(pos1< tr.getStart()) continue;
						if(pos1> tr.getEnd()) break;
						Element use = ctx.element("use");
						g.appendChild(use);
						use.setAttribute("href", "#strand"+(tr.isPositiveStrand()?"F":"R"));
						use.setAttribute("x",String.valueOf(pixX));
						use.setAttribute("y",String.valueOf(midY));
						}
					
					/* exons */
					for(final UcscTranscript.Exon exon: tr.getExons())
						{
						if(!exon.overlaps(ctx.loc)) continue;
						Element rect = ctx.element("rect");
						g.appendChild(rect);
						rect.setAttribute("class",getId()+"kgexon");
						
						rect.setAttribute("x",String.valueOf(ctx.pos2pixel(ctx.trimpos(exon.getStart()))));
						rect.setAttribute("y",String.valueOf(midY-exonHeight/2));
						rect.setAttribute("width",String.valueOf(
								ctx.pos2pixel(ctx.trimpos(exon.getEnd()))-ctx.pos2pixel(ctx.trimpos(exon.getStart()))));
						rect.setAttribute("height",String.valueOf(exonHeight));
						rect.appendChild(ctx.title(exon.getName()));
						}
					
					}
				
				// legend
				{
				final Set<String> geneNames = row.stream().map(KG->StringUtils.ifBlank(KG.getName2(), KG.getTranscriptId())).collect(Collectors.toCollection(LinkedHashSet::new));
				if(geneNames.size()<10) {
					final Element legend = ctx.element("text",String.join(" ", geneNames));
					legend.setAttribute("style", "text-anchor:end;font-size:10px;");
					legend.setAttribute("x", format(-5));
					legend.setAttribute("y",String.valueOf(ctx.y+10));
					g0.appendChild(legend);
					}
				}

				
				ctx.y += exonHeight;
				ctx.y += 2;
				}
			}
		}
	
	
	public static class VariantToStringBean implements BiFunction<VCFHeader, VariantContext, String>{
		@Override
		public String apply(VCFHeader t, VariantContext vc) {
			return vc.getContig()+":"+vc.getStart()+" "+vc.getType().name();
			}
		}
	
	public static class VcfTrack extends Track {
		private String path = null;
		private BiPredicate<VCFHeader,VariantContext> filter = (H,V)->true;
		private  BiFunction<VCFHeader, VariantContext, String> variantToString = new VariantToStringBean();
		private  BiFunction<VCFHeader, VariantContext, String> variantToStyle =  (H,V)->"fill:red;stroke:blue;";
		private  BiFunction<VCFHeader, VariantContext, String> variantToURL =  (H,V)->"http://en.wikipediag.org";
		public void setPath(String path) {
			this.path = path;
			}
		public String getPath() {
			return path;
			}
		
		@Override
		void paint(SVGContext ctx) {
			
			final Element g = ctx.element("g");
			ctx.tracksNode.appendChild(g);
			insertTitle(g,ctx);
			double featureHeight = 10;
			
			try(VCFReader reader = VCFReaderFactory.makeDefault().open(getPath(), true)) {
				final VCFHeader header= reader.getHeader();
				try(CloseableIterator<VariantContext> iter = reader.query(ctx.loc)) {
					while(iter.hasNext()) {
						final VariantContext vc = iter.next();
						if(!this.filter.test(header, vc)) continue;
						Element rect = ctx.element("rect");
						
						String url = variantToURL.apply(header, vc);
						if(StringUtils.isBlank(url)) {
							g.appendChild(rect);
							}
						else
							{
							Element a = ctx.element("a");
							a.setAttribute("href", url);
							g.appendChild(a);
							a.appendChild(rect);
							}
							
						
						rect.setAttribute("x",String.valueOf(ctx.pos2pixel(ctx.trimpos(vc.getStart()))));
						rect.setAttribute("y",String.valueOf(ctx.y));
						rect.setAttribute("width",String.valueOf(
								ctx.pos2pixel(ctx.trimpos(vc.getEnd()+1))-ctx.pos2pixel(ctx.trimpos(vc.getStart()))));
						rect.setAttribute("height",String.valueOf(featureHeight));
						String title = this.variantToString.apply(header,vc);
						if(!StringUtils.isBlank(title)) rect.appendChild(ctx.title(title));
						
						String style = variantToStyle.apply(header, vc);
						if(!StringUtils.isBlank(style)) rect.setAttribute("style", style);

						}
					}
				ctx.y += featureHeight+1;
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		}
	
	
	public static class BamCoverageTrack extends Track {
		private String path = null;
		private int featureHeight = 100;
		private Predicate<SAMRecord> samFilter = R->true;
		private int minCoverage = -1;
		public BamCoverageTrack() {
			}
		public void setPath(String path) {
			this.path = path;
			}
		public String getPath() {
			return path;
			}
		public void setMinCoverage(int minCoverage) {
			this.minCoverage = minCoverage;
			}
		public int getMinCoverage() {
			return minCoverage;
			}
		
		@Override
		public String getShortDesc() {
			return StringUtils.isBlank(super.shortDesc)?this.getPath():super.shortDesc;
			}
		
		public void setFeatureHeight(int featureHeight) {
			this.featureHeight = featureHeight;
			}
		public int getFeatureHeight() {
			return Math.max(10,featureHeight);
			}
		
		@Override
		void paint(SVGContext ctx) {
			final int width = (int)ctx.image_width;
			final double featureHeight = getFeatureHeight();
			String sampleName = getPath();
			final SamReaderFactory srf = SamReaderFactory.make();
			try(SamReader sr=srf.open(SamInputResource.of(getPath()))) {
				final SAMFileHeader header = sr.getFileHeader();
				sampleName = header.getReadGroups().stream().map(RG->RG.getSample()).filter(S->!StringUtils.isBlank(S)).findFirst().orElse(getPath());
				final double[] array = new double[width];
				Arrays.fill(array, -999);
				try(CloseableIterator<SAMRecord> iter = sr.queryOverlapping(ctx.loc.getContig(), ctx.loc.getStart(), ctx.loc.getEnd())) {
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						if(!samFilter.test(rec) || rec.getReadUnmappedFlag()) continue;
						int prev=-1;
						for(AlignmentBlock block : rec.getAlignmentBlocks()) {
							for(int x=0;x< block.getLength();++x) {
								int p1 = block.getReferenceStart()+x;
								if(p1 < ctx.loc.getStart()) continue;
								if(p1 > ctx.loc.getEnd()) break;
								int p2 = (int)ctx.pos2pixel(p1);
								if(prev==p2) continue;
								prev=p2;
								if(p2<0 || p2>array.length) continue;
								if(array[p2]<0) array[p2]=0;
								array[p2]++;
								}
							}
						}
					}
				int x=0;
				while(x<array.length) {
					if(array[x]<0) {
						array[x] = x>0?array[x-1]:0;
						}
					x++;
					}
				double maxcov = Math.max(1, Arrays.stream(array).max().orElse(1));
				if(getMinCoverage()>0 && maxcov < getMinCoverage()) maxcov = getMinCoverage();
				//defs
				{
				Element grad= ctx.element("linearGradient");
				ctx.defsNode.appendChild(grad);
				grad.setAttribute("id",getId()+"grad");
				grad.setAttribute("x1","0%");
				grad.setAttribute("y1","0%");
				grad.setAttribute("x2","0%");
				grad.setAttribute("y2","100%");
				//grad.setAttribute("gradientTransform","rotate(90)");
				
				Element stop = ctx.element("stop");
				grad.appendChild(stop);
				stop.setAttribute("offset", "0%");
				stop.setAttribute("stop-color",(maxcov>50?"red":maxcov>20?"green":"blue"));
				
				stop = ctx.element("stop");
				grad.appendChild(stop);
				stop.setAttribute("offset", "100%");
				stop.setAttribute("stop-color", "darkblue");
				}
				
				
				final Element g0 = ctx.element("g");
				ctx.tracksNode.appendChild(g0);
				insertTitle(g0,ctx);
				
				final Element legend = ctx.element("text",sampleName);
				legend.setAttribute("x", format(-ctx.left_margin));
				legend.setAttribute("y",String.valueOf(ctx.y+10));
				g0.appendChild(legend);
				
				final Element g = ctx.element("g");
				g0.appendChild(g);
				g.setAttribute("transform", "translate(0,"+ctx.y+")");
				ctx.tracksNode.appendChild(g);
				
				final Element path = ctx.element("polyline");
				path.setAttribute("style", "stroke:red;fill:url('#"+getId()+"grad')");
				g.appendChild(path);
				final List<Point2D> points = new ArrayList<>(array.length+2);
				points.add(new Point2D.Double(0,featureHeight));
				for(int i=0;i< array.length;i++) {
					points.add(new Point2D.Double(i, featureHeight - (array[i]/maxcov)*featureHeight));
					}
				points.add(new Point2D.Double(array.length-1,featureHeight));
				points.add(new Point2D.Double(0,featureHeight));
				x=1;
				while(x<points.size()) {
					if(x>0 && x+1 < points.size() &&  points.get(x-1).getY()==points.get(x).getY() && points.get(x).getY()==points.get(x+1).getY()) {
						points.remove(x);
						}
					else
						{
						++x;
						}
					}
								
				path.setAttribute("points", points.stream().
						map(P->format(P.getX())+","+format(P.getY())).
						collect(Collectors.joining(" "))
						);
				// yaxis
				ctx.styleNode.appendChild(ctx.text("."+getId()+"ticks{stroke:darkgray;stroke-width:1px;}\n"));
				ctx.styleNode.appendChild(ctx.text("."+getId()+"lbl{text-anchor:end;font-size:10px;}\n"));
				final int nTicks = 10;
				for(int i=0;i< nTicks;i++) {
					final double cov = (maxcov/nTicks)*i;
					final double v = featureHeight - (cov/maxcov) * featureHeight;
					final Element line = ctx.element("line");
					line.setAttribute("x1", format(0));
					line.setAttribute("y1", format(v));
					line.setAttribute("x2", format(-5));
					line.setAttribute("y2", format(v));
					line.setAttribute("class",getId()+"ticks");

					g.appendChild(line);
					
					final Element text = ctx.element("text",format(cov));
					text.setAttribute("class",getId()+"lbl");
					text.setAttribute("x", format(-5));
					text.setAttribute("y", format(v));
					g.appendChild(text);
					
				}
				
				ctx.y+=featureHeight;
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		
	}
	
	private ApplicationContext configuration = null;
	@Override
	public int doWork(List<String> args) {
		try {
			if(args.isEmpty()) {
				LOG.error("Spring config empty");
				return -1;
				}
			final SVGContext svgContext = new SVGContext();
			svgContext.image_width = this.image_width;
			svgContext.image_width = this.image_width;
			this.configuration =  new FileSystemXmlApplicationContext(args.toArray(new String[args.size()]));
			svgContext.referenceBean = this.configuration.getBean("reference", ReferenceBean.class);
			svgContext.referenceBean.open();
			
			
			
			svgContext.loc = IntervalParserFactory.newInstance(svgContext.referenceBean.getSequenceDictionary()).
					make().apply(this.intervalStr).
					get();

		
			
			final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			final DocumentBuilder db = dbf.newDocumentBuilder();
			final Document svgDom = db.newDocument();
			svgContext.svgDom = svgDom;
			final Element svgRoot = svgContext.element("svg");
			svgRoot.setAttribute("width", String.valueOf(1+this.image_width + svgContext.left_margin));
			svgDom.appendChild(svgRoot);
			
			svgContext.styleNode = svgContext.element("style");
			svgRoot.appendChild(svgContext.styleNode);

			
			svgContext.defsNode = svgContext.element("defs");
			svgRoot.appendChild(svgContext.defsNode);
			svgContext.tracksNode = svgContext.element("g");
			svgContext.tracksNode.setAttribute("transform", "translate("+svgContext.left_margin+",0)");
			svgRoot.appendChild(svgContext.tracksNode);
			final Element frame =  svgContext.element("rect");
			svgContext.tracksNode.appendChild(frame);
			
			frame.setAttribute("style", "fill:white;stroke:black;");
			frame.setAttribute("x", "0");
			frame.setAttribute("y", "0");
			frame.setAttribute("width", String.valueOf(svgContext.image_width+ svgContext.left_margin-1));
			

			for(Object ot:this.configuration.getBean("tracks",List.class)) {
				final Track track = Track.class.cast(ot);
				if(!track.isVisible()) continue;
				double prev_y = svgContext.y;
				track.paint(svgContext);
				if(svgContext.y>prev_y) {
					final Element frame2 =  svgContext.element("rect");
					frame2.setAttribute("style","stroke:darkgray;fill:none;");
					frame2.setAttribute("x", String.valueOf(-svgContext.left_margin));
					frame2.setAttribute("y",String.valueOf(prev_y));
					frame2.setAttribute("width", String.valueOf(svgContext.image_width+ svgContext.left_margin-1));
					frame2.setAttribute("height", String.valueOf(svgContext.y-prev_y));
					svgContext.tracksNode.appendChild(frame2);
					}
				}
			
			
			svgRoot.setAttribute("height", String.valueOf(svgContext.y));
			frame.setAttribute("height", String.valueOf(svgContext.y));
			
			TransformerFactory.newInstance().newTransformer().
				transform(new DOMSource(svgDom),
				this.outputPath==null?
				new StreamResult(stdout()):
				new StreamResult(outputPath.toFile())
				);
			svgContext.referenceBean.close();
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			configuration = null;
			}
		}
	
	public static void main(String[] args) {
		new GenomeToSvg().instanceMain(args);
	}
}
