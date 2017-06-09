package com.github.lindenb.jvarkit.tools.genscan;

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/**

BEGIN_DOC


END_DOC

*/
@Program(
		name="genscan2",
		description="Paint a Genome Scan picture from a Tab delimited file (CHROM/POS/VALUE1/VALUE2/....).",
		keywords={"chromosome","reference","chart","visualization"},
		generate_doc=false
		)
public class GenScan2 extends Launcher {
	
	private static final Logger LOG = Logger.build(GenScan2.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidxFile = null;
	@Parameter(names={"-r","--region"},description="One or more Specific region to observe. empty string is whole genome. Or the name of a chromosome. Or an interval.")
	private Set<String> regions = new HashSet<>();
	@Parameter(names={"-miny","--miny","--ymin","-ymin"},description="min Y value")
	private double minY=0;
	@Parameter(names={"-maxy","--maxy","--ymax","-ymax"},description="max Y value")
	private double maxY=100.0;
	@Parameter(names={"-W","--width"},description="width")
	private int WIDTH=1000;
	@Parameter(names={"-H","--height"},description="height")
	private int HEIGHT=700;
	@Parameter(names={"-dbc","--distancebetweencontigs"},description="number of pixels between contigs")
	private double distanceBetweenContigs = 1;
	
	private SAMSequenceDictionary dict=null;
	private List<DisplayRange> displayRanges = new ArrayList<>();
	private long genomeViewLength=0L;
	private final Insets insets = new Insets(100, 100, 100, 50);
	private final Hershey hershey = new Hershey();
	private final ColorUtils colorUtils = new ColorUtils();
	private enum ShapeType { circle,square}; 
	
	private class Style
		{
		Color penColor=Color.BLACK;
		Color paperColor=Color.ORANGE;
		float alpha=1.0f;
		double size=2;
		float strokeWidth=0.5f;
		ShapeType shapeType=ShapeType.circle;
		Style() {
			
			}
		Style(final Style cp) {
			this.penColor = cp.penColor;
			this.paperColor = cp.paperColor;
			this.alpha = cp.alpha;
			this.size = cp.size;
			this.strokeWidth = cp.strokeWidth;
			this.shapeType=cp.shapeType;
			}
		void update( String line) {
			if(line.startsWith("#!")) line= line.substring(2);
			for(final String token:line.split("[;]"))
				{
				if(token.trim().isEmpty()) continue;
				int colon=token.indexOf(':');
				if(colon==-1 |- colon==0) continue;
				final String key = token.substring(0,colon).toLowerCase().trim();
				final String value = token.substring(colon+1).trim();
				if(key.equals("color") || key.equals("stoke") || key.equals("stoke-color"))
					{
					this.penColor=(value.equals("null") || value.isEmpty()?null:colorUtils.parse(value));
					}
				else if(key.equals("fill") || key.equals("fill-color") || key.equals("background-color"))
					{
					this.paperColor=(value.equals("null") || value.isEmpty()?null:colorUtils.parse(value));
					}
				else if(key.equals("size") || key.equals("width") )
					{
					this.size=Double.parseDouble(value);
					}
				else if(key.equals("alpha") || key.equals("opacity") )
					{
					this.alpha=Float.parseFloat(value);
					}
				else if(key.equals("stroke-width") )
					{
					this.strokeWidth=Float.parseFloat(value);
					}
				else if(key.equals("transparency")  )
					{
					this.alpha= 1f - Float.parseFloat(value);
					}
				else if(key.equals("shape")  )
					{
					this.shapeType= ShapeType.valueOf(value);
					}
				else
					{
					LOG.error("unknown selector "+token);
					}
				}
			}
		}
	private final Stack<Style> styleStack = new Stack<>();

	
	private class Data implements Locatable
		{
		private double value=0;
		private String contig;
		private int pos;
		private Style style= GenScan2.this.styleStack.peek();
		public double getValue() { return this.value;}
		public Style getStyle(){ return this.style;}
		@Override
		public String getContig() {
			return contig;
			}
		@Override
		public int getStart() {
			return pos;
			}
		@Override
		public int getEnd() {
			return pos;
			}
		}

	
	private abstract class DisplayRange implements Locatable
		{
		double startX;
		double width=0;
		
		
		public abstract String getLabel();
		public abstract int length();
		public boolean contains(final Locatable data) {
			if(!data.getContig().equals(this.getContig())) return false;
			if(data.getEnd() < this.getStart()) return false;
			if(data.getStart() > this.getEnd()) return false;
			
			return true;
			}
		public int getY()
			{
			return GenScan2.this.insets.top;
			}
		public int getHeight() {
			return GenScan2.this.HEIGHT-(GenScan2.this.insets.top + GenScan2.this.insets.bottom);
			}
		public Rectangle2D.Double getBounds()
			{
			return new Rectangle2D.Double(
				this.startX,
				getY(),
				this.width,
				getHeight());
			}
		
		final Point2D.Double convertToPoint(final Data data)
			{
			final double x = this.startX + 
					((data.getStart()-this.getStart())/(double)(this.length()))*this.width;
			
			final double y = GenScan2.this.HEIGHT - 
					( GenScan2.this.insets.bottom + ((data.getValue() - GenScan2.this.minY)/( GenScan2.this.maxY- GenScan2.this.minY))*(double)(getHeight()));
			
			return new Point2D.Double(x,y);
			}
		}
	private class DisplayInterval extends DisplayRange
		{
		private final Interval interval;
		DisplayInterval(final String s)
			{
			final IntervalParser intervalParser = new IntervalParser(GenScan2.this.dict);
			this.interval = intervalParser.parse(s);
			if(this.interval==null) {
				throw new IllegalArgumentException("cannot parse interval "+s);
				}
			}
		@Override
		public String getLabel() {
			return this.getContig()+":"+this.getStart()+"-"+this.getEnd();
			}
		
		@Override
		public String getContig() {
			return this.interval.getContig();
			}
		@Override
		public int getStart() {
			return this.interval.getStart();
			}
		
		@Override
		public int getEnd() {
			return this.interval.getEnd();
			}
		
		@Override
		public int length() {
			return this.interval.length();
			}
		}
	
	private class DisplayContig extends DisplayRange
		{
		private final SAMSequenceRecord ssr;
		DisplayContig(final String contig) {
			this.ssr = GenScan2.this.dict.getSequence(contig);
			if(this.ssr==null) {
				throw new IllegalArgumentException("cannot find contig in dictionary "+contig);
				}
			}
		DisplayContig(final SAMSequenceRecord ssr) {
			this.ssr = ssr;
			if(this.ssr==null) {
				throw new IllegalArgumentException("SAMSequenceRecord==null");
				}
			}
		
		@Override
		public String getLabel() {
			return getContig();
			}
		
		@Override
		public String getContig() {
			return this.ssr.getSequenceName();
			}
		
		@Override
		public int getStart() {
			return 1;
			}
		
		@Override
		public int getEnd() {
			return this.ssr.getSequenceLength();
			}
		
		@Override
		public int length() {
			return this.ssr.getSequenceLength();
			}
		}
	
	
	private Data parseData(final String line) {
		String tokens[]=line.split("[\t]");
		Data data=new Data();
		data.contig=tokens[0];
		data.pos=Integer.parseInt(tokens[1]);
		data.value=Double.parseDouble(tokens[2]);
		
		return data;
	}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.faidxFile==null) {
			LOG.error("Reference missing");
			return -1;
			}
		if(this.maxY<=this.minY)
			{
			LOG.error("MaxY <= MinY");
			return -1;
			}
		
		BufferedReader r=null;
		try {
			final Rectangle2D.Double drawingRect = new Rectangle2D.Double(
					this.insets.left,
					this.insets.top,
					WIDTH-(this.insets.left+this.insets.right),
					HEIGHT-(this.insets.top+this.insets.bottom)
					);
			
			final Style styleBase = new Style();
			this.styleStack.add(styleBase);
			
			
			this.dict = SAMSequenceDictionaryExtractor.extractDictionary(this.faidxFile);
			for(final String rgnstr: this.regions) {
				if(rgnstr.trim().isEmpty()) continue;
				final DisplayRange contig;
				if(rgnstr.contains(":")) {
					contig = new DisplayInterval(rgnstr);
				} else
					{
					contig = new DisplayContig(rgnstr);
					}
				if(contig.length()==0) continue;
				this.displayRanges.add(contig);
				}
			// no region, add whole genome
			for(final SAMSequenceRecord ssr:this.dict.getSequences()) {
				this.displayRanges.add(new DisplayContig(ssr));
			}
			// sort on chrom/start
			Collections.sort(this.displayRanges, (DR1,DR2)->{
				final int tid1 = dict.getSequenceIndex(DR1.getContig());
				final int tid2 = dict.getSequenceIndex(DR2.getContig());
				int i=tid1-tid2;
				if(i!=0) return i;
				i = DR1.getStart() - DR2.getStart();
				if(i!=0) return i;
				i = DR1.getEnd() - DR2.getEnd();
				return i;
			});
			this.genomeViewLength = this.displayRanges.stream().mapToLong(DR->DR.length()).sum();
			final double adjustedImgWidth = drawingRect.getWidth() - this.distanceBetweenContigs*(this.displayRanges.size()-1);
			if(adjustedImgWidth<=0)
				{
				LOG.error("with such settings image width would be empty");
				return -1;
				}
			double x = this.insets.left;
			for(final DisplayRange displayRange: this.displayRanges) {
				displayRange.startX = x;
				displayRange.width = (displayRange.length()/(double)genomeViewLength)*adjustedImgWidth;
				x+=displayRange.width;
				x+=this.distanceBetweenContigs;
				}
			
			final BufferedImage img = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
			
			
			final Graphics2D g = img.createGraphics();
			g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g.setColor(ColorUtils.cornsilk);
			g.fillRect(0, 0, WIDTH, HEIGHT);
			
			//plot background of each contig
			for(int rangeIdx=0;rangeIdx< this.displayRanges.size();++rangeIdx) {
				final DisplayRange displayRange = this.displayRanges.get(rangeIdx);
				final Rectangle2D rec= displayRange.getBounds();
				g.setColor(rangeIdx%2==0?ColorUtils.antiquewhite:ColorUtils.navajowhite);
				g.fill(rec);
				}
			
			
			r = super.openBufferedReader(oneFileOrNull(args));
			String line;
			while((line=r.readLine())!=null) {
				if(line.trim().isEmpty()) continue;
				if(line.startsWith("#"))
					{
					if(line.equals("#!push"))
						{
						this.styleStack.push(new Style(this.styleStack.peek()));
						}
					else if(line.equals("#!pop"))
						{
						if(this.styleStack.size()==1) {
							LOG.error("Cannot pop bottom style");
							return -1;
							}
						this.styleStack.pop();
						}
					else if(line.startsWith("#!log"))
						{
						LOG.info(line);
						}
					else if(line.startsWith("#!"))
						{
						this.styleStack.peek().update(line);
						}
					continue;
					}
				final Data data = parseData(line);
				if(data==null) continue;
				
				
				if(data.getValue()< GenScan2.this.minY) {continue;}
				if(data.getValue()> GenScan2.this.maxY) {continue;}
				
				final Style style= data.getStyle();
				
				for(final DisplayRange dr: GenScan2.this.displayRanges) {
					if(!dr.contains(data)) continue;
					final Rectangle2D.Double bounds= dr.getBounds();

					final Point2D.Double point= dr.convertToPoint(data);
					
					final Shape oldClip = g.getClip();
					final Composite oldComposite = g.getComposite();
					g.setClip(bounds);
					//
					g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, style.alpha));
					final Shape dataShape;
					switch(style.shapeType)
						{
						case circle:
							dataShape = new Ellipse2D.Double(
									point.x-style.size/2.0,
									point.y-style.size/2.0,
									style.size,
									style.size
									);
							break;
						case square:
							dataShape = new Rectangle2D.Double(
									point.x-style.size/2.0,
									point.y-style.size/2.0,
									style.size,
									style.size
									);
							break;
						default:throw new IllegalArgumentException();
						}
					
					if(style.paperColor!=null)
						{
						g.setColor(style.paperColor);
						g.fill(dataShape);
						}
					if(style.penColor!=null && style.strokeWidth>0)
						{
						final Stroke oldStroke = g.getStroke();
						final Stroke stroke = new BasicStroke(
								style.strokeWidth,
								BasicStroke.CAP_BUTT,
								BasicStroke.JOIN_ROUND
								);
						g.setStroke(stroke);
						g.setColor(style.penColor);
						g.draw(dataShape);
						g.setStroke(oldStroke);
						}
					
					
					//
					g.setClip(oldClip);
					
					g.setComposite(oldComposite);
					}
					
				}
			
			
			//plot frame of each contig
			for(final DisplayRange displayRange : this.displayRanges) {
				final Rectangle2D rec= displayRange.getBounds();
				g.setColor(ColorUtils.aquamarine);
				g.draw(rec);
				final String label= displayRange.getLabel();
				int fontSize=insets.top-2;
				final double labelLength=Math.min( label.length()*fontSize,rec.getWidth());
				this.hershey.paint(g, label, new Rectangle2D.Double(
						rec.getCenterX()-labelLength/2.0,
						1,
						labelLength,
						(insets.top-2)
						));
				//plot ticks
				for(int i=0;i<= 10;i++)
					{
					double tx=drawingRect.x+ i*(drawingRect.getWidth()/10.0);
					double ty = drawingRect.getMaxY();
					g.setColor(ColorUtils.blueviolet);
					g.draw(new Line2D.Double(tx,ty, tx, ty+5));
					if(i>0) {
						g.translate(tx, ty+5);
						g.rotate(Math.PI/2.0);
						g.drawString(String.valueOf(i), 0, 0);
						g.rotate(-Math.PI/2.0);
						g.translate(-tx, -(ty+5));
						}
					}
				}

			//plot y axis
			
			for(int i=0;i<= 10;++i)
				{
				double v= this.minY+ ((this.maxY-this.minY)/10.0)*i;
				double y= drawingRect.getMaxY() - (drawingRect.getHeight()/10.0)*i;
				double font_size=12;
				
				g.setColor(Color.LIGHT_GRAY);
				g.draw(new Line2D.Double(
						insets.left-5, y,
						WIDTH-insets.right,
						y));
				
				g.setColor(Color.BLACK);
				this.hershey.paint(g, String.format("%8s",v), new Rectangle2D.Double(
						0,y-font_size,
						insets.left,
						font_size
						));
				}
			
			g.dispose();
			
			r.close();r=null;
			
			if(this.outputFile==null)
				{
				ImageIO.write(img, "PNG", stdout());
				}
			else
				{
				ImageIO.write(img,
						this.outputFile.getName().toLowerCase().endsWith(".png")?"PNG":"JPG",
						this.outputFile
						);
				}
			return 0;
			} 
		catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		finally {
			CloserUtil.close(r);
			this.dict =null;
			}
		}
	
	
	public static void main(String[] args) {
		new GenScan2().instanceMainWithExit(args);
	}

}
