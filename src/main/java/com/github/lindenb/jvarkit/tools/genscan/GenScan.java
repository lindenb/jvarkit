package com.github.lindenb.jvarkit.tools.genscan;

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics2D;
import java.awt.Insets;
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
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/**

BEGIN_DOC

## INPUT



Input consists in 3 fields delimited with a tabulation. An optional 4th column can be used to set a specific style for the point

```
(CHROM)<tab>(POS)<tab>(VALUE)(<tab>STYLE)?
```

### directives

* `#!push` push a copy of the current style on the stack
* `#!pop` pop the current style off the stack

## style




## Example

```
$ samtools depth in.bam |\
		java -jar dist/genscan.jar  --removeContigsSmallerThan 500000 --width 2000 -R  human_g1k_v37.fasta -o out.png
```	

```
(echo "#!track:S1;"; samtools depth S1.bam ;echo "#!track:S2;"; samtools depth S2.bam ) | java -jar dist/genscan.jar --track "S1,S2" --width 2000 --removeContigsSmallerThan 500000 -R  human_g1k_v37.fasta -o out.png
```


END_DOC

*/
@Program(
		name="genscan",
		description="Paint a Genome Scan picture from a Tab delimited file (CHROM/POS/VALUE1/VALUE2/....).",
		keywords={"chromosome","reference","chart","visualization"},
		generate_doc=true
		)
public class GenScan extends Launcher {
	
	private static final Logger LOG = Logger.build(GenScan.class).make();
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
	@Parameter(names={"-W","--width"},description="Image width")
	private int WIDTH=1000;
	@Parameter(names={"-H","--height"},description="Image height")
	private int HEIGHT=700;
	@Parameter(names={"-dbc","--distancebetweencontigs"},description="number of pixels between contigs")
	private double distanceBetweenContigs = 1;
	@Parameter(names={"-track","--track"},description="Add a track by specifying it's name")
	private Set<String> trackNames = new TreeSet<>();
	@Parameter(names={"-style","--style"},description="Default style")
	private String defaultStyleStr="";
	@Parameter(names={"--removeContigsSmallerThan"},description="When displaying a whole reference, don't use the configs having a length lower than this number. Useful to only display the main chromosomes.")
	private int removeContigsSmallerThan=0;

	/** SAM sequence dict associated to the data, will be used to get the whole contigs */
	private SAMSequenceDictionary dict=null;
	private List<DisplayRange> displayRanges = new ArrayList<>();
	private long genomeViewLength=0L;
	private final Insets insets = new Insets(100, 100, 100, 50);
	private final Hershey hershey = new Hershey();
	private final ColorUtils colorUtils = new ColorUtils();
	private enum ShapeType { circle,square}; 
	/** default track, only defined if trackNames isEmpty() */
	private Track defaultTrack = null;
	
	protected static final Color ALMOST_BLACK = new Color(20,20,20);
	protected static final Color ALMOST_WHITE = new Color(240,240,240);

	
	/** CSS - like selector */
	private enum Selector {
		penColor,
		paperColor,
		f_alpha,
		d_size,
		f_strokeWidth,
		shapeType,
		track_name
		}
	
	
	/** CSS style: just a map of selectors */ 
	private class Style
		{
		final Map<Selector,Object> selectors ;

		Style() {
			this.selectors = new HashMap<>();
			this.selectors.put(Selector.penColor, Color.BLACK);
			this.selectors.put(Selector.paperColor, Color.ORANGE);
			this.selectors.put(Selector.f_alpha, 1.0f);
			this.selectors.put(Selector.d_size, 1.0 );
			this.selectors.put(Selector.f_strokeWidth, 0.5f);
			this.selectors.put(Selector.shapeType, ShapeType.circle);
			}
		Style(final Style cp) {
			this.selectors = new HashMap<>(cp.selectors);
			}
		void update( final String line) {
			GenScan.this.parseSelectors(this.selectors,line);
			}
		public float getAlpha() {
			return Float.class.cast(this.selectors.get(Selector.f_alpha));
			}
		public float getStrokeWidth() {
			return Float.class.cast(this.selectors.get(Selector.f_alpha));
			}
		public Color getPenColor() {
			return Color.class.cast(this.selectors.get(Selector.penColor));
			}
		public Color getPaperColor() {
			return Color.class.cast(this.selectors.get(Selector.paperColor));
			}
		public ShapeType getShapeType() {
			return ShapeType.class.cast(this.selectors.get(Selector.shapeType));
			}
		public double getSize() {
			return Double.class.cast(this.selectors.get(Selector.d_size));
			}
		}
	/** stack of style when reading data */
	private final Stack<Style> styleStack = new Stack<>();

	/** an horizontal track */
	private class Track
		{
		double startY=0;
		double height=0.0;
		final String label;
		Track(final String label) {
			this.label=label;
			}
		final double convertToY(final Data data)
			{
			return convertToY(data.getValue());
			}
		final double convertToY(final double value)
			{
			return this.startY + 
					this.height - 
					((value-GenScan.this.minY)/(double)(GenScan.this.maxY-GenScan.this.minY))*this.height;
			}
		}
	
	private final Map<String,Track> name2track = new HashMap<>();
	
	
	private class Data implements Locatable
		{
		private double value=0;
		private String contig;
		private int pos;
		private Style style= GenScan.this.styleStack.peek();
		private Track track;
		
		
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

	
	/**
	 * a genomic region to display
	 */
	private abstract class DisplayRange implements Locatable
		{
		double startX;
		double width=0;
		
		/** label for this range */
		public abstract String getLabel();
		/** length of this interval in bp */
		public abstract int length();
		public boolean contains(final Locatable data) {
			if(!data.getContig().equals(this.getContig())) return false;
			if(data.getEnd() < this.getStart()) return false;
			if(data.getStart() > this.getEnd()) return false;
			return true;
			}
				
		final double convertToX(final Data data)
			{
			return convertToX(data.getStart());
			}
		final double convertToX(final double pos)
			{
			return this.startX + 
					((pos-this.getStart())/(double)(this.length()))*this.width;
			}
		@Override
		public String toString() {
			return getLabel();
			}
		}
	
	/** instance of display interval for a given Interval */
	private class DisplayInterval extends DisplayRange
		{
		private final Interval interval;
		DisplayInterval(final String intervalStr)
			{
			final IntervalParser intervalParser = new IntervalParser(GenScan.this.dict);
			this.interval = intervalParser.parse(intervalStr);
			if(this.interval==null) {
				throw new IllegalArgumentException("cannot parse interval \""+intervalStr+"\"");
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
	
	/** instance of display interval for a whole contig/interval */
	private class DisplayContig extends DisplayRange
		{
		private final SAMSequenceRecord ssr;
		DisplayContig(final String contig) {
			this.ssr = GenScan.this.dict.getSequence(contig);
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
	
	/** intersection track & display range */
	private class Pane
		{
		final Track track;
		final DisplayRange range;
		Pane(final DisplayRange range,final Track track) {
			this.range = range;
			this.track = track;
			}
		public double getX() { return this.range.startX;}
		public double getY() { return this.track.startY;}
		public double getWidth() { return this.range.width;}
		public double getHeight() { return this.track.height;}
		public Rectangle2D.Double getBounds() {
			return new Rectangle2D.Double(getX(),getY(),getWidth(),getHeight());
			}
		
		public Point2D.Double convertToPoint(final Data data) {
			return new Point2D.Double(
				this.range.convertToX(data),
				this.track.convertToY(data)
				);
			}
		
		}
	
	private Color parseColor(final String s) {
		if(s==null || s.trim().isEmpty()) return null;
		return this.colorUtils.parse(s);
	}
	
	private void parseSelectors(final Map<Selector,Object> selectors,String line) {
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
				selectors.put(Selector.penColor, parseColor(value));
				}
			else if(key.equals("fill") || key.equals("fill-color") || key.equals("background-color"))
				{
				selectors.put(Selector.paperColor, parseColor(value));
				}
			else if(key.equals("size") || key.equals("width") )
				{
				selectors.put(Selector.d_size, Double.parseDouble(value));
				}
			else if(key.equals("alpha") || key.equals("opacity") )
				{
				selectors.put(Selector.f_alpha, Float.parseFloat(value));
				}
			else if(key.equals("stroke-width") )
				{
				selectors.put(Selector.f_strokeWidth, Float.parseFloat(value));
				}
			else if(key.equals("transparency")  )
				{
				selectors.put(Selector.f_alpha, 1f-Float.parseFloat(value));
				}
			else if(key.equals("shape")  )
				{
				selectors.put(Selector.shapeType, ShapeType.valueOf(value));
				}
			else if(key.equals("track")  )
				{
				selectors.put(Selector.track_name, value);
				}
			else
				{
				LOG.error("unknown selector "+token);
				}
			}
		
		}
	
	
	
	private Data parseData(final String line) {
		final Pattern tab=Pattern.compile("[\t]");
		final String tokens[]=tab.split(line,4);
		if(tokens.length<3) 
			{
			LOG.warning("Not enought tokens in "+line);
			return null;
			}
		final Data data=new Data();
		data.contig=tokens[0];
		if(this.dict.getSequence(data.contig)==null)
			{
			LOG.warn("no such contig in dictionary "+line);
			return null;
			}
		try {
			data.pos=Integer.parseInt(tokens[1]);
		} catch (NumberFormatException e) {
			LOG.warn("bad genomic position in "+line);
			return null;
		}
		data.value = Double.parseDouble(tokens[2]);
		if(tokens.length>3)
			{
			data.style= new Style(data.style);
			data.style.update(tokens[3]);
			}
		final String trackName=String.class.cast(data.style.selectors.get(Selector.track_name));
		if(GenScan.this.defaultTrack!=null)
			{
			if(!(trackName==null || trackName.isEmpty()) )
				{
				LOG.warn("Track name "+trackName+" defined for "+line+" but no track with this name was defined by user.");
				return null;
				}
			data.track = GenScan.this.defaultTrack;
			}
		else
			{
			if(trackName==null || trackName.isEmpty()) 
				{
				LOG.warn("No Track name defined for data "+line);
				return null;
				}
			data.track = this.name2track.get(trackName);
			if(data.track==null)
				{
				LOG.warn("Track name "+trackName+" defined for "+line+" but no track with this name was defined by user.");
				return null;
				}
			}
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
		this.trackNames.stream().
			flatMap(S->Arrays.stream(S.split("[ ,;]"))).
			filter(S->!S.isEmpty()).
			collect(Collectors.toSet()).
				forEach(LABEL->{
				this.name2track.put(LABEL, new Track(LABEL));
				});
		/* no track specified, add a default one */
		if(this.name2track.isEmpty())
			{
			/** define default track */
			this.defaultTrack = new Track("");
			this.name2track.put(this.defaultTrack.label, this.defaultTrack);
			}
		
		final DecimalFormat niceIntFormat = new DecimalFormat("###,###");
		BufferedReader r=null;
		try {
			final Rectangle2D.Double drawingRect = new Rectangle2D.Double(
					this.insets.left,
					this.insets.top,
					WIDTH-(this.insets.left+this.insets.right),
					HEIGHT-(this.insets.top+this.insets.bottom)
					);
			
			final Style styleBase = new Style();
			styleBase.update(this.defaultStyleStr);
			this.styleStack.add(styleBase);
			
			final double adjustedHeight = drawingRect.getHeight() - this.distanceBetweenContigs*(this.name2track.size()-1);
			double y= drawingRect.y;
			for(final Track track:this.name2track.values() )
				{
				track.startY = y;
				track.height = (adjustedHeight/this.name2track.size());
				y += this.distanceBetweenContigs;
				y += track.height;
				}
			
			
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
			
			if(this.displayRanges.isEmpty())
				{
				// no region, add whole genome
				for(final SAMSequenceRecord ssr:this.dict.getSequences()) {
					if(ssr.getSequenceLength()<this.removeContigsSmallerThan)
						{
						LOG.warn("Ignoring "+ssr.getSequenceName()+" because length ="+ssr.getSequenceLength()+"  < "+ this.removeContigsSmallerThan);
						continue;
						}
					this.displayRanges.add(new DisplayContig(ssr));
					}
				if(this.displayRanges.isEmpty())
					{
					LOG.error("No range to display");
					return -1;
					}
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
			g.setColor(ALMOST_WHITE);
			g.fillRect(0, 0, WIDTH, HEIGHT);
			
			//plot background of each contig
			for(int rangeIdx=0;rangeIdx< this.displayRanges.size();++rangeIdx) {
				final DisplayRange displayRange = this.displayRanges.get(rangeIdx);
				int trackidx=0;
				for(final Track track: this.name2track.values()) {
					final Rectangle2D rec= new Pane(displayRange,track).getBounds();
					g.setColor(rangeIdx%2!=trackidx%2?
							ColorUtils.antiquewhite:
							ColorUtils.navajowhite
							);
					g.fill(rec);
					trackidx++;
					}
				}
			
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(this.dict);
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
				
				progress.watch(data.getContig(), data.getStart());
				
				if(data.getValue()< GenScan.this.minY) {continue;}
				if(data.getValue()> GenScan.this.maxY) {continue;}
				if(data.track==null) continue;
				if(data.style.getAlpha()<=0f) continue;

				final Style style= data.getStyle();
				
				for(final DisplayRange dr: GenScan.this.displayRanges) {
					if(!dr.contains(data)) continue;
					final Pane panel = new Pane(dr, data.track);
					
					final Rectangle2D.Double bounds= panel.getBounds();
					
					final Point2D.Double point= panel.convertToPoint(data);
					
					final Shape oldClip = g.getClip();
					final Composite oldComposite = g.getComposite();
					g.setClip(bounds);
					//
					g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, style.getAlpha()));
					final Shape dataShape;
					final double shape_size= style.getSize();
					switch(style.getShapeType())
						{
						case circle:
							dataShape = new Ellipse2D.Double(
									point.x-shape_size/2.0,
									point.y-shape_size/2.0,
									shape_size,
									shape_size
									);
							break;
						case square:
							dataShape = new Rectangle2D.Double(
									point.x-shape_size/2.0,
									point.y-shape_size/2.0,
									shape_size,
									shape_size
									);
							break;
						default:throw new IllegalArgumentException();
						}
					if(style.getPaperColor()!=null)
						{
						g.setColor(style.getPaperColor());
						g.fill(dataShape);
						}
					final float stroke_width = style.getStrokeWidth();
					if(style.getPenColor()!=null && stroke_width>0)
						{
						final Stroke oldStroke = g.getStroke();
						final Stroke stroke = new BasicStroke(
								stroke_width,
								BasicStroke.CAP_BUTT,
								BasicStroke.JOIN_ROUND
								);
						g.setStroke(stroke);
						g.setColor(style.getPenColor());
						g.draw(dataShape);
						g.setStroke(oldStroke);
						}
					
					
					//
					g.setClip(oldClip);
					
					g.setComposite(oldComposite);
					}
				}
			progress.finish();
			
			double lastTickX=0.0;
			//plot frame of each contig
			for(final DisplayRange displayRange : this.displayRanges) {
				for(final Track track:this.name2track.values())
					{
					final Rectangle2D rec= new Pane(displayRange,track).getBounds();
					g.setColor(ALMOST_BLACK);
					g.draw(rec);
					}
				final String label= displayRange.getLabel();
				int fontSize=insets.top-2;
				
				final double labelLength=Math.min( label.length()*fontSize,displayRange.width);
					if( labelLength > this.distanceBetweenContigs) {
					g.setColor(ALMOST_BLACK);
					this.hershey.paint(g, label, new Rectangle2D.Double(
							displayRange.startX + displayRange.width/2.0 -labelLength/2.0,
							1,
							labelLength,
							(insets.top-2)
							));
					}
				//plot ticks
				for(int i=0;i<= 10;i++)
					{
					double bp = displayRange.getStart()+ (displayRange.length()/10.0)*i; 
					double tx = displayRange.convertToX(bp);
					double ty = drawingRect.getMaxY();
					g.setColor(ALMOST_BLACK);
					g.draw(new Line2D.Double(tx,ty, tx, ty+5));
					if(tx - lastTickX>5) {
						g.translate(tx, ty+5);
						g.rotate(Math.PI/2.0);
						g.drawString(niceIntFormat.format((long)bp), 0, 0);
						g.rotate(-Math.PI/2.0);
						g.translate(-tx, -(ty+5));
						lastTickX=tx+5;
						}
					
					
					}
				}

			//plot y axis
			for(final Track track:this.name2track.values()) {
				for(int i=0;i<= 10;++i)
					{
					double v= this.minY+ ((this.maxY-this.minY)/10.0)*i;
					double ty= track.convertToY(v);
					double font_size=12;
					
					// horizontal line
					g.setColor(Color.LIGHT_GRAY);
					g.draw(new Line2D.Double(
							drawingRect.getX() - 5,
							ty,
							drawingRect.getMaxX(),
							ty));
					
					g.setColor(ALMOST_BLACK);
					this.hershey.paint(g, String.format("%8s",v), 
						new Rectangle2D.Double(
							0,
							ty-font_size,
							insets.left,
							font_size
							));
					}
				// plot track label if no default track
				if(this.defaultTrack==null)
					{
					g.setColor(ALMOST_BLACK);
					g.translate(WIDTH-insets.right+2, track.startY);
					g.rotate(Math.PI/2.0);
					g.drawString(track.label, 0, 0);
					g.rotate(-Math.PI/2.0);
					g.translate(-(WIDTH-insets.right+2), -(track.startY));
					}
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
	
	
	public static void main(final String[] args) {
		new GenScan().instanceMainWithExit(args);
	}

}
