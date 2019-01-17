/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.burden;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
/**
BEGIN_DOC

## Example

Create Exac/Gnomad from a VCF annotated with both databases:

```
java -jar casectrlcanvas.jar -caseAtt AC_NFE/AN_NFE -ctrlAtt gnomad_exome_AF_NFE -o exome_gnomad_vs_exac.png.NFE.png -title 'Case:ExacAC/AN_NFE Ctrl:gnomad_exome_AF_NFE' -opacity 0.2 input.vcf
```


Create Exac/my-controls from a VCF annotated with both databases:

```
java -jar casectrlcanvas.jar -p myped.ped -caseAtt AC_NFE/AN_NFE   -o out.png -opacity 0.2 input.vcf
```

## Note to self: create a mosaic of images:

```
montage -geometry 1000x1000+2+2 file1.png file2.png fileN.png out.png
```


## Note to self: Running in headless mode (no X11 available)

see [http://www.oracle.com/technetwork/articles/javase/headless-136834.html](http://www.oracle.com/technetwork/articles/javase/headless-136834.html)

```
java -Djava.awt.headless=true -jar casectrlcanvas.jar ....
```


END_DOC

 */
@Program(
		name="casectrlcanvas",
		description="draw a chart of case/control maf from a stream of X/Y values",
		keywords={"vcf","case","control","visualization","jfx","chart","maf"}
		)
public class CaseControlCanvas implements Consumer<Point2D> {
	private static final Logger LOG = Logger.build(CaseControlCanvas.class).make();
	private static final Color ALMOST_BLACK = new Color(20,20,20);
	private static final Color ALMOST_WHITE = new Color(240,240,240);
	
	
	public static interface PainterXY
		{
		public void paint(Graphics2D gc,double x,double y);
		}
	
	public static class DefaultPainterXY implements PainterXY
		{
		enum ShapeType { oval,square,cross}
		@Parameter(names={"--pen","--foreground"},description="pen color. "+ColorUtils.Converter.OPT_DESC,converter=ColorUtils.Converter.class)
		private Color pen = Color.ORANGE;
		@Parameter(names={"--paper","--background"},description="background color. "+ColorUtils.Converter.OPT_DESC,converter=ColorUtils.Converter.class)
		private Color paper = Color.RED;
		@Parameter(names={"-opacity","--opacity"},description="opacity [0-1]")
		private float opacity=0.6f;
		@Parameter(names={"--pointsize"},description="points width")
		private double size=10;
		@Parameter(names={"--pointshape"},description="Point Shape")
		private ShapeType shapeType=ShapeType.oval;
		
		DefaultPainterXY(final DefaultPainterXY cp)
			{
			this.shapeType=cp.shapeType;
			this.pen=cp.pen;
			this.paper=cp.paper;
			this.opacity=cp.opacity;
			this.size=cp.size;
			}
		
		DefaultPainterXY()
			{
			}

		public float getOpacity() {
			return opacity;
		}
		
		public double getSize() {
			return size;
		}
		protected Shape createShape(double x,double y)
			{
			switch(this.shapeType)
				{
				case cross : GeneralPath p=new GeneralPath();
					p.moveTo(x-size/2, y);
					p.lineTo(x+size/2, y);
					p.moveTo(x,y-size/2);
					p.lineTo(x,y+size/2);
					return p;
				case square: return new Rectangle2D.Double(x-size/2.0, y-size/2.0, size,size);
				default: return new Ellipse2D.Double(x-size/2.0, y-size/2.0, size,size);	
				}
			}
		@Override
		public void paint(final Graphics2D gc,double x,double y) {
			final Composite old=gc.getComposite();
			gc.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,this.opacity));
			final Shape shape=createShape(x, y);
			if(paper!=null && this.shapeType!=ShapeType.cross)
				{
				gc.setColor(paper);
				gc.fill(shape);
				}
			if(pen!=null)
				{
				gc.setColor(pen);
				gc.draw(shape);
				}
			gc.setComposite(old);
			}
		
		}
	
	
	private final BufferedImage canvas;
	private Graphics2D gc;
	private Insets insets;
	private final Config cfg;
	private PainterXY painterXY = new DefaultPainterXY();
	
	
	public CaseControlCanvas(String title,final Config cfg) {
		this.cfg=cfg;
		this.insets = new Insets(30, 30, 30, 10);
		this.canvas = new BufferedImage(
			cfg.width+(this.insets.left+this.insets.right),
			cfg.width+(this.insets.top+this.insets.bottom),
			BufferedImage.TYPE_INT_RGB);
		this.gc = this.canvas.createGraphics();
		this.gc.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		this.gc.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		this.gc.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);

		this.gc.setColor(ALMOST_WHITE);
		this.gc.fillRect(0,0,this.canvas.getWidth(),this.canvas.getHeight());
		
		this.gc.setColor(ALMOST_BLACK);
		Font font=new Font("Times", Font.PLAIN, 16);
		this.gc.setFont(font);
		if(title==null) title="Untitled";
		this.gc.drawString(title, this.insets.top+cfg.width/2 - title.length()/2*12,16);
		
		font=new Font("Times", Font.PLAIN, 7);
		this.gc.setFont(font);
		for(int i=0;i<= 10.0;i++) {
			this.gc.setColor(Color.LIGHT_GRAY);
			int x= (int)(this.insets.left+(cfg.width/10.0)*i);
			this.gc.drawLine( x, this.insets.top, x, this.insets.top+cfg.width);
			int y= (int)(this.insets.top+(cfg.width/10.0)*i);
			this.gc.drawLine(this.insets.left, y,this.insets.left+cfg.width,y);
			this.gc.setColor(ALMOST_BLACK);
			this.gc.drawString(String.format("%.1f",( i/10.0)),x, this.insets.top+cfg.width+8);
			this.gc.drawString(String.format("%.1f",(1.0 - i/10.0)),this.insets.left-11,y);
			}
		this.gc.drawLine(
				this.insets.left, 
				this.insets.top+cfg.width,
				this.insets.left+cfg.width,
				this.insets.top
				);
	
		this.gc.setColor(ALMOST_BLACK);
		this.gc.drawRect(
				this.insets.left,
				this.insets.top,
				cfg.width,
				cfg.width
				);
		
		font=new Font("Times", Font.PLAIN, 12);
		this.gc.setFont(font);
		this.gc.setColor(ALMOST_BLACK);
		this.gc.drawString("Cases", this.insets.left+cfg.width/2,this.insets.top+cfg.width+20);
		
		AffineTransform orig = this.gc.getTransform();
		this.gc.translate(7, 0);
		this.gc.rotate(Math.PI/2);
		this.gc.drawString("Controls", this.insets.top+cfg.width/2,0);
		this.gc.setTransform(orig);
		}
	
	public CaseControlCanvas()
		{
		this("",new Config());
		}
	
	public void setPainter(final PainterXY painterXY) {
		this.painterXY = painterXY;
		}
	
	@Override
	public void accept(final Point2D pt) {
		if(pt==null) throw new IllegalArgumentException("null");
		if(gc==null) throw new IllegalArgumentException("gc disposed");
		if(pt.getX()<0 || pt.getX()>1.0 || pt.getY()<0 || pt.getY()>1.0) return;
		
		double x = this.insets.left+(pt.getX()*this.cfg.width);
		double y = this.insets.top + this.cfg.width - (pt.getY()*this.cfg.width);
		this.painterXY.paint(this.gc, x, y);
		}
	
	public void saveAs(final File fileOrStdout) throws IOException{
            this.gc.dispose();
            this.gc=null;
            if(fileOrStdout!=null) {
	            String format="png";
	            if((fileOrStdout.getName().toLowerCase().endsWith(".jpg") || fileOrStdout.getName().toLowerCase().endsWith(".jpeg")))
	            	{
	            	format="jpg";
	            	}
            	ImageIO.write(this.canvas, format, fileOrStdout);
	            }
            else
            	{
            	ImageIO.write(this.canvas, "png", System.out);
            	}
		}
	
/** configuration for multiple canvas */
public static class Config
	{
	@Parameter(names={"--width"},description="Canvas width")
	private int width=700;

	}	
	
public static class Main extends Launcher
	{
	@Parameter(names={"-o","--out"},description="Image file name. "+OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-txt","--txt","--text","-tsv","--tsv"},description="Input is a tab delimited file. containg x=case and y=controls")
	private boolean inputIsText=false;
	@Parameter(names={"-tee","--tee"},description="Tee input to stdout, useful in linux pipelines to get intermediary results. Requires that -o 'file' is set.")
	private boolean teeStdout=false;
	@Parameter(names={"-xyAttribute","--xyAttribute"},description="When using 'tee', add this Attribute containing the MAF for case and control")
	private String outputXYAttribute="MAFCASECTRL";
	@Parameter(names={"-format","--format"},description="How to print doubles, printf-like precision format.")
	private String precisionFormat="%.5f";
	@Parameter(names={"-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION+" If not defined, I will try to extract the pedigree from the VCFheader.")
	private File pedigreeFile;
	@Parameter(names={"-caseAtt","--caseAttribute"},description="Do not calculate MAF for cases, but use this tag to get Controls' MAF. Notation 'AC/AN' will use two attributes")
	private String controlTag =null;
	@Parameter(names={"-ctrlAtt","--ctrlAttribute"},description="Do not calculate MAF for controls, but use this tag to get Cases' MAF. Notation 'AC/AN' will use two attributes")
	private String caseTag =null;
	@Parameter(names={"-title","--title"},description="Title")
	private String title ="";
	@ParametersDelegate
	private MafCalculator.Factory mafCalcFactory = new MafCalculator.Factory();
	@ParametersDelegate
	private DefaultPainterXY defaultPainter=new DefaultPainterXY();
	@ParametersDelegate
	private Config configuration = new Config();
	
	private static interface MafExtractor
		{
		public Double apply(final VariantContext ctx,Allele alt);
		}

	private static class AttributeMafExtractor implements MafExtractor
		{
		private boolean reportedError=false; 
		private final String attribute;
		public AttributeMafExtractor(final String tag) {this.attribute = tag;}
		public String getAttribute()
			{
			return attribute;
			}
		@Override
		public Double apply(final VariantContext ctx,final Allele alt) {
			final String att = this.getAttribute();
			if(att==null || att.isEmpty()) {
				if(!reportedError) {
					LOG.warn("Attribute missing for AttributeMafExtractor");
					reportedError=true;
					}
				return null;
				}
			if(!ctx.hasAttribute(att)) return null;
			
			int index=ctx.getAlleleIndex(alt);
			if( index<=0) {
				LOG.warn("Allele="+alt+"/attribute=["+att+"]strange.. asking for REF allele ?? index="+index);
				return null;//can't be REF==0
				}
			// shift -1 because REF=0 and we'ere looking at INFO.TYPE=A, so first ALT allele is index=0
			index--;
			
			final List<Object> L = ctx.getAttributeAsList(att);
			if(index>=L.size()) {
				return null;
				}
			final Object o = L.get(index);
			if(o==null || ".".equals(o)) return null;
			try {
				double f = Double.parseDouble(String.valueOf(o));
				if(f<0.0 || f>1.0) return null;
				
				return f;
				}
			catch(NumberFormatException err) {
				LOG.info("cannot parse MAF for attribute "+att);
				return null;
				}
			}
		}
	
	
	private static class AC_DIV_AN_MafExtractor implements MafExtractor
		{
		private final String acAttribute;
		private final String anAttribute;
		public AC_DIV_AN_MafExtractor(final String tag) {
			int slash=tag.indexOf("/");
			if(slash==-1) throw new IllegalArgumentException(tag);
			this.acAttribute=tag.substring(0,slash);
			this.anAttribute=tag.substring(slash+1);
			}
		@Override
		public Double apply(final VariantContext ctx,final Allele alt) {
			if(!ctx.hasAttribute(this.anAttribute))
				{
				return null;
				}
			if(!ctx.hasAttribute(this.acAttribute)){
				return null;
				}
		
			final double an;
			try 
				{
				an = ctx.getAttributeAsDouble(this.anAttribute, 0);
				}
			catch(NumberFormatException err)
				{
				return null;
				}
			if(an<=0) return null;
			
			
			int index=ctx.getAlleleIndex(alt);
			if( index<=0) {
				LOG.warn("Allele="+alt+"/attribute=["+this.acAttribute+"]strange.. asking for REF allele ?? index="+index);
				return null;//can't be REF==0
				}
			// shift -1 because REF=0 and we'ere looking at INFO.TYPE=A, so first ALT allele is index=0
			index--;
			
			final List<Object> L = ctx.getAttributeAsList(this.acAttribute);
			if(index>=L.size()) {
				LOG.error("boum?");
				return null;
				}
			final Object o = L.get(index);
			if(o==null || ".".equals(o)) return null;
			try {
				final double ac = Double.parseDouble(String.valueOf(o));
				double f= ac/(double)an;
				if(f<0.0 || f>1.0) return null;
				return f;
				}
			catch(NumberFormatException err) {
				LOG.info("cannot parse MAF for attribute "+this.acAttribute);
				return null;
				}
			}
		}

	
	private  class GenotypeMafExtractor implements MafExtractor
		{
		private final Set<Pedigree.Person> samples;
		GenotypeMafExtractor(final Set<Pedigree.Person> samples)
			{
			this.samples = samples;
			if(this.samples.isEmpty()) throw new IllegalArgumentException("samples is empty");
			}
		@Override
		public Double apply(final VariantContext ctx,final Allele alt) {
			final MafCalculator calc= mafCalcFactory.create(alt, ctx.getContig());
			for(final Pedigree.Person sample: this.samples) {
				calc.add(ctx.getGenotype(sample.getId()), sample.isMale());
				}
			return calc.isEmpty()?null:calc.getMaf();
			}
		}
	
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.teeStdout && this.outputFile==null) {
			LOG.error("cannot 'tee' while output file is stdout");
			return -1;
			}
		BufferedReader r=null;
		VCFIterator in = null;
		VariantContextWriter teeVariantWriter = null;
		try {
			final CaseControlCanvas instance=new  CaseControlCanvas(this.title,this.configuration);
			if(!this.inputIsText) {
				instance.setPainter(this.defaultPainter);
				in = super.openVCFIterator(oneFileOrNull(args));
				
				final VCFHeader header= in.getHeader();
				final Pedigree pedigree;
				
				if( this.pedigreeFile!=null) {
					pedigree = Pedigree.newParser().parse(this.pedigreeFile);
					}
				else
					{
					pedigree = Pedigree.newParser().parse(header);
					}
				if(this.caseTag==null && this.controlTag==null && (pedigree==null || pedigree.isEmpty()) ) {
					LOG.error("No pedigree defined , or it is empty");
					return -1;
					}
				final Set<Pedigree.Person> casepersons = pedigree.getPersons().
						stream().
						filter(F->F.isAffected() && header.getSampleNamesInOrder().indexOf(F.getId())!=-1).
						collect(Collectors.toSet());
						
				if(this.caseTag==null && casepersons.isEmpty()){
						LOG.error("No Affected individuals in pedigree/header");
						return -1;
						}
				
				final Set<Pedigree.Person> controlpersons = pedigree.getPersons().
						stream().
						filter(F->F.isUnaffected() && header.getSampleNamesInOrder().indexOf(F.getId())!=-1).
						collect(Collectors.toSet());
				
				if(this.controlTag==null && controlpersons.isEmpty()){
						LOG.error("No Unaffected individuals in pedigree/header");
						return -1;
						}
				final VCFInfoHeaderLine xyInfoHeaderLine;
				if( this.teeStdout)
					{
					teeVariantWriter = super.openVariantContextWriter(null);
					xyInfoHeaderLine=new VCFInfoHeaderLine(
							this.outputXYAttribute,
							VCFHeaderLineCount.A,
							VCFHeaderLineType.String,
							"MAF Cases|MAF Controls"
							);
					header.addMetaDataLine(xyInfoHeaderLine);
					teeVariantWriter.writeHeader(header);
					}
				else
					{
					xyInfoHeaderLine=null;
					}
				
				final MafExtractor casesMafExtractor = (
						this.caseTag==null?
							new GenotypeMafExtractor(casepersons):
							this.caseTag.contains("/")?
							new AC_DIV_AN_MafExtractor(this.caseTag):
							new AttributeMafExtractor(this.caseTag)
						);
				final MafExtractor controlMafExtractor = (
						this.controlTag==null?
								new GenotypeMafExtractor(controlpersons):
								this.controlTag.contains("/")?
								new AC_DIV_AN_MafExtractor(this.controlTag):
								new AttributeMafExtractor(this.controlTag)
							);
				
				final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
				while(in.hasNext())
					{
					final VariantContext ctx = progress.watch(in.next());
					final List<String> attributes=new ArrayList<>();
					final List<Allele> alts = ctx.getAlternateAlleles();
					for(int altidx=0; altidx < alts.size();++altidx)
						{
						final Allele alt= alts.get(altidx);
						final Double casex = casesMafExtractor.apply(ctx, alt);
						if(casex==null || casex<0.0 || casex>1.0) {
							attributes.add(".");
							continue;
							}
						final Double ctrly = controlMafExtractor.apply(ctx, alt);
						if(ctrly==null || ctrly<0.0 || ctrly>1.0) {
							attributes.add(".");
							continue;
							}
						attributes.add(String.format(precisionFormat,casex)+"|"+String.format(precisionFormat,ctrly));
						instance.accept( new Point2D.Double(casex,ctrly));
						}
					
					if(teeVariantWriter!=null)
						{
						if(attributes.isEmpty() )
							{
							teeVariantWriter.add(ctx);
							}
						else
							{
							teeVariantWriter.add(
								new VariantContextBuilder(ctx).
									attribute(xyInfoHeaderLine.getID(),attributes).
									make()
								);
							}
						}
					
					}
				
				progress.finish();
				if(teeVariantWriter!=null) {
					teeVariantWriter.close();
					teeVariantWriter=null;
					}
				in.close();in=null;
				}
			else
				{
				final ColorUtils colorUtils=new ColorUtils();
				final Pattern tab = Pattern.compile("[\t]");
				r = super.openBufferedReader(oneFileOrNull(args));
				String line;
				while((line=r.readLine())!=null) {
					if(this.teeStdout)
						{
						this.stdout().println(line);
						}
					if(line.trim().isEmpty() || line.startsWith("#")) continue;
					final String tokens[]=tab.split(line);
					if(tokens.length<2) {
						LOG.error("Bad line in "+line);
						return -1;
						}
					double values[]=new double[2];
					for(int i=0;i<2;i++)
						{
						try {
							values[i]=Double.parseDouble(tokens[i]);
							}
						catch(final NumberFormatException err) {
							LOG.error(err);
							return -1;
							}
						}
					 DefaultPainterXY painter=new DefaultPainterXY(this.defaultPainter);
					 if(tokens.length>2)
					 	{
						for(final String kv:tokens[2].split(";"))
							{
							if(kv.isEmpty()) continue;
							int eq=kv.indexOf("=");
							if(eq==-1) {
								LOG.warn("'=' missing in "+line);
								continue;
								}
							final String key = kv.substring(0,eq).trim().toLowerCase();
							final String value = kv.substring(eq+1).trim();
							if(key.equals("size")) {
								painter.size=Integer.parseInt(value);
								}
							else if(key.equals("shape"))
								{
								painter.shapeType=DefaultPainterXY.ShapeType.valueOf(value);
								}
							else if(key.equals("pen") || key.equals("color") || key.equals("foreground"))
								{
								painter.pen=colorUtils.parse(value);
								}
							else if(key.equals("paper")|| key.equals("background"))
								{
								painter.paper=colorUtils.parse(value);
								}
							else if(key.equals("opacity"))
								{
								painter.opacity=Float.parseFloat(value);
								}
							}
					 	}
					instance.setPainter(painter);
					final Point2D point = new Point2D.Double(values[0],values[1]);
					instance.accept(point);
					}
				r.close();
				r=null;
				}
			instance.saveAs(outputFile);
			return 0;
		} catch(Exception err) {
			LOG.error(err);
			return -1;
		} finally
		{
			CloserUtil.close(in);
			CloserUtil.close(teeVariantWriter);
			CloserUtil.close(r);
		}
	}
	
	public static void main(final String[] args) {
			new Main().instanceMainWithExit(args);
		}
	}
public static void main(String[] args) {
	Main.main(args);
	}
}
