package com.github.lindenb.jvarkit.tools.burden;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.function.Consumer;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;
import javax.script.CompiledScript;
import javax.script.ScriptException;
import javax.script.SimpleBindings;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.util.CloserUtil;

@Program(
		name="casectrlcanvas",
		description="draw a chart of case/control maf from a stream of X/Y values",
		keywords={"vcf","case","control","visualization","jfx","chart","maf"}
		)
public class CaseControlCanvas implements Consumer<Point2D> {
	private static final Logger LOG = Logger.build(CaseControlCanvas.class).make();
	private static final Color ALMOST_BLACK=new Color(20,20,20);
	private static final Color ALMOST_WHITE=new Color(240,240,240);
	
	
	public static interface PainterXY
		{
		public void paint(Graphics2D gc,double x,double y);
		}
	
	public static class DefaultPainterXY implements PainterXY
		{
		private Color pen = Color.ORANGE;
		private Color paper = Color.RED;
		private float opacity=0.6f;
		private double size=10;
		
		public float getOpacity() {
			return opacity;
		}
		
		public double getSize() {
			return size;
		}
		protected Shape createShape(double x,double y)
			{
			return new Ellipse2D.Double(x-size/2.0, y-size/2.0, size,size);	
			}
		@Override
		public void paint(final Graphics2D gc,double x,double y) {
			final Composite old=gc.getComposite();
			gc.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER));
			final Shape shape=createShape(x, y);
			if(paper!=null)
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
	
	
	public CaseControlCanvas(final Config cfg) {
		this.cfg=cfg;
		this.insets = new Insets(10, 10, 10, 10);
		this.canvas = new BufferedImage(
			cfg.width+(this.insets.left+this.insets.right),
			cfg.width+(this.insets.top+this.insets.bottom),
			BufferedImage.TYPE_INT_RGB);
		this.gc = this.canvas.createGraphics();
		this.gc.setColor(ALMOST_WHITE);
		this.gc.fillRect(0,0,this.canvas.getWidth(),this.canvas.getHeight());
		
		this.gc.setColor(ALMOST_BLACK);
		this.gc.drawRect(
				this.insets.left,
				this.insets.top,
				cfg.width,
				cfg.width
				);
		this.gc.setColor(Color.LIGHT_GRAY);
		for(int i=1;i< 10.0;i++) {
			int x= (int)(this.insets.left+(cfg.width/10.0)*i);
			this.gc.drawLine( x, this.insets.top, x, this.insets.top+cfg.width);
			int y= (int)(this.insets.top+(cfg.width/10.0)*i);
			this.gc.drawLine(this.insets.left, y,this.insets.left+cfg.width,y);
			}
		this.gc.drawLine(
				this.insets.left, 
				this.insets.top+cfg.width,
				this.insets.left+cfg.width,
				this.insets.top
				);

		this.gc.setColor(ALMOST_BLACK);
		this.gc.drawString("Cases",0,this.canvas.getHeight()-10);
		
		AffineTransform orig = this.gc.getTransform();
		this.gc.translate(5, 5);
		this.gc.rotate(Math.PI/2);
		this.gc.drawString("Controls",0,0);
		this.gc.setTransform(orig);
		}
	
	public CaseControlCanvas()
		{
		this(new Config());
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
	private int width=600;
	@Parameter(names={"--opacityexpr"},description="Javascript expression for opacity.")
	private String opacityexpr = "0.6";
	@Parameter(names={"--colorexpr"},description="Javascript expression for opacity.")
	private String colorexpr = "black";

	private final ColorUtils colorUtils=new ColorUtils();
	private CompiledScript colorScript =  null;
	private CompiledScript opacityScript = null;
	
	public Color getColor(final SimpleBindings bindings) 
		{
		if(colorScript==null) {
			colorScript = compileExpression(this.colorexpr);
			}
		final Object o;
		
		try {
			o= colorScript.eval(bindings);
		} catch (ScriptException e) {
			throw new JvarkitException.ScriptingError(e);
		}
		
		if(o==null)
			{
			throw new JvarkitException.ScriptingError("script returned null color");
			}
		if(o instanceof Color) {
			return (Color)o;
			}
		else
			{
			return colorUtils.parse(String.valueOf(o));
			}
		}
	public float Opacity(final SimpleBindings bindings) 
		{
		if(opacityScript==null) {
			opacityScript = compileExpression(this.opacityexpr);
			}
		final Object o;
		try {
			o= opacityScript.eval(bindings);
		} catch (ScriptException e) {
			throw new JvarkitException.ScriptingError(e);
			}
		if(o==null)
			{
			throw new JvarkitException.ScriptingError("script returned null opacity");
			}
		if(o instanceof Number) {
			return Number.class.cast(o).floatValue();
			}
		else
			{
			return Float.valueOf(String.valueOf(o));
			}
		}
	private CompiledScript compileExpression(final String jsExpr) {
		final javax.script.ScriptEngineManager manager = new javax.script.ScriptEngineManager();
		final javax.script.ScriptEngine engine = manager.getEngineByName("js");
		if(engine==null)
			{
			throw new JvarkitException.JavaScriptEngineNotFound();
			}
		final javax.script.Compilable compilingEngine = (javax.script.Compilable)engine;
		try {
			CompiledScript compiled = compilingEngine.compile(jsExpr);
			return compiled;
			}
		catch(final Exception err)
			{
			throw new RuntimeException(err);
			}
		}

	}	
	
public static class Main extends Launcher
	{
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;

	@ParametersDelegate
	private Config contig = new Config();
	
	private CaseControlCanvas instance=new  CaseControlCanvas();
	
	@Override
	public int doWork(final List<String> args) {
		BufferedReader r=null;
		try {
			final Pattern tab = Pattern.compile("[\t]");
			r = super.openBufferedReader(oneFileOrNull(args));
			String line;
			while((line=r.readLine())!=null) {
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
				if(values[0]<0 || values[0]>1.0) {
					LOG.error("Bad X: ignoring "+line);
					continue;
					}
				if(values[1]<0 || values[1]>1.0) {
					LOG.error("Bad Y:ignoring "+line);
					continue;
					}
				final Point2D point = new Point2D.Double(values[0],values[1]);
				instance.accept(point);
				}
			r.close();
			r=null;
			instance.saveAs(outputFile);
			return 0;
		} catch(Exception err) {
			LOG.error(err);
			return -1;
		} finally
		{
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
