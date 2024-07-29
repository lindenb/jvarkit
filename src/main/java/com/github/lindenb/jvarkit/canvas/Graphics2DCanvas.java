package com.github.lindenb.jvarkit.canvas;

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.nio.file.Path;

import javax.imageio.ImageIO;

import com.github.lindenb.jvarkit.util.FunctionalMap;

public class Graphics2DCanvas extends Canvas {
	private final int width;
	private final int height;
	private final Path outputFile;
	private final BufferedImage image;
	private final String outputFormat;
	private Graphics2D g2d;
	public Graphics2DCanvas(Path out,int width,int height,String format,FunctionalMap<String,Object> params) {
		this.width=width;
		this.height=height;
		this.outputFormat=format;
		this.image= new BufferedImage(width, height,format.equalsIgnoreCase("png")?BufferedImage.TYPE_INT_ARGB:BufferedImage.TYPE_INT_RGB);
		this.g2d = this.image.createGraphics();
		this.g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
		this.outputFile = out;
		}
	
	private BasicStroke getStroke() {
		return new BasicStroke(
				(float)getStrokeWidth(),
				getLineCap(),
				getLineJoin()
				);
		}
	

	
	@Override
	public Canvas text(String s, double x, double y, FunctionalMap<String, Object> fm) {
		this.states.push(this.states.peek().plus(fm));

		Color c= getFillColor();
		if(c!=null) {
			this.g2d.setColor(c);
			this.g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, (float)getFillOpacity()));
			this.g2d.drawString(s, (float)x, (float)y);
			}
		c= getStrokeColor();
		BasicStroke stroke = getStroke();
		if(stroke.getLineWidth() > 0f && c!=null) {
			this.g2d.setStroke(stroke);
			this.g2d.setColor(c);
			this.g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, (float)getStrokeOpacity()));
			this.g2d.drawString(s, (float)x, (float)y);
			}
		
		this.states.pop();
		return this;
		}
	
	@Override
	public Canvas shape(Shape shape, FunctionalMap<String, Object> fm) {
		this.states.push(this.states.peek().plus(fm));
		
		Color c= getFillColor();
		if(c!=null) {
			this.g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, (float)getFillOpacity()));
			this.g2d.setColor(c);
			this.g2d.fill(shape);
			}
		c= getStrokeColor();
		BasicStroke stroke = getStroke();
		if(stroke.getLineWidth() > 0f && c!=null) {
			this.g2d.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, (float)getStrokeOpacity()));
			this.g2d.setStroke(stroke);
			this.g2d.setColor(c);
			this.g2d.draw(shape);
			}

		this.states.pop();
		return this;
		}
	
	public BufferedImage getImage() {
		return image;
		}
	
	@Override
	public void close() throws IOException {
		if(this.outputFile==null)  {
			ImageIO.write(image,this.outputFormat, System.out);
			}
		else
			{
			ImageIO.write(image,this.outputFormat, this.outputFile.toFile());
			}
		}
	@Override
	public int getWidth() {
		return width;
		}
	@Override
	public int getHeight() {
		return height;
		}
	}
