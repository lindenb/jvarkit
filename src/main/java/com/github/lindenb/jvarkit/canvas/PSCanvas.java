package com.github.lindenb.jvarkit.canvas;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.PathIterator;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;


import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.FunctionalMap;




public class PSCanvas extends Canvas {
	private final int width;
	private final int height;
	private final PrintWriter w;
	public PSCanvas(Path out,int width,int height,boolean compressed,FunctionalMap<String,Object> params) throws IOException {
		super();
		this.width=width;
		this.height=height;
		
		
		if(out==null) {
			w = compressed ? 
					new PrintWriter(new GZIPOutputStream(System.out)):
					new PrintWriter(System.out)
					;
			}
		else
			{
			w = compressed ? 
					new PrintWriter(new GZIPOutputStream(Files.newOutputStream(out))):
					new PrintWriter(Files.newOutputStream(out))
					;
			}
		this.w.println("%PS");
		this.w.println("%%BoundingBox: 0 0 "+width+" "+height);
		}
	
	@Override
	public void close() throws IOException {
		this.w.println(" showpage");
		this.w.println("%EOF");
		this.w.flush();
		this.w.close();
		}
	String coord(double x,double y) {
		return String.valueOf(x)+" "+String.valueOf(getHeight()-y);
		}
	
	private String setrgbcolor(Color c) {
		if(c==null) return "";
		if(c.getRed()==c.getGreen() && c.getGreen()==c.getBlue()) {
			return String.valueOf(c.getRed()/255f)+" setgray";
			}
		
		return String.join(" ",
				String.valueOf(c.getRed()/255f),
				String.valueOf(c.getGreen()/255f),
				String.valueOf(c.getBlue()/255f),
				"setrgbcolor"
				);
		}
	
	private Canvas fillAndStroke() {
		Color c= getFillColor();
		if(c!=null) {
			w.print(" gsave");
			w.print(" ");
			w.print(setrgbcolor(c));
			w.print(" fill");
			w.print(" grestore");
			}
		
		double linewidth=getStrokeWidth();
		if(linewidth>0) {
			c= getStrokeColor();
			if(c!=null) {
				w.print(" ");
				w.print(String.valueOf(linewidth));
				w.print(" setlinewidth ");
				w.print(setrgbcolor(c));
				w.print(" stroke");
				}
			}
		
		return this;
		}
	
	private String escape(final CharSequence s) {
		final StringBuilder sb = new StringBuilder(s.length());
		for(int i=0;i< s.length();i++)
			{
			final char c = s.charAt(i);
			switch(c) {
				case '\n' : sb.append("\\n");break;
				case '\r' : sb.append("\\r");break;
				case '\t' : sb.append("\\t");break;
				case '\\' : sb.append("\\\\");break;
				case '\'' : sb.append("\\\'");break;
				case '\"' : sb.append("\\\"");break;
				case '(' : sb.append("\\(");break;
				case ')' : sb.append("\\)");break;
				default:sb.append(c);break;
				}
			}
		return sb.toString();
		}
	@Override
	public Canvas text(String text, double x, double y, FunctionalMap<String, Object> fm) {
		if(StringUtils.isBlank(text)) return this;
		fm = this.states.push(this.states.peek().plus(fm));
		w.print(" /"+ fm.getOrDefault(KEY_FONT_FAMILY, "Times-Roman")+" findfont");
		w.print(" "+ fm.getOrDefault(KEY_FONT_SIZE, "12")+" scalefont");
		w.print(" setfont newpath ");
		w.print(coord(x,y)+" moveto");
		w.print(" ("+escape(text)+") show");
		this.states.pop();
		return this;
		}
	
	@Override
	public Canvas comment(String s) {
		if(!StringUtils.isBlank(s)) {
			w.print("\n% ");
			w.print(s);
			w.println();
			}
		return this;
		}
	
	@Override
	public Canvas shape(Shape shape, FunctionalMap<String, Object> fm) {
		fm = this.states.push(this.states.peek().plus(fm));
		float coords[]=new float[6];
		w.append(" newpath");
		PathIterator iter = shape.getPathIterator(null);
		while(!iter.isDone()) {
			switch(iter.currentSegment(coords))
			{
			case PathIterator.SEG_MOVETO:
				{
				w.append(" ");
				w.append(coord(coords[0],coords[1]));
				w.append(" moveto");
				break;
				}
			case PathIterator.SEG_LINETO:
				{
				w.append(" ");
				w.append(coord(coords[0],coords[1]));
				w.append(" lineto");
				break;
				}
			case PathIterator.SEG_QUADTO:
				{
				//TODO
				break;
				}
			case PathIterator.SEG_CUBICTO:
				{
				w.append(" ");
				w.append(coord(coords[0],coords[1]));
				w.append(" ");
				w.append(coord(coords[2],coords[3]));
				w.append(" ");
				w.append(coord(coords[4],coords[5]));
				w.append(" curveto");
				break;
				}
			case PathIterator.SEG_CLOSE:
				{
				w.append(" closepath");
				break;
				}
			}
		
		iter.next();
		}
		fillAndStroke();
		this.states.pop();
		return this;
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
