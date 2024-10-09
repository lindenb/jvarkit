/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.canvas;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
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
		this.w.println("% generated with jvarkit. Author: Pierre Lindenbaum PhD.");
		this.w.println("%EOF");
		this.w.flush();
		this.w.close();
		}
	
	private String round(double x) {
		String s= String.format("%.2f",x);
		if(s.endsWith(".00")) s=s.substring(0,s.length()-3);
		return s;
	}
	
	private String coord(double x,double y) {
		return round(x)+" "+round(getHeight()-y);
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
	
	private Canvas fillAndStroke(Shape shape) {
		Color c= getFillColor();
		if(c!=null ) {
			w.print(" ");
			pathIterator(shape);
			w.print(" ");
			w.print(setrgbcolor(c));
			w.print(" fill");
			}
		
		final double linewidth=getStrokeWidth();
		if(linewidth>0) {
			c= getStrokeColor();
			if(c!=null) {
				pathIterator(shape);
				w.print(" ");
				w.print(round(linewidth));
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
	public Canvas polyline(final List<Point2D> points,FunctionalMap<String, Object> fm) {
		fm = this.states.push(this.states.peek().plus(fm));
		final double linewidth=getStrokeWidth();
		if(linewidth>0 && points.size()>1) {
			final Color c= getStrokeColor();
			if(c!=null) {
				w.print(" ");
				w.print(round(linewidth));
				w.print(" setlinewidth ");
				w.print(setrgbcolor(c));

				for(int i=0;i< points.size();++i) {
					final Point2D pt= points.get(i);
					w.append(" ");
					w.append(coord(pt.getX(),pt.getY()));
					w.append(" ");
					w.append(i==0?"moveto":"lineto");
					}
				w.print(" stroke");
				}
			}

		this.states.pop();
		return this;
		}
	
	private void pathIterator(final Shape shape) {
		float coords[]=new float[6];
		final PathIterator iter = shape.getPathIterator(null);
		w.append(" newpath");
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
	}
	
	@Override
	public Canvas shape(Shape shape, FunctionalMap<String, Object> fm) {
		fm = this.states.push(this.states.peek().plus(fm));
		fillAndStroke(shape);
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
