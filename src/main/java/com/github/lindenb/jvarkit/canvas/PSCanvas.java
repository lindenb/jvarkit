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
import java.awt.geom.Line2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
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
	
	private static double inch(double pixel) {
		return pixel;
		}
	
	public PSCanvas(
		final Path out,
		final int width,
		final int height,
		final boolean compressed,
		final FunctionalMap<String,Object> params) throws IOException {
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
		this.w.println("%%Creator: jvarkit");
		this.w.println("%%BoundingBox: 0 0 "+(int)inch(width)+" "+(int)inch(height));
		final String title = params.getOrDefault(KEY_TITLE, "").toString();
		if(!StringUtils.isBlank(title)) {
			this.w.println("%%Title: "+title);
			}
		this.w.println("/rect { 4 dict begin /h exch def /w exch def /y exch def /x exch def "+
			"newpath "+
			"x y moveto " +
			"w 0 rlineto " +
			"0 h neg rlineto " +
			"w neg 0 rlineto " +
			"0 h rlineto " +
			"closepath " +
			"end } def"
			);
		this.w.println("/square { 3 dict begin /w exch def /y exch def /x exch def "+
			"x y w w rect " +
			"end } def"
			);
		this.w.println("/ll { 4 dict begin /y2 exch def /x2 exch def /y1 exch def /x1 exch def "+
			"newpath "+
			"x1 y1 moveto " +
			"x2 y2 lineto " +
			"closepath " +
			"end } def"
			);
		}
	

	
	@Override
	public void close() throws IOException {
		this.w.println("\nshowpage");
		this.w.println("% generated with jvarkit. Author: Pierre Lindenbaum PhD.");
		this.w.println("%EOF");
		this.w.flush();
		this.w.close();
		}
	
	private String round(double x) {
		String s= String.format("%.3f",x);
		while(s.contains(".") && s.endsWith("0")) s=s.substring(0,s.length()-1);
		if(s.endsWith("."))  s=s.substring(0,s.length()-1);
		return s;
	}
	
	private String coord(double x,double y) {
		return round(inch(x))+" "+round(inch(getHeight()-y));
		}
	
	private String setrgbcolor(final Color c) {
		if(c==null) return "";
		if(c.getRed()==c.getGreen() && c.getGreen()==c.getBlue()) {
			return round(c.getRed()/255f)+" setgray";
			}
		
		return String.join(" ",
				round(c.getRed()/255f),
				round(c.getGreen()/255f),
				round(c.getBlue()/255f),
				"setrgbcolor"
				);
		}
	
	private Canvas fillAndStroke(Shape shape) {
		boolean gsave_set=false;
		Color c= getFillColor();
		if(c!=null ) {
			w.print(" gsave ");
			pathIterator(shape);
			w.print(" ");
			w.print(setrgbcolor(c));
			w.print(" fill");
			gsave_set = true;
			}
		
		final double linewidth=getStrokeWidth();
		if(linewidth>0) {
			c= getStrokeColor();
			if(c!=null) {
				if(!gsave_set) {
					w.print(" gsave ");
					pathIterator(shape);
					gsave_set = true;
					}
				w.print(" ");
				w.print(round(inch(linewidth)));
				w.print(" setlinewidth ");
				w.print(setrgbcolor(c));
				w.print(" stroke");
				}
			}
		if(gsave_set) {
			w.print(" grestore");
			}
		return this;
		}
	
	private String escape(final CharSequence s) {
		return StringUtils.escapePostscript(s);
		}
	@Override
	public Canvas text(String text, double x, double y, FunctionalMap<String, Object> fm) {
		if(StringUtils.isBlank(text)) return this;
		fm = this.states.push(this.states.peek().plus(fm));
		final Color c = getFillColor();
		if(c!=null) {
			w.print(" gsave ");
			w.print(setrgbcolor(c));
			w.print(" /"+ fm.getOrDefault(KEY_FONT_FAMILY, "Times-Roman")+" findfont");
			w.print(" "+ fm.getOrDefault(KEY_FONT_SIZE, "12")+" scalefont");
			w.print(" setfont newpath ");
			w.print(coord(x,y)+" moveto");
			w.print(" ("+escape(text)+") show");
			w.print(" grestore ");
			}
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
				w.print(" gsave ");
				w.print(round(inch(linewidth)));
				w.print(" setlinewidth ");
				w.print(setrgbcolor(c));

				for(int i=0;i< points.size();++i) {
					final Point2D pt= points.get(i);
					w.append(" ");
					w.append(coord(pt.getX(),pt.getY()));
					w.append(" ");
					w.append(i==0?"moveto":"lineto");
					}
				w.print(" stroke grestore ");
				}
			}

		this.states.pop();
		return this;
		}
	
	private void pathIterator(final Shape shape) {
		if(shape instanceof Rectangle2D) {
			final Rectangle2D rect = Rectangle2D.class.cast(shape);
			w.append(coord(rect.getX(),rect.getY()));
			w.append(" ");
			if(rect.getWidth()==rect.getHeight()) {
				w.append(round(inch(rect.getWidth())));
				w.append(" square ");
				}
			else
				{
				w.append(round(inch(rect.getWidth())));
				w.append(" ");
				w.append(round(inch(rect.getHeight())));
				w.append(" rect ");
				}
			return;
			}
		if(shape instanceof Line2D) {
			final Line2D line = Line2D.class.cast(shape);
			w.append(coord(line.getX1(),line.getY1()));
			w.append(" ");
			w.append(coord(line.getX2(),line.getY2()));
			w.append(" ll ");
			return;
			}
		
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
	public Canvas line(double x1, double y1, double x2, double y2, FunctionalMap<String, Object> fm) {
		return super.line(x1, y1, x2, y2, fm);
		}
	
	@Override
	public Canvas shape(final Shape shape, FunctionalMap<String, Object> fm) {
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
