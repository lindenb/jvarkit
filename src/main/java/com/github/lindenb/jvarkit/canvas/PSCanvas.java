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
package com.github.lindenb.jvarkit.canvas;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
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



/**
 * Implementation of Canvas for POSTSCRIPT
 *
 */
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
		
		// https://opensource.adobe.com/dc-acrobat-sdk-docs/library/pdfmark/pdfmark_Syntax.html
		this.w.println("%%BeginProlog");
		this.w.println("  /pdfmark where");
		this.w.println("% Is pdfmark already available?");
		this.w.println("");
		this.w.println("      { pop }");
		this.w.println("% Yes: do nothing (use that definition)");
		this.w.println("");
		this.w.println("      {");
		this.w.println("% No: define pdfmark as follows:");
		this.w.println("");
		this.w.println("      /globaldict where");
		this.w.println("% globaldict is preferred because");
		this.w.println("");
		this.w.println("          { pop globaldict }");
		this.w.println("% globaldict is always visible; else,");
		this.w.println("");
		this.w.println("          { userdict }");
		this.w.println("% use userdict otherwise.");
		this.w.println("");
		this.w.println("      ifelse");
		this.w.println("      /pdfmark /cleartomark load put");
		this.w.println("      }");
		this.w.println("% Define pdfmark to remove all objects");
		this.w.println("");
		this.w.println("  ifelse");
		this.w.println("% up to and including the mark object.");
		this.w.println("%%EndProlog");

		
		
		this.w.println("/rd { rlineto } def ");
		for(int i=1;i<=10;i++) {
			this.w.println("/slw"+i+" { "+i+" setlinewidth } def");
		}
		
		this.w.println("/rect { 4 dict begin /h exch def /w exch def /y exch def /x exch def "+
			"newpath "+
			"x y moveto " +
			"w 0 rd " +
			"0 h neg rd " +
			"w neg 0 rd " +
			"0 h rd " +
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
		this.w.println("/circle { 3 dict begin /r exch def /y exch def /x exch def "+
				"newpath "+
				"x y r 0 360 arc " +
				"closepath " +
				"end } def"
				);
		
		this.w.println("/cross { 3 dict begin /r exch def /y exch def /r exch def "+
				"newpath "+
				"x r 2.0 div add y moveto " +
				"x r 2.0 div sub y lineto " +
				"x y r 2.0 div add moveto " +
				"x y r 2.0 div sub lineto " +
				"closepath " +
				"end } def"
				);
		
		this.w.println("/hl { 1 dict begin  /dx exch def "+
				"dx 0 rlineto " +
				"end } def"
				);
		this.w.println("/vl { 1 dict begin  /dy exch def "+
				"0 dy neg rd " +
				"end } def"
				);
		this.w.println("/irgb { 3 dict begin /b exch def /g exch def /r exch def "+
				"r 255.0 div " +
				"g 255.0 div " +
				"b 255.0 div " +
				" setrgbcolor " +
				"end } def"
				);
		this.w.println("/igray { 1 dict begin /g exch def  "+
				" g g g irgb " +
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
	
	private String round(final double x) {
		String s= String.format("%.3f",x);
		while(s.contains(".") && s.endsWith("0")) s=s.substring(0,s.length()-1);
		if(s.endsWith("."))  s=s.substring(0,s.length()-1);
		return s;
	}
	
	private String coord(final double x,final double y) {
		return round(inch(x))+" "+round(inch(getHeight()-y));
		}
	
	private String setrgbcolor(final Color c) {
		if(c==null) return "";
		if(c.getRed()==c.getGreen() && c.getGreen()==c.getBlue()) {
			return String.valueOf(c.getRed())+" igray";
			}
		
		return String.join(" ",
				String.valueOf(c.getRed()),
				String.valueOf(c.getGreen()),
				String.valueOf(c.getBlue()),
				"irgb"
				);
		}
	
	private String setLineWidth(double linewidth) {
		for(int i=1;i<=10;i++) {
			if(i==linewidth) return "slw"+(int)linewidth;
			}
		return round(inch(linewidth))+" setlinewidth";
		}
	
	private Canvas fillAndStroke(final Shape shape) {
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
				w.print(setLineWidth(linewidth));
				w.print(" ");
				w.print(setrgbcolor(c));
				w.print(" stroke");
				}
			}
		if(gsave_set) {
			w.print(" grestore");
			}
		return this;
		}
	
	private static String escape(final CharSequence s) {
		return StringUtils.escapePostscript(s);
		}
	@Override
	public Canvas text(final String text, double x, double y, FunctionalMap<String, Object> fm) {
		if(StringUtils.isBlank(text)) return this;
		fm = this.states.push(this.states.peek().plus(fm));
		final Color c = getFillColor();
		if(c!=null) {
			final int fontSize = Integer.parseInt(String.valueOf(fm.getOrDefault(KEY_FONT_SIZE, "12")));
			w.print(" gsave ");
			w.print(setrgbcolor(c));
			w.print(" /"+ fm.getOrDefault(KEY_FONT_FAMILY, "Times-Roman")+" findfont");
			w.print(" "+ fontSize +" scalefont");
			w.print(" setfont newpath ");
			w.print(coord(x,y)+" moveto");
			w.print(" ("+escape(text)+") show");
			w.print(" grestore ");
			
			hyperlink(new Rectangle2D.Double(x,y-fontSize,text.length()*fontSize,fontSize));
			}
		this.states.pop();
		return this;
		}
	
	@Override
	public Canvas comment(final String s) {
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
				w.print(setLineWidth(linewidth));
				w.print("  ");
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
		else if(shape instanceof Line2D) {
			final Line2D line = Line2D.class.cast(shape);
			w.append(coord(line.getX1(),line.getY1()));
			w.append(" ");
			w.append(coord(line.getX2(),line.getY2()));
			w.append(" ll ");
			return;
			}
		else if(shape instanceof Ellipse2D) {
			final Ellipse2D ellipse = Ellipse2D.class.cast(shape);
			if(ellipse.getWidth()==ellipse.getHeight()) {
				w.append(coord(ellipse.getCenterX(),ellipse.getCenterY()));
				w.append(" ");
				w.append(round(inch(ellipse.getHeight()/2.0)));
				w.append(" circle ");
				return;
				}
			}
		else if(shape instanceof Canvas.CrossShape) {
			final Ellipse2D ellipse = Ellipse2D.class.cast(shape);
			if(ellipse.getWidth()==ellipse.getHeight()) {
				w.append(coord(ellipse.getCenterX(),ellipse.getCenterY()));
				w.append(" ");
				w.append(round(inch(ellipse.getHeight()/2.0)));
				w.append(" cross ");
				return;
				}
			}
		
		float prev_x=0;
		float prev_y=0;
		final float coords[]=new float[6];
		PathIterator iter = shape.getPathIterator(null);
		boolean has_close_path=false;
		while(!iter.isDone()) {
			if(iter.currentSegment(coords)== PathIterator.SEG_CLOSE)
				{
				has_close_path=true;
				break;
				}
			iter.next();
			}
		iter = shape.getPathIterator(null);
		if(has_close_path) w.append(" newpath");
		while(!iter.isDone()) {
			switch(iter.currentSegment(coords))
			{
			case PathIterator.SEG_MOVETO:
				{
				w.append(" ");
				w.append(coord(coords[0],coords[1]));
				w.append(" moveto");
				prev_x = coords[0];
				prev_y = coords[1];
				break;
				}
			case PathIterator.SEG_LINETO:
				{
				if(prev_y==coords[1]) {
					w.append(" ");
					w.append(round(inch(coords[0]-prev_x)));
					w.append(" hl");//horizontal line
					}
				else if(prev_x==coords[0]) {
					w.append(" ");
					w.append(round(inch(coords[1]-prev_y)));
					w.append(" vl");//vertical line
					}
				else
					{
					w.append(" ");
					w.append(round(inch(coords[0]-prev_x)));
					w.append(" ");
					w.append(round(inch(-(coords[1]-prev_y))));
					w.append(" rd");//rlineto
					}
				prev_x = coords[0];
				prev_y = coords[1];
				break;
				}
			case PathIterator.SEG_QUADTO:
				{
				throw new IllegalArgumentException("Quad curves are not supported in postscript canvas");
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
				prev_x = coords[4];
				prev_y = coords[5];
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
	public Canvas shape(final Shape shape, FunctionalMap<String, Object> fm) {
		fm = this.states.push(this.states.peek().plus(fm));
		fillAndStroke(shape);
		hyperlink(shape);
		this.states.pop();
		return this;
		}
	
	private void hyperlink(final Shape shape) {
		final Object o= this.states.peek().getOrDefault(KEY_HREF,"");
		final String href = o==null?"":o.toString();
		if(!StringUtils.isBlank(href)) {
			final Rectangle2D r=shape.getBounds2D();
			if(r.getWidth()>=1 && r.getHeight()>=1) {
				final Point2D.Double p1=new  Point2D.Double(r.getX(),r.getY());
				final Point2D.Double p2=new  Point2D.Double(r.getMaxX(),r.getMaxY());
				this.w.println(" [ /Rect ["+coord(p1.getX(),p1.getY()) +" "+coord(p2.getX(),p2.getY())+"]");
				this.w.println("/Action << /Subtype /URI /URI ("+ escape(href) +") >>");
				this.w.println("/Border [ 0 0 0 ]");
				this.w.println("/Subtype /Link");
				this.w.println("/ANN pdfmark ");
				}
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
