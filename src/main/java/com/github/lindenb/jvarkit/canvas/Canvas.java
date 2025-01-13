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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.Closeable;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Stack;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

/**
 * Abstract drawing class for drawing SVG, PNG, PS, etc..
 * @author lindenb
 *
 */
public abstract class Canvas implements Closeable {
	public static final String KEY_FILL = "fill";
	public static final String KEY_STROKE = "stroke";
	public static final String KEY_STROKE_WIDTH = "stroke-width";
	public static final String KEY_TITLE="title";
	public static final String KEY_FONT_FAMILY="font-family";
	public static final String KEY_FONT_SIZE="font-size";
	public static final String KEY_STROKE_LINE_JOIN="stroke-linejoin";
	public static final String KEY_STROKE_LINE_CAP="stroke-linecap";
	public static final String KEY_HREF="href";
	public static final String KEY_FILL_OPACITY="fill-opacity";
	public static final String KEY_STROKE_OPACITY="stroke-opacity";
	public static final String KEY_OPACITY="opacity";
	public static final String KEY_SVG_META = "-svg-meta";// instance of Consumer<XMLStreamWriter> 

	protected abstract class AbstractShape implements Shape {
		protected final GeneralPath gp = new GeneralPath();
		AbstractShape() {
			}
		@Override
		public Rectangle getBounds() {
			return gp.getBounds();
		}

		@Override
		public Rectangle2D getBounds2D() {
			return gp.getBounds2D();
		}

		@Override
		public boolean contains(double x, double y) {
			return gp.contains(x,y);
		}

		@Override
		public boolean contains(Point2D p) {
			return gp.contains(p);
		}

		@Override
		public boolean intersects(double x, double y, double w, double h) {
			return gp.intersects(x,y,w,h);
		}

		@Override
		public boolean intersects(Rectangle2D r) {
			return gp.intersects(r);		
			}

		@Override
		public boolean contains(double x, double y, double w, double h) {
			return gp.contains(x,y,w,h);
		}

		@Override
		public boolean contains(Rectangle2D r) {
			return gp.contains(r);
		}

		@Override
		public PathIterator getPathIterator(AffineTransform at) {
			return gp.getPathIterator(at);
		}

		@Override
		public PathIterator getPathIterator(AffineTransform at, double flatness) {
			return gp.getPathIterator(at,flatness);
		}
		
	}

	
	protected class CrossShape extends AbstractShape {
		CrossShape(double x,double y, double r) {
			gp.moveTo(x-r, y);
			gp.lineTo(x+r, y);
			gp.moveTo(x, y-r);
			gp.lineTo(x, y+r);
			}
	}
	
	
	protected Stack<FunctionalMap<String, Object>> states = new Stack<>();
	protected Hershey hershey = new Hershey();
	protected final ColorUtils colorUtils = new ColorUtils();
	protected Canvas() {
		states.add(new FunctionalMap<String, Object>().
				plus(KEY_FILL,null).
				plus(KEY_STROKE,Color.darkGray).
				plus(KEY_STROKE_WIDTH,1).
				plus(KEY_FILL_OPACITY,1.0).
				plus(KEY_STROKE_OPACITY,1.0).
				plus(KEY_FONT_SIZE,10)
				);
		}
	protected double getOpacity() {
		Object o =this.states.peek().getOrDefault(KEY_OPACITY,1.0);
		if(o!=null && o instanceof Number) {
			return Math.min(1.0, Math.max(0.0, Number.class.cast(o).doubleValue()));
			}
		throw new IllegalArgumentException("cannot convert "+o+" to double");
		}
	protected double getFillOpacity() {
		Object o =this.states.peek().getOrDefault(KEY_FILL_OPACITY,getOpacity());
		if(o!=null && o instanceof Number) {
			return Math.min(1.0, Math.max(0.0, Number.class.cast(o).doubleValue()));
			}
		throw new IllegalArgumentException("cannot convert "+o+" to double");
		}
	protected double getStrokeOpacity() {
		Object o =this.states.peek().getOrDefault(KEY_STROKE_OPACITY,getOpacity());
		if(o!=null && o instanceof Number) {
			return Math.min(1.0, Math.max(0.0, Number.class.cast(o).doubleValue()));
			}
		throw new IllegalArgumentException("cannot convert "+o+" to double");
		}

	protected int getLineCap() {
		Object o = this.states.peek().getOrDefault(KEY_STROKE_LINE_CAP,BasicStroke.CAP_ROUND);
		if(o!=null) {
			if(o instanceof Integer) return Integer.class.cast(o);
			if(o instanceof String) {
				String s= String.valueOf(o).toLowerCase();
				if(s.equals("butt")) return BasicStroke.CAP_BUTT;
				if(s.equals("round")) return BasicStroke.CAP_ROUND;
				if(s.equals("square")) return BasicStroke.CAP_SQUARE;
				}
			}
		return BasicStroke.CAP_ROUND;
		}
	
	protected int getLineJoin() {
		Object o =this.states.peek().getOrDefault(KEY_STROKE_LINE_JOIN,BasicStroke.JOIN_MITER);
		if(o!=null) {
			if(o instanceof Integer) return Integer.class.cast(o);
			if(o instanceof String) {
				String s= String.valueOf(o).toLowerCase();
				if(s.equals("bevel")) return BasicStroke.JOIN_BEVEL;
				if(s.equals("round")) return BasicStroke.JOIN_ROUND;
				if(s.equals("miter")) return BasicStroke.JOIN_MITER;
				}
			}
		return BasicStroke.JOIN_MITER;		
		}
	
	protected OptionalDouble toDouble(Object o) {
		if(o==null) return OptionalDouble.empty();
		if(o instanceof Number ) return OptionalDouble.of(Number.class.cast(o).doubleValue());
		if(o instanceof String) {
			if(o.equals("none")) return null;
			return OptionalDouble.of(Double.parseDouble(String.class.cast(o)));
			}
		return OptionalDouble.empty();
		}
	
	protected Color toColor(Object o) {
		if(o==null) return null;
		if(o instanceof Color ) return Color.class.cast(o);
		if(o instanceof String) {
			return colorUtils.parse(String.class.cast(o));
			}
		throw new IllegalArgumentException("cannot convert "+o+" to color");
		}
	protected Color getFillColor() {
		return toColor(this.states.peek().getOrDefault(KEY_FILL, null));
		}
	protected Color getStrokeColor() {
		return toColor(this.states.peek().getOrDefault(KEY_STROKE, null));
		}
	protected double getStrokeWidth() {
		return toDouble(this.states.peek().getOrDefault(KEY_STROKE_WIDTH, null)).orElse(0.0);
		}
		
	public abstract int getWidth();
	public abstract int getHeight();
	public Canvas begin(final FunctionalMap<String, Object> fm) {
		this.states.push(this.states.peek().plus(fm));
		return this;
		}
	public Canvas end() {
		this.states.pop();
		return this;
		}
	
	public Canvas comment(String s) {
		return this;
		}
	
	public abstract Canvas text(String s,double x,double y,FunctionalMap<String, Object> fm);
	
	public Canvas text(String s,double x,double y) {
		return text(s,x,y,FunctionalMap.make());
		}
	
	public Canvas polyline(List<Point2D> points,FunctionalMap<String, Object> fm) {
		return polyX(points,false,fm);
		}
	public Canvas polyline(List<Point2D> points) {
		return polyline(points,FunctionalMap.make());
		}
	
	public Canvas polygon(final List<Point2D> points,final FunctionalMap<String, Object> fm) {
		return polyX(points,true,fm);
		}
	
	public Canvas polygon(List<Point2D> points) {
		return polygon(points,FunctionalMap.make());
		}
	
	
	
	protected Canvas polyX(List<Point2D> points,boolean closed,FunctionalMap<String, Object> fm) {
		if(points.isEmpty()) return this;
		final GeneralPath gp=new GeneralPath();
		for(int i=0;i< points.size();i++) {
			final Point2D p = points.get(i);
			if(i==0) {
				gp.moveTo(p.getX(),p.getY());
				}
			else
				{
				gp.lineTo(p.getX(),p.getY());
				}
			}
		if(closed) gp.closePath();
		return shape(gp,closed?fm:fm.plus(KEY_FILL,"none"));
		}
	
	public Canvas ellipse(double cx, double cy, double rx, double ry,FunctionalMap<String, Object> fm) {
		return shape(new Ellipse2D.Double(cx-rx, cy-ry, rx*2, ry*2),fm);
		}
	public final Canvas ellipse(double cx, double cy, double rx,double ry) {
		return ellipse(cx,cy,rx,ry,FunctionalMap.make());
		}
	
	public Canvas circle(double cx, double cy, double r,FunctionalMap<String, Object> fm) {
		return ellipse(cx,cy,r,r,fm);
		}
	public  final Canvas circle(double cx, double cy, double r) {
		return circle(cx,cy,r,FunctionalMap.make());
		}
	public Canvas rect(double x, double y, double width,double height,FunctionalMap<String, Object> fm) {
		return shape(new Rectangle2D.Double(x, y, width, height),fm);
		}
	public  final Canvas rect(double x, double y, double width,double height) {
		return rect(x,y,width,height,FunctionalMap.make());
		}
	
	public Canvas line(double x1, double y1, double x2,double y2,FunctionalMap<String, Object> fm) {
		return shape(new Line2D.Double(x1, y1, x2, y2),fm.plus(KEY_FILL,null));	
		}
	
	public  final Canvas line(double x1, double y1, double x2,double y2) {
		return line(x1,y1,x2,y2,FunctionalMap.make());
		}
	
	public Canvas hersey(String s,double x, double y, double width,double height,FunctionalMap<String, Object> fm) {
		if(StringUtils.isBlank(s)) return this;
		return shape(hershey.toShape(s, x, y, width, height),fm);
		}
	
	public  final Canvas hersey(String s,double x, double y, double width,double height) {
		return hersey(s,x,y,width,height,FunctionalMap.make());
		}
	
	public abstract Canvas shape(Shape shape,FunctionalMap<String, Object> fm);
	
	public final Canvas shape(Shape shape) {
		return shape(shape,FunctionalMap.make());
		}
	
	
	public Canvas cross(double x, double y, double r,FunctionalMap<String, Object> fm) {
		return shape(new CrossShape(x, y, r),fm.plus(KEY_FILL,null));	
		}
	}
