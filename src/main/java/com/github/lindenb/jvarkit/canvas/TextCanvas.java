package com.github.lindenb.jvarkit.canvas;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.PathIterator;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import com.github.lindenb.jvarkit.util.FunctionalMap;

public class TextCanvas extends Canvas {
	private final int width;
	private final int height;
	private final PrintWriter w;
	private final char[][] drawingarea;

	public TextCanvas(final Path out,
			final int width,
			final int height,
			final FunctionalMap<String,Object> params) throws IOException {
		super();
		this.width=width;
		this.height=height;	
		this.drawingarea = new char[height][];
		for(int y=0;y< height;++y) {
			drawingarea[y]=new char[width];
			Arrays.fill(drawingarea[y], ' ');
			}
		if(out==null) {
			w = new PrintWriter(System.out);
			}
		else
			{
			w =new PrintWriter(Files.newOutputStream(out));
			}
		}
	
	
	
	private void set(int x,int y,char c) {
		if(x<0 || x>=this.width) return;
		if(y<0 || y>=this.height) return;
		drawingarea[y][x]=c;
		}
	
	private void bresenham( int x1, int y1, int x2, int y2, char color) {
	        int deltax = Math.abs(x2 - x1);
	        int deltay = Math.abs(y2 - y1);
	        int error = 0;
	        int y = y1;
	        for( int x=x1; x<x2; x++) {
	            set(x, y, color);
	            error = error + deltay;
	            if( 2*error >= deltax ) {
	                y = y + 1;
	                error=error - deltax;
	            	}
	        	}
	    	}



	@Override
	public void close() throws IOException {
		for(char[] s: this.drawingarea) {
		w.println(new String(s));
		}
		w.flush();
		w.close();
	}



	@Override
	public int getWidth() {
		return this.width;
	}



	@Override
	public int getHeight() {
		return this.height;
	}

	private char colorToChar(Color c)  {
		if(c.equals(Color.WHITE)) return ' ';
		return '#';
	}
	
	@Override
	public Canvas line(double x1, double y1, double x2, double y2, FunctionalMap<String, Object> fm) {
		fm = this.states.push(this.states.peek().plus(fm));
		Color col = toColor(fm.get(KEY_STROKE));
		if(col!=null) {
			char c = colorToChar(col);
			bresenham((int)x1, (int)y1, (int)x2, (int)y2,c);
			}
		this.states.pop();
		return this;
		}
	
	@Override
	public Canvas text(String s, double x, double y, FunctionalMap<String, Object> fm) {
		fm = this.states.push(this.states.peek().plus(fm));
		Color col = toColor(fm.get(KEY_FILL));
		if(col!=null) {
			for(int i=0;i< s.length();i++) {
				set((int)(x+i),(int)y,s.charAt(i));
				}
			}
		this.states.pop();
		return this;
		}



	@Override
	public Canvas shape(Shape shape, FunctionalMap<String, Object> fm) {
		fm = this.states.push(this.states.peek().plus(fm));
		Rectangle r1 = shape.getBounds();
		// fill
		Color col= toColor(fm.get(KEY_FILL));
		if(col!=null) {
			char c = colorToChar(col);
			for(int x=(int)r1.getMinX();x<=(int)r1.getMaxX();++x) {
				for(int y=(int)r1.getMinY();y<=(int)r1.getMaxY();++y) {
					if(!shape.contains(x, y)) continue;
					set(x,y,c);
					}
				}
			}
		// fill
		col= toColor(fm.get(KEY_STROKE));
		if(col!=null) {
			char c = colorToChar(col);
			float prev_x=0;
			float prev_y=0;
			final float coords[]=new float[6];
			PathIterator iter=shape.getPathIterator(null);
			while(!iter.isDone()) {
				switch(iter.currentSegment(coords))
				{
				case PathIterator.SEG_MOVETO:
					{
					prev_x = coords[0];
					prev_y = coords[1];
					break;
					}
				case PathIterator.SEG_LINETO:
					{
					bresenham((int)prev_x, (int)prev_y, (int)coords[0], (int)coords[1],c );
					prev_x = coords[0];
					prev_y = coords[1];
					break;
					}
				case PathIterator.SEG_QUADTO:
					{
					prev_x = coords[2];
					prev_y = coords[3];
					break;
					}
				case PathIterator.SEG_CUBICTO:
					{
					prev_x = coords[4];
					prev_y = coords[5];
					break;
					}
				case PathIterator.SEG_CLOSE:
					{
					break;
					}
				}
			iter.next();
			}
			}
		this.states.pop();
		return this;
	}
}
