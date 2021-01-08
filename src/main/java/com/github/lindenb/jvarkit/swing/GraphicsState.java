/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.swing;

import java.awt.Composite;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;

/**
 * save `Graphics2D` as a java resource
 *
 */
public class GraphicsState implements AutoCloseable
	{
	private final Graphics2D g;
	private final Font font;
	private final Composite composite;
	private final Stroke stroke;
	private final Paint paint;
	private final Shape shape;
	private final AffineTransform transform;
	private GraphicsState(final Graphics2D g) {
		this.g = g;
		this.font = g.getFont();
		this.composite = g.getComposite();
		this.stroke = g.getStroke();
		this.paint = g.getPaint();
		this.shape = g.getClip();
		this.transform = g.getTransform();
	}
	
	@Override
	public void close()
		{
		g.setFont(font);
		g.setPaint(paint);
		g.setComposite(composite);
		g.setStroke(stroke);
		g.setClip(shape);
		g.setTransform(transform);
		}
	
	public Composite getComposite() {
		return composite;
		}
	
	public static GraphicsState of(final Graphics2D g0) {
		return new GraphicsState(g0);
		}
	
	public Graphics2D getGraphics()  {
		return this.g;
		}
	}
