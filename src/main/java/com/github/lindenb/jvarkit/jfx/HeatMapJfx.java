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
package com.github.lindenb.jvarkit.jfx;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.layout.StackPane;
import javafx.scene.paint.Color;

public class HeatMapJfx extends StackPane {
	
private static class Pair
	{
	final String key1;
	final String key2;
	Pair(final String key1,final String key2)
		{
		if(key1.compareTo(key2)<=0)
			{
			this.key1=key1;
			this.key2=key2;
			}
		else
			{
			this.key1=key2;
			this.key2=key1;
			}
		}
	@Override
	public int hashCode() {
		return this.key1.hashCode()*31+this.key2.hashCode();
		}
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof Pair)) return false;
		final Pair o=Pair.class.cast(obj);
		return this.key1.equals(o.key1) &&
				 this.key2.equals(o.key2);
		}
	}
private final Set<String> keys = new TreeSet<>();
private final Map<Pair,Double> pair2vals = new HashMap<>();
private Double minValue = null;
private Double maxValue = null;
private Color lowColor = Color.BLUE;
private Color highColor = Color.RED;
private final Canvas canvas;
HeatMapJfx()
	{
	this.canvas=new Canvas()
		{
		
		@Override
		public boolean isResizable() {
			return true;
			}
		@Override
		public double prefWidth(double height) {
			return this.getWidth();
			}
		@Override
		public double prefHeight(double width) {
			return this.getHeight();
			}
		};
	this.canvas.widthProperty().bind( this.widthProperty());
	this.canvas.heightProperty().bind( this.heightProperty());
	// Redraw canvas when size changes.
	this.canvas.widthProperty().addListener(evt -> draw());
	this.canvas.heightProperty().addListener(evt -> draw());
	this.getChildren().add(this.canvas);
	}


public void put(final String key1,final String key2,double value)
	{
	this.keys.add(key1);
	this.keys.add(key2);
	this.pair2vals.put(new Pair(key1,key2),value);
	if(this.minValue==null || value < this.minValue.doubleValue())
		{
		this.minValue = value;
		}
	if(this.maxValue==null || value > this.maxValue.doubleValue())
		{
		this.maxValue = value;
		}
	}
private void draw() {
	final double width =  this.canvas.getWidth();
	final double height =  this.canvas.getHeight();
	final GraphicsContext gc = this.canvas.getGraphicsContext2D();
	gc.clearRect(0, 0, width, height);
	if(this.minValue==null || this.maxValue==null || this.keys.isEmpty()) return;
	final int labelWith=100;
	final int scaleWith=200;
	double sqareSize= Math.min(width, height);
	if(sqareSize<=0) return;
	final List<String> keyList = new ArrayList<>(this.keys);
	final double unitSize = sqareSize/keyList.size();
	if(unitSize<=0) return;
	gc.translate(0, 0);
	for(int x=0;x< keyList.size();++x)
		{
		for(int y=0;y< keyList.size();++y)
			{
			final Pair p = new Pair(keyList.get(x), keyList.get(y));
			final Double val = this.pair2vals.get(p);
			if(val==null) continue;
			float t;
			if(this.maxValue.equals(this.minValue)) {
				t=0f;
				}
			else
				{
				t=(float)((val-this.minValue.doubleValue())/(this.maxValue-this.minValue));
				}
			final Color c = this.lowColor.interpolate(this.highColor, t);
			gc.setFill(c);
			gc.fillRect(x*unitSize, y*unitSize, unitSize, unitSize);
			}
		}
	//draw grid
	gc.setStroke(Color.LIGHTGREY);
	for(int z=0;z<= keyList.size();++z)
		{
		gc.strokeLine(0,z*unitSize, sqareSize,z*unitSize);
		gc.strokeLine(z*unitSize,0,z*unitSize, unitSize);
		}
	gc.translate(0, 0);
	
	//LinearGradient lg1 = new LinearGradient(125, 0, 225, 0, false, CycleMethod.NO_CYCLE, this.lowColor,this.highColor);
	}
}
