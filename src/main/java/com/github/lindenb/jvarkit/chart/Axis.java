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
package com.github.lindenb.jvarkit.chart;

public abstract class Axis<T> {
	public static enum Type {STRING,NUMBER};
	private String label="";
	private double tickLabelRotation = 0.0;
	
	protected Axis(final String label) {
		this.label = label;
	}
	
	public String getLabel() {
		return this.label;
		}

	public void setLabel(final String label) {
		this.label = label;
		}

	@Override
	public String toString() {
		return getClass().getSimpleName()+":"+getLabel();
		}
	
	public void setTickLabelRotation(double tickLabelRotation) {
		this.tickLabelRotation = tickLabelRotation;
	}
	
	public double getTickLabelRotation() {
		return tickLabelRotation;
	}
	public abstract Type getType(); 
	public final boolean isNumber() { return getType().equals(Type.NUMBER);}
	public final boolean isString() { return getType().equals(Type.STRING);}
	}
