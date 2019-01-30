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

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public abstract class XYChart<X, Y> extends Chart {
public static class Data<X,Y>
	extends AbstractData<X,Y>
	{
	public Data(final X x,final Y y) {
		super(x,y);
		}
	public X getXValue() { return this.getX();}
	public Y getYValue() { return this.getY();}
	}
public static class Series<X,Y>
	{
	private String name;
	private final List<Data<X,Y>> data;
	public Series(final String name,final List<Data<X,Y>> data) {
		this.name = name;
		this.data = new ArrayList<>(data);
		}
	public Series() {
		this(Collections.emptyList());
		}
	public Series(final List<Data<X,Y>> data) {
		this("",new ArrayList<>(data));
		}
	public Series(final String name) {
		this(name,Collections.emptyList());
		}
	public void setName(final String name) {
		this.name = name;
		}
	public String getName() {
		return this.name;
		}
	public List<Data<X,Y>> getData() {
		return this.data;
		}
	@Override
	public String toString() {
		return "Series name="+this.name+" N="+this.data.size();
		}
	}
private Axis<X> xAxis;
private Axis<Y> yAxis;
private List<XYChart.Series<X,Y>> data;
private boolean verticalGridLinesVisible​ = false;

protected XYChart(final Axis<X> xAxis,final Axis<Y> yAxis, final List<XYChart.Series<X,Y>> data) {
	this.xAxis = xAxis;
	this.yAxis = yAxis;
	this.data = new ArrayList<>(data);
	}
protected XYChart(final Axis<X> xAxis,final Axis<Y> yAxis) {
	this(xAxis,yAxis,Collections.emptyList());
	}

public Axis<X> getXAxis() {
	return this.xAxis;
	}

public Axis<Y> getYAxis() {
	return this.yAxis;
	}

public List<XYChart.Series<X, Y>> getData() {
	return data;
	}



@Override
public void update() {
	if(getXAxis().isString() && getYAxis().isNumber())
		{
		final Set<String> set =  getData().stream().
				flatMap(S->S.getData().stream()).
				map(D->String.class.cast(D.getXValue())).
				collect(Collectors.toCollection(LinkedHashSet::new));
		final CategoryAxis catx=CategoryAxis.class.cast(getXAxis());
		catx.getCategories().clear();
		catx.getCategories().addAll(set);
		
		
		
		}
	}

public void setVerticalGridLinesVisible(boolean verticalGridLinesVisible​) {
	this.verticalGridLinesVisible​ = verticalGridLinesVisible​;
}
public boolean isVerticalGridLinesVisible​() {
	return verticalGridLinesVisible​;
}

@Override
public String toString() {
	return getClass().getSimpleName()+" xaxis:"+getXAxis()+" yaxis:"+getYAxis();
	}
}
