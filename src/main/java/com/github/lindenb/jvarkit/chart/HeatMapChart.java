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

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;


public class HeatMapChart<Z extends Number> extends Chart {
public static class Data
	extends AbstractData<String,String>
	{
	public Data(final String x,final String y) {
		super(x,y);
		}
	@Override
	public int hashCode() {
		return getX().hashCode()*31 + getY().hashCode();
		}
	public boolean equals(final Object o) {
		if(o==this) return true;
		if(o==null || !(o instanceof Data)) return false;
		final Data d=Data.class.cast(o);
		return	this.getX().equals(d.getX()) &&
				this.getY().equals(d.getY());
		}
	}
private class HeatCategoryAxis extends CategoryAxis
	{
	private boolean dirty_flag=true;
	private final Function<Data,String> extractor;
	HeatCategoryAxis(final Function<Data,String> extractor) {
		this.extractor = extractor;
		}
	HeatCategoryAxis update() {
		if(dirty_flag) {
			dirty_flag=false;
			getCategories().clear();
			getCategories().addAll(HeatMapChart.this.getData().keySet().stream().map(this.extractor).collect(Collectors.toCollection(()->new TreeSet<>())));
			}
		return this;
		}
	}



private final Map<Data,Z> data = new HashMap<>();
private final HeatCategoryAxis xAxis = new HeatCategoryAxis(D->D.getX());
private final HeatCategoryAxis yAxis = new HeatCategoryAxis(D->D.getY());



	
public HeatMapChart() {
	}

public Optional<Z> get(final String x,final String y) {
	return get(new Data(x,y));
	}

public Optional<Z> get(final Data key) {
	return Optional.ofNullable(this.getData().get(key));
	}

public boolean contains(final Data key) {
	return this.data.containsKey(key);
	}
public boolean contains(final String x,final String y) {
	return this.contains(new Data(x,y));
	}

public void put(final Data key,final Z z) {
	this.data.put(key, z);
	}
public void put(final String x,final String y,final Z z) {
	this.xAxis.dirty_flag=true;
	this.yAxis.dirty_flag=true;
	this.put(new Data(x,y),z);
	}

public CategoryAxis getXAxis() {
	return this.xAxis.update();
	}

public CategoryAxis getYAxis() {
	return this.yAxis.update();
	}


public boolean isEmpty() {
	return this.getData().isEmpty();
	}

public Map<Data,Z>  getData() {
	return data;
	}

@Override
public void update() {
	this.xAxis.update();
	this.yAxis.update();
	}


@Override
public String toString() {
	return getClass().getSimpleName();
	}
}
