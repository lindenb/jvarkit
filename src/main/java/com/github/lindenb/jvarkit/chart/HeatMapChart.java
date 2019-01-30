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
import java.util.Set;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.chart.PieChart.Data;

public class HeatMapChart<X, Y,Z> extends Chart {
public static class Key<X,Y>
	extends AbstractData<X,Y>
	{
	public Key(final X x,final Y y) {
		super(x,y);
		}
	public X getXValue() { return this.getX();}
	public Y getYValue() { return this.getY();}
	@Override
	public int hashCode() {
		return getX().hashCode()*31 + getY().hashCode();
		}
	public boolean equals(final Object o) {
		if(o==this) return true;
		if(o==null || !(o instanceof Data)) return false;
		final Key<?,?> d=Key.class.cast(o);
		return	this.getX().equals(d.getX()) &&
				this.getY().equals(d.getY());
		}
	}

private final Map<Key<X,Y>,Z> data = new HashMap<>();

public HeatMapChart() {
	}

public Optional<Z> get(final X x,final Y y) {
	return get(new Key<>(x,y));
	}


public Optional<Z> get(final Key<X,Y> key) {
	return Optional.ofNullable(this.getData().get(key));
	}

public boolean contains(final Key<X,Y> key) {
	return this.data.containsKey(key);
	}
public boolean contains(final X x,final Y y) {
	return this.contains(new Key<>(x,y));
	}

public void put(final Key<X,Y> key,final Z z) {
	this.data.put(key, z);
	}
public void put(final X x,final Y y,final Z z) {
	this.put(new Key<>(x,y),z);
	}

public Set<X> getXValues() {
	return getData().keySet().stream().map(K->K.getX()).collect(Collectors.toSet());
	}
public Set<Y> getYValues() {
	return getData().keySet().stream().map(K->K.getY()).collect(Collectors.toSet());
	}

public boolean isEmpty() {
	return this.getData().isEmpty();
	}

public Map<Key<X,Y>,Z>  getData() {
	return data;
	}

@Override
public void update() {
	}


@Override
public String toString() {
	return getClass().getSimpleName();
	}
}
