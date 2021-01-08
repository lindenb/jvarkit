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
package com.github.lindenb.jvarkit.iterator;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.Locatable;

public class SlidingWindowIterator<V extends Locatable> extends AbstractIterator<Map.Entry<? extends Locatable, List<V>>> {

private static class Window<V> implements Locatable
	{
	final SimpleInterval interval;
	final List<V> variants = new ArrayList<>();
	
	
	Window(final SimpleInterval  r) {
		this.interval = r;
		}
	@Override
	public String getContig() {
		return interval.getContig();
		}
	@Override
	public int getStart() {
		return interval.getStart();
		}
	@Override
	public int getEnd() {
		return interval.getEnd();
		}
	}
	
	
	
private final Iterator<V> delegate;
private final int window_size;
private final int window_shift;
private final List<Window<V>> buffer= new ArrayList<>();
private final Map<SimpleInterval,Window<V>> interval2win=new HashMap<>();
private final LinkedList<Window<V>> to_be_released = new LinkedList<>();
private String prevContig=null;

public SlidingWindowIterator(final Iterator<V> delegate,int window_size,int window_shift) {
	this.delegate = delegate;
	this.window_size = window_size;
	this.window_shift = window_shift;
	if(this.delegate==null) throw new IllegalArgumentException("null delegate iterator");
	if(this.window_size <=0 ) throw new IllegalArgumentException("window_size <=0 :" + window_size);
	if(this.window_shift <=0 ) throw new IllegalArgumentException("window_shift <=0 :" + window_shift);
	}

@Override
protected Entry<? extends Locatable, List<V>> advance() {
		for(;;) {
			if(!to_be_released.isEmpty())
				{
				final Window<V> w = to_be_released.pollFirst();
				return new AbstractMap.SimpleEntry<>(w.interval,w.variants);
				}
						
			final V ctx = this.delegate.hasNext()?this.delegate.next():null;
			// new contig?
			if(ctx==null || !ctx.getContig().equals(prevContig)) {
				to_be_released.addAll(this.buffer);
				if(ctx==null && to_be_released.isEmpty()) return null;//EOF
				interval2win.clear();
				buffer.clear();
				}
			
			if(ctx==null) continue;// because to_be_released might be not empty
			this.prevContig = ctx.getContig();
			
			//remove previous windows
			int i=0;
			while(i< buffer.size()) {
				final Window<V> w = buffer.get(i);
				if(w.getEnd() < ctx.getStart()) {
					to_be_released.add(w);
					buffer.remove(i);
					interval2win.remove(w.interval);
					}
				else
					{
					i++;
					}
				}
			
			prevContig=ctx.getContig();
			int x1 = ctx.getStart() -   ctx.getStart()%this.window_shift;
			while(x1-this.window_shift+this.window_size>=ctx.getStart()) {
				x1 -= this.window_shift;
				}
			for(;;) {
				final SimpleInterval r= new SimpleInterval(
						ctx.getContig(),
						Math.max(1,x1),
						Math.max(1,x1+this.window_size)
						);
				
				if(r.getStart()>ctx.getEnd()) break;
				if(r.overlaps(ctx)) {
					Window<V> w = interval2win.get(r);
					if(w==null) {
						w= new Window<V>(r);
						interval2win.put(r, w);
						buffer.add(w);
						}
					w.variants.add(ctx);
					}
				x1+=this.window_shift;
				}
			
			}
	}

	@Override
	public String toString() {
		return this.getClass().getName()+"(size:"+this.window_size+";shift:"+this.window_shift+")";
	}
}
