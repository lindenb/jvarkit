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
package com.github.lindenb.jvarkit.samtools.util;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.function.BiPredicate;

import htsjdk.samtools.util.Locatable;

/** pileup locatable for graphical visualizations */
	public class Pileup<T extends Locatable> implements Iterable<List<T>> {
	private final List<List<T>> rows = new ArrayList<>();
	private String prevContig=null;
	private int prev_start=0;
	private final BiPredicate<T, T> noCollisionTest;
	
	/** BiPredicate returns true if the two item DO NOT collidate */
	public Pileup(final BiPredicate<T, T> noCollisionTest) {
		this.noCollisionTest = noCollisionTest;
		}
	
	public Pileup() {
		this((left,right)->left.getEnd()+1 < right.getStart());
		}
	
	/** add a new item to the pileup , return the item */
	public T add(final T item) {
		if(prevContig==null) {
			this.prevContig = item.getContig();
			this.prev_start = item.getStart();
			}
		else {
			if(!item.getContig().equals(this.prevContig)) {
				throw new IllegalArgumentException("Cannot pileup with different contigs: "+this.prevContig+" "+item.getContig());
				}
			if(item.getStart()< this.prev_start) {
				throw new IllegalArgumentException("Cannot pileup with unordered element: "+this.prev_start+" > "+item.getStart());
				}
			}
		
		int y=0;
		while(y< this.rows.size()) {
			final List<T> row = this.rows.get(y);
			final T last = row.get(row.size()-1);
			if(overlap(last,item)) {
				y++;
				continue;
				}
			row.add(item);
			break;
			}
		if(y==this.rows.size()) {
			final List<T> row = new ArrayList<>();
			row.add(item);
			this.rows.add(row);
			}
		return item;
		}
	
	/** test wether two items can be added to the same row */
	protected boolean overlap(final T left,final T right) {
		if(left.overlaps(right)) return true;
		return !this.noCollisionTest.test(left, right);
		}
	
	public boolean isEmpty() {
		return this.getRows().isEmpty();
		}
	
	public int getRowCount() {
		return this.getRows().size();
		}
	
	public List<List<T>> getRows() {
		return this.rows;
		}
	public List<T> getRow(int y) {
		return this.getRows().get(y);
		}
	
	public void clear() {
		this.prevContig = null;
		this.prev_start = 0;
		this.rows.clear();
		}
	
	@Override
	public Iterator<List<T>> iterator()
		{
		return getRows().iterator();
		}
	
	@Override
	public String toString()
		{
		return "Pileup("+getRowCount()+")";
		}
	}
