
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
package com.github.lindenb.jvarkit.util.iterator;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;


public class EqualRangeIterator<T>
	extends AbstractIterator<List<T>>
	implements CloseableIterator<List<T>>
	{	
	private final Iterator<T> delegate;
	private final Comparator<T> comparator;
	private T lastPeek = null;
	private boolean closed = false;
	public EqualRangeIterator(
			final Iterator<T> delegate,
			final Comparator<T> comparator)
		{
		this.delegate = delegate;
		this.comparator = comparator;
		}
	
	/** advance until the current item in the iterator is equal to target.
	 *  result is never null.
	 */
	public List<T> next(final T target)
		{
		if(this.closed) return Collections.emptyList();
		if( target == null) throw new NullPointerException("target is null");
		final List<T> buffer=new ArrayList<>();
		for(;;)
			{
			if(this.lastPeek!=null) {
				final int diff = this.comparator.compare(this.lastPeek, target);
				if(diff > 0) 
					{
					return buffer;
					}
				else if(diff==0)
					{
					buffer.add(this.lastPeek);
					this.lastPeek = null;
					}
				else
					{
					this.lastPeek = null;
					}
				}
			if(!this.delegate.hasNext()) {
				close();
				return buffer;
				}
			final T item = this.delegate.next();
			this.lastPeek =  item;
			}
		}
	
	@Override
	protected List<T> advance() {
		if(this.closed) return null;
		final List<T> buffer=new ArrayList<>();
		
		if(this.lastPeek!=null) {
			buffer.add(this.lastPeek);
			this.lastPeek = null;
		}
		while(this.delegate.hasNext()) {
			final T curr = this.delegate.next();
			if(buffer.isEmpty()) {
				buffer.add(curr);
				}
			else
				{
				final int d= this.comparator.compare(buffer.get(0), curr);
				if(d==0 ) {
					buffer.add(curr);
					}
				else if(d < 0) {
					this.lastPeek = curr;
					break;
					}
				else {
					throw new IllegalStateException("got "+curr +" after "+ buffer.get(0));
					}
				}
			}
		return buffer.isEmpty()?null:buffer;
		}
	
	@Override
	public void close() {
		CloserUtil.close(this.delegate);
		this.closed = true;
		}
	}
