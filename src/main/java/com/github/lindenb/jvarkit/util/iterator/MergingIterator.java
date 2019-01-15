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

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.PeekableIterator;

public class MergingIterator<T> 
	extends AbstractIterator<T>
	implements CloseableIterator<T>
	{
	private final List<PeekableIterator<T>> buffer;
	private final Comparator<T> comparator;
	private T lastForChecking=null;
	
	public MergingIterator( final Comparator<T> comparator,final List<? extends Iterator<T>> delegates)
		{
		this.comparator = Objects.requireNonNull(comparator, "comparator is null");
		this.buffer = new ArrayList<>(delegates.stream().map(I->new PeekableIterator<>(I)).collect(Collectors.toList()));
		}	
	
	@Override
	protected T advance() {
		T smallest= null;
		int smallest_index=-1;
		int i=0;
		while(i< this.buffer.size())
			{
			final PeekableIterator<T> delegate = this.buffer.get(i);
			if(!delegate.hasNext())
				{
				CloserUtil.close(delegate);
				this.buffer.remove(i);
				}
			else
				{
				final T item = delegate.peek();
				if(smallest==null || this.comparator.compare(item, smallest)<0)
					{
					smallest = item;
					smallest_index = i;
					}
				i++;
				}
			}
		if(smallest_index!=-1)
			{
			this.buffer.get(smallest_index).next();//consumme
			if(this.lastForChecking!=null &&  this.comparator.compare(smallest, lastForChecking)<0)
				{
				throw new IllegalStateException("Data are not ordered... got "+ 
						smallest+" after "+lastForChecking +" comparator(curr,previous) returns: "+
						this.comparator.compare(smallest, lastForChecking)
						);
				}
			lastForChecking = smallest;	
			return smallest;
			}
		return null;
		}
	
	@Override
	public void close() {
		CloserUtil.close(this.buffer);
		}
	}
