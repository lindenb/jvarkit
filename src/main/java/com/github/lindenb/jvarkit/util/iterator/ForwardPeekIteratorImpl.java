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
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;

public class ForwardPeekIteratorImpl<E> 
	extends  AbstractIterator<E>
	implements ForwardPeekIterator<E>
	{
	private final Iterator<E> delegate;
	private final List<E> stack = new ArrayList<>();
	public ForwardPeekIteratorImpl(final Iterator<E> delegate) {
		this.delegate = delegate;
		}
	
	private void fillTo(final int index) {
		while(this.stack.size()<=index && this.delegate.hasNext()) {			
			final E N= this.delegate.next();
			if(N==null) throw new NullPointerException("delegate returned null"); 
			this.stack.add(N);
			}
		}
	@Override
	public boolean hasNext() {
		if(!this.stack.isEmpty()) return true;
		return this.delegate.hasNext();
		}
	
	@Override
	public E next() {
		if(!this.stack.isEmpty())
			{
			return this.stack.remove(0);
			}
	
		if(!this.delegate.hasNext()) throw new NoSuchElementException();
		return this.delegate.next();
		}
	@Override
	protected E advance() {
		if(!this.stack.isEmpty())
			{
			return this.stack.remove(0);
			}
		return this.delegate.hasNext()?
				this.delegate.next():
				null;
		}
	
	@Override
	public E remove(int index) {
		fillTo(index);
		if(index < this.stack.size()) {
			return this.stack.remove(index);
			}
		else
			{
			throw new NoSuchElementException("cannot remove "+index+"-th item");
			}
		}
	
	@Override
	public E peek(final int index) {
		fillTo(index);
		return (index<this.stack.size()?this.stack.get(index):null);
		}
	@Override
	public void close() {
		this.stack.clear();
		CloserUtil.close(this.delegate);
		}
	}
