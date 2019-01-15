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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.util;

import java.util.AbstractList;
import java.util.Arrays;
import java.util.List;

/**
 * Instance of java.util.List with fast push_back, pop_front capabilities
 * not all methods currently implemented
 * @author lindenb
 *
 */
public class BufferedList<T>
extends AbstractList<T>
	{
	private int start=0;
	private int length=0;
	private Object container[]=null;
	private void resize()
		{
		
		if( this.container.length>10 && this.start>=this.container.length/2)
			{
			System.arraycopy(
					this.container, start,
					this.container, 0,
					this.length);
			//help GC collectory, setting data to null
			Arrays.fill(this.container, this.length,this.container.length,null);
			this.start=0;
			}
		}

	
	public BufferedList()
		{
		this(100);
		}
	public BufferedList(int capacity)
		{
		this.container=new Object[capacity];
		}
	
	public BufferedList(List<T> L)
		{
		this.container=new Object[L.size()];
		for(int i=0;i< L.size();++i)
			{
			this.container[i]=L.get(i);
			}
		this.length=L.size();
		}
	
	public T getFirst()
		{
		return this.get(0);
		}
	
	public T getLast()
		{
		return this.get(size()-1);
		}	
	
	@Override
	@SuppressWarnings("unchecked")
	public T remove(int index)
		{
		if(index<0 || index>=size()) throw new IndexOutOfBoundsException();
		if(index==0) return removeFirst();
		if(index+1==this.size()) return removeLast();
		T old= (T)this.container[this.start+index];
		//shift left
		System.arraycopy(this.container,
				this.start+index+1,//Src
				this.container,
				this.start+index,//dest
				(this.length-(index+1))
				);
		this.container[this.start+this.length-1]=null;
		this.length--;
		return old;
		}
	
	@SuppressWarnings("unchecked")
	public final T removeFirst()
		{
		if(isEmpty()) throw new IllegalStateException();
		T old= (T)this.container[this.start];
		this.container[this.start]=null;
		this.start++;
		this.length--;
		resize();
		return old;
		}
	@SuppressWarnings("unchecked")
	public final  T removeLast()
		{
		if(isEmpty()) throw new IllegalStateException();
		T old= (T)this.container[this.start+this.length-1];
		this.container[this.start+this.length-1]=null;
		this.length--;
		return old;
		}
	
	@Override
	public boolean add(T e)
		{
		if(this.start+this.length>=this.container.length)
			{
			int new_length = 1+this.length+this.length/2;
			Object[] copy=new Object[new_length];
			System.arraycopy(this.container,this.start,copy, 0, this.length);
			this.start=0;
			this.container=copy;
			}
		this.container[this.start+this.length]=e;
		this.length++;
		return true;
		}
	
	@SuppressWarnings("unchecked")
	@Override
	public T get(int index) {
		if(index<0 || index>=size())
			{
			throw new IndexOutOfBoundsException("0<"+index+"<="+size());
			}
		return (T)this.container[this.start+index];
		}

	@Override
	public int size() {
		return this.length;
		}
	
	@Override
	public void clear() {
		this.start=0;
		this.length=0;
		Arrays.fill(this.container, null);
		}
	

	}
