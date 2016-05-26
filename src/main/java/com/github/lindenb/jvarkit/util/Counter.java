/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Counter<T>
	{
	private Map<T,Long> object2count=new HashMap<T,Long>();
	private long total=0L;
	
	public Counter()
		{
		}
	
	public void initializeIfNotExists(final T key)
		{
		initializeIfNotExists(key,0L);
		}
	
	public void initializeIfNotExists(final T key,long initialValue)
		{
		if(!this.object2count.containsKey(key))
			{
			if(initialValue<0) throw new IllegalArgumentException("n<0 :"+initialValue);
			this.object2count.put(key,initialValue);
			this.total+=initialValue;
			}
		}
	
	/** increase by 1 returns the new count */
	public long incr(final T object)
		{
		return incr(object,1L);
		}
	/**  increase by n, returns the new count */
	public long incr(final T object,long n)
		{
		if(n<=0) throw new IllegalArgumentException("n<=0 :"+n);
		if(object==null) throw new IllegalArgumentException("null argument in "+getClass());
		Long count=this.object2count.get(object);
		if(count==null) count=0L;
		//increase
		count+=n;
		this.object2count.put(object,count);
		this.total+=n;
		return count;
		}
	
	public void putAll(Counter<T> other)
		{
		if(this==other)  throw new IllegalArgumentException("cannot put to self");
		for(T k: other.keySet())
			{
			this.incr(k,other.count(k));
			}
		}
	
	public long getTotal()
		{
		return total;
		}
	
	/** count number of times object was seen. returns 0 if object never seen */
	public long count(final T object)
		{
		Long count=this.object2count.get(object);
		return count==null?0L:count;
		}
	
	public Set<T> keySet()
		{
		return this.object2count.keySet();
		}
	
	public T getMostFrequent()
		{
		T key=null;
		for(T o:this.object2count.keySet())
			{
			if(key==null || count(key)< count(o))
				{
				key=o;
				}
			}
		return key;
		}
	
	public List<T> keySetDecreasing()
		{
		List<T> L=new ArrayList<T>(this.object2count.keySet());
		Collections.sort(L, new Comparator<T>()
			{
			@Override
			public int compare(T o1, T o2)
				{
				long n= count(o2)-count(o1);
				return (n<0L?-1:n>0L?1:0);
				}
			});
		return L;
		}
	
	public List<T> keySetIncreasing()
		{
		List<T> L=new ArrayList<T>(this.object2count.keySet());
		Collections.sort(L, new Comparator<T>()
			{
			@Override
			public int compare(T o1, T o2)
				{
				long n= count(o1)-count(o2);
				return (n<0L?-1:n>0L?1:0);
				}
			});
		return L;
		}
	
	/** return the number of categories */
	public int getCountCategories()
		{
		return this.object2count.size();
		}
	
	public boolean isEmpty()
		{
		return this.object2count.isEmpty();
		}
	
	@Override
	public String toString() {
		return "Counter "+this.getTotal();
		}
	}
