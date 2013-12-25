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
	
	public void incr(final T object)
		{
		incr(object,1L);
		}
	public void incr(final T object,long n)
		{
		if(n<=0) throw new IllegalArgumentException("n<=0 :"+n);
		if(object==null) throw new IllegalArgumentException("null argument in "+getClass());
		Long count=this.object2count.get(object);
		if(count==null) count=0L;
		this.object2count.put(object, count+n);
		this.total+=n;
		}
	public long getTotal()
		{
		return total;
		}
	public long count(final T object)
		{
		Long count=this.object2count.get(object);
		return count==null?0L:count;
		}
	
	public Set<T> keySet()
		{
		return this.object2count.keySet();
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
