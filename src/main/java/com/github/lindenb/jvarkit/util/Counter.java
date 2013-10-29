package com.github.lindenb.jvarkit.util;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Counter<T>
	{
	private Map<T,Long> object2count=new HashMap<T,Long>();
	private long total=0L;
	
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
	@Override
	public String toString() {
		return "Counter "+this.getTotal();
		}
	}
