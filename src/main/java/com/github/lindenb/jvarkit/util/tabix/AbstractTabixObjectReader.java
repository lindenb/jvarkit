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

*/
package com.github.lindenb.jvarkit.util.tabix;

import java.io.Closeable;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;

import htsjdk.samtools.util.Interval;


public abstract class AbstractTabixObjectReader<T> implements Closeable
	//,Iterable<VariantContext> NO, not a true iterator
	{
    protected final TabixFileReader tabix;
    protected final String uri;
  
    
    protected AbstractTabixObjectReader(final String uri) throws IOException
    	{
    	this.uri=uri;
    	this.tabix=new TabixFileReader(uri);
    	}
    
    public String getURI()
    	{
    	return this.uri;
    	}
    
    public Set<String> getChromosomes()
    	{
    	return this.tabix.getChromosomes();
    	}
    
    public Iterator<T> iterator(final String chrom)
		{
    	return iterator(tabix.iterator(chrom));
		}

    public Iterator<T> iterator(final String chrom,int start)
		{
    	return iterator(tabix.iterator(chrom+":"+start));
		}
    public Iterator<T> iterator(final String chrom,int start,int end)
    	{
    	return iterator(tabix.iterator(chrom+":"+start+"-"+end));
    	}
    
    public Iterator<T> iterator(final Interval interval)
		{
		return iterator(tabix.iterator(interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd()));
		}
    
    /** returns a list of items in the interval */
    public List<T> getItemsInInterval(final Interval interval)
		{
    	final Iterator<T> t = this.iterator(interval);
    	if(!t.hasNext()) return Collections.emptyList();
    	final List<T> items= new ArrayList<>();
    	while(t.hasNext())
    		{
    		items.add(t.next());
    		}
    	return items;
		}

    
    protected  abstract Iterator<T> iterator(Iterator<String> delegate);
    
    @Override
    public void close()
    	{
    	this.tabix.close();
    	}
    
    protected abstract class AbstractMyIterator
    	implements Iterator<T>
    	{
    	protected final Iterator<String> delegate;
    	protected AbstractMyIterator(final Iterator<String> delegate)
    		{
    		this.delegate=delegate;
    		}
    	@Override
    	public boolean hasNext()
    		{
    		return this.delegate.hasNext();
    		}
    	@Override
    	public abstract T next();
    	@Override
    	public void remove() {
    		throw new UnsupportedOperationException();
    		}
    	}	
    
    @Override
    public String toString() {
    	return getURI();
    	}
	}
