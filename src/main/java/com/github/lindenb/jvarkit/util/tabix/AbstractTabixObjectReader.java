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

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.Locatable;


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
    
    public Iterator<T> iterator(final Locatable interval)
		{
		return iterator(interval.getContig() , interval.getStart() , interval.getEnd());
		}

    public Iterator<T> iterator(final String contig,int start,int end)
		{
		return iterator(tabix.iterator(contig+":"+start+"-"+end));
		}
    
    /** returns a list of items in the interval */
    public List<T> getItemsInInterval(final Locatable interval)
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
    	extends AbstractIterator<T>
    	{
    	protected final Iterator<String> delegate;
    	protected AbstractMyIterator(final Iterator<String> delegate)
    		{
    		this.delegate=delegate;
    		}
    	
    	/** convert line to item. item is ignored if it returns null */
    	protected abstract T convert(String line) ;
    	
    	@Override
    	protected T advance() {
    		while(this.delegate.hasNext()) {
    			final String line  = this.delegate.next();
    			final T item = convert(line);
    			if(item==null) continue;
    			return item;
    			}
    		return null;
    		}

    	}	
    
    @Override
    public String toString() {
    	return getURI();
    	}
	}
