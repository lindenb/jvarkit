package com.github.lindenb.jvarkit.util.tabix;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;

import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;


public abstract class AbstractTabixObjectReader<T> implements Closeable
	//,Iterable<VariantContext> NO, not a true iterator
	{
    protected TabixFileReader tabix=null;
    protected String uri;
  
    
    protected AbstractTabixObjectReader(String uri) throws IOException
    	{
    	this.uri=uri;
    	this.tabix=new TabixFileReader(uri);
    	}
    
    public String getURI()
    	{
    	return this.uri;
    	}
    
    public Iterator<T> iterator(String chrom)
		{
    	return iterator(tabix.iterator(chrom));
		}

    public Iterator<T> iterator(String chrom,int start)
		{
    	return iterator(tabix.iterator(chrom+":"+start));
		}
    public Iterator<T> iterator(String chrom,int start,int end)
    	{
    	return iterator(tabix.iterator(chrom+":"+start+"-"+end));
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
    	protected Iterator<String> delegate;
    	protected AbstractMyIterator(Iterator<String> delegate)
    		{
    		this.delegate=delegate;
    		}
    	@Override
    	public boolean hasNext()
    		{
    		return delegate.hasNext();
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
