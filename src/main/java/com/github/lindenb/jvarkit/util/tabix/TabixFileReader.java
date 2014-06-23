package com.github.lindenb.jvarkit.util.tabix;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.Set;
import java.util.logging.Logger;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.tribble.readers.TabixReader;

/**
 * Safe wrapper around org.broad.tribble.readers.TabixReader (won't return a null iterator )
 * @author lindenb
 *
 */
public class TabixFileReader implements Closeable
	//,Iterable<VariantContext> NO, not a true iterator
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	private TabixReader tabix=null;
    private String uri;
    
    /** return true if 'f' is a file, path ends with '.gz' and there is an associated .tbi file */
    public static final boolean isValidTabixFile(File f)
    	{
    	return VCFUtils.isTabixVcfFile(f);
    	}
    
    public TabixFileReader(String uri) throws IOException
    	{
    	this.uri=uri;
    	this.tabix=new TabixReader(uri);
    	}
    
    public TabixReader getTabix() {
		return tabix;
		}
    
    public Set<String> getChromosomes()
    	{
    	return Collections.unmodifiableSet(getTabix().mChr2tid.keySet());
    	}
    
    public String getURI()
    	{
    	return this.uri;
    	}
    
    public String readLine() throws IOException
    	{
    	if(isClosed()) return null;
    	return this.tabix.readLine();
    	}
    
    protected int[] _parseReg(String rgn)
    	{
    	if(isClosed()) return null;
    	int parseReg[]=this.tabix.parseReg(rgn);
    	if(parseReg==null || parseReg.length!=3 ||
				parseReg[0]==-1 || parseReg[1]>parseReg[2])
			{
    		if(parseReg[0]==-1)
    			{
    			LOG.warning("unknown chromosome in \""+rgn+"\". Available are: "+getChromosomes());
    			}
    		LOG.warning("cannot parse region "+rgn);
			return null;
			}
    	return parseReg;
    	}
    
    public Iterator<String> iterator(String chrom)
		{
    	if(isClosed()) return Collections.emptyIterator();
    	return iterator(_parseReg(chrom));
		}

    public Iterator<String> iterator(String chrom,int start)
		{
    	if(isClosed()) return Collections.emptyIterator();
    	return iterator(_parseReg(chrom+":"+start));
		}
    public Iterator<String> iterator(String chrom,int start,int end)
    	{
    	if(isClosed()) return Collections.emptyIterator();
    	return iterator(_parseReg(chrom+":"+start+"-"+end));
    	}
    private  Iterator<String> iterator(int parseReg[])
    	{
    	if(isClosed()) return Collections.emptyIterator();
    	if(parseReg==null)
    		{
			return Collections.emptyIterator();
			}
		TabixReader.Iterator titer=this.tabix.query(parseReg[0], parseReg[1],parseReg[2]);
		if(titer==null)
			{
			return Collections.emptyIterator();
			}
		return new MyIterator(titer);
		}
    
    @Override
    public void close()
    	{
    	if(tabix!=null) this.tabix.close();
    	tabix=null;
    	}
    
    public boolean isClosed()
    	{
    	return tabix==null;
    	}
    
    private class MyIterator
    	implements Iterator<String>
    	{
    	TabixReader.Iterator delegate;
    	String _next=null;
    	MyIterator(TabixReader.Iterator delegate)
    		{
    		this.delegate=delegate;
    		}
    	@Override
    	public boolean hasNext()
    		{
    		if(isClosed()) return false;
    		if(delegate==null) return false;
    		if(_next==null)
    			{
    			try
    				{
    				_next=delegate.next();
    				}
    			catch(IOException err)
    				{
    				throw new PicardException("Tabix",err);
    				}
    			if(_next==null) delegate=null;
    			}
    		return _next!=null;
    		}
    	@Override
    	public String next() {
    		if(!hasNext()) throw new IllegalStateException(
    				"no next iterator: closed:"+isClosed()+ " delegate.null:"+(delegate==null) 
    				);
    		String s=_next;
    		_next=null;
    		return s;
    		}
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
