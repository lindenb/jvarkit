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
package com.github.lindenb.jvarkit.util.tabix;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.Set;

import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;

/**
 * Safe wrapper around org.broad.tribble.readers.TabixReader (won't return a null iterator )
 * @author lindenb
 *
 */
public class TabixFileReader implements Closeable
	//,Iterable<VariantContext> NO, not a true iterator
	{
	private static final Logger LOG=Logger.build(TabixFileReader.class).make();
	private TabixReader tabix=null;
    private final String uri;
    
    /** return true if 'f' is a file, path ends with '.gz' and there is an associated .tbi file */
    public static final boolean isValidTabixFile(final File f)
		 	{
			return VCFUtils.isTabixVcfFile(f);
			}
    
    public TabixFileReader(final String uri) throws IOException
    	{
    	this.uri=uri;
    	this.tabix=new TabixReader(uri);
    	}
    
    public TabixFileReader(final String uri,final String idxFn) throws IOException
		{
		this.uri=uri;
		this.tabix=new TabixReader(uri,idxFn);
		}

    
    /** returns the underlying TabixReader */
    public TabixReader getTabix() {
		return tabix;
		}
    
    /** return contigs associated to the tabix file*/
    public Set<String> getChromosomes()
    	{
    	return getTabix().getChromosomes();
    	}
    
    /** return URL of the TABIX file */
    public String getURI()
    	{
    	return this.uri;
    	}
    
    public String readLine() throws IOException
    	{
    	if(isClosed()) return null;
    	return this.tabix.readLine();
    	}
    
    protected int[] _parseReg(final  String rgn)
    	{
    	if(isClosed()) return null;
    	final int parseReg[]=this.tabix.parseReg(rgn);
    	if(parseReg==null || parseReg.length!=3 ||
				parseReg[0]==-1 || parseReg[1]>parseReg[2])
			{
    		if(parseReg[0]==-1)
    			{
    			//LOG.warning("unknown chromosome in \""+rgn+"\". Available are: "+getChromosomes());
    			}
    		//LOG.warning("cannot parse region "+rgn);
			return null;
			}
    	return parseReg;
    	}
    
    public Iterator<String> iterator(final String chrom)
		{
    	if(isClosed()) return Collections.emptyIterator();
    	return iterator(_parseReg(chrom));
		}

    public Iterator<String> iterator(final String chrom,final int start)
		{
    	if(isClosed()) return Collections.emptyIterator();
    	return iterator(_parseReg(chrom+":"+start));
		}
    public Iterator<String> iterator(final String chrom,final int start ,final  int end)
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
		final TabixReader.Iterator titer=this.tabix.query(parseReg[0], parseReg[1],parseReg[2]);
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
    	extends AbstractIterator<String>
    	{
    	TabixReader.Iterator delegate;
    	MyIterator(final TabixReader.Iterator delegate)
    		{
    		this.delegate=delegate;
    		}
    	
    	@Override
    	protected String advance()
    		{
    		try
    			{
    			if(isClosed() || delegate==null ) return null;
    			final String s= delegate.next();
    			if(s==null ) delegate=null;
    			return s;
    			}
    		catch(final IOException err)
    			{
    			throw new RuntimeIOException(err);
    			}	
    		}
    	}	
    
    @Override
    public String toString() {
    	return getURI();
    	}
	}
