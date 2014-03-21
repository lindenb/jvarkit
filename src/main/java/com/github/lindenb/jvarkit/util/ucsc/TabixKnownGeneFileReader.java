package com.github.lindenb.jvarkit.util.ucsc;

import java.io.IOException;
import java.util.Iterator;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;


public class TabixKnownGeneFileReader extends AbstractTabixObjectReader<KnownGene>
	{
    public TabixKnownGeneFileReader(String uri) throws IOException
    	{
    	super(uri);
    	}
    @Override
    protected  Iterator<KnownGene> iterator(Iterator<String> delegate)
		{
		return new MyIterator(delegate);
		}
    
    private class MyIterator
    	extends AbstractMyIterator
    	{
    	private Pattern tab=Pattern.compile("[\t]");
    	MyIterator(Iterator<String> delegate)
    		{
    		super(delegate);
    		}
    	
    	@Override
    	public KnownGene next() {
    		return new KnownGene(this.tab.split(delegate.next()));
    		}
    	}
	}
