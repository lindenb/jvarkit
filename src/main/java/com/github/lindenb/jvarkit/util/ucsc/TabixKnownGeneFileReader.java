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
package com.github.lindenb.jvarkit.util.ucsc;

import java.io.IOException;
import java.util.Iterator;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;


public class TabixKnownGeneFileReader extends AbstractTabixObjectReader<KnownGene>
	{
    public TabixKnownGeneFileReader(final String uri) throws IOException
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
    	private final CharSplitter tab=CharSplitter.TAB;
    	MyIterator(final Iterator<String> delegate)
    		{
    		super(delegate);
    		}
    	
    	@Override
    	public KnownGene next() {
    		return new KnownGene(this.tab.split(delegate.next()));
    		}
    	}
	}

