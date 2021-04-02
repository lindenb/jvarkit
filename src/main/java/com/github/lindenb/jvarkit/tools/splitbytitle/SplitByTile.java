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
*/
package com.github.lindenb.jvarkit.tools.splitbytitle;


/**
BEGIN_DOC


END_DOC
 */

import java.util.Collections;
import java.util.Set;

import com.github.lindenb.jvarkit.jcommander.AbstractBamSplitter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMRecord;

@Program(name="splitbytile",description="Split Bam By tile",
		keywords={"sam","bam"},
		modificationDate="20210401",
		creationDate="20130406"
		)
public class SplitByTile  extends AbstractBamSplitter<Integer> {
    private static final Logger LOG = Logger.build(SplitByTile.class).make();
    
    
    @Override
    protected Logger getLogger() {
    	return LOG;
    	}
   
    @Override
    protected Set<Integer> createKeys(final SAMRecord rec) {
    	final CharSplitter  colon= CharSplitter.COLON;
     	final String tokens[]=colon.split(rec.getReadName(),6);
		if(tokens.length<5)
				{
	    		LOG.error("Cannot get the 6th field in "+rec.getReadName());
	    		return null;
				}
		int tile=-1;
		try {
			tile=Integer.parseInt(tokens[4]);
			}
		catch(final Throwable e)
			{
			tile=-1;
			}
		if(tile<0)
			{
			LOG.error("Bad tile in read: "+rec.getReadName());
			return null;
			}
		return Collections.singleton(tile);
    	}
	
    public static void main(final String[] argv)
		{
	    new SplitByTile().instanceMainWithExit(argv);
		}

	}
