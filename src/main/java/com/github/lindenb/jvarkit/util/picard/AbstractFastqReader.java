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
package com.github.lindenb.jvarkit.util.picard;


import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.fastq.FastqRecord;
import java.util.NoSuchElementException;

import com.github.lindenb.jvarkit.util.log.Logger;

import java.io.*;


/**
 * the original picard FastqReader doesn't allow empty lines...
 */
public abstract class AbstractFastqReader
	implements FastqReader
	{
	private static final Logger LOG=Logger.build(AbstractFastqReader.class).make();
	private File fastqFile=null;
    private ValidationStringency validationStringency=ValidationStringency.STRICT;
    private FastqRecord nextRecord=null;
    protected String seqHeader=null;
    
    protected AbstractFastqReader(final File fastqFile)
    	{
    	this.nextRecord=null;
    	this.fastqFile=fastqFile;
    	}

    
    public void setValidationStringency(final ValidationStringency validationStringency) {
		this.validationStringency = validationStringency;
		}
    
    public ValidationStringency getValidationStringency()
    	{
		return validationStringency;
		}
    

    
    protected void throw_error(String msg)
    	{
    	switch(getValidationStringency())
    		{
    		case LENIENT: LOG.warning(msg);break;
    		case STRICT: throw new RuntimeException(msg);
    	    default:break;
    		}
    	}
    
    protected abstract FastqRecord readNextRecord();
    
    private boolean _has_next_called=false;
    

    public boolean hasNext()
    	{
    	if(!_has_next_called)
    		{
    		this._has_next_called=true;
    		this.nextRecord= readNextRecord();
    		}
    	return nextRecord != null;
    	}

    public FastqRecord next()
    	{
        if (!hasNext())
        	{
            throw new NoSuchElementException("next() called when !hasNext()");
        	}
        final FastqRecord rec = nextRecord;
        nextRecord = null;
        _has_next_called=false;
        seqHeader=null;
        return rec;
    	}

    @Override
    public void remove() { throw new UnsupportedOperationException("Unsupported operation"); }


    /**
     * @return Name of FASTQ being read, or null if not known.
     */
    public File getFile() { return fastqFile ; }

   

    protected void checkLine(final String line, final String kind)
    	{
        if (line == null)
        	{
            throw new RuntimeException(error("File is too short - missing "+kind+" line"));
        	}
        if (StringUtil.isBlank(line))
        	{
            throw_error(error("Missing "+kind));
        	}
    	}
    
    protected abstract String location();
    

    protected String error(final String msg)
    	{
        return msg + " at line "+location()+" in fastq \""+getAbsolutePath()+"\"."+
        		( this.seqHeader==null?"": " Read name:"+ this.seqHeader)
        		;
    	}

    protected String getAbsolutePath()
    	{
        return (fastqFile == null ?"":fastqFile.getAbsolutePath());
    	}
	}
