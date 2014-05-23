package com.github.lindenb.jvarkit.util.picard;

import com.github.lindenb.jvarkit.util.picard.PicardException;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.fastq.FastqRecord;
import java.util.NoSuchElementException;
import java.util.logging.Logger;
import java.io.*;


/**
 * the original picard FastqReader doesn't allow empty lines...
 */
public abstract class AbstractFastqReader
	implements FastqReader
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	private File fastqFile=null;
    private ValidationStringency validationStringency=ValidationStringency.STRICT;
    private FastqRecord nextRecord=null;
    protected String seqHeader=null;
    
    protected AbstractFastqReader(File fastqFile)
    	{
    	this.nextRecord=null;
    	this.fastqFile=fastqFile;
    	}

    
    public void setValidationStringency( ValidationStringency validationStringency) {
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
    		case STRICT: throw new PicardException(msg);
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
            throw new PicardException(error("File is too short - missing "+kind+" line"));
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
