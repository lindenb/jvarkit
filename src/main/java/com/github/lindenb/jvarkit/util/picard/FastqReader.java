package com.github.lindenb.jvarkit.util.picard;

import net.sf.picard.PicardException;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.RuntimeIOException;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.fastq.FastqConstants;
import net.sf.picard.fastq.FastqRecord;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.logging.Logger;
import java.io.*;

import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.LineReaderUtil;


import com.github.lindenb.jvarkit.io.IOUtils;

/**
 * the original picard FastqReader doesn't allow empty lines...
 */
public class FastqReader implements Iterator<FastqRecord>, Closeable
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
    final private File fastqFile;
    private FastqRecord nextRecord;
    private LineReader lineReader;
    private ValidationStringency validationStringency=ValidationStringency.STRICT;
    private long nLines=0;
    private String seqHeader=null;
    
    public FastqReader(final File file)
    	{
        try {
        	this.fastqFile = file;
        	this.lineReader=LineReaderUtil.fromBufferedStream(IOUtils.openFileForReading(file));
        	}
        catch (IOException ioe) 
        	{
            throw new RuntimeIOException(ioe);
        	}
    	}
    
    public FastqReader(InputStream in)
		{
		this.fastqFile = null;
		this.lineReader=LineReaderUtil.fromBufferedStream(in);
		}

    
    public void setValidationStringency( ValidationStringency validationStringency) {
		this.validationStringency = validationStringency;
		}
    
    public ValidationStringency getValidationStringency()
    	{
		return validationStringency;
		}
    

    private String readLine() throws IOException
    	{	
    	return this.lineReader.readLine();
    	}
    
    private void throw_error(String msg)
    	{
    	switch(getValidationStringency())
    		{
    		case LENIENT: LOG.warning(msg);break;
    		case STRICT: throw new PicardException(msg);
    	    default:break;
    		}
    	}
    
    private FastqRecord readNextRecord() {
        try {

            // Read sequence header
        	this.seqHeader = this.readLine();
            if (this.seqHeader == null) return null ;
            
            ++nLines;
            
            if (StringUtil.isBlank(this.seqHeader))
            	{
                throw new PicardException(error("Missing sequence header"));
            	}
            
            if (!this.seqHeader.startsWith(FastqConstants.SEQUENCE_HEADER))
            	{
                throw new PicardException(error("Sequence header must start with "+ FastqConstants.SEQUENCE_HEADER));
            	}

            // Read sequence line
            final String seqLine = this.readLine();
            ++nLines;
            checkLine(seqLine,"sequence line");

            // Read quality header
            final String qualHeader = this.readLine();
            ++nLines;
            
            checkLine(qualHeader,"quality header");
            if (!qualHeader.startsWith(FastqConstants.QUALITY_HEADER))
            	{
                throw new PicardException(error("Quality header must start with "+ FastqConstants.QUALITY_HEADER+": "+qualHeader));
            	}

            // Read quality line
            final String qualLine = this.readLine();
            ++nLines;
            checkLine(qualLine,"quality line");

            // Check sequence and quality lines are same length
            if (seqLine.length() != qualLine.length()) {
                throw new PicardException(error("Sequence and quality line must be the same length"));
            }

            final FastqRecord frec = new FastqRecord(seqHeader.substring(1, seqHeader.length()), seqLine,
                    qualHeader.substring(1, qualHeader.length()), qualLine);
            this.seqHeader=null;
            return frec ;

        } catch (IOException e) {
            throw new PicardException(String.format("Error reading fastq '%s'", getAbsolutePath()), e);
        }
    }
    
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
        if (!hasNext()) {
            throw new NoSuchElementException("next() called when !hasNext()");
        }
        final FastqRecord rec = nextRecord;
        nextRecord = null;
        _has_next_called=false;
        return rec;
    }

    @Override
    public void remove() { throw new UnsupportedOperationException("Unsupported operation"); }


    public long getLineNumber() { return this.nLines ; }


    /**
     * @return Name of FASTQ being read, or null if not known.
     */
    public File getFile() { return fastqFile ; }

    @Override
    public void close()
    	{
        CloserUtil.close(this.lineReader);
    	}

    private void checkLine(final String line, final String kind) {
        if (line == null) {
            throw new PicardException(error("File is too short - missing "+kind+" line"));
        }
        if (StringUtil.isBlank(line)) {
            throw_error(error("Missing "+kind));
        }
    }

    private String error(final String msg) {
        return msg + " at line "+nLines+" in fastq \""+getAbsolutePath()+"\"."+
        		( this.seqHeader==null?"": " Read name:"+ this.seqHeader)
        		;
    }

    private String getAbsolutePath() {
        if (fastqFile == null) return "";
        else return fastqFile.getAbsolutePath();
    }
}
