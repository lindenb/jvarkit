package com.github.lindenb.jvarkit.util.picard;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqRecord;

import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.LineReaderUtil;


import com.github.lindenb.jvarkit.io.IOUtils;

/**
 * the original picard FastqReader doesn't allow empty lines...
 */
public class FourLinesFastqReader
	extends AbstractFastqReader
	{
	//private static final java.util.logging.Logger LOG=java.util.logging.Logger.getLogger("jvarkit");
    private LineReader lineReader;
    private long nLines=0;
   
    
    public FourLinesFastqReader(final File file) throws IOException
    	{
    	super(file);
    	this.lineReader=LineReaderUtil.fromBufferedStream(IOUtils.openFileForReading(file));
    	}
    
    public FourLinesFastqReader(InputStream in)
		{
		super(null);
		this.lineReader=LineReaderUtil.fromBufferedStream(in);
		}
    

    private String readLine() throws IOException
    	{	
    	return this.lineReader.readLine();
    	}
    
    @Override
    protected FastqRecord readNextRecord()
    	{
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
    
    public long getLineNumber() { return this.nLines ; }
    
    @Override
    protected String location()
    	{
    	return String.valueOf(getLineNumber());
    	}

    @Override
    public void close()
    	{
        CloserUtil.close(this.lineReader);
    	}



	}
