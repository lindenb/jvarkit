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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.tribble.readers.LineReader;

import com.github.lindenb.jvarkit.io.IOUtils;

/**
 * the original picard FastqReader didn't allow empty lines... I created that file.
 * Use the new htjk {@link FastqReader}
 */
public class FourLinesFastqReader
	extends AbstractFastqReader
	{
    private final LineReader lineReader;
    private long nLines=0;
   
    
    public FourLinesFastqReader(final File file) throws IOException
    	{
    	super(file);
    	this.lineReader= IOUtils.openFileForLineReader(file);
    	}
    
    public FourLinesFastqReader(final InputStream in)
		{
		super(null);
		try {
			this.lineReader= IOUtils.openStreamForLineReader(in);
		} catch (final IOException e) {
			throw new RuntimeIOException(e);
			}
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
                throw new RuntimeException(error("Missing sequence header"));
            	}
            
            if (!this.seqHeader.startsWith(FastqConstants.SEQUENCE_HEADER))
            	{
                throw new RuntimeException(error("Sequence header must start with "+ FastqConstants.SEQUENCE_HEADER));
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
                throw new RuntimeException(error("Quality header must start with "+ FastqConstants.QUALITY_HEADER+": "+qualHeader));
            	}

            // Read quality line
            final String qualLine = this.readLine();
            ++nLines;
            checkLine(qualLine,"quality line");

            // Check sequence and quality lines are same length
            if (seqLine.length() != qualLine.length()) {
                throw new RuntimeException(error("Sequence and quality line must be the same length"));
            }

            final FastqRecord frec = new FastqRecord(seqHeader.substring(1, seqHeader.length()), seqLine,
                    qualHeader.substring(1, qualHeader.length()), qualLine);
            this.seqHeader=null;
            return frec ;

        } catch (IOException e) {
            throw new RuntimeException(String.format("Error reading fastq '%s'", getAbsolutePath()), e);
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
