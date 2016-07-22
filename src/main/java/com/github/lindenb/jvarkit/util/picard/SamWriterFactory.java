/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
import java.io.OutputStream;
import java.util.logging.Logger;
import java.util.zip.Deflater;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMTextWriter;
import htsjdk.samtools.util.BlockCompressedOutputStream;

/** should replace with htsjk now */
@Deprecated
public abstract class SamWriterFactory
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	private boolean binary=false;
	private int compressionLevel= BlockCompressedOutputStream.getDefaultCompressionLevel();
	private SamWriterFactory()
		{
		}
	

	public SamWriterFactory setBinary(boolean binary)
		{
		this.binary=binary;
		return this;
		}
	
	public boolean isBinary()
		{
		return binary;
		}
	
	public SamWriterFactory setUncompressed()
		{
		return setCompressionLevel(Deflater.NO_COMPRESSION);
		}
	
	public SamWriterFactory setBestCompressed()
		{
		return setCompressionLevel(Deflater.BEST_COMPRESSION);
		}
	
	public SamWriterFactory setCompressionLevel(int compressionLevel)
		{
		this.compressionLevel =
				Math.min( Math.max(
						Deflater.NO_COMPRESSION,compressionLevel),
						Deflater.BEST_COMPRESSION);
		return this;
		}
	
	public int getCompressionLevel()
		{
		return compressionLevel;
		}
	
	public static SamWriterFactory newInstance()
		{
		return new SamWriterFactory()
			{
			
			};
		}
	
	public SAMFileWriter make(SAMFileHeader header)
		{
		return make(header, System.out);
		}
	
	public SAMFileWriter make(SAMFileHeader header,OutputStream os)
		{
        if(!isBinary())
        	{
        	LOG.info("opening SAM file to stream");
            SAMTextWriter w=new SAMTextWriter(os);
            w.setHeader(header);
            return w;
        	}
        else
        	{
        	SAMFileWriterFactory swf=new SAMFileWriterFactory();
        	LOG.info("opening BAM file to stream");
        	return swf.makeBAMWriter(header, false, os);
        	}
        	
		}
	public SAMFileWriter make(SAMFileHeader header,File outputFile)
		{
		final String filename = outputFile.getName().toLowerCase();
        boolean sam_output=true;
        
		if (filename.endsWith(".bam"))
        	{
        	sam_output=false;
        	}
		else if (filename.endsWith(".sam"))
        	{
        	sam_output=true;
        	}
		else if(isBinary())
        	{
        	sam_output=false;
        	}
		
        if(sam_output)
        	{
        	LOG.info("opening SAM file to "+outputFile);
            SAMTextWriter w=new SAMTextWriter(outputFile);
            w.setHeader(header);
            return w;
        	}
        else
        	{
        	SAMFileWriterFactory swf=new SAMFileWriterFactory();
        	swf.setCreateIndex(false);
        	swf.setCreateMd5File(false);
        	LOG.info("opening BAM file to "+outputFile);
        	return swf.makeBAMWriter(header, false,outputFile, getCompressionLevel());
        	}
		}
	
	@Override
	protected SamWriterFactory clone()
		{
		SamWriterFactory swf=new SamWriterFactory() {};
		swf.setBinary( this.isBinary() );
		swf.setCompressionLevel( this.getCompressionLevel() );
		return swf;
		}
	
	@Override
	public String toString()
		{
		return "SamWriterFactory " +
				" compression:"+getCompressionLevel()+
				" binary:"+isBinary();
		}
	
	}
