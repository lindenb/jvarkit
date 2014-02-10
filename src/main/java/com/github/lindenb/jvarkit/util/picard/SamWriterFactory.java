package com.github.lindenb.jvarkit.util.picard;

import java.io.File;
import java.io.OutputStream;
import java.util.logging.Logger;
import java.util.zip.Deflater;

import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMTextWriter;
import net.sf.samtools.util.BlockCompressedOutputStream;

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
        	LOG.info("opening BAM file to stream");
        	BAMFileWriter w= new BAMFileWriter(os,null);
        	w.setHeader(header);
        	return w;
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
        	LOG.info("opening BAM file to "+outputFile);
        	BAMFileWriter w= new BAMFileWriter(outputFile,getCompressionLevel());
        	w.setHeader(header);
        	return w;
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
