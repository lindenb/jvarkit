package com.github.lindenb.jvarkit.util.picard;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;


import net.sf.samtools.Defaults;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.BlockCompressedOutputStream;

public class IOUtils {
	public static InputStream openFileForReading(File file) throws IOException
		{
		net.sf.picard.io.IoUtil.assertFileIsReadable(file);
		if(file.getName().endsWith(".gz"))
			{
			return new BlockCompressedInputStream(file);
			}
		else
			{
			return new FileInputStream(file);
			}
		}
	
	public static BufferedReader openFileForBufferedReading(File file) throws IOException
		{
		return  new BufferedReader(new InputStreamReader(openFileForReading(file), Charset.forName("UTF-8")));
		}

    public static BufferedWriter openFileForBufferedWriting(final File file)  throws IOException
    	{
        return new BufferedWriter(new OutputStreamWriter(openFileForWriting(file)), Defaults.BUFFER_SIZE);
    	}

    public static OutputStream openFileForWriting(final File file) throws IOException
    	{
        if (file.getName().endsWith(".gz"))
        	{
            return new BlockCompressedOutputStream(file);
        	}
        else
        	{
            return new FileOutputStream(file);
        	}         
    	}

}
