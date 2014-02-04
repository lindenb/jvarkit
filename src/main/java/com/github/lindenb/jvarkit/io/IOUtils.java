package com.github.lindenb.jvarkit.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PushbackInputStream;
import java.net.URL;
import java.nio.charset.Charset;
import java.util.zip.GZIPInputStream;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.LineReaderUtil;


import net.sf.samtools.Defaults;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.BlockCompressedOutputStream;

public class IOUtils {
	
	public static void copyTo(InputStream in,File f) throws IOException
		{
		FileOutputStream fous=new FileOutputStream(f);
		copyTo(in,fous);
		fous.flush();
		fous.close();
		}
	public static void copyTo(InputStream in,OutputStream out) throws IOException
		{
		byte buffer[]=new byte[2048];
		int nRead;
		while((nRead=in.read(buffer))!=-1)
			{
			out.write(buffer,0,nRead);
			}
		out.flush();
		}

	
	
	public static boolean isRemoteURI(String uri)
		{	
	    return  uri.startsWith("http://") ||
				uri.startsWith("https://") ||
				uri.startsWith("ftp://")
				;
		}
	private static InputStream tryBGZIP(InputStream in) throws IOException
		{
		byte buffer[]=new byte[2048];
		
		PushbackInputStream push_back=new PushbackInputStream(in,buffer.length+10);
		int nReads=push_back.read(buffer);
		push_back.unread(buffer, 0, nReads);
		
		try
			{
			BlockCompressedInputStream bgz=new BlockCompressedInputStream(new ByteArrayInputStream(buffer, 0, nReads));
			bgz.read();
			bgz.close();
			return new BlockCompressedInputStream(push_back);
			}
		catch(Exception err)
			{
			return new GZIPInputStream(push_back);
			}
		}
	
	public static InputStream openURIForReading(String uri) throws IOException
		{
		if(isRemoteURI(uri))
			{
			URL url=new URL(uri);
			InputStream in=url.openStream();
			//do we have .... azdazpdoazkd.vcf.gz?param=1&param=2
			int question=uri.indexOf('?');
			if(question!=-1) uri=uri.substring(0, question);
			if(uri.endsWith(".gz"))
				{
				return tryBGZIP(in);
				}
			return in;
			}
		if(uri.startsWith("file://"))
			{
			uri=uri.substring(7);
			}
		return openFileForReading(new File(uri));
		}
	
	public static BufferedReader openURIForBufferedReading(String uri) throws IOException
		{
		return  new BufferedReader(new InputStreamReader(openURIForReading(uri), Charset.forName("UTF-8")));
		}


	
	@SuppressWarnings("resource")
	public static InputStream openFileForReading(File file) throws IOException
		{
		net.sf.picard.io.IoUtil.assertFileIsReadable(file);
		InputStream in= new FileInputStream(file);
		if(file.getName().endsWith(".gz"))
			{
			in=tryBGZIP(in);
			}
		return in;
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
        if (file.getName().endsWith(".vcf.gz"))
        	{
            return new BlockCompressedOutputStream(file);
        	}
        else
        	{
            return new FileOutputStream(file);
        	}         
    	}
    
    public static LineReader openFileForLineReader(File file) throws IOException
		{
		return  LineReaderUtil.fromBufferedStream(openFileForReading(file));
		}
    /** @return a LineIterator that should be closed with CloserUtils */
    public static LineIterator openFileForLineIterator(File file) throws IOException
  		{
  		return  new LineIteratorImpl(openFileForLineReader(file));
  		}
    
    public static LineReader openStdinForLineReader() throws IOException
		{
		return  LineReaderUtil.fromBufferedStream(System.in);
		}
    /** @return a LineIterator that should be closed with CloserUtils */
    public static LineIterator openStdinForLineIterator() throws IOException
  		{
  		return  new LineIteratorImpl(openStdinForLineReader());
  		}

    
    public static LineReader openURIForLineReader(String uri) throws IOException
		{
		return  LineReaderUtil.fromBufferedStream(openURIForReading(uri));
		}
    /** @return a LineIterator that should be closed with CloserUtils */
    public static LineIterator openURIForLineIterator(String uri) throws IOException
  		{
  		return  new LineIteratorImpl(openURIForLineReader(uri));
  		}

    
    
}
