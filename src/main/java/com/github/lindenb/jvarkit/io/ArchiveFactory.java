/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.io;

import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;


public abstract class ArchiveFactory
	implements Closeable{

	private ArchiveFactory()
		{
		}
	
	public abstract OutputStream openOuputStream(final String filename) throws IOException;
	public PrintWriter openWriter(final String filename) throws IOException
		{
		return new PrintWriter(openOuputStream(filename), true);
		}

	
	public static ArchiveFactory open(final File f)  throws IOException
		{
		if( f == null ) throw new IllegalArgumentException("Cannot open(null)");
		if(f.getName().toLowerCase().endsWith(".zip"))
			{
			return new ZipInstance(f);
			}
		else
			{
			return new FileInstance(f);
			}
		}
	
	
	private static class ZipInstance
		extends ArchiveFactory
		{
		FileOutputStream fout;
		ZipOutputStream zout;
		
		ZipInstance(final File f) throws IOException
			{
			fout=new FileOutputStream(f);
			zout=new ZipOutputStream(fout);
			}
		
		@Override
		public OutputStream openOuputStream(String filename) throws IOException
			{
			ZipOS os =new ZipOS(filename);
			return os;
			}
		
		@Override
		public void close() throws IOException
			{
			if(zout!=null)
				{
				zout.finish();
				zout.flush();
				fout.flush();
				zout.close();
				fout.close();
				zout=null;
				fout=null;
				}
			}
		
		private class ZipOS extends OutputStream
			{
			ZipEntry ze;
			File tmp;
			OutputStream out;
			
			ZipOS(String filename) throws IOException
				{
				while(filename.startsWith("/")) filename=filename.substring(1);
				this.ze=new ZipEntry(filename);
				
				tmp=File.createTempFile("tmp", ".zipentry");
				tmp.deleteOnExit();
				out=new FileOutputStream(tmp);
				}
			
			@Override
			public void write(int b) throws IOException {
				if(out!=null) out.write(b);
				}
			
			@Override
			public void write(byte[] b) throws IOException {
				if(out!=null)  out.write(b);
				}
			
			@Override
			public void write(byte[] b, int off, int len) throws IOException {
				if(out!=null) out.write(b, off, len);
				}
			
			@Override
			public void flush() throws IOException
				{
				if(out!=null) out.flush();
				}
			
			@Override
			public void close() throws IOException
				{
				if(out!=null)
					{
					out.flush();
					out.close();
					
					if(ZipInstance.this.zout!=null)
						{
						zout.putNextEntry(this.ze);
						FileInputStream fin=new FileInputStream(tmp);
						IOUtils.copyTo(fin,ZipInstance.this.zout );
						fin.close();
						ZipInstance.this.zout.flush();
						ZipInstance.this.zout.closeEntry();
						}
					ze=null;
					out=null;
					tmp.delete();
					tmp=null;
					}
				}
			}
		
		}
	
	private static class FileInstance
	extends ArchiveFactory
		{
		private final File baseDir;
		
		FileInstance(final File baseDir) throws IOException
			{
			this.baseDir=baseDir;
			if(baseDir.exists() && !baseDir.isDirectory())
				{
				throw new IOException("Not a directory:"+baseDir);
				}
				
			}
		
		@Override
		public OutputStream openOuputStream(String filename) throws IOException
			{
			while(filename.startsWith("/")) filename=filename.substring(1);
			File f=new File(baseDir, filename);
			if(f.getParentFile()!=null)
				{
				f.getParentFile().mkdirs();
				}
			return new FileOutputStream(f);
			}
		@Override
		public void close() throws IOException
			{
			
			}
		}
}
