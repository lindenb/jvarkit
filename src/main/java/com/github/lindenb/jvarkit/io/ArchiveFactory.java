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
package com.github.lindenb.jvarkit.io;

import java.io.Closeable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;


public interface ArchiveFactory extends Closeable{

	public static final String OPT_DESC="An existing directory or a filename ending with the '.zip' suffix.";
	
	public abstract OutputStream openOuputStream(final String filename) throws IOException;
	
	/** copy whole file into the archive */
	public default void copyTo(final Path externalFile,final String filename) throws IOException {
		try ( OutputStream os = openOuputStream(filename))
			{
			IOUtils.copyTo(externalFile, os);
			os.flush();
			}
		}
	
	public default PrintWriter openWriter(final String filename) throws IOException
		{
		return new PrintWriter(openOuputStream(filename), true);
		}

	
	/** open a new ArchiveFactory, if filename ends with '.zip' it will be a zip instance
	 * otherwise it will be a FileInstance */
	public static ArchiveFactory open(final Path f)  throws IOException
		{
		if( f == null ) throw new IllegalArgumentException("Cannot open(null)");
		if(f.getFileName().toString().toLowerCase().endsWith(".zip"))
			{
			return new ZipInstance(f);
			}
		else
			{
			return new FileInstance(f);
			}
		}

	
	/** open a new ArchiveFactory, if filename ends with '.zip' it will be a zip instance
	 * otherwise it will be a FileInstance */
	public static ArchiveFactory open(final File f)  throws IOException
		{
		return open(f==null?null:f.toPath());
		}
	
	
	static class ZipInstance
		implements ArchiveFactory
		{
		OutputStream fout;
		ZipOutputStream zout;
		
		ZipInstance(final Path f) throws IOException
			{
			fout= Files.newOutputStream(f);
			zout=new ZipOutputStream(fout);
			}
		
		@Override
		public OutputStream openOuputStream(final String filename) throws IOException
			{
			final ZipOS os =new ZipOS(filename);
			return os;
			}
		
		@Override
		public void close() throws IOException
			{
			if(this.zout!=null)
				{
				this.zout.finish();
				this.zout.flush();
				this.fout.flush();
				this.zout.close();
				this.fout.close();
				this.zout=null;
				this.fout=null;
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
				
				this.tmp=File.createTempFile("tmp", ".zipentry");
				this.tmp.deleteOnExit();
				this.out=new FileOutputStream(this.tmp);
				}
			
			@Override
			public void write(final int b) throws IOException {
				if(out!=null) out.write(b);
				}
			
			@Override
			public void write(final byte[] b) throws IOException {
				if(out!=null)  out.write(b);
				}
			
			@Override
			public void write(final byte[] b, final int off, final int len) throws IOException {
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
						IOUtils.copyTo(this.tmp,ZipInstance.this.zout );
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
	
	static class FileInstance implements ArchiveFactory
		{
		private final Path baseDir;
		
		FileInstance(final Path baseDir) throws IOException
			{
			this.baseDir=baseDir;
			if(Files.exists(baseDir) && !Files.isDirectory(baseDir))
				{
				throw new IOException("Not a directory:"+baseDir);
				}				
			}
		
		@Override
		public OutputStream openOuputStream(String filename) throws IOException
			{
			while(filename.startsWith(File.separator)) filename=filename.substring(1);
			final Path f= this.baseDir.resolve(filename);
			if(f.getParent()!=null)
				{
				Files.createDirectories(f.getParent());
				}
			return Files.newOutputStream(f);
			}
		@Override
		public void close() throws IOException
			{
			}
		}
}
