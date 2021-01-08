/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream;


public interface ArchiveFactory extends Closeable{

	public static final String OPT_DESC="An existing directory or a filename ending with the '.zip' or '.tar' or '.tar.gz' suffix.";
	
	/** return true is archive is instance of ZipInstance or Tar*/
	public boolean isTarOrZipArchive();
	
	public abstract OutputStream openOuputStream(final String filename) throws IOException;
	
	/** copy whole file into the archive */
	public default void copyTo(final Path externalFile,final String filename) throws IOException {
		try ( OutputStream os = openOuputStream(filename))
			{
			IOUtils.copyTo(externalFile, os);
			os.flush();
			}
		}
	
	/** open a writer to this archive */
	public default PrintWriter openWriter(final String filename) throws IOException
		{
		return new PrintWriter(openOuputStream(filename), true);
		}

	/** set compression level for zip archives */
	public void setCompressionLevel(int level);
	 
	/** open a new ArchiveFactory, if filename ends with '.zip' it will be a zip instance
	 * otherwise it will be a FileInstance */
	public static ArchiveFactory open(final Path f)  throws IOException
		{
		if( f == null ) throw new IllegalArgumentException("Cannot open(null)");
		final String fn = f.getFileName().toString().toLowerCase();
		if(fn.endsWith(".zip"))
			{
			return new ZipInstance(f);
			}
		else if(fn.endsWith(".tar") || fn.endsWith(".tar.gz") )
			{
			return new TarInstance(f);
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
	
	
	static abstract class AbstractOutputStream extends OutputStream
		{
		protected OutputStream out;
		
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
		public final boolean isTarOrZipArchive() { return true;}

		@Override
		public void setCompressionLevel(int level) {
			this.zout.setLevel(Math.max(Deflater.NO_COMPRESSION, Math.min(Deflater.BEST_COMPRESSION, level)));
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
		
		private class ZipOS extends AbstractOutputStream
			{
			ZipEntry ze;
			File tmp;
			
			ZipOS(String filename) throws IOException
				{
				while(filename.startsWith("/")) filename=filename.substring(1);
				this.ze=new ZipEntry(filename);
				
				this.tmp=File.createTempFile("tmp", ".zipentry");
				this.tmp.deleteOnExit();
				super.out=new FileOutputStream(this.tmp);
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
	
	static class TarInstance
	implements ArchiveFactory
		{
		OutputStream fout;
		GZIPOutputStream gzout = null;
		TarArchiveOutputStream tarout;
		
		TarInstance(final Path f) throws IOException
			{
			fout= Files.newOutputStream(f);
			if(f.getFileName().toString().toLowerCase().endsWith(".gz")) {
				gzout = new GZIPOutputStream(fout);
				tarout=new TarArchiveOutputStream(gzout);
				}
			else
				{
				tarout=new TarArchiveOutputStream(fout);
				}
			}
		
		@Override
		public boolean isTarOrZipArchive() {
			return true;
			}
		@Override
		public void setCompressionLevel(int level) {
			//ignore
			}
		
		
		@Override
		public OutputStream openOuputStream(final String filename) throws IOException
			{
			final TarOs os =new TarOs(filename);
			return os;
			}
		
		@Override
		public void close() throws IOException
			{
			this.tarout.close();
			if(this.gzout!=null)
				{
				this.gzout.finish();
				this.gzout.flush();
				}
			if(this.fout!=null) {
				this.fout.flush();
				this.fout.close();
				}
			}
		
		private class TarOs extends AbstractOutputStream
			{
			final String entryName;
			File tmp;
			
			TarOs(String filename) throws IOException
				{
				while(filename.startsWith("/")) filename=filename.substring(1);
				this.entryName = filename;
				this.tmp=File.createTempFile("tmp", ".tarentry");
				this.tmp.deleteOnExit();
				super.out=new FileOutputStream(this.tmp);
				}
			
			
			@Override
			public void close() throws IOException
				{
				if(out!=null)
					{
					out.flush();
					out.close();

					if(TarInstance.this.tarout!=null)
						{
						final TarArchiveEntry tarEntry  = new TarArchiveEntry(this.tmp, this.entryName);
						TarInstance.this.tarout.putArchiveEntry(tarEntry);
												
						IOUtils.copyTo(this.tmp,TarInstance.this.tarout );
						TarInstance.this.tarout.closeArchiveEntry();
						}
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
		public void setCompressionLevel(int level) {
			// do nothing
			}
		
		@Override
		public boolean isTarOrZipArchive() {
			return false;
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
