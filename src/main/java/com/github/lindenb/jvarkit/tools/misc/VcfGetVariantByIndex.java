/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

public class VcfGetVariantByIndex extends AbstractVcfGetVariantByIndex
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfGetVariantByIndex.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractVcfGetVariantByIndex.AbstractVcfGetVariantByIndexCommand
		{		
		private static final String STANDARD_EXTENSION=".ith";
		private static abstract class IndexFile
		implements Closeable
		{
		protected File vcfFile;
		protected File indexFile;
		protected VCFUtils.CodecAndHeader cah;
		private RandomAccessFile rafile;
		private long count=0L;
		public IndexFile(File vcfFile) throws IOException
			{
			this.vcfFile=vcfFile;
			this.indexFile = new File(
					vcfFile.getParentFile(),
					vcfFile.getName()+STANDARD_EXTENSION
					);
			}
		public VCFHeader getHeader() {
			return cah.header;
			}
		public AbstractVCFCodec getCodec()
			{
			return cah.codec;
			}
		
		public abstract int read() throws IOException;
		
		protected  String readLine() throws IOException
			{
	        StringBuffer buf = new StringBuffer();
	        int c;
	        while ((c = this.read()) >= 0 && c != '\n')
	            buf.append((char) c);
	        if (c < 0) return null;
	        return buf.toString();
			}
		
		public void open() throws IOException
			{
			if(this.indexFile.exists())
				{
				LOG.info("Index exists reading "+this.indexFile);
				openForReading();
				}
			else
				{
				LOG.info("Writing index for "+this.indexFile);
				openForBuilding();
				}
			}
		
		protected abstract void openVcfFile() throws IOException;
		protected abstract long getFilePointer() throws IOException;
		private void openForBuilding() throws IOException
			{
			try
				{
				this.rafile=new RandomAccessFile(this.indexFile,"rw");
				openVcfFile();
				
				String line=null;
				List<String> headerLines=new ArrayList<>();
				while((line=readLine())!=null)
					{
					headerLines.add(line);
					if(line.startsWith("#CHROM"))  break;
					}
				VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(headerLines);
				this.cah=cah;
				this.count=0L;
				long virtualPtr=this.getFilePointer();
				while((line=readLine())!=null)
					{
					this.rafile.writeLong(virtualPtr);
					this.count++;
					virtualPtr=this.getFilePointer();
					}
				}
			catch(Exception err)
				{
				throw new IOException(err);
				}
			}
		
		private void openForReading() throws IOException
			{
			try
				{
				this.rafile=new RandomAccessFile(this.indexFile,"r");
				openVcfFile();
				String line=null;
				List<String> headerLines=new ArrayList<>();
				while((line=readLine())!=null)
					{
					headerLines.add(line);
					if(line.startsWith("#CHROM"))  break;
					}
				VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(headerLines);
				this.cah=cah;
				
				this.count = (this.indexFile.length()/8L);
				}
			catch(Exception err)
				{
				throw new IOException(err);
				}
			}
		
		public long getVirtualPtr(long index) throws IOException
			{
			if(index<0 || index>=this.count) throw new IndexOutOfBoundsException();
			this.rafile.seek(index*8);//8 =sizeof(long)
			return this.rafile.readLong();
			}
		
		public abstract String getLine(long index)   throws IOException;

		
		public long size()
			{
			return count;
			}
		@Override
		public void close() throws IOException {
			CloserUtil.close(this.rafile);
			this.rafile=null;
			this.cah=null;
			this.count=0;
			}
		}
	
	private static class BGZIndexFile extends IndexFile
		{
		private BlockCompressedInputStream bgzin=null;
		BGZIndexFile(File f) throws IOException
			{
			super(f);
			}
		@Override
		protected void openVcfFile() throws IOException {
			this.bgzin = new BlockCompressedInputStream(this.vcfFile);
			}
		
		@Override
		public void close() throws IOException {
			super.close();
			CloserUtil.close(this.bgzin);
			this.bgzin=null;
			}
		@Override
		public  int read() throws IOException
			{
			return this.bgzin.read();
			}
		@Override
		public String getLine(long index)   throws IOException
			{
			long offset = getVirtualPtr(index);
			this.bgzin.seek(offset);
			return readLine();
			}
		@Override
		protected long getFilePointer()  throws IOException {
			return this.bgzin.getFilePointer();
			}

		}
	private static class RandomAccessIndexFile extends IndexFile
		{
		RandomAccessFile vcfraf=null;
		RandomAccessIndexFile(File f) throws IOException
			{
			super(f);
			}
		
		@Override
		protected void openVcfFile() throws IOException {
			this.vcfraf=new RandomAccessFile(this.vcfFile, "r");
			}
		
		@Override
		public void close() throws IOException {
			super.close();
			CloserUtil.close(this.vcfraf);
			this.vcfraf=null;
			}
		@Override
		public  int read() throws IOException
			{
			return this.vcfraf.read();
			}
		@Override
		public String getLine(long index)   throws IOException
			{
			long offset = getVirtualPtr(index);
			this.vcfraf.seek(offset);
			return readLine();
			}
		@Override
		protected long getFilePointer() throws IOException {
			return this.vcfraf.getFilePointer();
			}

		}

	@Override
	public Collection<Throwable> initializeKnime() {
		if(this.fileListOfIndexes==null)
			{
			return wrapException("undefined list of indexes");
			}
		return super.initializeKnime();
		}
	@SuppressWarnings("resource")
	@Override
	public Collection<Throwable> call(String inpuName) throws Exception
		{
			if(inpuName==null)
				{
				return wrapException("Expected only one vcf file on input");
				}
			File vcfFile=new File(inpuName);
			VariantContextWriter w=null;
			IndexFile indexFile=null;
			BufferedReader r=null;
			String line;
			try {
				LOG.info("Opening "+vcfFile);
				
				if(vcfFile.getName().endsWith(".vcf.gz"))
					{
					indexFile = new BGZIndexFile(vcfFile);
					}
				else if(vcfFile.getName().endsWith(".vcf"))
					{
					indexFile = new RandomAccessIndexFile(vcfFile);
					}
				else
					{
					return wrapException("Not a .vcf or .vcf.gz file: "+vcfFile);
					}
	
				
				indexFile.open();
				w = super.openVariantContextWriter();
				w.writeHeader(indexFile.getHeader());
				r=IOUtils.openFileForBufferedReading(fileListOfIndexes);
				while((line=r.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					long ith;
					try {
						ith=Long.parseLong(line);
						} 
					catch (Exception e) {
						return wrapException("Bad index in "+line+" ignoring");
						}
					ith--;//0-based index
					if(ith<0 || ith>=indexFile.size())
						{
						return wrapException("Index out of bound in "+line+" ignoring");
						}
					String varStr = indexFile.getLine(ith);
					w.add(indexFile.getCodec().decode(varStr));
					}
				return RETURN_OK;
				} 
			catch (Exception e)
				{
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(indexFile);
				CloserUtil.close(w);
				CloserUtil.close(r);
				}
			}
	
		
		}
	
	public static void main(String[] args) throws IOException
		{
		new VcfGetVariantByIndex().instanceMain(args);
		}
	}
