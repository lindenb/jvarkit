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


*/
package com.github.lindenb.jvarkit.tools.vcflist;

import java.io.Closeable;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.AbstractList;
import java.util.List;

import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Interface used to access a VCF by offset. Static method are used to create a VCF index if needed.
 * @author lindenb
 *
 */
public interface VcfList extends List<VariantContext>,Closeable {
	public static final String INDEX_EXTENSION =".offsets";
	static final byte MAGIC[]="vcfindex.0.1".getBytes();

	public VCFHeader getHeader();

static class VcfFileList extends AbstractList<VariantContext>
	implements VcfList
	{
	private final File vcfFile;
	private final VCFHeader header;
	private RandomAccessFile indexio;
	private final BlockCompressedInputStream bgzfin;
	private final RandomAccessFile vcfrandom;
	private final VCFCodec codec = new VCFCodec();
	private final int _size;
	private VcfFileList(final File vcf) throws IOException {
		this.vcfFile = vcf;
		final File indexFile=new File(this.vcfFile.getParentFile(),this.vcfFile.getName()+INDEX_EXTENSION);
		IOUtil.assertFileIsReadable(indexFile);
		IOUtil.assertFileIsReadable(this.vcfFile);
		if(indexFile.lastModified()< this.vcfFile.lastModified()) {
			System.err.println("[WARNING] index "+indexFile+" is older than vcf file "+this.vcfFile);
			}
		try (final VCFFileReader r=new VCFFileReader(this.vcfFile, false)){
			this.header = r.getFileHeader();
			}
		this.codec.readHeader(VCFUtils.convertVCFHeaderToLineIterator(header));
		if(vcf.getName().endsWith(".gz"))
			{
			this.bgzfin = new  BlockCompressedInputStream(vcf);
			this.vcfrandom = null;
			}
		else
			{
			this.vcfrandom = new RandomAccessFile(vcf, "r");
			this.bgzfin = null;
			}
		long fileLength = indexFile.length();
		if(fileLength< MAGIC.length) {
			close();
			throw new IOException("index file doesn't contain magic header " + indexFile);
			}
		fileLength-=MAGIC.length;
		if(fileLength%Long.BYTES!=0) {
			close();
			throw new IOException("bad index file  " + indexFile);
			}
		this.indexio = new RandomAccessFile(indexFile, "r");
		this._size=(int)(fileLength/Long.BYTES);
		}
	@Override
	public VCFHeader getHeader() {
		return this.header;
		}
	
	@Override
	public VariantContext get(final int index) {
		if(index<0 || index>=this.size()) throw new IndexOutOfBoundsException("0<"+index+"<"+size() +" in "+vcfFile);
		try {
			this.indexio.seek((long)MAGIC.length+ (long)index*(long)Long.BYTES);
			long offset = this.indexio.readLong();
			final String line;
			if(this.bgzfin!=null) {
				this.bgzfin.seek(offset);
				line = this.bgzfin.readLine();
				}
			else
				{
				this.vcfrandom.seek(offset);
				line = this.vcfrandom.readLine();
				}
			return this.codec.decode(line);
			}
		catch(final IOException err)
			{
			throw new RuntimeIOException(err);
			}
		}
	@Override
	public int size() {
		return this._size;
		}
	@Override
	public void close() throws IOException {
		CloserUtil.close(this.bgzfin);
		CloserUtil.close(this.vcfrandom);
		CloserUtil.close(this.indexio);
		}
	}

public static File getIndexFile(final File vcf) {
	return new File(vcf.getParentFile(),vcf.getName()+INDEX_EXTENSION);
}

public static VcfList open(final File vcfFile) throws IOException
	{
	return new VcfFileList(vcfFile);
	}

public static File indexVcfFile(final File vcfFile) throws IOException
	{
	IOUtil.assertFileIsReadable(vcfFile);
	final File indexFile = VcfList.getIndexFile(vcfFile);
	DataOutputStream daos = null;
	BlockCompressedInputStream bgzin = null;
	AsciiLineReader ascii = null;

	try {
		daos = new DataOutputStream(new FileOutputStream(indexFile));
		daos.write(MAGIC);
		if(vcfFile.getName().endsWith(".gz")) {
			bgzin = new BlockCompressedInputStream(vcfFile);
			ascii = null;			
			}
		else
			{
			bgzin  = null;
			ascii  = new AsciiLineReader(new FileInputStream(vcfFile));
			}
		for(;;)
			{
			final long offset = (ascii==null?bgzin.getPosition():ascii.getPosition());
			final String line = (ascii==null?bgzin.readLine():ascii.readLine());
			if(line==null) break;
			if(line.startsWith("#")) continue;
			daos.writeLong(offset);
			}
		daos.flush();
		daos.close();
		return indexFile;
		}
	catch(final IOException err)
		{
		indexFile.delete();
		throw err;
		}
	finally
		{
		CloserUtil.close(ascii);
		CloserUtil.close(bgzin);
		CloserUtil.close(daos);
		}
	}

}
