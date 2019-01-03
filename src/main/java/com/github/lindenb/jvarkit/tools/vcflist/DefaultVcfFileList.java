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
package com.github.lindenb.jvarkit.tools.vcflist;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.AbstractList;
import java.util.Arrays;

import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/** default implementation of a VcfList */
class DefaultVcfFileList extends AbstractList<VariantContext>
	implements VcfList
	{
	private static final Logger LOG=Logger.build(DefaultVcfFileList.class).make();

	private final File vcfFile;
	private final VCFHeader header;
	private RandomAccessFile indexio;
	private final BlockCompressedInputStream bgzfin;
	private final RandomAccessFile vcfrandom;
	private final VCFCodec codec = new VCFCodec();
	private final int _size;
	private int last_list_index = -1;
	
	DefaultVcfFileList(final File vcf) throws IOException {
		this(vcf,VcfOffsetsIndexFactory.getDefaultIndexFile(vcf));
		}
	
	DefaultVcfFileList(final File vcf,final File indexFile) throws IOException {
		this.vcfFile = vcf;
		IOUtil.assertFileIsReadable(indexFile);
		IOUtil.assertFileIsReadable(this.vcfFile);
		if(indexFile.lastModified()< this.vcfFile.lastModified()) {
			LOG.warn("index "+indexFile+" is older than vcf file "+this.vcfFile);
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
		if(fileLength< VcfOffsetsIndexFactory.MAGIC.length) {
			close();
			throw new IOException("index file doesn't contain magic header " + indexFile);
			}
		fileLength-= VcfOffsetsIndexFactory.MAGIC.length;
		if(fileLength%Long.BYTES!=0) {
			close();
			throw new IOException("bad index file  " + indexFile);
			}
		this.indexio = new RandomAccessFile(indexFile, "r");
		final byte magic[]=new byte[VcfOffsetsIndexFactory.MAGIC.length];
		this.indexio.readFully(magic);
		if(!Arrays.equals(magic, VcfOffsetsIndexFactory.MAGIC)) {
			close();
			throw new IOException("bad index file (magic)  " + indexFile);
			}
		
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
			final String line;
			if(this.last_list_index==-1 || this.last_list_index+1!=index)
				{
				this.indexio.seek((long)VcfOffsetsIndexFactory.MAGIC.length+ (long)index*(long)Long.BYTES);
				final long offset = this.indexio.readLong();
				
				if(this.bgzfin!=null) {
					this.bgzfin.seek(offset);
					line = this.bgzfin.readLine();
					}
				else
					{
					this.vcfrandom.seek(offset);
					line = this.vcfrandom.readLine();
					}
				}
			else
				{
				if(this.bgzfin!=null) {
					line = this.bgzfin.readLine();
					}
				else
					{
					line = this.vcfrandom.readLine();
					}
				}
			this.last_list_index = index;
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
	@Override
	public String toString() {
		return "VcfList: "+this.vcfFile;
		}
	}
