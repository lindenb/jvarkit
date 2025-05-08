/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.bgen;

import java.io.FileOutputStream;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.BitSet;
import java.util.List;
import java.util.zip.DeflaterOutputStream;

import org.apache.commons.compress.compressors.zstandard.ZstdCompressorOutputStream;
import org.apache.commons.compress.utils.IOUtils;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.RuntimeIOException;


public class BGenWriter extends BGenUtils implements AutoCloseable {
	private static final Logger LOG = Logger.build(BGenWriter.class).make();

	private enum State {expect_write_header,expect_variant,expect_genotypes};
	private final BinaryCodec binaryCodec;
	private final OutputStream out;
	private State state= State.expect_write_header;
	private boolean debug_flag=true;
	private Compression compression=Compression.e_ZlibCompression;
	private Layout layout = Layout.e_Layout2;
	private ByteBuffer buffer1 = new ByteBuffer();
	private ByteBuffer buffer2 = new ByteBuffer();
	private long n_variants=0L;
	private long layout1_expect_n_samples;
	/** class counting the number of bytes written. Useful to skip bytes after the header */
	private static class ByteCountOutputStream extends FilterOutputStream {
		long count=0L;
		ByteCountOutputStream(OutputStream delegate) {
			super(delegate);
			}
		@Override
		public void write(int b) throws IOException {
			super.write(b);
			count++;
			}
		@Override
		public void write(byte[] b, int off, int len) throws IOException {
			super.write(b, off, len);
			count+=len;
			}
		
		@Override
		public void close() throws IOException {
			//do not close !
			}
		}
	
	public BGenWriter(Path out) throws IOException{
		this(Files.newOutputStream(out));
		}
	
	public BGenWriter(OutputStream out) {
		this.out = out;
		this.binaryCodec = new BinaryCodec(this.out);
		}
	public void setCompression(Compression compression) {
		assertState(State.expect_write_header);
		this.compression = compression;
		}
	public void setLayout(Layout layout) {
		assertState(State.expect_write_header);
		this.layout = layout;
		}
	
	public void writeHeader(List<String> samples) throws IOException{
		_writeHeader(samples,samples.size());
		}
	public void writeHeader(int n_samples) throws IOException{
		_writeHeader(null,n_samples);
		}
	private void _writeHeader(final List<String> samples,int n_samples) throws IOException {
		assertState(State.expect_write_header);
		buffer1.reset();
		this.layout1_expect_n_samples = n_samples;
		BinaryCodec bc2 = new BinaryCodec(buffer1);
		
		/* num variants , empty for now */
		writeNBytes(bc2, Integer.BYTES);
	      
	      /* num samples */
		bc2.writeUInt(n_samples);

	      
	     /* magic */
		bc2.writeBytes(BGEN_MAGIC);
	     
		//TODO
		writeNBytes(bc2, 00);//TODO
	      
	      
	     BitSet bitSet=new BitSet(4*8);
	     //TODO set flags
	     bc2.writeBytes(bitSet.toByteArray());
	     
	     
	     /* samples */
	     if(samples!=null) {
	    	 bc2.writeUInt(n_samples);
	    	 for(int i=0;i< n_samples;++i) {
	      	     writeStringUInt16(bc2,samples.get(i));
	             }      
	     	}
	     
	    this.binaryCodec.writeUInt(this.buffer1.size()+Integer.BYTES);
	    this.state=State.expect_variant;
		}
	private void assertState(State st) {
		if(!this.state.equals(st)) throw new IllegalStateException("expected "+st+" but got "+this.state);
		}
	
	private boolean isWritingToFile() {
		return this.out instanceof FileOutputStream;
		}
	
	public void writeVariant(String contig,int pos,String variantId,String rsId,List<String> alleles) throws IOException{
		assertState(State.expect_variant);
		if(layout.equals(Layout.e_Layout1)) {
			this.binaryCodec.writeUInt(layout1_expect_n_samples);
			}
		writeStringUInt16(this.binaryCodec, variantId);
		writeStringUInt16(this.binaryCodec, rsId);
		writeStringUInt16(this.binaryCodec, contig);
		this.binaryCodec.writeUInt(pos);
		
		if(layout.equals(Layout.e_Layout1)) {
			if(alleles.size()!=2) throw new IllegalArgumentException("expected 2 alleles");
			}
		else
			{
			this.binaryCodec.writeUShort(alleles.size());
			}
		for(String a : alleles) {
			writeStringUInt32(this.binaryCodec, a);
			}
		this.state=State.expect_genotypes;
		}
	
	public void writeGenotypes(List<double[]> probs) throws IOException {
		assertState(State.expect_variant);
		buffer1.reset();
		buffer2.reset();
		
		byte[] compressed;
		switch(this.compression) {
			case e_NoCompression : compressed = buffer1.toByteArray(); break;
			case e_ZlibCompression:
				try(java.util.zip.DeflaterOutputStream compressor=new DeflaterOutputStream(this.buffer2)) {
					IOUtils.copy(buffer1.toByteArrayInputStream(), compressor);
					}
				compressed = buffer2.toByteArray();
				break;
			case e_ZstdCompression:
				try(ZstdCompressorOutputStream compressor=new ZstdCompressorOutputStream(this.buffer2)) {
					IOUtils.copy(buffer1.toByteArrayInputStream(), compressor);
					}
				compressed = buffer2.toByteArray();
				break;
			default: throw new IllegalStateException();
			}
		binaryCodec.writeBytes(compressed);
		this.state=State.expect_variant;
		}
	
	@Override
	public void close()  {
		if(isWritingToFile()) {
			try {
				FileOutputStream fos=FileOutputStream.class.cast(this.out);
				fos.getChannel().position(122);
				this.binaryCodec.writeUInt(n_variants);
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		try {this.out.flush();}
		catch(IOException err) {}
		try {this.out.close();}
		catch(IOException err) {}
		this.binaryCodec.close();
		}
}
