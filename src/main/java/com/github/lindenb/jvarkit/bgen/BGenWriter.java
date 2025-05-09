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
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
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
	/** output compression */
	private Compression compression=Compression.e_ZlibCompression;
	/** output layout */
	private Layout layout = Layout.e_Layout2;
	/** temporary buffer1 */
	private ByteBuffer buffer1 = new ByteBuffer();
	/** temporary buffer2 */
	private ByteBuffer buffer2 = new ByteBuffer();
	/** number of variants written so far */
	private long n_variants=0L;
	/** genotype to export for each variant */
	private ExportGenotype[] export_genotypes = null;
	/** theorical number of output variants */
	private int theoritical_num_variants=0;
	/** force output to be anonymous */
	private boolean set_anonymous_flag=false;
	/** number of bits for layout2 */
	private int nbits=8;
	
	/** class holding a genotype data for each sample */
	private static class ExportGenotype {
		boolean missing;
		boolean phased;
		int ploidy;
		double[] probs;
		
	}
	
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
	
	public BGenWriter(final Path out) throws IOException{
		this(Files.newOutputStream(out));
		}
	
	public BGenWriter(final OutputStream out) {
		this.out = Objects.requireNonNull(out);
		this.binaryCodec = new BinaryCodec(this.out);
		if(isDebugging()) {
			LOG.debug("the classe of outputstream is "+out.getClass());
			}
		}
	
	private boolean isDebugging() {
		return this.debug_flag;
	}
	
	public void setNBits(int nbits) {
		this.nbits = nbits;
		}
	
	public void setCompression(final Compression compression) {
		assertState(State.expect_write_header);
		this.compression = compression;
		}
	public void setLayout(final Layout layout) {
		assertState(State.expect_write_header);
		this.layout = layout;
		}
	/** set theorical number of variants, only useful when writing to a stream rather to a file */
	public void setNumVariants(int theoritical_num_variants) {
		assertState(State.expect_write_header);
		this.theoritical_num_variants = theoritical_num_variants;
		}
	/** do not export names */
	public void setAnonymousSamples(boolean set_anonymous_flag) {
		assertState(State.expect_write_header);
		this.set_anonymous_flag = set_anonymous_flag;
		}
	
	public void writeHeader(final List<String> samples) throws IOException{
		_writeHeader(samples,samples.size());
		}
	public void writeHeader(final int n_samples) throws IOException{
		_writeHeader(null,n_samples);
		}
	private void _writeHeader(final List<String> samplesOrNull,final int n_samples) throws IOException {
		assertState(State.expect_write_header);
		if(samplesOrNull!=null && samplesOrNull.size()!=n_samples) {
			throw new IllegalArgumentException();
			}
		buffer1.reset();
		
		
		this.export_genotypes = new ExportGenotype[n_samples];
		
		
		try(final BinaryCodec bc2 = new BinaryCodec(this.buffer1)) {	
			/* size of the header, no free area, so it's 5 bytes */
			bc2.writeUInt(5*Integer.BYTES);
			/* num variants  */
			bc2.writeUInt(this.theoritical_num_variants);
		      /* num samples */
			bc2.writeUInt(n_samples);
		     /* magic */
			bc2.writeBytes(BGEN_MAGIC);
		      
		      
		     final BitSet bitSet=new BitSet(Integer.BYTES*8);
		     // bit 0 and 1 are compression flags
		     switch(this.compression) {
		     	case e_NoCompression:break;
		     	case e_ZlibCompression: bitSet.set(0); break;
		     	case e_ZstdCompression: bitSet.set(1); break;
		     	default : throw new IllegalStateException();
		     	}
		     switch(this.layout) {
		     	case e_Layout1: bitSet.set(2); break;
		     	case e_Layout2: bitSet.set(3); break;
		     	default : throw new IllegalStateException();
		     	}
		     if(!(samplesOrNull==null || this.set_anonymous_flag)) {
		    	 bitSet.set(31);
		     	} 
		    // to byteArray doesn't work, as it ignores the non-set bits, so we need padding
		     final byte[] raw = bitSet.toByteArray();
		     final byte[] flags = new byte[Integer.BYTES];
		     System.arraycopy(raw, 0, flags, 0,raw.length);
		     bc2.writeBytes(flags);
		     
		     
		     /* samples */
		     if(!(samplesOrNull==null || this.set_anonymous_flag)) {
		    	 final Set<String> seen = new HashSet<>(n_samples);
		    	 bc2.writeUInt(n_samples);
		    	 for(int i=0;i< n_samples;++i) {
		    		 final String sn = samplesOrNull.get(i);
		    		 if(seen.contains(sn)) {
		    			throw new IllegalArgumentException("duplicate sample in header : "+sn); 
		    		 	}
		    		 seen.add(sn);
		      	     writeStringUInt16(bc2,sn);
		             }      
		     	}
		     
		    final byte[] full_header = this.buffer1.toByteArray();
		    
		    //write the OFFSET value
		    this.binaryCodec.writeUInt(full_header.length);
		    //write the header + samples data
		    this.binaryCodec.writeBytes(full_header);
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
	
	public void writeVariant(final BGenVariant ctx) throws IOException{
		writeVariant(ctx.getContig(),ctx.getPosition(),ctx.getId(),ctx.getRsId(),ctx.getAlleles());
		}
	
	public void writeVariant(String contig,int pos,String variantId,String rsId,List<String> alleles) throws IOException{
		assertState(State.expect_variant);
		if(layout.equals(Layout.e_Layout1)) {
			this.binaryCodec.writeUInt(export_genotypes.length);//number of samples
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
	
	public void setGenotype(int sample_index,boolean phased, int ploidy,boolean missing,double[] probs) {
		assertState(State.expect_genotypes);
		if(export_genotypes[sample_index]!=null) {
			throw new IllegalStateException("Variant index["+sample_index+"] was already set");
			}
		final ExportGenotype gt = new ExportGenotype();
		gt.phased = phased;
		gt.ploidy= ploidy;
		gt.missing=missing;
		gt.probs = Arrays.copyOf(probs, probs.length);
		}
	
	public void writeGenotypes() throws IOException {
		assertState(State.expect_genotypes);
		for(int i=0;i< export_genotypes.length;++i) {
			if(export_genotypes[i]==null) {
				throw new IllegalStateException("Variant index["+i+"] was not set");
			}
		}
		byte[] uncompressed_genotypes;
		switch(this.layout)
			{
			case e_Layout1: uncompressed_genotypes = writeGenotypeLayout1(); break;
			case e_Layout2:uncompressed_genotypes = writeGenotypeLayout2(); break;
			default: throw new IllegalStateException();
			}
				
		
		byte[] compressed;
		switch(this.compression) {
			case e_NoCompression : compressed = uncompressed_genotypes; break;
			case e_ZlibCompression:
				try(java.util.zip.DeflaterOutputStream compressor=new DeflaterOutputStream(this.buffer2)) {
					compressor.write(uncompressed_genotypes);
					compressor.finish();
					compressor.flush();
					}
				compressed = buffer2.toByteArray();
				break;
			case e_ZstdCompression:
				try(ZstdCompressorOutputStream compressor=new ZstdCompressorOutputStream(this.buffer2)) {
					compressor.write(uncompressed_genotypes);
					compressor.flush();
					}
				compressed = buffer2.toByteArray();
				break;
			default: throw new IllegalStateException();
			}
		// The total length C of the compressed genotype probability data for this variant. S
		binaryCodec.writeUInt(compressed.length);
		// Genotype probability data for the SNP for each of the N individuals in the cohort 
		binaryCodec.writeBytes(compressed);
		
		//reset genotypes
		Arrays.fill(this.export_genotypes, null);
		this.state=State.expect_variant;
		}
	
	private byte[] writeGenotypeLayout1() {
		this.buffer1.reset();
		try(BinaryCodec bc2= new BinaryCodec(this.buffer1)) {
			for(int i=0;i< export_genotypes.length;++i) {
				final ExportGenotype gt= export_genotypes[i];
				if(gt==null) throw new IllegalStateException();
				if(gt.ploidy!=2) throw new IllegalArgumentException("Ploidy="+gt.ploidy+" for sample["+i+"] but layout "+layout+" only supports ploidy=2");
				if(gt.missing) {
					//write 0 0 0 (3*2 bytes)
					writeNBytes(bc2, 6,(byte)0);
					}
				else
					{
					if(gt.probs==null) throw new IllegalArgumentException("probs==null for sample["+i+"]");
					if(gt.probs.length!=3)  throw new IllegalArgumentException("probs.lenth="+gt.probs.length +"("+Arrays.toString(gt.probs) +") for sample["+i+"] but layout "+layout+" only supports 3 values AA AB BB");
					double sum=0;
					for(int x=0;x<3;++x) {
						double f=gt.probs[x];
						sum+=f;
						if(f<0 || f>1)  throw new IllegalArgumentException("0<=probs["+x+"]="+f+"<1  in("+Arrays.toString(gt.probs) +") for sample["+i+"]");
						int val=(int)(f*USHRT_MAX);
						bc2.writeUShort(val);
						}
					if(sum>1)  throw new IllegalArgumentException("sum of probs>1  in("+Arrays.toString(gt.probs) +") for sample["+i+"] ");
					}
				}
			}
		return buffer1.toByteArray();
		}
	
	private byte[] writeGenotypeLayout2() {
		this.buffer1.reset();
		int min_ploidy = Arrays.stream(export_genotypes).mapToInt(E->E.ploidy).min().orElse(0);
		int max_ploidy = Arrays.stream(export_genotypes).mapToInt(E->E.ploidy).max().orElse(0);
		
		byte[] meta_bits = new byte[export_genotypes.length];
		for(int i=0;i< export_genotypes.length;++i) {
			
			}
		for(int i=0;i< export_genotypes.length;++i) {
			
			}
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void close()  {
		if(isWritingToFile()) {
			try {
				if(isDebugging()) {
					LOG.debug("rewind to write n_variants:"+n_variants);
					}
				final FileOutputStream fos=FileOutputStream.class.cast(this.out);
				// rewind.
				// first integer is 'offset'
				// second integer is 'header size'
				// then we have the n-variants integer
				fos.getChannel().position(2*Integer.BYTES);
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
