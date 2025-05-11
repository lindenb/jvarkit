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
import java.io.InputStream;
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

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.RuntimeIOException;


public class BGenWriter extends BGenUtils implements AutoCloseable {
	private static final Logger LOG = Logger.of(BGenWriter.class).setDebug();

	private enum State {expect_write_header,expect_variant,expect_genotypes};
	private final BinaryCodec binaryCodec;
	private final OutputStream out;
	private State state= State.expect_write_header;
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
	/** number phased status for layout2 */
	private boolean variant_is_phased;
	/** last variant set when writing */
	private LastVariant lastVariant = null;
	
	private static class LastVariant {
		int n_alleles;
		}
	
	/** class holding a genotype data for each sample */
	private static class ExportGenotype {
		boolean missing;
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
		if(LOG.isDebug()) {
			LOG.debug("the classe of outputstream is "+out.getClass());
			}
		}
	
	
	
	public void setNBits(int nbits) {
		this.nbits = nbits;
		}
	public void setPhased(boolean variant_is_phased) {
		this.variant_is_phased = variant_is_phased;
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
	
	long ftell() throws IOException {
		if(isWritingToFile()) {
			final FileOutputStream fos=FileOutputStream.class.cast(this.out);
			return fos.getChannel().position();
		}
		return -1L;
	}
	
	private void _writeHeader(final List<String> samplesOrNull,final int n_samples) throws IOException {
		if(LOG.isDebug()) {
      	    LOG.debug("Writing BGEN writer");
            }
		
		assertState(State.expect_write_header);
		if(samplesOrNull!=null && samplesOrNull.size()!=n_samples) {
			throw new IllegalArgumentException();
			}
		buffer1.reset();
		
		
		this.export_genotypes = new ExportGenotype[n_samples];
		
		
		try(final BinaryCodec bc2 = new BinaryCodec(this.buffer1)) {	
			/* header block size */
			bc2.writeUInt(5*Integer.BYTES);
			/* num variants  */
			bc2.writeUInt(this.theoritical_num_variants);
		      /* num samples */
			bc2.writeUInt(n_samples);
		     /* magic */
			bc2.writeBytes(BGEN_MAGIC);
		      
		      
		     final BitSet bitSet=new BitSet(Integer.BYTES*8 /* number of bits=32 */);
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
		     bc2.writeBytes(toByteArray(bitSet,Integer.BYTES));
		     if(LOG.isDebug()) {
	      	    LOG.debug("bitSet :"+bitSet);
	            }
		     
		     /* samples */
		     if(!(samplesOrNull==null || this.set_anonymous_flag)) {
		    	 final ByteBuffer samplesBuffer=new ByteBuffer();
		    	 try(BinaryCodec bc3=new BinaryCodec(samplesBuffer)) {
		    		 final Set<String> seen = new HashSet<>(n_samples);
		    		 // number of samples
			    	 bc3.writeUInt(n_samples);
			    	 for(int i=0;i< n_samples;++i) {
			    		 final String sn = samplesOrNull.get(i);
			    		 if(seen.contains(sn)) {
			    			throw new IllegalArgumentException("duplicate sample in header : "+sn); 
			    		 	}
			    		 seen.add(sn);
			      	     writeStringUInt16(bc3,sn);
			             }
			    	 }
		    	 /* An unsigned integer LSI indicating the length in bytes of the sample identifier block. */
		    	 if(LOG.isDebug()) {
		      	    LOG.debug("sample_block_size :"+samplesBuffer.size());
		            }
		    	 bc2.writeUInt(samplesBuffer.size()+Integer.BYTES/* count itself */ );
		    	 /* the sample block itself */
		    	 samplesBuffer.writeTo(bc2.getOutputStream());
		     	}
		     else
		     	{
		    	 if(LOG.isDebug()) {
		    		 LOG.debug("anonymous samples");
		    		}
		     	}
		     
		    final byte[] full_header = this.buffer1.toByteArray();
		    //write the  offset
		    this.binaryCodec.writeUInt(full_header.length+Integer.BYTES);
		    //write the header + samples data
		    this.binaryCodec.writeBytes(full_header);
			}
	     
	    //this.binaryCodec.writeUInt(this.buffer1.size()+Integer.BYTES);
	    if(LOG.isDebug() && isWritingToFile()) {
	    	LOG.debug("end writing header. Now file offset is "+ftell());
	    	}
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
		if(StringUtils.isBlank(contig)) throw new IllegalArgumentException("blank chromosome");
		if(pos<1) throw new IllegalArgumentException("bad position:"+pos);
		
		writeStringUInt16(this.binaryCodec, StringUtils.ifBlank(variantId,""));
		writeStringUInt16(this.binaryCodec, StringUtils.ifBlank(rsId,""));
		writeStringUInt16(this.binaryCodec, contig);
		this.binaryCodec.writeUInt(pos);
		
		if(layout.equals(Layout.e_Layout1)) {
			if(alleles.size()!=2) throw new IllegalArgumentException("expected 2 alleles for layout 1 but got "+String.join(",", alleles));
			}
		else
			{
			this.binaryCodec.writeUShort(alleles.size());
			}
		for(String a : alleles) {
			writeStringUInt32(this.binaryCodec, a);
			}
		this.state=State.expect_genotypes;
		this.lastVariant = new LastVariant();
		this.lastVariant.n_alleles=alleles.size();
		}
	
	public void setGenotype(int sample_index,int ploidy,boolean missing,double[] probs) {
		assertState(State.expect_genotypes);
		if(LOG.isDebug()) {
			LOG.debug("setting genotype["+sample_index+"]");
			}
		if(sample_index<0 || sample_index>=export_genotypes.length) {
			throw new ArrayIndexOutOfBoundsException("0<="+sample_index+"<"+export_genotypes.length);
			}
		if(export_genotypes[sample_index]!=null) {
			throw new IllegalStateException("Variant index["+sample_index+"] was already set");
			}
		if(!missing) {
			final double sum = Arrays.stream(probs).sum();
			if(sum<0 || sum>1) throw new ArrayIndexOutOfBoundsException("sample["+sample_index+"] got sum of probs 0<="+sum+"<=1.0 : "+Arrays.toString(probs));
			}
		final ExportGenotype gt = new ExportGenotype();
		gt.ploidy= ploidy;
		gt.missing=missing;
		gt.probs = Arrays.copyOf(probs, probs.length);
		export_genotypes[sample_index]=gt;
		}
	
	public void writeGenotypes() throws IOException {
		assertState(State.expect_genotypes);
		for(int i=0;i< export_genotypes.length;++i) {
			if(export_genotypes[i]==null) {
				throw new IllegalStateException("Variant index["+i+"] was not set");
			}
		}
		buffer1.reset();
		buffer2.reset();
		switch(this.layout)
			{
			case e_Layout1: writeGenotypeLayout1(); break;
			case e_Layout2: writeGenotypeLayout2(); break;
			default: throw new IllegalStateException();
			}
				
		
		ByteBuffer compressed;
		switch(this.compression) {
			case e_NoCompression : compressed = buffer1; break;
			case e_ZlibCompression:
				try(java.util.zip.DeflaterOutputStream compressor=new DeflaterOutputStream(this.buffer2)) {
					try(InputStream raw_in=buffer1.toByteArrayInputStream()) {
						IOUtils.copy(raw_in, compressor);
						}
					compressor.finish();
					compressor.flush();
					}
				compressed = buffer2;
				break;
			case e_ZstdCompression:
				try(ZstdCompressorOutputStream compressor=new ZstdCompressorOutputStream(this.buffer2)) {
					try(InputStream raw_in=buffer1.toByteArrayInputStream()) {
						IOUtils.copy(raw_in, compressor);
						}
					compressor.flush();
					}
				compressed = buffer2;
				break;
			default: throw new IllegalStateException();
			}
		
		// The total length C of the compressed genotype probability data for this variant. S
		binaryCodec.writeUInt(compressed.size());
		// Genotype probability data for the SNP for each of the N individuals in the cohort 
		try(InputStream compressed_in=compressed.toByteArrayInputStream()) {
			IOUtils.copy(compressed_in, binaryCodec.getOutputStream());
			}
		
		//reset genotypes
		Arrays.fill(this.export_genotypes, null);
		this.state=State.expect_variant;
		}
	
	private void writeGenotypeLayout1() {
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
						int val=(int)(f*BinaryCodec.MAX_USHORT);
						bc2.writeUShort(val);
						}
					if(sum>1)  throw new IllegalArgumentException("sum of probs>1  in("+Arrays.toString(gt.probs) +") for sample["+i+"] ");
					}
				}
			}
		}
	
	private ByteBuffer writeGenotypeLayout2() throws IOException {
		this.buffer1.reset();
		try(BinaryCodec bc2= new BinaryCodec(this.buffer1)) {
			/* The number of individuals for which probability data is stored. */
			bc2.writeUInt(export_genotypes.length);
			/* The number of alleles, encoded as an unsigned 16-bit integer. This must equal K as defined in the variant identifying data block. */
			bc2.writeUShort(this.lastVariant.n_alleles);
			
			/* The minimum ploidy Pmin of samples in the row */
			final int min_ploidy = Arrays.stream(export_genotypes).mapToInt(E->E.ploidy).min().orElse(0);
			if(min_ploidy<0 || min_ploidy>63) throw new IllegalArgumentException("min_ploidy 0<="+min_ploidy+"<=63");
			bc2.writeByte((byte)min_ploidy);
			/* The maximum ploidy Pmax of samples in the row */
			final  int max_ploidy = Arrays.stream(export_genotypes).mapToInt(E->E.ploidy).max().orElse(0);
			if(max_ploidy<0 || max_ploidy>63) throw new IllegalArgumentException("v 0<="+max_ploidy+"<=63");
			bc2.writeByte((byte)max_ploidy);
			
			/* A list of N bytes, where the nth byte is an unsigned integer representing the ploidy and missingness of the nth sample.  */
			byte[] meta_bits = new byte[export_genotypes.length];
			for(int i=0;i< export_genotypes.length;++i) {
				final ExportGenotype gt=export_genotypes[i];
				if(gt.ploidy<1 || gt.ploidy>63) throw new IllegalArgumentException("bad ploidy in gt["+i+"]="+gt.ploidy);
				byte meta=0;
				meta |= (byte) (gt.ploidy & 0b111111); 
				if(gt.missing) meta |=  0b10000000; 
				meta_bits[i]=meta;
				}
			bc2.writeBytes(meta_bits);
			
			/* Flag, denoted Phased indicating what is stored in the row */
			bc2.writeByte(variant_is_phased?1:0);
			
			/* Unsigned integer B representing the number of bits used to store each probabilit */
			if(this.nbits<0 || this.nbits>32) throw new IllegalArgumentException("nbits 0<="+this.nbits+"<=32");
			bc2.writeByte(this.nbits);
			
			
			
			BitWriter bitWriter = new BitWriter(bc2.getOutputStream());
			BitNumWriter numWriter=new BitNumWriter(bitWriter, this.nbits);
			for(int i=0;i< export_genotypes.length;++i) {
				final ExportGenotype gt=export_genotypes[i];
				final int n_expected_probs = calculateTotalCombinationsForLayout2(
						this.variant_is_phased,
						this.lastVariant.n_alleles,
						gt.ploidy
						);
				if(n_expected_probs!=gt.probs.length) {
					throw new IllegalArgumentException("for ploidy="+gt.ploidy+" phased:"+this.variant_is_phased+" n-alleles:"+this.lastVariant.n_alleles+
							" it was expected "+n_expected_probs+" but got "+
							gt.probs.length + "("+Arrays.toString(gt.probs)+")");
					}
				
				for(int j=0;j +1 /* do not store last allele */<n_expected_probs;++j) {
					
					
					
					if(gt.missing) {
						numWriter.write(0);
						}
					else
						{
						numWriter.write(gt.probs[i]);
						}
					}
				}
			bitWriter.flush();
			}
		return buffer1;
		}
	
	@Override
	public void close()  {
		if(isWritingToFile()) {
			try {
				if(LOG.isDebug()) {
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
