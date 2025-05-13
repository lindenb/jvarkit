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

import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;

import htsjdk.io.HtsPath;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.seekablestream.SeekableBufferedStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.LocationAware;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.BinaryFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import org.apache.commons.compress.compressors.zstandard.ZstdCompressorInputStream;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.MathUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.log.Logger;


public class BGenCodec extends AbstractFeatureCodec<BGenVariant,SeekableStream> {
	private static final Logger LOG = Logger.of(BGenCodec.class).setDebug();

	private enum State {expect_variant_def, expect_genotypes };
	private BGenHeader header;
	private long snps_offset ;
	private State state;
	private VariantImpl lastVariant=null;
	private int layout1_expect_n_samples=-1;
	private final BGenUtils.ByteBuffer buffer1 = new BGenUtils.ByteBuffer();
	/** auxiliary indexed VCF file that will be used as an coordinate index for the VCF */
	private final BinaryCodec binaryCodec = new BinaryCodec();
	private final boolean skip_genotypes;

	/** class counting the number of bytes read. Useful to skip bytes after the header */
	private static class ByteCountInputStream extends FilterInputStream {
		long count=0L;
		ByteCountInputStream(InputStream delegate) {
			super(delegate);
			}
		@Override
		public int read() throws IOException {
			int c= super.read();
			if(c!=-1) count++;
			return c;
			}
		@Override
		public int read(byte[] b, int off, int len) throws IOException {
			int n=super.read(b, off, len);
			if(n!=-1) count+=n;
			return n;
			}
		@Override
		public long skip(long n) throws IOException {
			n= super.skip(n);
			count+=n;
			return n;
			}
		@Override
		public void close() throws IOException {
			//do not close !
			}
		}
	
	
	
	private static abstract class AbstractBenVariant implements BGenVariant {
		protected AbstractBenVariant() {}
		
		@Override
		public int hashCode() {
			return Objects.hash(getContig(), getPosition(), getRsId(), getId(),this.getAlleles());
		}
		
		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (!(obj instanceof BGenVariant))
				return false;
			final BGenVariant other = (BGenVariant) obj;
			return contigsMatch(other) &&
					getPosition() == other.getPosition() &&
					Objects.equals(getAlleles(), other.getAlleles());
			}
		protected void toString(final StringBuilder builder) {
			builder.append("chrom=");
			builder.append(getContig());
			builder.append(", pos=");
			builder.append(getPosition());
			builder.append(", variant_id=");
			builder.append(getId());
			builder.append(", rs_id=");
			builder.append(getRsId());
			builder.append(", alleles=");
			builder.append(String.join(",", getAlleles()));
			}
		
		
		}	
	
	private class VariantImpl extends AbstractBenVariant {
		String contig;
		int position;
		String variantId;
		String rsId;
		List<String> alleles;
		long offset;
		
		@Override
		public int getBitsPerProb() {
			return -1;
			}
		
		@Override
		public String getContig() {
			return contig;
			}
		@Override
		public int getPosition() {
			return position;
			}
		@Override
		public List<String> getAlleles() {
			return alleles;
			}
		@Override
		public String getId() {
			return variantId;
			}
		@Override
		public String getRsId() {
			return rsId;
			}
		@Override
		public long getOffset() {
			return this.offset;
			}
		@Override
		public final boolean hasGenotypes() {
			return false;
			}
		@Override
		public final List<BGenGenotype> getGenotypes() {
			throw new IllegalStateException("this variant was read without genotype");
			}
		@Override
		public final BGenGenotype getGenotype(final String sn) {
			throw new IllegalStateException("this variant was read without genotype");
			}
		
		@Override
		public boolean isPhased() {
			throw new IllegalStateException("this variant was read without genotype");
			}
		
		@Override
		public String toString() {
			final StringBuilder builder = new StringBuilder();
			builder.append("VariantImpl [=");
			toString(builder);
			builder.append("]");
			return builder.toString();
			}
		
		
		}
	
	private abstract class AbstractVariantAndGenotypes extends AbstractBenVariant {
		private final VariantImpl delegate;
		
		protected abstract class AbstractGenotypelayoutXX implements BGenGenotype {
			final int sample_index;
			protected AbstractGenotypelayoutXX(int sample_index) {
				this.sample_index=sample_index;
				}
			@Override
			public String getSample() {
				return BGenCodec.this.getHeader().getSamples().get(this.sample_index);
				}
			
			@Override
			public String toString() {
				final StringBuilder sb=new  StringBuilder("Genotype(");
				sb.append("layout=").append(BGenCodec.this.getHeader().getLayout().name()).append(",");
				sb.append("phased=").append(isPhased()).append(",");
				sb.append("ploidy=").append(getPloidy()).append(",");
				sb.append("n-alleles=").append(AbstractVariantAndGenotypes.this.getNAlleles()).append(",");
				sb.append("freq=").append(Arrays.toString(getProbs())).append(",");
				sb.append("sample").append(getSample());
				sb.append(")");
				return sb.toString();
				}
		}
		
		AbstractVariantAndGenotypes(final VariantImpl delegate) {
			this.delegate = delegate;
			}
		@Override
		public String getContig() {
			return this.delegate.getContig();
			}
		@Override
		public int getPosition() {
			return this.delegate.getPosition();
			}
		@Override
		public List<String> getAlleles() {
			return this.delegate.getAlleles();
			}
		@Override
		public String getId() {
			return this.delegate.getId();
			}
		@Override
		public String getRsId() {
			return this.delegate.getRsId();
			}
		@Override
		public long getOffset() {
			return this.delegate.getOffset();
			}
		@Override
		public final boolean hasGenotypes() {
			return true;
			}
		@Override
		public BGenGenotype getGenotype(int i) {
			return getGenotypes().get(i);
			}
		@Override
		public final BGenGenotype getGenotype(final String sn) {
			final int i = BGenCodec.this.getHeader().getSampleIndex(sn);
			if(i==-1) throw new IllegalArgumentException("no such sample "+sn);
			return getGenotype(i);
			}
		
		
		}
	
	private final class VariantAndGenotypesLayout1 extends AbstractVariantAndGenotypes {
		private final List<BGenGenotype> genotypes;
		private final double[] probs;
		private class GenotypeLayout1 extends  AbstractGenotypelayoutXX {
			GenotypeLayout1(int sample_index){
				super(sample_index);
				}
			@Override
			public final int getPloidy() {
				return 2;
				}
			
			@Override
			public boolean isPhased() {
				return false;
				}
			@Override
			public boolean isMissing() {
				return get(0)==0 && get(1)==0 && get(2)==0;
				}
			double get(int idx) {
				return probs[super.sample_index*3+idx];
				}
			@Override
			public double[] getProbs() {
				return new double[] {get(0),get(1),get(2)};
				}
			
			
			}
		
		VariantAndGenotypesLayout1(final VariantImpl delegate,final double[] probs) {
			super(delegate);
			this.probs = probs;
			this.genotypes = new AbstractList<BGenGenotype>() {
				@Override
				public int size() {
					return probs.length;
					}
				@Override
				public BGenGenotype get(int index) {
					return new GenotypeLayout1(index);
					}
				};
			}
		
		@Override
		public final boolean isPhased() {
			return false;
			}
		@Override
		public final int getBitsPerProb() {
			return 16;
			}
		@Override
		public final List<BGenGenotype> getGenotypes() {
			return this.genotypes;
			}

		@Override
		public String toString() {
			final StringBuilder builder = new StringBuilder();
			builder.append("VariantAndGenotypesLayout1 [=");
			toString(builder);
			builder.append("]");
			return builder.toString();
			}
		}
	
	private final class VariantAndGenotypesLayout2 extends AbstractVariantAndGenotypes {
		private final List<BGenGenotype> genotypes;
		private short phased;
		private int nbits;
		private class GenotypeLayout2 extends AbstractGenotypelayoutXX {
			byte metaByte;
			double[] probs;
			GenotypeLayout2(int sample_index){
				super(sample_index);
				}
			
			
			void setProbs(List<Double> v) {
				this.probs = new double[v.size()];
				for(int i=0;i< v.size();i++) {
					this.probs[i] = v.get(i).doubleValue();
				}
			}
			
			/* Missingness is encoded by the most significant bit;  */
			@Override
			public boolean isMissing() {
				return ((metaByte &  0b10000000) >>> 7)==1;
			}
			/*  Ploidy (possible values 0-63) is encoded in the least significant 6 bits of this value. */
			public int getPloidy() {
				return (metaByte &  0b00111111);
				}
			
			@Override
			public boolean isPhased() {
				return VariantAndGenotypesLayout2.this.phased==1;
				}
			@Override
			public double[] getProbs() {
				return this.probs;
				}
			
			}
		
		VariantAndGenotypesLayout2(final VariantImpl delegate,final int n_samples) {
			super(delegate);
			this.genotypes = new ArrayList<>(n_samples);
			for(int i=0;i< n_samples;i++) {
				this.genotypes.add(new GenotypeLayout2(i));
				}
			}
		@Override
		public final boolean isPhased() {
			return this.phased==1;
			}
		@Override
		public int getBitsPerProb() {
			return this.nbits;
			}
		
		private GenotypeLayout2 at(int idx) {
			return GenotypeLayout2.class.cast(this.genotypes.get(idx));
		}
		
		@Override
		public List<BGenGenotype> getGenotypes() {
			return genotypes;
			}
		
		
		
		@Override
		public String toString() {
			final StringBuilder builder = new StringBuilder();
			builder.append("VariantAndGenotypesLayout2 [=");
			toString(builder);
			builder.append("]");
			return builder.toString();
			}
		}
	
	
	@Override
	public FeatureCodecHeader readHeader(SeekableStream source) throws IOException {
		this.binaryCodec.setInputStream(source);
		try(ByteCountInputStream bc=new ByteCountInputStream(this.binaryCodec.getInputStream())) {
			final BinaryCodec codec1=new BinaryCodec(bc);
			this.snps_offset =  codec1.readUInt();
			this.header = BGenHeader.of(codec1);
			final long skip = this.snps_offset + 4 - bc.count;
			if(LOG.isDebug()) {
				LOG.info("skip bytes: this.snps_offset="+this.snps_offset+" + 4 - header.count="+bc.count+"="+skip );
				}
			this.binaryCodec.getInputStream().skipNBytes(skip);
			if(LOG.isDebug() && isSeekableStream()) {
				LOG.info("end reading header and now file offset is "+ftell(source) );
				}
			}
		this.state= State.expect_variant_def;
		return new FeatureCodecHeader(this.header, source.position());
		}
	
	
	
	public BGenCodec(boolean skip_genotypes) throws IOException {
		super(BGenVariant.class);
		this.skip_genotypes  = skip_genotypes;
		}
	
	public BGenCodec() throws IOException {
		this(false);
		}
	
	public BGenHeader getHeader() {
		return header;
		}
	
	public BGenUtils.Layout getLayout() {
		return getHeader().getLayout();
		}
	
	public BGenUtils.Compression getCompression() {
		return getHeader().getCompression();
	}
	
	public long getSnpsOffset() {
		return snps_offset;
		}
	
	/** change the reading pointer offset
	 * 
	 * @param position the physical offset
	 * @throws IOException the inutstream is not an instance of SeekableStream
	 */
	public void fseek(PositionalBufferedStream source , long position)  throws IOException {
		assertState(State.expect_variant_def);
		this.binaryCodec.setInputStream(source);
		if(position<=0) throw new IllegalArgumentException("seek.pos<=0: "+position);
		if(!(isSeekableStream())) {
			throw new IOException("cannot do random-access with this king of input "+source.getClass());
			}
		SeekableStream.class.cast(source).seek(position);
		}
	
	private boolean isSeekableStream() {
		return this.binaryCodec.getInputStream() instanceof SeekableStream;
	}
	
	private long ftell(SeekableStream source) throws IOException {
		if(isSeekableStream()) {
			return SeekableStream.class.cast(source).position();
			}
		return -1L;
		}
	
	private List<String> readNAlleles(final int num_alleles) throws IOException {
	      final String[] alleles = new String[num_alleles];
	      for(int k=0;k< num_alleles;++k) {
	            alleles[k] = BGenUtils.readStringUInt32(binaryCodec);
	            }
	      return Arrays.asList(alleles);
	      }
	
	@Override
	public BGenVariant decodeLoc(final SeekableStream source) throws IOException {
		final BGenVariant ctx=readVariant(source);
		skipGenotypes(source);
		return ctx ;		
		}
	
	@Override
	public BGenVariant decode(final SeekableStream source) throws IOException {
		BGenVariant ctx=readVariant(source);
		if(this.skip_genotypes) {
			skipGenotypes(source);
			}
		else {
			ctx = readGenotypes(source);
			}
		return ctx ;
		}
	
	/** read the next variant or returns null if end of file */
	public BGenVariant readVariant(SeekableStream source) throws IOException {
		assertState(State.expect_variant_def);
		this.binaryCodec.setInputStream(source);
		final VariantImpl ctx = new VariantImpl();
		ctx.offset = ftell(source);
		
		if(header.getLayout().equals(BGenUtils.Layout.LAYOUT_1)) {
			try {
				this.layout1_expect_n_samples = BGenUtils.longToUnsignedInt(this.binaryCodec.readUInt());
				}
			catch(Throwable eof) {
				if(LOG.isDebug()) {
					LOG.debug("not more variant layout 1 "+eof.getClass());
					eof.printStackTrace();
					}
				if(eof instanceof RuntimeEOFException) {
					
					return null;
					}
				//cannot read the n-samples it's EOF for layout 1?
				throw eof;
				}
			 }
		try {
			ctx.variantId = BGenUtils.readStringUInt16(this.binaryCodec);
			}
		catch(Throwable eof) {
			if(LOG.isDebug()) {
				LOG.debug("not more variant "+header.getLayout()+" "+eof.getClass());
				eof.printStackTrace();
				}
			if(!header.getLayout().equals(BGenUtils.Layout.LAYOUT_1)) {
				if(eof instanceof RuntimeEOFException) return null;
				}
			throw eof;
			}
	    ctx.rsId =  BGenUtils.readStringUInt16(this.binaryCodec);
	    ctx.contig =  BGenUtils.readStringUInt16(this.binaryCodec);
	     if(StringUtil.isBlank( ctx.contig)) throw new IOException("empty chromosome");
	    ctx.position = BGenUtils.longToUnsignedInt(this.binaryCodec.readUInt());
	    
	     switch(header.getLayout()) {
	     	case LAYOUT_1: ctx.alleles = readNAlleles(2);break;
	     	case LAYOUT_2: 
	     		final int  num_alleles  = this.binaryCodec.readUShort(); 
	     		ctx.alleles =  readNAlleles(num_alleles);
	     		break;
	     	default: throw new IllegalStateException(header.getLayout().name());
	     }
	    this.state = State.expect_genotypes;
	    this.lastVariant = ctx;
	    return ctx;
		}
	
	public void skipGenotypes(SeekableStream source) throws IOException{
		assertState(State.expect_genotypes);
		this.binaryCodec.setInputStream(source);
		final long skip_n =  this.binaryCodec.readUInt();
		if(LOG.isDebug()) {
			LOG.debug("skipping "+skip_n+" bytes");
			}
		this.binaryCodec.getInputStream().skipNBytes(skip_n);
	    this.state = State.expect_variant_def;
		}
	

	
	@SuppressWarnings("resource")
	public BGenVariant readGenotypes(SeekableStream source)  throws IOException{
		assertState(State.expect_genotypes);
		this.binaryCodec.setInputStream(source);

		
		if(LOG.isDebug()) {
			LOG.debug("reading genotypes. offset= "+ftell(source));
			}
		
		int genotype_block_size = BGenUtils.longToUnsignedInt( this.binaryCodec.readUInt());
		if(LOG.isDebug()) {
			LOG.debug("expect compressed size "+genotype_block_size);
			}
		
		/* The total length D of the probability data after uncompression.
		 * If CompressedSNPBlocks = 0, this field is omitted and the total length of the probability data */
		if(!getHeader().getCompression().equals(BGenUtils.Compression.NONE)) {
			final int uncompressed_data_size_D = BGenUtils.longToUnsignedInt( this.binaryCodec.readUInt());
			if(LOG.isDebug()) {
				LOG.debug("expect uncompressed size "+uncompressed_data_size_D+" now offset="+ftell(source));
				}
			genotype_block_size-= Integer.BYTES;
			}
		
		this.buffer1.reset();
		this.buffer1.copyTo(this.binaryCodec.getInputStream(),genotype_block_size);
		if(this.buffer1.size()!=genotype_block_size) {
			throw new IllegalStateException("bad size");
			}
		if(LOG.isDebug()) {
			LOG.debug("after reading genotypes. offset is now = "+ftell(source));
			}
		InputStream uncompressedStream;
		
		
		
		switch(getHeader().getCompression()) {
			case NONE: uncompressedStream = this.buffer1.toByteArrayInputStream();break;
			case ZLIB: uncompressedStream = zlibDecompress(this.buffer1.toByteArrayInputStream()); break;
			case ZSTD: uncompressedStream = zstdDecompress(this.buffer1.toByteArrayInputStream()); break;
	     	default: throw new IllegalStateException(header.getCompression().name());
			}
		final BinaryCodec codec2=new BinaryCodec(uncompressedStream);
		AbstractVariantAndGenotypes vcg;
		switch(getHeader().getLayout()) {
			case LAYOUT_1: vcg=decodeLayout1(codec2); break;
			case LAYOUT_2: vcg=decodeLayout2(codec2); break;
	     	default: throw new IllegalStateException(header.getLayout().name());
			}
		uncompressedStream.close();
		codec2.close();
		
		if(LOG.isDebug()) {
			LOG.debug("End reading genotypes. offset= "+ftell(source));
			}
		
		this.state = State.expect_variant_def;
		return vcg;
		}
	
	private InputStream zlibDecompress(InputStream in) throws IOException {
		return new java.util.zip.InflaterInputStream(in);
		}
	private InputStream zstdDecompress(InputStream in)  throws IOException  {
		return new ZstdCompressorInputStream(in);
		}
	
	private void assertState(State st) {
		if(!this.state.equals(st)) throw new IllegalStateException("expected "+st+" but got "+this.state);
		}
	
	
	
	private VariantAndGenotypesLayout1 decodeLayout1(final BinaryCodec bc) throws IOException {
		final double[] genotypes =new double[ this.layout1_expect_n_samples*3];
		for(int i=0;i< this.layout1_expect_n_samples;++i) {
			for(int k=0;k< 3 ;++k) {
				genotypes[i*3+k] = (float)bc.readUShort() / BinaryCodec.MAX_USHORT;
				}
			}
		return new VariantAndGenotypesLayout1(this.lastVariant,genotypes);
		}
	
	private VariantAndGenotypesLayout2 decodeLayout2(final BinaryCodec codec2) throws IOException {
        /* The number of individuals for which probability data is stored. 
         * This must equal N as defined in the header block.  */
		final int n_samples = BGenUtils.longToUnsignedInt( codec2.readUInt());
        if(n_samples != this.getHeader().getNSamples()) throw new IOException("expected "+this.getHeader().getNSamples()+" but got  " +n_samples +")");
        
        /* The number of alleles, encoded as an unsigned 16-bit integer. 
         * This must equal K as defined in the variant identifying data block. */
        final int n_alleles = codec2.readUShort();
      	if(n_alleles != this.lastVariant.getNAlleles()) throw new IOException("expected "+ this.lastVariant.getNAlleles() +" but got  " +n_alleles +")");
       
      	final VariantAndGenotypesLayout2 vcg = new VariantAndGenotypesLayout2(this.lastVariant,n_samples);
      	/* 1 	The minimum ploidy Pmin of samples in the row. Values between 0 and 63 are allowed. */
      	final int min_ploidy  = (int)codec2.readUByte();
        if(min_ploidy  <0 || min_ploidy >63) throw new IOException("bad ctx.min_ploidy  ("+ min_ploidy +")");
        
        /* The maximum ploidy Pmax of samples in the row. Values between 0 and 63 are allowed. */
        final int max_ploidy  = (int)codec2.readUByte();
        if(max_ploidy  <0 || max_ploidy >63) throw new IOException("bad ctx.max_ploidy  ("+ max_ploidy +")");
        
       // https://github.com/GenomicsNX/systemsgenetics/blob/db00ede6b79bad79f8a5638a45d2ca23c6f23eb4/Genotype-IO/src/main/java/org/molgenis/genotype/bgen/BgenGenotypeData.java#L686C4-L686C34
        
        
        /* A list of N bytes, where the nth byte is an unsigned integer representing the ploidy and
         *  missingness of the nth sample. */
        for(int x=0;x < n_samples ;x++) {
        	vcg.at(x).metaByte =  codec2.readByte();
        	}
        /* Flag, denoted Phased indicating what is stored in the row. */
        vcg.phased = codec2.readUByte();
	      
        /* Unsigned integer B representing the number of bits used to store each probability in this row. This must be between 1 and 32 inclusive. */
       	vcg.nbits=  (int) codec2.readUByte();
        if(vcg.nbits  <0 || vcg.nbits > 32) throw new IOException("bad nbits  ("+ vcg.nbits +")");
        if(LOG.isDebug()) {
   			LOG.debug("nbits="+vcg.nbits);
   			}
        
        final BGenUtils.BitReader bitReader = new BGenUtils.BitReader(codec2.getInputStream());
        final BGenUtils.BitNumReader numReader = new BGenUtils.BitNumReader(bitReader, vcg.nbits);
        final List<Double> values= new ArrayList<>();
        if(vcg.phased==1) {
        	if(LOG.isDebug()) {
       			LOG.debug("phased==1");
       			}
    	   /* read haplotypes */
       		for(int i=0;i< n_samples;++i) {
       			values.clear();
       			VariantAndGenotypesLayout2.GenotypeLayout2 gt= vcg.at(i);
       			if(LOG.isDebug()) {
           			LOG.debug("read phased genotype["+i+"] ploidy :"+gt.getPloidy()+" missing:"+gt.isMissing()+" n-alleles:"+this.lastVariant.getNAlleles());
           			new BGenUtils.AlleleCombinations(
           					BGenUtils.Layout.LAYOUT_2,
           					this.lastVariant.getAlleles(), 
           					gt.getPloidy(),
           					true
           					).debug();
       				}
       			for(int py=0;py < gt.getPloidy();++py) {
       				double sum = 0;
	       			for(int ax=0;ax < this.lastVariant.getNAlleles() -1 ;++ax) {
       					double v = numReader.next();
       					if(LOG.isDebug()) {
       						LOG.debug("next value is "+v);
       						}
       					if(gt.isMissing()) v=0;
       					values.add(v);
       					sum += v;
       					}
	       			//set value for the last allele
	       			if(gt.isMissing()) {
	       				values.add(0.0);
	       				}
	       			else
	       				{
	       				if(sum>1.0) throw new IllegalStateException("sum>1");
	       				values.add(1.0-sum);
	       				}
	       			gt.setProbs(values);
       				}
       			}
        	
        	
       	} else  if(vcg.phased==0) {
       		if(LOG.isDebug()) {
       			LOG.debug("phased==0");
       			}
       		/* read genotypes */
       		for(int i=0;i< n_samples;++i) {
       			values.clear();
       			
       			final VariantAndGenotypesLayout2.GenotypeLayout2 gt= vcg.at(i);
       			if(LOG.isDebug()) {
           			LOG.debug("read unphased genotype["+i+"] ploidy :"+gt.getPloidy()+" missing:"+gt.isMissing()+" n-alleles:"+this.lastVariant.getNAlleles());
           			new BGenUtils.AlleleCombinations(
           					BGenUtils.Layout.LAYOUT_2,
           					this.lastVariant.getAlleles(), 
           					gt.getPloidy(),
           					false
           					).debug();
       				}
       			final BigInteger bi=MathUtils.combination( gt.getPloidy()+ this.lastVariant.getNAlleles()-1,this.lastVariant.getNAlleles()-1);// whyyyy ?
       			if(LOG.isDebug()) {
           			LOG.debug("combination="+bi);
           			}
       			final int n_combinaison_as_int= bi.intValueExact();
       			if(n_combinaison_as_int!= BGenUtils.calculateTotalCombinationsForLayout2(
       						false,this.lastVariant.getNAlleles(), gt.getPloidy()) 
       					) {
       				throw new IllegalStateException();
       				}
       			double sum = 0;
       			for(int j=0;j < n_combinaison_as_int -1 ;++j) {
   					double v = numReader.next();
   					if(LOG.isDebug()) {
   						LOG.debug("next value is "+v);
   						}
   					if(gt.isMissing()) v=0;
   					values.add(v);
   					sum += v;
   					}
       			//set value for the last allele
       			if(gt.isMissing()) {
       				values.add(0.0);
       				}
       			else
       				{
       				if(sum>1.0) throw new IllegalStateException("sum>1");
       				values.add(1.0-sum);
       				}
       			gt.setProbs(values);
       			}
       	} else
       	{
       		throw new IOException("bad phased  ("+ vcg.phased +")");
       	}
        
        return vcg;
		}
	
	/*
	private static List<int[]> buildPhasedIndexes(int n_alleles,int n_ploidy) {
		final List<int[]> container = new ArrayList<>(n_alleles*n_ploidy);
		buildPhasedIndexes(n_alleles,container,new int[n_ploidy],0);
		return container;
		}
	private static void buildPhasedIndexes(
			int n_alleles,
			List<int[]> container,
			int[] buffer,
			int ploidy_index
			) {
			if(ploidy_index==buffer.length) {
				container.add(Arrays.copyOf(buffer, buffer.length));
				}
			else
				{
				for(int a=0;a<n_alleles;++a) {
					int[] buffer2 =Arrays.copyOf(buffer, buffer.length);
					buffer2[ploidy_index]=a;
					buildPhasedIndexes(n_alleles,container,buffer2,ploidy_index+1);
					}
				}
			}*/
	@Override
	public void close(SeekableStream source) {
		try {source.close();}
		catch(Throwable err) {}
		this.binaryCodec.close();
		}
	@Override
	public boolean canDecode(final String path) {
		if(path==null || !path.endsWith(BGenUtils.FILE_SUFFIX)) return false;
		SeekableStream in = null;
		try {
			in = SeekableStreamFactory.getInstance().getBufferedStream(SeekableStreamFactory.getInstance().getStreamFor(Objects.requireNonNull(path,"null input")),Defaults.NON_ZERO_BUFFER_SIZE);
			try(SeekableStream bi=makeSourceFromStream(in)) {
				readHeader(bi);
				}
			return true;
			}
		catch(Throwable err) {
			return false;
			}
		finally {
			if(in!=null) try {
				in.close();
				}
			catch(IOException err) {
				LOG.warn(err);
				}
			}
		}
	
	@Override
    public SeekableStream makeSourceFromStream(final InputStream is) {
        if (is instanceof SeekableBufferedStream)
            return SeekableBufferedStream.class.cast(is);
        else if (is instanceof SeekableStream)
        	return SeekableStreamFactory.getInstance().getBufferedStream(SeekableStream.class.cast(is));
        throw new IllegalArgumentException();
    }

	@Override
	public boolean isDone(SeekableStream source) {
		try {
			return source.eof();
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
    /** {@link PositionalBufferedStream} is already {@link LocationAware}. */
    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        throw new IllegalArgumentException();
    }
	
	
	@Override
	public String toString() {
		return "BGenCodec()";
		}
}
