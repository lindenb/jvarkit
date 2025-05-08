package com.github.lindenb.jvarkit.bgen;

import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.nio.file.Path;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

import org.apache.commons.compress.compressors.zstandard.ZstdCompressorInputStream;

import com.github.lindenb.jvarkit.math.MathUtils;
import com.github.lindenb.jvarkit.util.log.Logger;


public class BGenReader extends BGenUtils implements AutoCloseable {
	private static final Logger LOG = Logger.build(BGenReader.class).make();

	private enum State {expect_variant_def, expect_genotypes };
	private final BGenHeader header;
	private final long snps_offset ;
	private final BinaryCodec binaryCodec;
	private final InputStream in;
	private State state;
	private VariantImpl lastVariant=null;
	private int layout1_expect_n_samples=-1;
	private final ByteBuffer buffer1 = new ByteBuffer();
	private boolean debug_flag=true;

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
	
	
	public interface Genotype {
		public String getSample();
		public boolean isMissing();
		public boolean isPhased();
		public double[] getProbs();
		public int getPloidy();
		}
	
	
	
	public static abstract class Variant implements Locatable {
		protected Variant() {}
		
		public final int getNAlleles() {
			return this.getAlleles().size();
			}
		public final String getAllele(int idx) {
			return this.getAlleles().get(idx);
			}
		
		
		public abstract int getPosition();
		@Override
		public final int getStart() {
			return getPosition();
			}
		@Override
		public final int getEnd() {
			return getPosition()+ getAlleles().stream().mapToInt(A->A.length()-1).max().orElse(0);
			}
		public abstract List<String> getAlleles();
		public abstract String getId();
		public abstract String getRsId();
		public abstract long getOffset();
		public abstract boolean hasGenotypes();
		public abstract List<? extends Genotype> getGenotypes();
		public Genotype getGenotype(int i) {
			return getGenotypes().get(i);
			}
		public int getNGenotypes() {
			return getGenotypes().size();
			}
		public abstract Genotype getGenotype(String sn);
		
		@Override
		public int hashCode() {
			return Objects.hash(getContig(), getPosition(), getRsId(), getId(),this.getAlleles());
		}
		
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (!(obj instanceof Variant))
				return false;
			final Variant other = (Variant) obj;
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
	
	private class VariantImpl extends Variant {
		String contig;
		int position;
		String variantId;
		String rsId;
		List<String> alleles;
		long offset;
		
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
		public final List<Genotype> getGenotypes() {
			throw new IllegalStateException("this variant was read without genotype");
			}
		@Override
		public final Genotype getGenotype(final String sn) {
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
	
	private abstract class AbstractVariantAndGenotypes extends Variant {
		private final VariantImpl delegate;
		
		protected abstract class AbstractGenotypelayoutXX implements Genotype {
			final int sample_index;
			protected AbstractGenotypelayoutXX(int sample_index) {
				this.sample_index=sample_index;
				}
			@Override
			public String getSample() {
				return BGenReader.this.getHeader().getSamples().get(this.sample_index);
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
		public Genotype getGenotype(int i) {
			return getGenotypes().get(i);
			}
		@Override
		public final Genotype getGenotype(final String sn) {
			int i = BGenReader.this.getHeader().getSampleIndex(sn);
			if(i==-1) throw new IllegalArgumentException("no such sample "+sn);
			return getGenotype(i);
			}
		
		
		}
	
	private final class VariantAndGenotypesLayout1 extends AbstractVariantAndGenotypes {
		private final List<Genotype> genotypes;
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
			this.genotypes = new AbstractList<Genotype>() {
				@Override
				public int size() {
					return probs.length;
					}
				@Override
				public Genotype get(int index) {
					return new GenotypeLayout1(index);
					}
				};
			}
		
		@Override
		public final List<? extends Genotype> getGenotypes() {
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
		private final List<GenotypeLayout2> genotypes;
		private short phased;
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
		public List<? extends Genotype> getGenotypes() {
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
	
	
	public BGenReader(final InputStream in) throws IOException {
		this.in =  Objects.requireNonNull(in,"null input");
		try(ByteCountInputStream bc=new ByteCountInputStream(in)) {
			final BinaryCodec codec1=new BinaryCodec(bc);
			this.snps_offset =  codec1.readUInt();
			this.header = BGenHeader.of(codec1);
			final long skip = this.snps_offset + 4 - bc.count;
			in.skipNBytes(skip);
			}
		this.binaryCodec=new BinaryCodec(this.in);

		this.state= State.expect_variant_def;
		}
	
	public BGenReader(final String path) throws IOException {
		this(SeekableStreamFactory.getInstance().getBufferedStream(SeekableStreamFactory.getInstance().getStreamFor(Objects.requireNonNull(path,"null input")),Defaults.NON_ZERO_BUFFER_SIZE));
		}
	
	public BGenReader(final Path path) throws IOException {
		this(Objects.requireNonNull(path,"null input").toString());
		}
	
	public BGenHeader getHeader() {
		return header;
		}
	
	public long getSnpsOffset() {
		return snps_offset;
		}
	
	private long ftell() throws IOException {
		if(this.in instanceof SeekableStream) {
			return SeekableStream.class.cast(this.in).position();
			}
		return -1L;
		}
	
	private class Iter extends AbstractIterator<Variant> {
		private final boolean skipGenotypes;
		Iter(boolean skipGenotypes) {
			this.skipGenotypes=skipGenotypes;
			}
		
		@Override
		protected Variant advance() {
			try {
				final Variant v= BGenReader.this.readVariant();
				if(v==null) {
					return null;
					}
				else if(this.skipGenotypes) {
					BGenReader.this.skipGenotypes();
					return v;
					}
				else
					{
					return BGenReader.this.readGenotypes();
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		};
	
	public AbstractIterator<Variant> iterator(boolean skipGenotypes) {
		return new Iter(skipGenotypes);
		}
	
	public Iterator<Variant> iterator() {
		return iterator(false);
		}
	
	private List<String> readNAlleles(final int num_alleles) throws IOException {
	      final String[] alleles = new String[num_alleles];
	      for(int k=0;k< num_alleles;++k) {
	            alleles[k] = readStringUInt32(binaryCodec);
	            }
	      return Arrays.asList(alleles);
	      } 
	
	/** read the next variant or returns null if end of file */
	public Variant readVariant() throws IOException {
		assertState(State.expect_variant_def);
		
		final VariantImpl ctx = new VariantImpl();
		ctx.offset = ftell();
		
		if(header.getLayout().equals(BGenHeader.Layout.e_Layout1)) {
			try {
				this.layout1_expect_n_samples = longToUnsignedInt(this.binaryCodec.readUInt());
				}
			catch(Throwable eof) {
				if(isDebugging()) {
					LOG.debug("not more variant "+eof.getClass());
					}
				if(eof instanceof RuntimeEOFException) return null;
				//cannot read the n-samples it's EOF for layout 1?
				throw eof;
				}
			 }
		try {
			ctx.variantId = readStringUInt16(this.binaryCodec);
			}
		catch(Throwable eof) {
			if(isDebugging()) {
				LOG.debug("not more variant "+eof.getClass());
				}
			if(!header.getLayout().equals(BGenHeader.Layout.e_Layout1)) {
				if(eof instanceof RuntimeEOFException) return null;
				}
			throw eof;
			}
	    ctx.rsId =  readStringUInt16(this.binaryCodec);
	    ctx.contig =  readStringUInt16(this.binaryCodec);
	     if(StringUtil.isBlank( ctx.contig)) throw new IOException("empty chromosome");
	    ctx.position = longToUnsignedInt(this.binaryCodec.readUInt());
	    
	     switch(header.getLayout()) {
	     	case e_Layout1: ctx.alleles = readNAlleles(2);break;
	     	case e_Layout2: 
	     		final int  num_alleles  = this.binaryCodec.readUShort(); 
	     		ctx.alleles =  readNAlleles(num_alleles);
	     		break;
	     	default: throw new IllegalStateException(header.getLayout().name());
	     }
	    this.state = State.expect_genotypes;
	    this.lastVariant = ctx;
	    return ctx;
		}
	
	public void skipGenotypes() throws IOException{
		assertState(State.expect_genotypes);
		final long skip_n =  this.binaryCodec.readUInt();
		this.binaryCodec.getInputStream().skipNBytes(skip_n);
	    this.state = State.expect_variant_def;
		}
	
	@SuppressWarnings("resource")
	public Variant readGenotypes()  throws IOException{
		assertState(State.expect_genotypes);
		final int compressed_data_size_C =  longToUnsignedInt( this.binaryCodec.readUInt());
		/* The total length D of the probability data after uncompression.
		 * If CompressedSNPBlocks = 0, this field is omitted and the total length of the probability data */
		if(!getHeader().getCompression().equals(BGenHeader.Compression.e_NoCompression)) {
			/* uncompressed_data_size_D */ longToUnsignedInt( this.binaryCodec.readUInt());
			}
		
		this.buffer1.reset();
		this.buffer1.copyTo(this.in,compressed_data_size_C);
		InputStream uncompressedStream;
		
		
		
		switch(getHeader().getCompression()) {
			case e_NoCompression: uncompressedStream = this.buffer1.toByteArrayInputStream();break;
			case e_ZlibCompression: uncompressedStream = zlibDecompress(this.buffer1.toByteArrayInputStream()); break;
			case e_ZstdCompression: uncompressedStream = zstdDecompress(this.buffer1.toByteArrayInputStream()); break;
	     	default: throw new IllegalStateException(header.getCompression().name());
			}
		final BinaryCodec codec2=new BinaryCodec(uncompressedStream);
		AbstractVariantAndGenotypes vcg;
		switch(getHeader().getLayout()) {
			case e_Layout1: vcg=decodeLayout1(codec2); break;
			case e_Layout2: vcg=decodeLayout2(codec2); break;
	     	default: throw new IllegalStateException(header.getLayout().name());
			}
		uncompressedStream.close();
		codec2.close();
		
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
	
	private boolean isDebugging() {
		return this.debug_flag;
	}
	
	private VariantAndGenotypesLayout1 decodeLayout1(final BinaryCodec bc) throws IOException {
		final double[] genotypes =new double[ this.layout1_expect_n_samples*3];
		for(int i=0;i< this.layout1_expect_n_samples;++i) {
			for(int k=0;k< 3 ;++k) {
				genotypes[i*3+k] = (float)bc.readUShort() / BGenUtils.USHRT_MAX;
				}
			}
		return new VariantAndGenotypesLayout1(this.lastVariant,genotypes);
		}
	
	private VariantAndGenotypesLayout2 decodeLayout2(final BinaryCodec codec2) throws IOException {
        /* The number of individuals for which probability data is stored. 
         * This must equal N as defined in the header block.  */
		final int n_samples = longToUnsignedInt( codec2.readUInt());
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
        	vcg.genotypes.get(x).metaByte =  codec2.readByte();
        	}
        /* Flag, denoted Phased indicating what is stored in the row. */
        vcg.phased = codec2.readUByte();
	      
        /* Unsigned integer B representing the number of bits used to store each probability in this row. This must be between 1 and 32 inclusive. */
       	final int nbits=  (int) codec2.readUByte();
        if(nbits  <0 || nbits > 32) throw new IOException("bad nbits  ("+ nbits +")");
        if(isDebugging()) {
   			LOG.debug("nbits="+nbits);
   			}
        
        final BitReader bitReader = new BitReader(codec2.getInputStream());
        final BitNumReader numReader = new BitNumReader(bitReader, nbits);
        final List<Double> values= new ArrayList<>();
        if(vcg.phased==1) {
        	if(isDebugging()) {
       			LOG.debug("phased==1");
       			}
    	   /* read haplotypes */
       		for(int i=0;i< n_samples;++i) {
       			values.clear();
       			VariantAndGenotypesLayout2.GenotypeLayout2 gt= vcg.genotypes.get(i);
       			if(isDebugging()) {
           			LOG.debug("read phased genotype["+i+"] ploidy :"+gt.getPloidy()+" missing:"+gt.isMissing()+" n-alleles:"+this.lastVariant.getNAlleles());
           			}
       			for(int py=0;py < gt.getPloidy();++py) {
       				double sum = 0;
	       			for(int ax=0;ax < this.lastVariant.getNAlleles() -1 ;++ax) {
       					double v = numReader.next();
       					if(isDebugging()) {
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
       		if(isDebugging()) {
       			LOG.debug("phased==0");
       			}
       		/* read genotypes */
       		for(int i=0;i< n_samples;++i) {
       			values.clear();
       			
       			final VariantAndGenotypesLayout2.GenotypeLayout2 gt= vcg.genotypes.get(i);
       			if(isDebugging()) {
           			LOG.debug("read unphased genotype["+i+"] ploidy :"+gt.getPloidy()+" missing:"+gt.isMissing()+" n-alleles:"+this.lastVariant.getNAlleles());
           			}
       			final BigInteger bi=MathUtils.combination( gt.getPloidy()+ this.lastVariant.getNAlleles()-1,this.lastVariant.getNAlleles()-1);// whyyyy ?
       			if(isDebugging()) {
           			LOG.debug("combination="+bi);
           			}
       			final int n_combinaison_as_int= bi.intValueExact();
       			double sum = 0;
       			for(int j=0;j < n_combinaison_as_int -1 ;++j) {
   					double v = numReader.next();
   					if(isDebugging()) {
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
	
	
	@Override
	public void close()  {
		try {this.in.close();}
		catch(IOException err) {}
		this.binaryCodec.close();
		}
}
