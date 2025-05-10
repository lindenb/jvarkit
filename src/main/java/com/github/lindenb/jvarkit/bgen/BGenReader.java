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

import java.io.Closeable;
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
import java.util.function.Consumer;
import java.util.function.ToIntFunction;

import htsjdk.io.HtsPath;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import org.apache.commons.compress.compressors.zstandard.ZstdCompressorInputStream;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.MathUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
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
	/** auxiliary indexed VCF file that will be used as an coordinate index for the VCF */
	private Path filenameOrNull = null;
	private VCFFileReader vcfIndexFile = null;

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
		public boolean equals(Object obj) {
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
				return BGenReader.this.getHeader().getSamples().get(this.sample_index);
				}
			
			protected int[] getAlleleIndexesForMissingGT() {
				final int[] array = new int[getPloidy()];
				Arrays.fill(array, -1);
				return array;
				}
			
			@Override
			public String toString() {
				final StringBuilder sb=new  StringBuilder("Genotype(");
				sb.append("layout=").append(BGenReader.this.getHeader().getLayout().name()).append(",");
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
			final int i = BGenReader.this.getHeader().getSampleIndex(sn);
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
			
			
			@Override
			public int[] getAllelesIndexes(ToIntFunction<double[]> fun) {
				if(isMissing()) return getAlleleIndexesForMissingGT();
				final int i = fun.applyAsInt(getProbs());
				if(i<0)  return getAlleleIndexesForMissingGT();
				switch(i) {
					case 0: return new int[] {0,0};
					case 1: return new int[] {0,1};
					case 2: return new int[] {1,1};
					default: throw new IllegalArgumentException("invalid index "+i);
					}
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
			
			@Override
			public int[] getAllelesIndexes(ToIntFunction<double[]> fun) {
				if(isMissing()) return getAlleleIndexesForMissingGT();
				final int i = fun.applyAsInt(getProbs());
				if(i<0)  return getAlleleIndexesForMissingGT();
				//simple case
				if(!isPhased() && isDiploid() && VariantAndGenotypesLayout2.this.getNAlleles()==2) {
					switch(i) {
						case 0: return new int[] {0,0};
						case 1: return new int[] {0,1};
						case 2: return new int[] {1,1};
						default: throw new IllegalArgumentException("invalid index "+i+"/"+getProbs());
						}
					}
				if(isPhased()) {
					
					}
				return null;
				}
			private int[] findAlleleIndexesPhased(int wanted_index) {
				AlleleCombinations ac= new AlleleCombinations(VariantAndGenotypesLayout2.this.getNAlleles(), getPloidy(), isPhased());
				return ac.findAlleleIndexesByIndex(wanted_index);
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
		if(this.in instanceof SeekableStream) {
			String src=SeekableStream.class.cast(this.in).getSource();
			if(!StringUtils.isBlank(src)) {
				final HtsPath source = new HtsPath(src);
				if(source.isPath()) {
					this.filenameOrNull = source.toPath();
					}
				}
			}
		}
	
	public BGenReader(final Path path) throws IOException {
		this(Objects.requireNonNull(path,"null input").toString());
		this.filenameOrNull = path;
		}
	
	public BGenHeader getHeader() {
		return header;
		}
	
	public long getSnpsOffset() {
		return snps_offset;
		}
	
	/** change the reading pointer offset
	 * 
	 * @param position the physical offset
	 * @throws IOException the inutstream is not an instance of SeekableStream
	 */
	public void fseek(long position)  throws IOException {
		assertState(State.expect_variant_def);;
		if(position<=0) throw new IllegalArgumentException("seek.pos<=0: "+position);
		if(!(this.in instanceof SeekableStream)) {
			throw new IOException("cannot do random-access with this king of input "+this.in.getClass());
			}
		SeekableStream.class.cast(this.in).seek(position);
		}
	
	private long ftell() throws IOException {
		if(this.in instanceof SeekableStream) {
			return SeekableStream.class.cast(this.in).position();
			}
		return -1L;
		}
	
	private class Iter extends AbstractIterator<BGenVariant> {
		private final boolean skipGenotypes;
		Iter(boolean skipGenotypes) {
			this.skipGenotypes=skipGenotypes;
			}
		
		@Override
		protected BGenVariant advance() {
			try {
				final BGenVariant v= BGenReader.this.readVariant();
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
	
	public AbstractIterator<BGenVariant> iterator(boolean skipGenotypes) {
		return new Iter(skipGenotypes);
		}
	
	public Iterator<BGenVariant> iterator() {
		return iterator(false);
		}
	
	private class QueryIterator extends AbstractIterator<BGenVariant> implements CloseableIterator<BGenVariant> {
		private Locatable query;
		private CloseableIterator<VariantContext> delegate;
		private boolean skipGenotypes;
		@Override
		protected BGenVariant advance() {
			try {
				for(;;) {
					if(delegate==null) return null;
					if(!delegate.hasNext()) {
						close();
						return null;
						}
					final VariantContext ctx = this.delegate.next();
					if(!ctx.overlaps(this.query)) continue;
					if(!ctx.hasAttribute(OFFSET_format_header_line.getID())) continue;
					final long offset= Long.parseLong(ctx.getAttributeAsString(OFFSET_format_header_line.getID(), "-1"));
					if(offset<=0L) throw new IllegalArgumentException("bad offset in "+ctx);
					BGenReader.this.fseek(offset);
					BGenVariant bv = BGenReader.this.readVariant();
					if(skipGenotypes) {
						BGenReader.this.skipGenotypes();
						}
					else
						{
						bv = BGenReader.this.readGenotypes();
						}
					return bv;
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close()  {
			if(this.delegate!=null) this.delegate.close();
			this.delegate=null;
			}
		}
	
	public CloseableIterator<BGenVariant> query(final Locatable loc, boolean skipGenotypes) {
		loadVcfIndex();
		final QueryIterator iter=new QueryIterator();
		iter.query=new SimpleInterval(loc);
		iter.delegate = Objects.requireNonNull(this.vcfIndexFile).query(loc);
		iter.skipGenotypes = skipGenotypes;
		return iter;
		}
	
	public boolean hasIndex() {
		if(this.vcfIndexFile!=null) return true;
		if(this.filenameOrNull==null)return false;
		final Path vcfPath = getVcfIndexName(this.filenameOrNull);
		if(!Files.exists(vcfPath)) return false;
		final String vcffilename = vcfPath.toAbsolutePath().toString();
		final  Path tbi = vcfPath.getFileSystem().getPath(ParsingUtils.appendToPath(vcffilename, FileExtensions.TABIX_INDEX));
		if(!Files.exists(tbi)) return false;
		try(VCFFileReader r=new VCFFileReader(vcfPath,false)) {
			final VCFHeader hder = r.getHeader();
			if(!hder.hasInfoLine(OFFSET_format_header_line.getID())) return false;
			}
		catch(Throwable err) {
			return false;
			}
		return true;
		}
	
	private void loadVcfIndex() {
		if(this.vcfIndexFile!=null) return;
		if(this.filenameOrNull==null) throw new RuntimeIOException("Cannot load a VCF index as the bgen is not a file (but probably a stream) ");
		final Path vcfPath = getVcfIndexName(this.filenameOrNull);
		if(!Files.exists(vcfPath))  throw new RuntimeIOException("cannot find associated VCF file "+vcfPath);
		this.vcfIndexFile=new VCFFileReader(vcfPath,true);
		final VCFHeader hder = this.vcfIndexFile.getHeader();
		if(!hder.hasInfoLine(OFFSET_format_header_line.getID())) {
			this.vcfIndexFile.close();
			this.vcfIndexFile=null;
			throw new RuntimeIOException("cannot use "+vcfPath+" as index because the header is missing ##INFO/"+OFFSET_format_header_line.getID());
			}
		}
	
	/** get the name of the VCF index for the bgen without checking it exists */
	public static Path getVcfIndexName(Path bgenPath) {
		final String filename = Objects.requireNonNull(bgenPath).toAbsolutePath().toString();
		final  String vcf = ParsingUtils.appendToPath(filename, VCF_INDEX_SUFFIX);
		return bgenPath.getFileSystem().getPath(vcf);
		}
	
	private List<String> readNAlleles(final int num_alleles) throws IOException {
	      final String[] alleles = new String[num_alleles];
	      for(int k=0;k< num_alleles;++k) {
	            alleles[k] = readStringUInt32(binaryCodec);
	            }
	      return Arrays.asList(alleles);
	      } 
	
	/** read the next variant or returns null if end of file */
	public BGenVariant readVariant() throws IOException {
		assertState(State.expect_variant_def);
		
		final VariantImpl ctx = new VariantImpl();
		ctx.offset = ftell();
		
		if(header.getLayout().equals(Layout.e_Layout1)) {
			try {
				this.layout1_expect_n_samples = longToUnsignedInt(this.binaryCodec.readUInt());
				}
			catch(Throwable eof) {
				if(isDebugging()) {
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
			ctx.variantId = readStringUInt16(this.binaryCodec);
			}
		catch(Throwable eof) {
			if(isDebugging()) {
				LOG.debug("not more variant "+header.getLayout()+" "+eof.getClass());
				eof.printStackTrace();
				}
			if(!header.getLayout().equals(Layout.e_Layout1)) {
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
		if(isDebugging()) {
			LOG.debug("skipping "+skip_n+" bytes");
			}
		this.binaryCodec.getInputStream().skipNBytes(skip_n);
	    this.state = State.expect_variant_def;
		}
	
	@SuppressWarnings("resource")
	public BGenVariant readGenotypes()  throws IOException{
		assertState(State.expect_genotypes);
		
		if(isDebugging()) {
			LOG.debug("reading genotypes. offset= "+ftell());
			}
		
		int genotype_block_size =  longToUnsignedInt( this.binaryCodec.readUInt());
		if(isDebugging()) {
			LOG.debug("expect compressed size "+genotype_block_size);
			}
		
		/* The total length D of the probability data after uncompression.
		 * If CompressedSNPBlocks = 0, this field is omitted and the total length of the probability data */
		if(!getHeader().getCompression().equals(Compression.e_NoCompression)) {
			final int uncompressed_data_size_D =longToUnsignedInt( this.binaryCodec.readUInt());
			if(isDebugging()) {
				LOG.debug("expect uncompressed size "+uncompressed_data_size_D+" now offset="+ftell());
				}
			genotype_block_size-= Integer.BYTES;
			}
		
		this.buffer1.reset();
		this.buffer1.copyTo(this.in,genotype_block_size);
		if(this.buffer1.size()!=genotype_block_size) {
			throw new IllegalStateException("bad size");
			}
		if(isDebugging()) {
			LOG.debug("after reading genotypes. offset is now = "+ftell());
			}
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
		
		if(isDebugging()) {
			LOG.debug("End reading genotypes. offset= "+ftell());
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
	
	private boolean isDebugging() {
		return this.debug_flag;
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
        	vcg.at(x).metaByte =  codec2.readByte();
        	}
        /* Flag, denoted Phased indicating what is stored in the row. */
        vcg.phased = codec2.readUByte();
	      
        /* Unsigned integer B representing the number of bits used to store each probability in this row. This must be between 1 and 32 inclusive. */
       	vcg.nbits=  (int) codec2.readUByte();
        if(vcg.nbits  <0 || vcg.nbits > 32) throw new IOException("bad nbits  ("+ vcg.nbits +")");
        if(isDebugging()) {
   			LOG.debug("nbits="+vcg.nbits);
   			}
        
        final BitReader bitReader = new BitReader(codec2.getInputStream());
        final BitNumReader numReader = new BitNumReader(bitReader, vcg.nbits);
        final List<Double> values= new ArrayList<>();
        if(vcg.phased==1) {
        	if(isDebugging()) {
       			LOG.debug("phased==1");
       			}
    	   /* read haplotypes */
       		for(int i=0;i< n_samples;++i) {
       			values.clear();
       			VariantAndGenotypesLayout2.GenotypeLayout2 gt= vcg.at(i);
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
       			
       			final VariantAndGenotypesLayout2.GenotypeLayout2 gt= vcg.at(i);
       			if(isDebugging()) {
           			LOG.debug("read unphased genotype["+i+"] ploidy :"+gt.getPloidy()+" missing:"+gt.isMissing()+" n-alleles:"+this.lastVariant.getNAlleles());
           			}
       			final BigInteger bi=MathUtils.combination( gt.getPloidy()+ this.lastVariant.getNAlleles()-1,this.lastVariant.getNAlleles()-1);// whyyyy ?
       			if(isDebugging()) {
           			LOG.debug("combination="+bi);
           			}
       			final int n_combinaison_as_int= bi.intValueExact();
       			if(n_combinaison_as_int!=new AlleleCombinations(this.lastVariant.getNAlleles(),gt.getPloidy(),false).calculateTotalCombinations()) {
       				throw new IllegalStateException();
       				}
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
	
	private static class AlleleCombinations {
		private final int n_alleles;
		private final int m_ploidy;
		private final boolean phased;
		
		private static class FindByIndex implements Consumer<int[]> {
			final int wanted_index;
			int curr=0;
			int[] indexes = null;
			FindByIndex(final int wanted_index) {
				this.wanted_index = wanted_index;
				}
			@Override
			public void accept(int[] t) {
				if(curr==this.wanted_index) {
					this.indexes=Arrays.copyOf(t, t.length);
					}
				curr++;
				}
			}
		
		AlleleCombinations(int n_alleles,int m_ploidy,boolean phased) {
			this.n_alleles=n_alleles;
			this.m_ploidy=m_ploidy;
			this.phased = phased;
			}
		
		 private int calculateTotalCombinations() {
			int n = n_alleles;
		    int k = m_ploidy;
			final BigInteger bi;
			if(phased) {
				bi = BigInteger.valueOf(n_alleles).pow(m_ploidy);
				}
			else
				{
				bi = MathUtils.factorial(n + k - 1).divide(MathUtils.factorial(k).multiply(MathUtils.factorial(n - 1)));
				}
			return bi.intValueExact();
			}
		
		  int[] findAlleleIndexesByIndex(int i) {
			  final int[] buffer=new int[this.m_ploidy];
			  FindByIndex andThen = new FindByIndex(i);
			  if(phased) {
		        	visitPhased(andThen, buffer, 0);
		        	}
		        else
		        	{
		        	visitUnphased(andThen, buffer, 0, 0);
		        	}
			  return andThen.indexes;
		  	}
		
		 
		  public List<int[]> getAllGenotypesPhased() {
	        final List<int[]> container = new ArrayList<>(calculateTotalCombinations());
        	final int[] buffer=new int[this.m_ploidy];
        	final Consumer<int[]> anThen =A->container.add(Arrays.copyOf(A, A.length));
	        if(phased) {
	        	visitPhased(anThen, buffer, 0);
	        	}
	        else
	        	{
	        	visitUnphased(anThen, buffer, 0, 0);
	        	}
	        return container;
		  	}
		  /*
		     alleles: 0,1,2
		     list:
		             0/0
		             0/1
		             0/2
		             1/1
		             1/2
		             2/2
		     */
	    private void visitUnphased(Consumer<int[]> anThen, int[] current, int start,int array_index) {	    	
	        if (array_index == m_ploidy) {
	        	anThen.accept(current);
	        	}
	        else
		        {
		        for (int i = start; i < n_alleles; i++) {
		            current[array_index]=i;
		            visitUnphased(anThen, current, i,array_index+1); // Allow repetition by not incrementing 'i'
		        	}
		        }
		    }
		  /*
	     alleles: 0,1,2
	     list:
	             0|0
	             0|1
	             0|2
	             1|0
	             1|1
	             1|2
	             2|0
	             2|1
	             2|2
	     */
	    private void visitPhased(Consumer<int[]> anThen, int[] current, int array_index) {
	        if (array_index == m_ploidy) {
	        	anThen.accept(current);
	        }
	        else {
	        for (int i = 0; i < n_alleles; i++) {
		            current[array_index] = i;
		            visitPhased(anThen, current, array_index + 1); // Move to the next depth
		        	}
		        }
	    	}
		}
	
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
			}
	
	@Override
	public void close()  {
		if(this.vcfIndexFile!=null) {
			vcfIndexFile.close();
			vcfIndexFile=null;
			}
		
		try {this.in.close();}
		catch(IOException err) {}
		this.binaryCodec.close();
		}
}
