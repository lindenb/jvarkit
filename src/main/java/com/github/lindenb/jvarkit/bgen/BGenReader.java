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

import java.io.IOException;
import java.nio.file.Path;
import java.util.Objects;


import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;


public class BGenReader extends BGenUtils implements AutoCloseable {
	private static final Logger LOG = Logger.of(BGenReader.class).setDebug();
	//private final AbstractFeatureReader<BGenVariant, SeekableStream> featureReader;
	/** auxiliary indexed VCF file that will be used as an coordinate index for the VCF */
	private final BGenCodec codec = new BGenCodec(false);
	private final BGenUtils.RandomAccessSeekableStream seekableStream;
	private VCFReader vcfIndexFileReader = null;
	private final String sourceName;
	

    public BGenReader(final String path) throws IOException {
		this.seekableStream = new RandomAccessSeekableStream(path);
		this.sourceName=path;
		this.codec.readActualHeader(this.seekableStream);
		}
	
 
    
	public BGenReader(final Path path) throws IOException {
		this(Objects.requireNonNull(path,"null input").toString());
		}
	
	public void attachVcfReaderAsIndex( VCFReader vcfIndexFileReaderOrNull) {
		this.vcfIndexFileReader = vcfIndexFileReaderOrNull;
		}
	
	/** rename samples using a plink sample file */
	public void bindPLinkSampleFile(final Path path) throws IOException {
		this.codec.getHeader().bindPLinkSampleFile(path);
		}
	
	public BGenHeader getHeader() {
		return this.codec.getHeader();
		}
	
	public BGenUtils.Layout getLayout() {
		return getHeader().getLayout();
		}
	
	public BGenUtils.Compression getCompression() {
		return getHeader().getCompression();
	}
	
	public long getSnpsOffset() {
		return this.codec.getSnpsOffset();
		}
	
	/** reset reader to first variant position */
	public void rewind() throws IOException {
		this.codec.rewind(this.seekableStream);
	}
	 
	/** change the reading pointer offset
	 * 
	 * @param position the physical offset
	 * @throws IOException the inutstream is not an instance of SeekableStream
	 */
	public void fseek(long position)  throws IOException {
		if(position<=0) throw new IllegalArgumentException("seek.pos<=0: "+position);
		this.seekableStream.fseek(position);
		}
	

	
	
	
	public BGenIterator iterator(boolean skipGenotypes) {
		return new Iter(skipGenotypes);
		}
	
	public BGenIterator iterator() {
		return iterator(false);
		}
	
	
	
	private abstract class BaseIterator extends AbstractIterator<BGenVariant> implements BGenIterator {
		protected final boolean skipGenotypes;
		BaseIterator(boolean skipGenotypes) {
			this.skipGenotypes = skipGenotypes;
			}
		
		@Override
		public void seekTo(long file_offset) throws IOException {
			BGenReader.this.fseek(file_offset);
			}
		@Override
		public BGenHeader getHeader() {
			return BGenReader.this.getHeader();
			}
		@Override
		public BGenVariant readVariant() throws IOException {
			return BGenReader.this.readVariant();
			}
		@Override
		public void skipGenotypes() throws IOException {
			BGenReader.this.skipGenotypes();
			}
		@Override
		public BGenVariant readGenotypes() throws IOException {
			return BGenReader.this.readGenotypes();
			}
		@Override
		public long getSnpsOffset() {
			return BGenReader.this.codec.getSnpsOffset();
			}
		}
	
	private class Iter extends  BaseIterator {
		Iter(boolean skipGenotypes) {
			super(skipGenotypes);
			}
	
		@Override
		protected BGenVariant advance() {
			try {
				final BGenVariant v= this.readVariant();
				if(v==null) {
					return null;
					}
				else if(super.skipGenotypes) {
					this.skipGenotypes();
					return v;
					}
				else
					{
					return this.readGenotypes();
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		
		@Override
		public void close() {
			
			}
		};

	private class OffsetIterator extends BaseIterator {
		private final CloseableIterator<Long> offsetIterator;
		OffsetIterator(final CloseableIterator<Long> offsetIterator,boolean skipGenotypes) throws IOException{
			super(skipGenotypes);
			this.offsetIterator= offsetIterator;
			}
		@Override
		protected BGenVariant advance() {
			try {
				for(;;) {
					if(!offsetIterator.hasNext()) {
						close();
						return null;
						}
					final Long offset= offsetIterator.next();
					if(offset==null) continue;
					if(offset.longValue()<=0L) throw new IllegalArgumentException("bad offset: "+offset);
					BGenReader.this.fseek(offset.longValue());
					BGenVariant bv = readVariant();
					if(skipGenotypes) {
						skipGenotypes();
						}
					else
						{
						bv = readGenotypes();
						}
					return bv;
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() {
			this.offsetIterator.close();
			}
		}
		
	private class OverlapIterator extends BaseIterator {
		final Locatable query;
		OverlapIterator(final Locatable query,boolean skipGenotypes) {
			super(skipGenotypes);
			this.query=query;
			}
		@Override
		protected BGenVariant advance() {
			try {
				for(;;) {
					BGenVariant bv = readVariant();
					if(bv==null) return null;
					if(!bv.overlaps(this.query)) {
						skipGenotypes();
						continue;
						}
					if(skipGenotypes) {
						skipGenotypes();
						}
					else
						{
						bv = readGenotypes();
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
			}
		}
	
	private static class VcfOffsetIterator extends AbstractIterator<Long> implements CloseableIterator<Long>{
		private CloseableIterator<VariantContext> delegate;
		VcfOffsetIterator(CloseableIterator<VariantContext> delegate) {
			this.delegate = delegate;
			}
		@Override
		protected Long advance() {
			for(;;) {
				if(!delegate.hasNext()) return null;
				final VariantContext ctx = this.delegate.next();
				if(!ctx.hasAttribute(OFFSET_format_header_line.getID())) continue;
				final long offset= Long.parseLong(ctx.getAttributeAsString(OFFSET_format_header_line.getID(), "-1"));
				if(offset<=0L) throw new IllegalArgumentException("bad offset in "+ctx);
				return offset;
				}
			}
		@Override
		public void close() {
			delegate.close();
			}
		}
	
	private class QueryIterator extends BaseIterator {
		private CloseableIterator<VariantContext> delegate;
		private final Locatable query;
		QueryIterator(final Locatable query,CloseableIterator<VariantContext> delegate,boolean skipGenotypes) {
			super(skipGenotypes);
			this.delegate=delegate;
			this.query = query;
			}

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
					BGenVariant bv = readVariant();
					if(skipGenotypes) {
						skipGenotypes();
						}
					else
						{
						bv = readGenotypes();
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
	
	public CloseableIterator<BGenVariant> query(final String loc,int start,int end) throws IOException {
		return query(new SimpleInterval(loc,start,end),false);
		}
	
	/**
	 * if the BGen is VCF-indexed, use the VCF index to scan for the variant overlapping loc
	 * if the BGen is not VCF-indexed, call rewind and scan all the vcf for an overlap
	 * @param loc
	 * @param skipGenotypes
	 * @return
	 * @throws IOException
	 */
	public BGenIterator query(final Locatable loc, boolean skipGenotypes) throws IOException  {
		if(hasIndex()) {
			return new OffsetIterator(
					new VcfOffsetIterator(Objects.requireNonNull(this.vcfIndexFileReader).query(loc)),
					skipGenotypes
					);
			}
		else
			{
			if(LOG.isDebug()) {
				LOG.debug("index=false, rewing and scan all for overlap");
				}
			rewind();
			return new OverlapIterator(loc,skipGenotypes);
			}
		}
	
	public boolean hasIndex() {
		if(this.vcfIndexFileReader==null) return false;
		if(!this.vcfIndexFileReader.isQueryable()) return false;
		final VCFHeader hder = vcfIndexFileReader.getHeader();
		if(!hder.hasInfoLine(OFFSET_format_header_line.getID())) return false;
		return true;
		}
	
	
	/** read the next variant or returns null if end of file */
	public BGenVariant readVariant() throws IOException {
		return this.codec.readVariant(this.seekableStream);
	}
	
	public void skipGenotypes() throws IOException{
		this.codec.skipGenotypes(this.seekableStream);
		}
	
	public BGenVariant readGenotypes()  throws IOException{
		return this.codec.readGenotypes(this.seekableStream);
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
	public void close()  {
		try {this.codec.close(this.seekableStream);}
		catch(Throwable err) {}
		}
	
	@Override
	public String toString() {
		return "BGenReader(file:"+this.sourceName+")";
		}
}
