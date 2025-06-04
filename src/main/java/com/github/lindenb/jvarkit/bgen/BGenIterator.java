package com.github.lindenb.jvarkit.bgen;

import java.io.IOException;
import java.io.InputStream;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;

public interface  BGenIterator extends CloseableIterator<BGenVariant> {

	
	public static BGenIterator of(InputStream in,boolean skipGenotypes) throws IOException {
		return new BGenIteratorImpl(in,skipGenotypes);
		}
	public static BGenIterator of(InputStream in) throws IOException {
		return of(in,false);
		}
	public BGenHeader getHeader();
	public void seekTo(long file_offset) throws IOException;
	public BGenVariant readVariant() throws IOException;
	public BGenVariant readGenotypes() throws IOException;
	public void skipGenotypes() throws IOException;
	public long getSnpsOffset();
	
	
		static class BGenIteratorImpl extends AbstractIterator<BGenVariant> implements BGenIterator {
			private final BGenCodec codec = new BGenCodec(false);
			private final BGenUtils.RandomAccessPositionalStream in;
			private final boolean skipGenotypes;
			
			
			
			public BGenIteratorImpl(final InputStream in,boolean skipGenotypes) throws IOException {
				this.in = new BGenUtils.RandomAccessPositionalStream(in);
				this.codec.readActualHeader(this.in);
				this.skipGenotypes=skipGenotypes;
				}
	
			@Override
			public BGenHeader getHeader() {
				return this.codec.getHeader();
				}
			
			@Override
			protected BGenVariant advance() {
				try {
					final BGenVariant v=readVariant();
					if(v==null) return null;
					if(this.skipGenotypes) {
						this.skipGenotypes();
						return v;
					}
					else
					{
						return readGenotypes();
					}
				}
				catch(final IOException err) {
					throw new RuntimeIOException(err);
				}
			}
	
			@Override
			public long getSnpsOffset() {
				return codec.getSnpsOffset();
				}
			
			@Override
			public void seekTo(long file_offset) throws IOException  {
				this.in.fseek(file_offset);
			}
			@Override
			public BGenVariant readVariant() throws IOException {
				if(this.in.eof()) return null;
				return this.codec.readVariant(this.in);	
			}
			@Override
			public BGenVariant readGenotypes() throws IOException  {
				return this.codec.readGenotypes(this.in);	
			}
			@Override
			public void skipGenotypes() throws IOException  {
				this.codec.skipGenotypes(this.in);
			}
	
			@Override
			public void close() {
				codec.close(in);
			}
		}
}
