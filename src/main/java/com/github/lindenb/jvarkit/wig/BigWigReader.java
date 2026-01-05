/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.wig;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;

import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.samtools.util.ExtendedLocatable;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;


/** 
 * BigWig fle reader
 */
public class BigWigReader implements AutoCloseable {
	public static final String OPT_DESC="Path/URI to bigwig file.";
	private final String biwWigFile;
	private BBFileReader bbFileReader=null;
	private ContigNameConverter contigNameConverter = null;

	public static interface WigItem extends ExtendedLocatable {
		public float getValue();
	}
	
	private static class WigItemImpl extends SimpleInterval implements WigItem {
		final float value;
		WigItemImpl(String userContig,final org.broad.igv.bbfile.WigItem delegate) {
			super(userContig,delegate.getStartBase()+1 /* OK, checked */,delegate.getEndBase());
			this.value = delegate.getWigValue();
			}
		@Override
		public float getValue() {
			return this.value;
			}
		@Override
		public String toString() {
			return super.toString()+"="+getValue();
			}
		}
	
	public BigWigReader(final Path path) throws IOException {
		this(Objects.requireNonNull(path).toString());
		}
	
	public BigWigReader(final String biwWigFile) throws IOException {
		this.biwWigFile = Objects.requireNonNull(biwWigFile);
		try {
			this.bbFileReader= new BBFileReader(this.biwWigFile);
			}
		catch(final IOException err)
			{
			throw new IOException("Cannot open "+this.biwWigFile,err);
			}
		if(!this.bbFileReader.isBigWigFile())
			{
			this.bbFileReader.close();
			throw new IOException(this.biwWigFile+" is not a bigWIG file.");
			}
		this.contigNameConverter = ContigNameConverter.fromContigSet(new HashSet<>(this.bbFileReader.getChromosomeNames()));
		}
	
	public Set<String> getChromosomeNames() {
		return new HashSet<>(this.bbFileReader.getChromosomeNames());
	}
	
	public CloseableIterator<WigItem> query(final Locatable locatable)
		{
		if(this.bbFileReader==null) return AbstractCloseableIterator.empty();
		final String resolvedContig = this.contigNameConverter.apply(locatable.getContig());
		if(StringUtils.isBlank(resolvedContig))  return AbstractCloseableIterator.empty();
		return new WigItemIterator(
				this.bbFileReader.getBigWigIterator(
				resolvedContig,
				locatable.getStart()-1,
				resolvedContig,
				locatable.getEnd(),
				false
				),locatable.getContig());
		}
	
	@Override
	public void close() {
		try { if(this.bbFileReader!=null) this.bbFileReader.close();
		} catch(final Throwable err){}
		this.bbFileReader=null;
		}
	
	public String getURI() {
		return biwWigFile;
		}
	
	@Override
	public String toString() {
		return "BigWig("+getURI()+")";
		}
	
	/** wraps a BigWigIterator */
	private static class WigItemIterator
		extends AbstractIterator<WigItem>
		implements CloseableIterator<WigItem>
		{
		private BigWigIterator delegate;
		private final String userContig;
		WigItemIterator(final BigWigIterator delegate,final String userContig ){
			this.delegate  = delegate;
			this.userContig = userContig;
			}
		@Override
		protected WigItem advance() {
			if(this.delegate==null) return null;
			if(!this.delegate.hasNext()) return null;
			return new WigItemImpl(this.userContig,this.delegate.next());
			}
		@Override
		public void close() {
			this.delegate=null;
			}
		}
	}
