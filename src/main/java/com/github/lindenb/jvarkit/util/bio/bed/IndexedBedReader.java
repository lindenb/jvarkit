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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.util.bio.bed;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.LineIterator;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;


/**
 * @author lindenb
 * 
 * Tabix or Tribble BED reader.
 */
public class IndexedBedReader
	implements Closeable
	{
	private static final Logger LOG=Logger.build(IndexedBedReader.class).make();

		
	private Object source;
	private AbstractIndexReader reader;
	
	public IndexedBedReader(final File bedFile) throws IOException
		{
		init(bedFile);
		}
	
	
	private void init(final File bedFile) throws IOException
		{
		this.source=bedFile;
    	if(bedFile==null) throw new NullPointerException("bed file==null");
    	IOUtil.assertFileIsReadable(bedFile);
    	if(bedFile.getName().endsWith(".gz"))
    		{
    		this.reader = new TabixReader(bedFile.getPath());
    		}
    	else
    		{
    		this.reader = new TribbleReader(bedFile);
    		}
		}
	
	public IndexedBedReader(final String onlyTabixOrLocal) throws IOException
		{
		if(onlyTabixOrLocal==null) throw new NullPointerException("bed file==null");
		
		if(IOUtil.isUrl(onlyTabixOrLocal))
			{
			this.reader = new TabixReader(onlyTabixOrLocal);
			this.source=onlyTabixOrLocal;
			}
		else
			{
			init(new File(onlyTabixOrLocal));
			}
		
		}

	
	private void checkOpen()
		{
		if(this.reader==null)
				throw new IllegalStateException("bed reader is closed "+getSource());
		}
	
	/** string of File */
	public Object getSource()
		{
		return this.source;
		}
	
	public CloseableIterator<BedLine>
		iterator(final String chrom,final int start,final int end)
		throws IOException
		{
		checkOpen();
		return this.reader.query(chrom, start, end);
		}
	
	/** return distinct contigs in this bed */
	public Set<String> getContigs() {
		return new LinkedHashSet<>(this.reader.getContigs());
		}
	
	public List<BedLine> getLines(final String chrom,final int start,final int end,final Predicate<BedLine> predicate) throws IOException
		{
		checkOpen();
		CloseableIterator<BedLine> iter=null;
		try
			{
			final List<BedLine> L= new ArrayList<>();
			iter = iterator(chrom,start,end);
			while(iter.hasNext())
				{
				final BedLine bed = iter.next();
				if( predicate!=null && !predicate.test(bed)) continue;
				L.add(bed);
				}
			return L;
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	
	public List<BedLine> getLines(final String chrom,final int start,final int end) throws IOException
		{
		return getLines(chrom,start,end,BED->true);
		}
	
	@Override
	public void close() throws IOException
		{
		CloserUtil.close(reader);
		this.reader=null;
		}
	
	private interface AbstractIndexReader
		{
		public Collection<String> getContigs();
		public CloseableIterator<BedLine> query(String chrom,int start,int end) throws IOException;
		public void close();
		}
	
	private class TribbleReader
		implements AbstractIndexReader
		{
		private Index tribbleIndex=null;
	    private AbstractFeatureReader<BedLine, LineIterator> reader;
	    private BedLineCodec bedCodec =new BedLineCodec();

		TribbleReader(final File bedFile)  throws IOException
			{
			final File indexFile = Tribble.indexFile(bedFile);
			if(indexFile.exists())
	    		{
				LOG.info("loading tribble index in memory for "+ bedFile);
				this.tribbleIndex=IndexFactory.loadIndex(indexFile.getPath());
	    		}
			else
			 	{
				LOG.info("creating index from file "+bedFile+" indexFile:"+indexFile);
				this.tribbleIndex=IndexFactory.createLinearIndex(bedFile, this.bedCodec);
			 	}
			this.reader = AbstractFeatureReader.getFeatureReader(
				bedFile.getPath(),
				this.bedCodec,
				this.tribbleIndex
				);
			}
		@Override
		public CloseableIterator<BedLine> query(
				final String chrom,
				final int start,
				final int end) throws IOException {
			return this.reader.query(chrom, start, end);
			}
		@Override
		public void close() {
			CloserUtil.close(reader);
			this.reader=null;
			this.tribbleIndex=null;
			}
		@Override
		public List<String> getContigs() {
			return this.reader.getSequenceNames();
			}
		}
	private class TabixReader
		extends AbstractTabixObjectReader<BedLine>
		implements AbstractIndexReader
		{
		TabixReader(final String bedFile) throws IOException
			{
			super(bedFile);
			}
		
		@Override
		public Collection<String> getContigs() {
			return super.getChromosomes();
		}
		
		@Override
		public CloseableIterator<BedLine> query(final String chrom, final int start,final int end)
				throws IOException {
			return (CloseableIterator<BedLine>)this.iterator(chrom, start, end);
			}
		
		@Override
		protected CloseableIterator<BedLine> iterator(final Iterator<String> delegate) {
			return new MyIterator(delegate);
			}
		
	    private class MyIterator
    	extends AbstractMyIterator
    	implements CloseableIterator<BedLine>
	    	{
	    	private final Pattern tab=Pattern.compile("[\t]");
	    	MyIterator(final Iterator<String> delegate)
	    		{
	    		super(delegate);
	    		}
	    	
	    	@Override
	    	public BedLine next() {
	    		final String tokens[]=this.tab.split(delegate.next());
	    		return new BedLine(tokens);
	    		}
	    	@Override
	    	public void close() {
	    		//nothing to do
	    		}
	    	}
		}
	}
