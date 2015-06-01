/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.LineIterator;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;


/**
 * @author lindenb
 *
 */
public class IndexedBedReader
	implements Closeable
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");

	
	public static class BedLine
		implements Feature
		{
		private String tokens[];
		
		BedLine(String tokens[])
			{
			this.tokens=tokens;
			}
		@Override
		@Deprecated
		public String getChr() {
			return tokens[0];
			}
		
		@Override
		public String getContig() {
			return tokens[0];
			}
		@Override
		public int getStart() {
			return Integer.parseInt(tokens[1]);
			}
		
		@Override
		public int getEnd() {
			return (tokens.length<3 ?getStart(): Integer.parseInt(tokens[2]));
			}

		public String get(int index)
			{
			return (index<tokens.length?tokens[index]:null);
			}
		
		public int getColumnCount()
			{
			return tokens.length;
			}
		
		}
	
	private File source;
	private AbstractIndexReader reader;
	
	public IndexedBedReader(File bedFile) throws IOException
		{
		this.source=bedFile;
    	if(bedFile==null) throw new NullPointerException("bed file==null");
    	if(!bedFile.isFile())  throw new IOException("bed is not a file "+bedFile);
    	if(!bedFile.canRead())  throw new IOException("cannot read "+bedFile);
    	if(bedFile.getName().endsWith(".gz"))
    		{
    		this.reader = new TabixReader(bedFile);
    		}
    	else
    		{
    		this.reader = new TribbleReader(bedFile);
    		}
		}
	
	private void checkOpen()
		{
		if(this.reader==null)
				throw new IllegalStateException("bed reader is closed "+getSource());
		}
	
	public File getSource()
		{
		return this.source;
		}
	

	public CloseableIterator<BedLine>
		iterator(String chrom,int start,int end)
		throws IOException
		{
		checkOpen();
		return this.reader.query(chrom, start, end);
		}
	
	public List<BedLine> getLines(String chrom,int start,int end) throws IOException
		{
		checkOpen();
		CloseableIterator<BedLine> iter=null;
		try
			{
			List<BedLine> L= new ArrayList<>();
			iter = iterator(chrom,start,end);
			while(iter.hasNext())
				{
				L.add(iter.next());
				}
			return L;
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	
	
	@Override
	public void close() throws IOException
		{
		CloserUtil.close(reader);
		this.reader=null;
		}
	
	private interface AbstractIndexReader
		{
		public abstract CloseableIterator<BedLine> query(String chrom,int start,int end) throws IOException;
		public abstract void close();
		}
	
	private class TribbleReader
		implements AbstractIndexReader
		{
		private Index tribbleIndex=null;
	    private AbstractFeatureReader<BedLine, LineIterator> reader;
	    private BedLineCodec bedCodec =new BedLineCodec();

		TribbleReader(File bedFile)  throws IOException
			{
			File indexFile=Tribble.indexFile(bedFile);
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
		public CloseableIterator<BedLine> query(String chrom, int start,
				int end) throws IOException {
			return this.reader.query(chrom, start, end);
			}
		@Override
		public void close() {
			CloserUtil.close(reader);
			this.reader=null;
			this.tribbleIndex=null;
			}
		}
	private class TabixReader
		extends AbstractTabixObjectReader<BedLine>
		implements AbstractIndexReader
		{
		TabixReader(File bedFile) throws IOException
			{
			super(bedFile.getPath());
			}
		
		@Override
		public CloseableIterator<BedLine> query(String chrom, int start, int end)
				throws IOException {
			return (CloseableIterator<BedLine>)this.iterator(chrom, start, end);
			}
		
		@Override
		protected CloseableIterator<BedLine> iterator(Iterator<String> delegate) {
			return new MyIterator(delegate);
			}
		
	    private class MyIterator
    	extends AbstractMyIterator
    	implements CloseableIterator<BedLine>
	    	{
	    	private Pattern tab=Pattern.compile("[\t]");
	    	MyIterator(Iterator<String> delegate)
	    		{
	    		super(delegate);
	    		}
	    	
	    	@Override
	    	public BedLine next() {
	    		String tokens[]=this.tab.split(delegate.next());
	    		return new BedLine(tokens);
	    		}
	    	@Override
	    	public void close() {
	    		//nothing to do
	    		}
	    	}
		}
	private static class BedLineCodec
	extends AsciiFeatureCodec<BedLine>
		{
		private Pattern tab=Pattern.compile("[\t]");
		public BedLineCodec() {
			super(BedLine.class);
			}
		
		@Override
		public BedLine decode(String line) {
			
			if (line.trim().isEmpty()) {
	            return null;
	        	}

	        if (line.startsWith("#") || line.startsWith("track") || line.startsWith("browser")) {
	            return null;
	        	}

	        String[] tokens = tab.split(line);
	        if(tokens.length<2) return null;
	       
	        return new BedLine(tokens);
	        }
		
		@Override
		public Object readActualHeader(LineIterator reader) {
			return null;
			}
		
		}

}
