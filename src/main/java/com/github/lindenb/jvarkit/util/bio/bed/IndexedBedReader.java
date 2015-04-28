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
import htsjdk.tribble.readers.LineIterator;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * @author lindenb
 *
 */
public class IndexedBedReader
	implements Closeable
	{
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
	
	public static class BedLine
		implements Feature
		{
		private String tokens[];
		
		BedLine(String tokens[])
			{
			this.tokens=tokens;
			}
		@Override
		public String getChr() {
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
	
	private BedLineCodec bedCodec=null;
    private AbstractFeatureReader<BedLine, LineIterator> reader;
	private Object source;
	
	private IndexedBedReader(Object bedFile,int ignore) throws IOException
		{
		this.source=bedFile;
		this.bedCodec=new BedLineCodec();
		}
	
	public IndexedBedReader(File bedFile) throws IOException
		{
		this(bedFile,0);
    	if(bedFile==null) throw new NullPointerException("bed file==null");
    	if(!bedFile.isFile())  throw new IOException("bed is not a file "+bedFile);
    	if(!bedFile.canRead())  throw new IOException("cannot read "+bedFile);
    	this.reader=AbstractFeatureReader.getFeatureReader(bedFile.getPath(), bedCodec,false);
 
		}
	public IndexedBedReader(String path) throws IOException
		{
		this(path,0);
    	this.reader=AbstractFeatureReader.getFeatureReader(path, bedCodec,false);
		}
	private void checkOpen()
		{
		if(this.reader==null)
				throw new IllegalStateException("bed reader is closed "+getSource());
		}
	
	public Object getSource()
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
		reader=null;
		}
	

	
}
