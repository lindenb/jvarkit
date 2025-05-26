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
package com.github.lindenb.jvarkit.wib;

import java.io.IOException;
import java.util.HashSet;
import java.util.function.UnaryOperator;

import com.github.lindenb.jvarkit.bed.BedInterval;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.readers.LineIterator;
/**
 * class reading UCSC wib files
 */
public class WibReader implements AutoCloseable {
	public static final String TABIX_DESC="A wib associated indexed tabix file.e.g: wget \"https://hgdownload.soe.ucsc.edu/gbdb/hg38/multiz100way/phastCons100way.wib\" && wget -O - \"https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/phastCons100way.txt.gz\" | gunzip -c | bgzip > phastCons100way.txt.gz && tabix --force -0 -s 2 -e 3 -e 4 phastCons100way.txt.gz\n";
	public static final String WIB_DESC="The wib file itself e.g: wget 'https://hgdownload.soe.ucsc.edu/gbdb/hg38/multiz100way/phastCons100way.wib' .";

	  static public interface WigItem extends Locatable, BedInterval {
		  public double getValue();
	  }
	  
	  static private class WigItemImpl extends SimpleInterval implements WigItem {
			final double value;
			WigItemImpl(String C,int S,int E,double value) {
				super(C,S,E);
				this.value=value;
			}
			@Override
			public double getValue() {
				return value;
				}
			@Override
			public String toString() {
				return super.toString()+" "+getValue();
				}
	  	}
	
	  static class WiggleForWib extends SimpleInterval implements Feature {
          int span;
          int count;
		  long offset;
          String wibFile;
          float lowerLimit;
          float dataRange;

		  WiggleForWib(final String[] tokens) {
			super(tokens[1],Integer.parseInt(tokens[2])+1,Integer.parseInt(tokens[3]));
			// name = tokens[4]
			this.span = Integer.parseInt(tokens[5]);
			this.count = Integer.parseInt(tokens[6]);
			this.offset = Long.parseLong(tokens[7]);
			this.wibFile = tokens[8];
			this.lowerLimit = Float.parseFloat(tokens[9]);
			this.dataRange = Float.parseFloat(tokens[10]);
			}

	  	}
	
	 static class WiggleForWibCodec extends AsciiFeatureCodec<WiggleForWib> {
		 WiggleForWibCodec() {
			 super(WiggleForWib.class);
		 	}
		@Override
		public WiggleForWib decode(final String s) {
			return new WiggleForWib(CharSplitter.TAB.split(s));
			}
		@Override
		public boolean canDecode(String fname) {
			return fname.endsWith(".bed.gz") ||  fname.endsWith(".txt.gz");
			}
		
		@Override
		public Object readActualHeader(LineIterator li) {
			return null;
			}
	 	}

	 /** simple path converter, exepect only one wibPath in the tabix file */
	private static class SimplePathConverter implements UnaryOperator<String> {
		private final String resolvedPath;
		private String wibName = null;
		SimplePathConverter( final String resolvedPath) {
			this.resolvedPath=resolvedPath;
			}
		@Override
		public String apply(String wibName) {
			if(this.wibName==null) {
				this.wibName = wibName;
				}
			else if(!this.wibName.equals(wibName)) {
				throw new IllegalArgumentException("expected only one wib path in the file but found '"+wibName+"' and '"+this.wibName+"'.");
				}
			return this.resolvedPath;
			}
		} 
	private final  AbstractFeatureReader<WiggleForWib, LineIterator> tabixReader;
	private final UnaryOperator<String> pathConverter;
	private ISeekableStreamFactory seekableStreamFactory;
	private SeekableStream wibSeekableStream=null;
	private String wibResolvedName = null;
	private final UnaryOperator<String> ctgConverter;
	
	public WibReader(final String tabixPath) throws IOException {
		this(tabixPath,X->X, SeekableStreamFactory.getInstance());
		}
	
	public WibReader(final String tabixPath, final String wibPath) throws IOException {
		this(tabixPath,new SimplePathConverter(wibPath), SeekableStreamFactory.getInstance());
		}
	public WibReader(final String tabixPath, final UnaryOperator<String> pathConverter,ISeekableStreamFactory seekableStreamFactory) throws IOException {	
		final WiggleForWibCodec codec = new WiggleForWibCodec();
		this.seekableStreamFactory = seekableStreamFactory;
		this.pathConverter = pathConverter ;
		this.tabixReader = AbstractFeatureReader.getFeatureReader(tabixPath,codec);
		this.ctgConverter = ContigNameConverter.fromContigSet(new HashSet<>(this.tabixReader.getSequenceNames()));
		this.wibSeekableStream = null;
		}
	
	/** return Item in the specified range */
	public CloseableIterator<WigItem> query(final Locatable loc) throws IOException {
		final String ctg = this.ctgConverter.apply(loc.getContig());
		if(StringUtils.isBlank(ctg)) return AbstractCloseableIterator.empty();
		final CloseableTribbleIterator<WiggleForWib> iter0=this.tabixReader.query(ctg, loc.getStart(), loc.getEnd());
		return new MyIter2(iter0,loc.getContig());
		}
	
	private void reopenSeekableStream(final String rawWibPath) throws IOException {
		final String resolvedPath = this.pathConverter.apply(rawWibPath);
		if(StringUtils.isBlank(resolvedPath)) throw new IllegalArgumentException("resolved path is blank for "+ rawWibPath);
		//open stream if needed
		if(this.wibResolvedName==null || !this.wibResolvedName.equals(resolvedPath)) {
			this.wibResolvedName = resolvedPath;
			if(this.wibSeekableStream!=null) this.wibSeekableStream.close();
			this.wibSeekableStream = this.seekableStreamFactory.getStreamFor(this.wibResolvedName);
		}
	}
	
	private class MyIter2 extends AbstractCloseableIterator<WigItem>  {
		private CloseableTribbleIterator<WiggleForWib> iter0;
		private WiggleForWib currentWiggleForWib =null;
		private final String userContig;
		byte[] buffer;
		int buffer_index =0;
		int currentGenomicPos1=0;
		MyIter2(CloseableTribbleIterator<WiggleForWib> iter0,final String userContig) throws IOException {
			this.iter0=iter0;
			this.userContig=userContig;
			}
		
		@Override
		protected WigItem advance() {
			try {
				for(;;) {
					if(this.currentWiggleForWib==null) {
						if(this.iter0==null || !this.iter0.hasNext()) {
							close();
							return null;
							}
						this.currentWiggleForWib = iter0.next();
						reopenSeekableStream(this.currentWiggleForWib.wibFile);
		                WibReader.this.wibSeekableStream.seek(this.currentWiggleForWib.offset);
		                this.buffer = new byte[this.currentWiggleForWib.count];
		                WibReader.this.wibSeekableStream.read(this.buffer);
		                this.buffer_index=0;
		                this.currentGenomicPos1=this.currentWiggleForWib.getStart();
						}
					
					if( this.buffer_index>=this.buffer.length) {
						this.currentWiggleForWib = null;
						this.buffer = null;
						continue;
						}
					final int wigStart = this.currentGenomicPos1;
	                final int b = (this.buffer[this.buffer_index]  & 0xff);
	                int len=1;
	                while(this.buffer_index+len < this.buffer.length && (this.buffer[this.buffer_index+len] & 0xff ) ==b) {
	                	len++;
	                	}
					this.buffer_index+=len;
					this.currentGenomicPos1+=len;
					
					if(b < 128) {
						double value = this.currentWiggleForWib.lowerLimit+(this.currentWiggleForWib.dataRange*((double)b/127.0));
						return new WigItemImpl(
								this.userContig,
								wigStart,
								wigStart+len-1,
								value
								);
						}		
					}
				} 
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() {
			this.buffer = null;
			currentWiggleForWib=null;
			if(this.iter0!=null) iter0.close();
			iter0=null;
			}
		
	}
	
	@Override
	public void close() {
		try {this.tabixReader.close(); } catch(Exception err) {}
		if(this.wibSeekableStream!=null) try {this.wibSeekableStream.close();} catch(Exception err) {}
		}
	}
