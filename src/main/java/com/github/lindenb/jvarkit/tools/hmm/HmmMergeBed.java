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

*/
package com.github.lindenb.jvarkit.tools.hmm;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

@Program(name="hmmmergebed",
	description="Pour Julien. Merge hmm bed",
	keywords={"bed","merge"},
	creationDate="20191108",
	modificationDate="20191108"
	)
public class HmmMergeBed extends Launcher {
	private static final Logger LOG = Logger.build(HmmMergeBed.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference","--dict"},description=DICTIONARY_SOURCE)
	private final Path dictPath = null;
	
	private static class Base implements Locatable
		{
		final BedLine bedline;
		final int pos;
		Base(final BedLine bedline,int pos) {
			this.bedline = bedline;
			this.pos = pos;
			}
		@Override
		public String getContig() {
			return bedline.getContig();
			}
		public int getPos() {
			return pos;
			}
		@Override
		public int getStart() {
			return getPos();
			}
		@Override
		public int getEnd() {
			return getPos();
			}
		@Override
		public String toString() {
			return getContig()+":"+getPos()+" in "+this.bedline;
			}
		}
	
	private static class BedLineToBaseIterator extends AbstractIterator<Base> implements CloseableIterator<Base> {
		final BedLineIterator delegate;
		BedLine bedline = null;
		int pos;
		BedLineToBaseIterator(final BedLineIterator delegate) {
			this.delegate = delegate;
			}
		@Override
		protected Base advance() {
			for(;;) {
				if(this.bedline==null || this.pos>this.bedline.getEnd()) {
					if(!this.delegate.hasNext()) {
						this.bedline=null;
						return null;
						}
					this.bedline = this.delegate.next();
					this.pos = this.bedline.getStart();
					continue;
					}
				final Base b = new Base(this.bedline,this.pos);
				this.pos++;
				return b;
				}
			}
		@Override
		public void close() {
			this.delegate.close();
			this.bedline=null;
			}
		
		}
	private static class BedLineIterator extends AbstractIterator<BedLine> implements Closeable {
		final BedLineCodec codec = new BedLineCodec();
		private BufferedReader br ;
		private BedLineIterator(final Path path) throws IOException{
			this.br = IOUtils.openPathForBufferedReading(path);
			}
		@Override
		protected BedLine advance() {
			try { 
				for(;;) {
					if(this.br==null) return null;
					final String line = this.br.readLine();
					if(line==null) {
						this.br.close();
						return null;
						}
					if(StringUtil.isBlank(line) || BedLine.isBedHeader(line)) {
						continue;
						}
					final BedLine bed = this.codec.decode(line);
					if(bed==null) continue;
					return bed;
					}
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() {
			if(this.br!=null) {
				try {this.br.close();} catch(final IOException err) {LOG.warn(err);}
				this.br =null;
				}
			}
		}
	@Override
	public int doWork(final List<String> args) {
		PrintWriter pw = null;
		try {
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.dictPath);
			final Comparator<Base> locCompare = new ContigDictComparator(dict).createLocatableComparator();
			final List<Path> paths = IOUtils.unrollPaths(args);
			if(paths.isEmpty()) {
				LOG.error("empty input");
				return -1;
				}
			final List<CloseableIterator<Base>> baseiterators=new ArrayList<>(paths.size());
			
			
			for(final Path bedPath:paths) {
				final BedLineIterator iter1 = new BedLineIterator(bedPath);
				final BedLineToBaseIterator iter2 = new BedLineToBaseIterator(iter1);
				baseiterators.add(iter2);
				}
			final MergingIterator<Base> mergingIter = new MergingIterator<>(locCompare, baseiterators);
			final EqualRangeIterator<Base> equalIter = new EqualRangeIterator<>(mergingIter, locCompare);
			
			
			pw = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			
			while(equalIter.hasNext()) {
				final List<Base> base = equalIter.next();
				}
			equalIter.close();
			mergingIter.close();
			
			pw.flush();
			pw.close();
			
			for(CloseableIterator<?> iter:baseiterators) iter.close();
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		
		}
	
	
	}
