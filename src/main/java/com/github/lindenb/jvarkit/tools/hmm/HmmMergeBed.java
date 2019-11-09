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

/**
BEGIN_DOC

## input

input is a set of sorted BED file , or a file with the '.list' suffix containing the paths to the bed file. 

END_DOC
*/
@Program(name="hmmmergebed",
	description="For @154sns. Merge hmm bed",
	keywords={"bed","merge"},
	creationDate="20191108",
	modificationDate="20191109"
	)
public class HmmMergeBed extends Launcher {
	private static final Logger LOG = Logger.build(HmmMergeBed.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference","--dict"},description="dictionary if data are sorted on a REF order." +DICTIONARY_SOURCE)
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
		private BedLine prevRec = null;
		private BufferedReader br ;
		private final Path path;
		private BedLineIterator(final Path path) throws IOException {
			this.path = path;
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
					if(prevRec!=null && prevRec.overlaps(bed)) {
						throw new IOException("got two bed overlapping "+bed+" "+prevRec+" in "+this.path);
						}
					prevRec=bed;
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
	
	private class BaseMerge implements Locatable
		{
		String contig;
		int start;
		int end;
		String opcode;
		
		@Override
		public String getContig() {
			return contig;
			}
		@Override
		public int getStart() {
			return start;
			}
		@Override
		public int getEnd() {
			return end;
			}
		}
	private class BaseMergerIterator extends AbstractIterator<BaseMerge> implements CloseableIterator<BaseMerge> {
		private final List<CloseableIterator<Base>> baseiterators;
		private final MergingIterator<Base> mergingIter;
		private final EqualRangeIterator<Base> equalIter;
		BaseMergerIterator(final List<CloseableIterator<Base>> baseiterators,final Comparator<Base> locCompare) {
			this.baseiterators = baseiterators ;
			this.mergingIter = new MergingIterator<>(locCompare, baseiterators);
			this.equalIter = new EqualRangeIterator<>(mergingIter, locCompare);
 			}
		String getOpcode(final  List<Base> s) {
			return "*";
			}
		
		@Override
		protected BaseMerge advance() {
			if(!equalIter.hasNext()) return null;
			final List<Base> L = equalIter.next();
			final Base first = L.get(0);
			final BaseMerge bm = new BaseMerge();
			bm.contig = first.getContig();
			bm.start = first.getPos();
			bm.end = first.getPos();
			bm.opcode = getOpcode(L);
			while(equalIter.hasNext()) {
				final List<Base> L2 = equalIter.peek();
				final Base second = L2.get(0);
				if(!second.getContig().equals(bm.contig)) break;
				if(bm.end+1 != second.getPos()) break;
				if(!bm.opcode.equals(getOpcode(L2))) break;
				equalIter.next();//consumme
				bm.end = second.getPos();
				}
			return bm;
			}
		@Override
		public void close() {
			this.equalIter.close();
			this.mergingIter.close();
			for(CloseableIterator<?> iter:baseiterators) iter.close();
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter pw = null;
		try {
			final Comparator<Base> locCompare;
			
			if(this.dictPath!=null) {
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.dictPath);
				locCompare = new ContigDictComparator(dict).createLocatableComparator();
				}
			else
				{
				locCompare = (A,B)->{
					int i= A.getContig().compareTo(B.getContig());
					if(i!=0) return i;
					i = Integer.compare(A.getStart(),B.getStart());
					if(i!=0) return i;
					return Integer.compare(A.getEnd(),B.getEnd());
					};
				}
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
			
			
			pw = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			
			final BaseMergerIterator bmiter = new BaseMergerIterator(baseiterators, locCompare);
			
			while(bmiter.hasNext()) {
				final BaseMerge bm = bmiter.next();
				pw.print(bm.contig);
				pw.print('\t');
				pw.print(bm.start-1);
				pw.print('\t');
				pw.print(bm.end);
				pw.print('\t');
				pw.print(bm.opcode);
				pw.println();
				}
			bmiter.close();
			
			pw.flush();
			pw.close();
			
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		
		}
	
	
	}
