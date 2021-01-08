/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC

## input

input is a set of sorted BED file , or a file with the '.list' suffix containing the paths to the bed file. 

## Example:

```
$ head -n 3  ~/DC10* 
==> DC1074_UCSC.bed.txt <==
chr22	0	64000	6_Heterochrom_lowsignal
chr22	64000	65400	1_Polycomb_repressed
chr22	65400	119600	6_Heterochrom_lowsignal

==> DC1076_UCSC.bed.txt <==
chr22	0	119600	6_Heterochrom_lowsignal
chr22	119600	120000	3_Weak_promoter
chr22	120000	120200	4_Active_promoter

==> DC1077_UCSC.bed.txt <==
chr22	0	118400	7_Heterochrom_lowsignal
chr22	118400	119600	1_Polycomb_repressed
chr22	119600	120200	2_Poised_promoter

==> DC1082_UCSC.bed.txt <==
chr22	0	119600	6_Heterochrom_lowsignal
chr22	119600	120200	3_Weak_promoter
chr22	120200	122200	6_Heterochrom_lowsignal


$ java -jar dist/hmmmergebed.jar -n 3 ~/DC10*.txt 
chr22	0	64000	6_Heterochrom_lowsignal
#chr22	64000	65400	1_Polycomb_repressed:1;6_Heterochrom_lowsignal:2;7_Heterochrom_lowsignal:1
chr22	65400	119600	6_Heterochrom_lowsignal
chr22	119600	120000	3_Weak_promoter
#chr22	120000	120200	2_Poised_promoter:1;3_Weak_promoter:1;4_Active_promoter:2
#chr22	120200	120400	1_Polycomb_repressed:1;3_Weak_promoter:1;4_Active_promoter:1;6_Heterochrom_lowsignal:1
chr22	120400	122000	6_Heterochrom_lowsignal
chr22	122000	122600	3_Weak_promoter
#chr22	122600	122800	3_Weak_promoter:2;6_Heterochrom_lowsignal:2
chr22	122800	130000	6_Heterochrom_lowsignal
(...)

```

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
	@Parameter(names={"-R","--reference","--dict"},description="optional dictionary if data are sorted on a REF order." +DICTIONARY_SOURCE)
	private Path dictPath = null;
	@Parameter(names={"-n","--num"},description="dictionary if data are sorted on a REF order." , required=true)
	private int numTreshold = -1;
	@Parameter(names={"-u"},description="Hide invalid lines")
	private boolean hide_invalid = false;

	/** one base in a BED record */
	private static class Base implements Locatable
		{
		final BedLine bedline;
		final int pos;
		Base(final BedLine bedline,int pos) {
			this.bedline = bedline;
			this.pos = pos;
			}
		String getOpcode() {
			return this.bedline.get(3);
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
	
	/** iterator over all the bases of a bed record */
	private static class BedLineToBaseIterator 
		extends AbstractCloseableIterator<Base> {
		final BedLineIterator delegate;
		/* current bed line */
		BedLine bedline = null;
		/* curren pos in bedline */
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
	
	private static class BedLineIterator extends AbstractCloseableIterator<BedLine> {
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
					if(bed.getColumnCount()<4) throw new IOException("Column $4 missing in "+bed+" in "+this.path);
					if(StringUtils.isBlank(bed.get(3)))  throw new IOException("Column $4 empty in "+bed+" in "+this.path);
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
		final String contig;
		final int start;
		int end;
		final String opcode ;
		final boolean valid_flag;
		
		BaseMerge(final List<Base> bases) {
			final Base first = bases.get(0);
			this.contig = first.getContig();
			this.start = first.getPos();
			this.end  = start;
			final Counter<String> count = new  Counter<>();
			for(final Base b:bases) {
				count.incr(b.getOpcode());
				}
			
			final String freq = count.getMostFrequent();
			final long num =  count.count(freq);
			if( num>= HmmMergeBed.this.numTreshold ) {
				//ambiguity on best: check other with the sam occurence
				if(count.keySet().
						stream().
						filter(K->count.count(K)==num).
						count()==1L) {
					this.opcode= freq;
					valid_flag= true;
					return;
					}
				}
			valid_flag = false;
			this.opcode= count.
				stream().
				sorted((A,B)->A.getKey().compareTo(B.getKey())).
				map(KV->KV.getKey()+":"+KV.getValue()).
				collect(Collectors.joining(";"));
			}
		
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
		boolean isValid() {
			
			return this.valid_flag;
			}
		
		String getOpcode() {
			return this.opcode;
			}
		}
	
	private class BaseMergerIterator extends AbstractCloseableIterator<BaseMerge>  {
		private final List<CloseableIterator<Base>> baseiterators;
		private final MergingIterator<Base> mergingIter;
		private final EqualRangeIterator<Base> equalIter;
		BaseMergerIterator(final List<CloseableIterator<Base>> baseiterators,final Comparator<Base> locCompare) {
			this.baseiterators = baseiterators ;
			this.mergingIter = new MergingIterator<>(locCompare, baseiterators);
			this.equalIter = new EqualRangeIterator<>(mergingIter, locCompare);
 			}
		
		
		
		@Override
		protected BaseMerge advance() {
			if(!this.equalIter.hasNext()) return null;
			final List<Base> L = equalIter.next();
			final BaseMerge bm = new BaseMerge(L);
			
			while(this.equalIter.hasNext()) {
				final List<Base> L2 = equalIter.peek();
				final BaseMerge bm2 = new BaseMerge(L2);
				
				if(!bm2.getContig().equals(bm.contig)) break;
				if(bm.end+1 != bm2.getStart()) break;
				if(!bm.getOpcode().equals(bm2.getOpcode())) break;
				this.equalIter.next();//consumme
				bm.end = bm2.getStart();
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
			if(this.numTreshold<1) {
				LOG.error("bad treshold.");
				return -1;
			}
			
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
			if(this.numTreshold>paths.size()) {
				LOG.error("bad treshold (> number of files).");
				return -1;
			}
			
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
				if(!bm.isValid()) {
					if(this.hide_invalid) continue;
					pw.print("#");
					}
				pw.print(bm.contig);
				pw.print('\t');
				pw.print(bm.start-1);
				pw.print('\t');
				pw.print(bm.end);
				pw.print('\t');
				pw.print(bm.getOpcode());
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
	
	public static void main(final String[] args) {
		new HmmMergeBed().instanceMainWithExit(args);
		}
	}
