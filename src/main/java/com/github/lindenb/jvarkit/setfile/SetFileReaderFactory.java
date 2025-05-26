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
package com.github.lindenb.jvarkit.setfile;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.locatable.SimplePosition;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;

/* 'SetFile' Reader for rvtest */
public class SetFileReaderFactory {
	private final Function<String, String> converter;
	private final Comparator<Locatable> comparator ;
	private boolean skipUnknownChromosome = false;
	private boolean skipEmptyRecord = false;
	private boolean mergeOverlappingRecords = false;
	private boolean allowOverlappingRecords = false;
	public SetFileReaderFactory() {
		this(null);
	}
	public SetFileReaderFactory(final SAMSequenceDictionary dict) {
		if(dict==null) {
			this.converter = S->S;
			this.comparator = (A,B)->{
				int i = A.getContig().compareTo(B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				return Integer.compare(A.getEnd(), B.getEnd());
				};
			}
		else
			{
			this.converter = ContigNameConverter.fromOneDictionary(dict);
			this.comparator= new ContigDictComparator(dict).createLocatableComparator();
			}
		}
	
	public SetFileReaderFactory setAllowOverlappingRecords(boolean allowOverlappingRecords) {
		this.allowOverlappingRecords = allowOverlappingRecords;
		return this;
		}
	
	public SetFileReaderFactory setMergeOverlappingRecords(boolean mergeOverlappingRecords) {
		this.mergeOverlappingRecords = mergeOverlappingRecords;
		return this;
		}
		
	public SetFileReaderFactory setSkipEmptyRecords(boolean skipEmptyRecord) {
		this.skipEmptyRecord = skipEmptyRecord;
		return this;
		}
	
	public SetFileReaderFactory setSkipUnknownChromosomes(boolean skipUnknownChromosome) {
		this.skipUnknownChromosome = skipUnknownChromosome;
		return this;
		}
	
	public CloseableIterator<SetFileRecord> open(final BufferedReader br) throws IOException {
		return new MyIter(br);
		}
	public CloseableIterator<SetFileRecord> open(final Path path) throws IOException {
		return open(IOUtils.openPathForBufferedReading(path));
		}
	

	
	private class MyIter extends AbstractCloseableIterator<SetFileRecord> {
		private final BufferedReader br;
		MyIter(final BufferedReader br) {
			this.br = br;
			}
		@Override
		protected SetFileRecord advance() {
			try {
				for(;;) {
					final String line  = this.br.readLine();
					if(line==null) return null;
					if(StringUtils.isBlank(line)) continue;
					if(line.startsWith("#")) continue;
					String[] tokens = line.split("[\\s]+", 2);
					if(StringUtils.isBlank(tokens[0])) throw new IOException("empty name ("+tokens[0]+") in "+line);
					if(tokens.length!=2) throw new IOException("only name in "+line);
					String name = tokens[0];
					tokens = CharSplitter.COMMA.split(tokens[1]);
					final List<Locatable> L = new ArrayList<>(tokens.length);
					for(int i=0;i< tokens.length;i++) {
						String s = tokens[i].trim();
						if(StringUtils.isBlank(s)) continue;
						final int colon = s.indexOf(':');
						if(colon<=0) throw new IOException("bad chromosome for item["+i+"] in line "+line );
						final String contig = converter.apply(s.substring(0,colon));
						if(StringUtils.isBlank(contig)) {
							if(skipUnknownChromosome) continue;
							throw new IOException("undefined chromosome for item["+i+"] in line"+line);
						}
						final int  hyphen = s.indexOf('-', colon+1);
						if(hyphen==-1) {
							final int pos = Integer.parseInt(s.substring(colon+1));
							L.add(new SimplePosition(contig, pos));
							}
						else {
							final int start = Integer.parseInt(s.substring(colon+1,hyphen));
							final int end = Integer.parseInt(s.substring(hyphen+1));
							if(start>end) throw new IOException("bad start>end for item["+i+"] in line");
							L.add(new SimpleInterval(contig, start, end));
							}
						}
					if(L.isEmpty()) {
						if(skipEmptyRecord)continue;
						throw new IOException("all items skipped in line "+line);
						}
					Collections.sort(L, comparator);
					int i=0;
					while(i+1 < L.size()) {
						final Locatable xi = L.get(i  );
						final Locatable xj = L.get(i+1);
						if(xi.overlaps(xj)) {
							if(!allowOverlappingRecords) {
								throw new IOException("two items overlapping "+xi+" and "+ xj+"in line "+line);
								}
							if(mergeOverlappingRecords) {
								L.set(i, new SimpleInterval(
										xi.getContig(),
										Math.min(xi.getStart(),xj.getStart()),
										Math.max(xi.getEnd(),xj.getEnd())
										));
								L.remove(i+1);
								}
						} else {
							i++;
						}
					}
					
					
					return SetFileRecord.create(name, Collections.unmodifiableList(L));
					}
				} 
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() {
			try {br.close();}
			catch(final IOException err) {}
			}
		}
}
