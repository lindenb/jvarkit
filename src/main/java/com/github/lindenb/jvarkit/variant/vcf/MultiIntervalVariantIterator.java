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
package com.github.lindenb.jvarkit.variant.vcf;

import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;

import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
 * creates a Multiple-interval iterator list
 */
public class MultiIntervalVariantIterator {
	
public static VCFIterator query(final VCFFileReader reader,List<? extends Locatable> intervals0) {
	final SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
	final List<Locatable> intervals   = LocatableUtils.mergeIntervals(intervals0);
	if(dict!=null && dict.size()>0) {
		intervals.removeIf(R->dict.getSequence(R.getContig())==null);
		Collections.sort(intervals,ContigDictComparator.createLocatableComparator(dict));
		}
	if(intervals.isEmpty()) {
		final VCFHeader h0 = reader.getFileHeader();
		return new VCFIterator() {
			@Override
			public VCFHeader getHeader() {return h0;};
			@Override
			public boolean hasNext() {
				return false;
				}
			@Override
			public VariantContext next() {
				throw new NoSuchElementException();
				}
			@Override
			public VariantContext peek() {
				return null;
				}
			public void close() {};
			};
		}
	return new Iter0(reader,intervals);
	}

	private static class Iter0 extends AbstractIterator<VariantContext>
		implements VCFIterator {
		private final VCFFileReader reader;
		private  CloseableIterator<VariantContext> delegate=null;
		private  final List<Locatable> intervals;
		private  int index=-1;
		Iter0(VCFFileReader reader,List<Locatable> intervals) {
			this.reader = reader;
			this.intervals = intervals;
			}
		@Override
		public VCFHeader getHeader() {
			return reader.getHeader();
			}
		
		@Override
		protected VariantContext advance() {
				for(;;) {
					if(delegate==null)  {
						if(index+1>= intervals.size()) return null;
						index++;
						delegate = this.reader.query(intervals.get(index));
						}
					if(!delegate.hasNext()) {
						delegate.close();
						delegate=null;
						continue;
						}
					final VariantContext curr= delegate.next();
					boolean overlap_prev=false;
					for(int j=index-1;j>=0;--j) {
						final Locatable prev = this.intervals.get(j);
						if(!prev.contigsMatch(intervals.get(index))) break;
						if(prev.overlaps(curr)) {
							overlap_prev = true;
							break;
							}
						}
					if(overlap_prev) continue;
						
					return curr;
					}
			}
		@Override
		public void close() {
			if(this.delegate!=null) {
				this.delegate.close();
				this.delegate=null;
				}
			index = intervals.size();
			}
		}
}
