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
package com.github.lindenb.jvarkit.samtools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

/**
 * Compute coverage over a segment in a SAM file
 */
public class CoverageFactory
	{
	public enum ScaleType {AVERAGE,MEDIAN};
	private int mappingQuality = 1;
	private Predicate<SAMRecord> samRecordFilter = R->SAMRecordDefaultFilter.accept(R);
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	private boolean useClip=false;
	public CoverageFactory() {
		}
	/** set min mapping quality */
	public CoverageFactory setMappingQuality(int mappingQuality) {
		this.mappingQuality = mappingQuality;
		return this;
		}
	/** set partition */
	public CoverageFactory setPartition(final SAMRecordPartition partition) {
		this.partition = partition;
		return this;
		}
	/** set use clipping */
	public CoverageFactory setUseClipping(final boolean useClip) {
		this.useClip = useClip;
		return this;
		}
	/** set function to accept read */
	public CoverageFactory setRecordFilter(final Predicate<SAMRecord> accept) {
		this.samRecordFilter = accept;
		return this;
		}

	 static interface BaseCoverage  {
		public OptionalDouble getMedian();
		public OptionalInt getMax();
		public IntStream stream();
		}
	
	private static interface MultiCoverage extends BaseCoverage,Iterable<SimpleCoverage> {
		public List<SimpleCoverage> getItems();
		@Override
		public default Iterator<SimpleCoverage> iterator() {
			return getItems().iterator();
			}
		}
	
	public static interface SimpleCoverage extends BaseCoverage,Locatable {
		public int[] toIntArray();
		
		/** reduce coverage to this loc. Loc **must** be contained in this coverage */
		public SimpleCoverage getSubCoverage(final Locatable loc);
		
		/** I should have called it toIntArray */
		@Deprecated
		public default int[] toByteArray() {
			return toIntArray();
			}
		public int get(int array_index);
		public OptionalDouble getMedian(final Locatable loc);
		public OptionalDouble getAverage(final Locatable loc);
		public OptionalInt getMax(final Locatable loc);
		public OptionalDouble getMedian();
		public OptionalDouble getAverage();
		public OptionalInt getMax();
		public default double[] scale(final ScaleType scaleType,int length) {
			switch(scaleType) {
				case MEDIAN: return scaleMedian(length);
				case AVERAGE: return scaleAverage(length);
				default: throw new IllegalArgumentException();
				}
			}
		public double[] scaleMedian(int length);
		public double[] scaleAverage(int length);
		}
	
	private static class MultiCoverageImpl implements MultiCoverage {
		private final IntervalTreeMap<SimpleCoverage> treeMap = new IntervalTreeMap<>();
		private final List<SimpleCoverage> items = new ArrayList<>();
		@Override
		public List<SimpleCoverage> getItems() {
			return items;
			}
		@Override
		public OptionalDouble getMedian() {
			final DiscreteMedian<Integer> med = new DiscreteMedian<>();
			for(SimpleCoverage cov: items) {
				final int[] array = cov.toByteArray();
				for(int i=0;i< array.length;i++) {
					med.add(array[i]);
					}
				}
			return med.getMedian();
			}
		@Override
		public IntStream stream() {
			return items.stream().flatMapToInt(L->Arrays.stream(L.toByteArray()));
			}
		@Override
		public OptionalInt getMax() {
			return stream().max();
			}
		}
	
	private static class SimpleCoverageImpl implements SimpleCoverage {
		final String contig;
		final int start1;
		final int[] coverage;
		SimpleCoverageImpl(final Locatable loc) {
			this.contig = loc.getContig();
			this.start1 = loc.getStart();
			this.coverage = new int[loc.getLengthOnReference()];
			Arrays.fill(this.coverage, 0);
			}
		
		public SimpleCoverage getSubCoverage(final Locatable loc) {
			if(loc==null) throw new IllegalArgumentException();
			if(!this.contains(loc)) throw new IllegalArgumentException(this.toString()+" doesn't contain "+loc.toString());
			final SimpleCoverageImpl another= new SimpleCoverageImpl(loc);
			System.arraycopy(this.coverage, loc.getStart() - this.getStart(), another.coverage,0,another.coverage.length);
			return another;
			}
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof SimpleCoverageImpl)) return false;
			final SimpleCoverageImpl other = SimpleCoverageImpl.class.cast(obj);
			return this.contigsMatch(other) &&
					this.getStart() == other.getStart() &&
					this.getEnd() == other.getEnd()
					;
			}
		@Override
		public int hashCode()
			{
			return (contig.hashCode() * 31 + Integer.hashCode(this.start1) ) * 31 + Integer.hashCode(this.coverage.length);
			}
		
		@Override
		public double[] scaleMedian(int length) {
			final double[] array = new double[length];
			final DiscreteMedian<Integer> med = new DiscreteMedian<>();
			for(int i=0;i< array.length;i++) {
				int idx0 = (int)(((i+0)/(double)length)*this.coverage.length);
				int idx1 = (int)(((i+1)/(double)length)*this.coverage.length);
				if(idx0==idx1) idx1++;
				while(idx0<idx1 && idx0 < this.coverage.length) {
					med.add(this.coverage[idx0]);
					idx0++;
					}
				array[i]=med.getMedian().orElse(0.0);
				med.clear();
				}
			return array;
			}
		@Override
		public double[] scaleAverage(int length) {
			final double[] array = new double[length];
			for(int i=0;i< array.length;i++) {
				int idx0 = (int)(((i+0)/(double)length)*this.coverage.length);
				int idx1 = (int)(((i+1)/(double)length)*this.coverage.length);
				if(idx0==idx1) idx1++;
				int N=0;
				double t=0.0;
				while(idx0<idx1 && idx0 < this.coverage.length) {
					t+=this.coverage[idx0];
					idx0++;
					N++;
					}
				array[i]=(N==0?0:t/N);
				}
			return array;
			}
		
		@Override
		public String getContig() {
			return contig;
			}
		@Override
		public int getStart() {
			return this.start1;
			}
		@Override
		public int getEnd() {
			return this.start1 + this.coverage.length -1;
			}
		@Override
		public OptionalDouble getMedian(final Locatable loc) {
			if(!this.contains(loc)) throw new IllegalArgumentException(loc.toString()+" is not contained in "+this.toString());
			final DiscreteMedian<Integer> med = new DiscreteMedian<>();
			final int i0 = loc.getStart() - this.getStart();
			final int L = loc.getLengthOnReference();
			for(int i=0;i< L;i++) {
				med.add(this.coverage[i0 + i]);
				}
			return med.getMedian(); 
			}
		
		@Override
		public OptionalDouble getAverage(final Locatable loc) {
			if(!this.contains(loc)) throw new IllegalArgumentException(loc.toString()+" is not contained in "+this.toString());
			long count=0L;
			double sum=0.0;
			final int i0 = loc.getStart() - this.getStart();
			final int L = loc.getLengthOnReference();
			for(int i=0;i< L;i++) {
				sum+= this.coverage[i0 + i];
				count++;
				}
			return count==0L?OptionalDouble.empty():OptionalDouble.of(sum/count);
			}
		
		@Override
		public OptionalInt getMax(final Locatable loc) {
			if(!this.contains(loc)) throw new IllegalArgumentException(loc.toString()+" is not contained in "+this.toString());
			return Arrays.stream(this.coverage,
					loc.getStart() - this.getStart(),
					1+ (loc.getEnd() - this.getStart())
					).max();
			}

		@Override
		public OptionalDouble getMedian()  {
			return getMedian(this);
			}
		@Override
		public OptionalDouble getAverage()  {
			return getAverage(this);
			}

		@Override
		public OptionalInt getMax() {
			return getMax(this);
			}

		 @Override
		public int get(int i) {
			return this.coverage[i];
			}
		@Override
		public IntStream stream() {
			return Arrays.stream(this.coverage);
			}
		@Override
		public int[] toIntArray() {
			return this.coverage;
			}
		@Override
		public String toString() {
			return getContig()+":"+getStart()+"-"+getEnd();
			}
		}
			
	
	private MultiCoverage getMultiCoverage(final SamReader reader,final Collection<Locatable> locCollection,final String sample) {
		if(reader==null) throw new IllegalArgumentException("reader==null");
		if(!reader.hasIndex()) throw new IllegalArgumentException("SamReader is not indexed. "+reader.getResourceDescription());
		if(locCollection==null || locCollection.isEmpty()) throw new IllegalArgumentException("no interval provided");
		final MultiCoverageImpl multi = new MultiCoverageImpl();
		
		for(final Locatable loc: locCollection) {
			final SimpleCoverage cov = getSimpleCoverage(reader,loc,sample);
			multi.items.add(cov);
			final Interval rgn = new Interval(cov);
			multi.treeMap.put(rgn, cov);
			}
		return multi;
		}	
	
	public SimpleCoverage getSimpleCoverage(final SamReader reader,final Locatable loc,final String sample) {
		return getSimpleCoverage(reader,Collections.singletonList(loc),sample);
	}
	
	/** 
	 * 
	 * @param reader 
	 * @param locCollection collection of items. must not be empty, all position must belong to the sampel contig
	 * @param sample
	 * @return
	 */
	public SimpleCoverage getSimpleCoverage(final SamReader reader,final Collection<? extends Locatable> locCollection,final String sample) {
		if(reader==null) throw new IllegalArgumentException("reader==null");
		if(!reader.hasIndex()) throw new IllegalArgumentException("SamReader is not indexed. "+reader.getResourceDescription());
		if(locCollection.isEmpty()) throw new IllegalArgumentException("Empty Query");
		final SAMFileHeader header= reader.getFileHeader();
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
		final Set<String> contigs = locCollection.stream().map(L->L.getContig()).collect(Collectors.toSet());
		if(contigs.size()!=1)  throw new IllegalArgumentException("Expected only one chromosome but got:"+String.join(",", contigs));
		final String contig = contigs.iterator().next();
		
		final int chromStart = locCollection.stream().mapToInt(L->L.getStart()).min().orElse(0);
		final int chromEnd = locCollection.stream().mapToInt(L->L.getEnd()).max().orElse(0);
						
		
		final Locatable loc = new SimpleInterval(contig,chromStart,chromEnd);
		final SimpleCoverageImpl cov = new SimpleCoverageImpl(loc);
		final int tid = dict.getSequenceIndex(contig);
		if(tid<0) {
			return cov;
		}
		
		if(!StringUtils.isBlank(sample) &&
				header.getReadGroups().stream().
				map(RG->this.partition.apply(RG)).
				noneMatch(S->sample.equals(S))) {
			return cov;
			}

		int idx=0;
		QueryInterval[] array = new QueryInterval[locCollection.size()];
		for(final Locatable item: locCollection) {
			array[idx] = new QueryInterval(tid, item.getStart(), item.getEnd());
			idx++;
		}
		array  = QueryInterval.optimizeIntervals(array);
		
		try(CloseableIterator<SAMRecord> iter = reader.query(array, false) ) {
			while(iter.hasNext()) {
				final SAMRecord rec = iter.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.getMappingQuality()< this.mappingQuality) continue;
				if(!this.samRecordFilter.test(rec)) continue;
				final Cigar cigar = rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue; 
				if(!StringUtils.isBlank(sample) && !SAMRecordPartition.any.equals(this.partition)) {
					if(!sample.equals(this.partition.getPartion(rec))) continue;
					}
				
				int maxEnd = loc.getEnd();
				if(rec.getReadPairedFlag() &&
					!rec.getMateUnmappedFlag() &&
					rec.getReferenceIndex().equals(rec.getMateReferenceIndex()))
					{
					int mateStart;
					
					if(this.useClip && SAMUtils.getMateCigar(rec)!=null) {
						mateStart = SAMUtils.getMateUnclippedStart(rec);
						}
					else
						{
						mateStart = rec.getMateAlignmentStart();
						}
					if( rec.getStart() < mateStart &&
						rec.getEnd() > mateStart)
						{		
						maxEnd = Math.min(maxEnd, mateStart );
						}
					}
				
				if(!this.useClip) {
					for(AlignmentBlock block:rec.getAlignmentBlocks()) {
						if(block.getReferenceStart() > maxEnd ) break;
						for(int t=0;t< block.getLength();t++) {
							final int pos1 = block.getReferenceStart() + t;
							if(pos1 < loc.getStart()) continue;
							if(pos1 > maxEnd ) break;
							cov.coverage[pos1-loc.getStart()]++;
							}
						}
					}
				else
					{
					int ref1= rec.getUnclippedStart();
					for(final CigarElement ce:cigar) {
						if(ref1> maxEnd ) break;
						final CigarOperator op =ce.getOperator();
						switch(op) {
							case P: break;
							case N: case D: ref1+=ce.getLength(); break;
							case I: break;
							case H: case S: case X: case EQ: case M:
								{
								for(int i=0;i< ce.getLength();i++) {
									final int pos1 = ref1+i;
									if(pos1 < loc.getStart()) continue;
									if(pos1 > maxEnd ) break;
									cov.coverage[pos1-loc.getStart()]++;
									}
								ref1+=ce.getLength();
								break;
								}
							default: throw new IllegalStateException(op.name());
							}
						}
					}
				}
			}
		return cov;
		}
	}
