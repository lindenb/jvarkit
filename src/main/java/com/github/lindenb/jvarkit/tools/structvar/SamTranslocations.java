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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.GenomicIndexUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
/**

BEGIN_DOC

## History

* 2018-09-18 :  rewriting
* 2017-12-13 :  refactoring for balanced translocation.

END_DOC
*/
@Program(name="samtranslocations",
	description="Explore balanced translocations between two chromosomes using discordant paired-end reads.",
	keywords={"sam","bam"},
	generate_doc=false
	)
public class SamTranslocations extends Launcher {
	private static final Logger LOG = Logger.build(SamTranslocations.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-d","--distance"},description="Max distance between two read to test if they both end at the same ~ position.")
	private int fuzzy_distance = 100;
	@Parameter(names={"-t","--trans"},description="Allow same-contig variants")
	private boolean allowSameContig=false;
	@Parameter(names={"-m","--min"},description="Min number of events")
	private int min_number_of_events=0;
	@Parameter(names={"-B","--bed"},description="Optional BED file. SV should overlap this bed.")
	private File bedFile = null;

	private class BinMap
		{
		private final Map<Integer, Set<BreakPointDigest>> bin2set = new HashMap<>(GenomicIndexUtil.MAX_BINS);
		
		private int bin(final BreakPointDigest bp)
			{
			return GenomicIndexUtil.regionToBin(
					start1(bp),
					end1(bp)
					);
			}
		
		private int start1(final BreakPointDigest bp)
			{
			return Math.max(1,(bp.start1-1)-fuzzy_distance);
			}
		private int end1(final BreakPointDigest bp) {
			return (bp.end1)+fuzzy_distance;
			}
		
		List<BreakPointDigest> toList() {
			return this.bin2set.
					values().
					stream().
					flatMap(S->S.stream().filter(P->P.count>=min_number_of_events)).
					sorted((A,B)->{
						int d = Integer.compare(A.tid1,B.tid1);
						if(d!=0) return d;
						d = Integer.compare(A.tid2,B.tid2);
						if(d!=0) return d;
						d = Integer.compare(A.start1,B.start1);
						if(d!=0) return d;
						return Integer.compare(A.start2,B.start2);
						}).
					collect(Collectors.toList());
			}
		
		void put(BreakPointDigest bp)
			{
			boolean done=false;
			while(!done)
				{
				done = true;
				final int b= bin(bp);
				if(b>=GenomicIndexUtil.MAX_BINS) throw new IllegalStateException();
				final BitSet bitset = GenomicIndexUtil.regionToBins(start1(bp), end1(bp));
				for(int x=0;x<GenomicIndexUtil.MAX_BINS;++x)
					{
					if(!bitset.get(x)) continue;
					final Set<BreakPointDigest> set  = this.bin2set.get(b);
					if(set==null) continue;
					BreakPointDigest overlapper = null;
					for(BreakPointDigest other : set) {
						if(bp.tid1!=other.tid1) throw new IllegalStateException();
						if(bp.tid2!=other.tid2) continue;
						if( overlap( bp.start1  ,bp.end1 , other.start1, other.end1 ) &&
							overlap( bp.start2  ,bp.end2 , other.start2, other.end2 ))
							{
							overlapper = other;
							break;
							}
						}
					if(overlapper!=null)
						{
						//System.err.println("merging "+bp+" "+overlapper);
						set.remove(overlapper);
						//
						bp.start1 = Math.min(bp.start1, overlapper.start1);
						bp.end1 = Math.max(bp.end1, overlapper.end1);
						//
						bp.start2 = Math.min(bp.start2, overlapper.start2);
						bp.end2 = Math.max(bp.end2, overlapper.end2);
						//
						bp.count += overlapper.count;
						done=false;
						break;
						}
					}
				}
			final int b= bin(bp);
			Set<BreakPointDigest> set  = this.bin2set.get(b);
			if(set==null) {
				set = new HashSet<>();
				this.bin2set.put(b,set);
				}
			set.add(bp);
			}
		
		}

	private static class BreakPoint
		{
		final int tid1;
		final int pos1;
		final int tid2;
		final int pos2;
		BreakPoint(final int tid1,int pos1,int tid2,int pos2) {
			if(tid1==tid2)
				{
				this.tid1 = tid1;
				this.tid2 = tid2;
				this.pos1 = pos1<pos2?pos1:pos2;
				this.pos2 = pos1<pos2?pos2:pos1;
				}
			else if(tid1<tid2)
				{
				this.tid1 = tid1;
				this.pos1 = pos1;
				this.tid2 = tid2;
				this.pos2 = pos2;
				}
			else
				{
				this.tid1 = tid2;
				this.pos1 = pos2;
				this.tid2 = tid1;
				this.pos2 = pos1;
				}
			}
		}
	private static long ID_GENERATOR = 0L;
	private static class BreakPointDigest
		{
		final long id = (++ID_GENERATOR);
		final int tid1;
		int start1 = 0;
		int end1 = 0;
		
		final int tid2;
		int start2 = 0;
		int end2 = 0;

		
		int count=0;
		BreakPointDigest(final int tid1,final int tid2) {
			this.tid1 = tid1;
			this.tid2 = tid2;
			}
		@Override
		public int hashCode() {
			return Long.hashCode(this.id);
			}
		@Override
		public boolean equals(final Object obj) {
			return this.id == BreakPointDigest.class.cast(obj).id;
			}
		@Override
		public String toString() {
			return tid1+":"+start1+"-"+end1+" -> "+tid2+":"+start2+"-"+end2;
			}
		}
	
	private boolean overlap(final int start, final int end, final int start2, final int end2) {
		return CoordMath.overlaps(
			start-this.fuzzy_distance, end+this.fuzzy_distance,
			start2-this.fuzzy_distance, end2+this.fuzzy_distance
			);
	}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.fuzzy_distance<=0) {
			LOG.error("fuzzy_distance <=0 ("+fuzzy_distance+")");
			return -1;
		}
		PrintWriter pw = null;
		SamReader samReader = null;
		SAMRecordIterator iter = null;
		final short SA = SAMTag.SA.getBinaryTag();
		try {
			String inputName = oneFileOrNull(args);
			samReader = super.openSamReader(inputName);
			final SAMFileHeader header = samReader.getFileHeader();
			
			final String sampleName = header.
					getReadGroups().
					stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtil.isBlank(S)).
					findFirst().
					orElse(StringUtil.isBlank(inputName)?"<BAM>":inputName);
			
			final SAMSequenceDictionary refDict = header.getSequenceDictionary();
			if(refDict==null ||refDict.isEmpty()) {
				LOG.error("Not enough contigs in sequence dictionary.");
				return -1;
			}
			final boolean allowedContigs[] = new boolean[refDict.size()];
			for(int tid=0;tid < refDict.size();++tid)
				{
				final String refName = refDict.getSequence(tid).getSequenceName();
				allowedContigs[tid]= true;
				if(!refName.matches("(chr)?([1-9][0-9]*|X|Y)"))
					{
					allowedContigs[tid]= false;
					}
				}
			final IntervalTreeMap<Boolean> bedIntervals;
			if(this.bedFile!=null)
				{
				java.io.BufferedReader r=null;
				try
					{
					final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(refDict);
					bedIntervals = new IntervalTreeMap<Boolean>();
					r=com.github.lindenb.jvarkit.io.IOUtils.openFileForBufferedReading(this.bedFile);
					final BedLineCodec bedCodec=new BedLineCodec();
					r.lines().
						filter(line->!(line.startsWith("#") ||  com.github.lindenb.jvarkit.util.bio.bed.BedLine.isBedHeader(line) ||  line.isEmpty())).
						map(line->bedCodec.decode(line)).
						filter(B->B!=null).
						map(B->B.toInterval()).
						filter(L->L.getStart()<L.getEnd()).
						forEach(B->{
							final String c = contigNameConverter.apply(B.getContig());
							if(c==null) return;
							int tid = refDict.getSequenceIndex(c);
							if(tid==-1 || !allowedContigs[tid]) return;
							bedIntervals.put(new Interval(c,B.getStart(),B.getEnd()),true);							
							});						
					}
				finally
					{
					htsjdk.samtools.util.CloserUtil.close(r);
					}
				}
			else
				{
				bedIntervals=null;
				}
			iter = samReader.iterator();
			
			final List<BreakPoint> breakPoints = new ArrayList<>(10_000_000);
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.
					newInstance().
					dictionary(refDict).
					logger(LOG).
					build();
			
			final Consumer<BreakPoint> addBreakPoint = (BP)->{				
				breakPoints.add(BP);
				if(breakPoints.size()%100_000==0) {
					LOG.warn("Breakpoints:"+breakPoints.size());
					}
			};
			
			while(iter.hasNext()) {
				final SAMRecord rec = progress.apply(iter.next());
				if(rec.getReadUnmappedFlag())  continue;
				if(rec.getDuplicateReadFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				final int tid1 = rec.getReferenceIndex();
				if(!allowedContigs[tid1]) continue;
				if(rec.getReadPairedFlag() && 
					!rec.getMateUnmappedFlag()
					)
					{
					final int tid2 = rec.getMateReferenceIndex();
					if(allowedContigs[tid2] && (this.allowSameContig || tid1!=tid2))
						{
						addBreakPoint.accept(new BreakPoint(
								tid1,
								rec.getReadNegativeStrandFlag()?rec.getEnd():rec.getStart(),
								tid2,
								rec.getMateAlignmentStart()
								));
						}
					}
				if(rec.getAttribute(SA)!=null) {
					for(final SAMRecord sup:SAMUtils.getOtherCanonicalAlignments(rec))
						{
						final int tid2 = sup.getReferenceIndex();
						if(!allowedContigs[tid2]) continue;
						if(!this.allowSameContig && tid1==tid2) continue;
						addBreakPoint.accept(new BreakPoint(
								tid1,
								rec.getReadNegativeStrandFlag()?rec.getEnd():rec.getStart(),
								tid2,
								sup.getAlignmentStart()
								));
						
						}
					}
				}
			progress.close();
			samReader.close();
			
			LOG.info("Done. Computing output");
			pw = openFileOrStdoutAsPrintWriter(this.outputFile);
			
			pw.print("#chrom1");
			pw.print('\t');
			pw.print("chrom1-start");
			pw.print('\t');
			pw.print("chrom1-end");
			pw.print('\t');
			pw.print("chrom2");
			pw.print('\t');
			pw.print("chrom2-start");
			pw.print('\t');
			pw.print("chrom2-end");
			pw.print('\t');
			pw.print("sample");
			pw.print('\t');
			pw.print("count");
			pw.println();

			
			for(int tid1=0;tid1 < refDict.size();++tid1)
				{
				if(!allowedContigs[tid1]) continue;
				final String ctg1 = refDict.getSequence(tid1).getSequenceName();
				for(int tid2=tid1 + (this.allowSameContig?0:1);tid2< refDict.size();++tid2)
					{
					if(!allowedContigs[tid2]) continue;
					final String ctg2 = refDict.getSequence(tid2).getSequenceName();
					final BinMap binMap = new BinMap();
					
					
					
					int x=0;
					while(x<breakPoints.size())
						{
						final BreakPoint bp = breakPoints.get(x);
						if(bp.tid1==tid1 && bp.tid2==tid2)
							{
							breakPoints.remove(x);
							final BreakPointDigest dg = new BreakPointDigest(tid1, tid2);
							dg.start1 = bp.pos1;
							dg.end1 = bp.pos1;
							dg.start2 = bp.pos2;
							dg.end2 = bp.pos2;
							dg.count = 1;
							binMap.put(dg);
							}
						else
							{
							++x;
							}
						}
					for(final BreakPointDigest dg:binMap.toList()) {
						
						if(bedIntervals!=null)
							{
							if(!(
								bedIntervals.containsOverlapping(new Interval(ctg1,dg.start1,dg.end1))) || 
								bedIntervals.containsOverlapping(new Interval(ctg2,dg.start2,dg.end2)))
								{
								continue;
								}
							}
						
						pw.print(ctg1);
						pw.print("\t");
						pw.print(dg.start1);
						pw.print("\t");
						pw.print(dg.end1);
						pw.print("\t");
						pw.print(ctg2);
						pw.print("\t");
						pw.print(dg.start2);
						pw.print("\t");
						pw.print(dg.end2);
						pw.print("\t");
						pw.print(sampleName);
						pw.print("\t");
						pw.print(dg.count);
						pw.println();
						}
					}
				}
			
			pw.flush();
			pw.close();pw=null;
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(samReader);
			CloserUtil.close(pw);
			}
		}
	
	public static void main(final String[] args) {
		new SamTranslocations().instanceMainWithExit(args);

	}

}
