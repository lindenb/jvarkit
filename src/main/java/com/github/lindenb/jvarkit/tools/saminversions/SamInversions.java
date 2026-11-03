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
package com.github.lindenb.jvarkit.tools.saminversions;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.function.Predicate;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.BAMUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;

/**
BEGIN_DOC


END_DOC

**/

@Program(name="saminversions",
	description="Scan inversions in SAM using split and FF/RR reads on same strand.",
	keywords={"sam","bam","sv","inversion"},
	modificationDate="20261103",
	creationDate="20200609",
	generate_doc=false
	)
public class SamInversions extends Launcher
	{
	private static final Logger LOG = Logger.of(SamInversions.class);
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=CRAM_INDEXED_REFENCE)
	private List<Path> faidx_list = new ArrayList<>();
	@Parameter(names={"--distance"},description="search for breakend (BND) around 'x' bases.",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int max_distance = 100 ;
	@Parameter(names={"-B","--bed"},description="Optional BED file. SV should overlap this bed.")
	private Path bedFile = null;
	@Parameter(names={"--exclude-bed","--exclude"},description="skip reads pverlapping that bed file")
	private Path excludeBedFile = null;
	@Parameter(names={"--mapq"},description="min mapping quality.")
	private int min_mapq = 30;
	@Parameter(names={"--control","--controls"},description="BAM/CRAM file that can be used as a control. I discordand reads are found in any control, the candidate inversion is rejected. A file with the '.list' suffix is interpreted as a file containing the path to the bam files.")
	private Path controlFiles = null;
	
	private class LookForwardIterator extends AbstractCloseableIterator<SAMRecord> {
		private final List<SAMRecord> buffer = new ArrayList<>(100_000);
		final CloseableIterator<SAMRecord> delegate;
		LookForwardIterator(final CloseableIterator<SAMRecord> delegate) {
			this.delegate = delegate;
			}
		@Override
		public SAMRecord peek() {
			return at(0);
			}
		
		
		SAMRecord at(int i) {
			while(i>=this.buffer.size() && delegate.hasNext()) {
				final SAMRecord rec = delegate.next();
				if(!SAMRecordDefaultFilter.accept(rec)) continue;
				this.buffer.add(rec);
				}
			return i<this.buffer.size()?buffer.get(i):null;
			}
		@Override
		protected SAMRecord advance() {
			if(!buffer.isEmpty()) {
				return buffer.remove(0);
				}
			while(delegate.hasNext()) {
				final SAMRecord rec = delegate.next();
				if(!SAMRecordDefaultFilter.accept(rec)) continue;
				return rec;
				}
			return null;
			}
		@Override
		public void close() {
			delegate.close();
			this.buffer.clear();
			}
		}
	
	
	private static class Event {
		int start;
		int end;
		
	}
	
	private Path findReferenceForBam(final Path bamPath) {
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(bamPath);
		for(final Path ref: this.faidx_list) {
			final SAMSequenceDictionary dict2 = SequenceDictionaryUtils.extractRequired(ref);
			if(SequenceUtil.areSequenceDictionariesEqual(dict, dict2)) {
				return ref;
				}
			}
		throw new IllegalArgumentException("Cannot find REFerence for "+bamPath+".");
		}
	
	/** get avg depth in the interval */
	private double meanDepth(final List<SAMRecord> records,int start,int end) {
		double total=0.0;
		for(SAMRecord rec:records) {
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.isSecondaryOrSupplementary()) continue;
			if(rec.getReadFailsVendorQualityCheckFlag()) continue;
			if(rec.getDuplicateReadFlag()) continue;
			if(rec.getMappingQuality()<= this.min_mapq) continue;
			if(!CoordMath.overlaps(start, end, rec.getStart(), rec.getEnd())) continue;
			for(AlignmentBlock ab:rec.getAlignmentBlocks()) {
				for(int x=0;x<ab.getLength();x++) {
					int pos1 = ab.getReferenceStart()+x;
					if(pos1<start) continue;
					if(pos1>end) break;
					total++;
					}
				}
			}
		return total/CoordMath.getLength(start, end);
		}
	
	/** get avg MAPQ in the interval */
	private OptionalDouble meanMAPQ(final List<SAMRecord> records,int start,int end) {
		double total=0.0;
		int n=0;
		for(SAMRecord rec:records) {
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.isSecondaryOrSupplementary()) continue;
			if(rec.getReadFailsVendorQualityCheckFlag()) continue;
			if(rec.getDuplicateReadFlag()) continue;
			if(!CoordMath.overlaps(start, end, rec.getStart(), rec.getEnd())) continue;
			total+=rec.getMappingQuality();
			n++;
			}
		
		return n==0L?OptionalDouble.empty():OptionalDouble.of(total/n);
		}
	
	
	private boolean discordantFoundInControl(final List<SamReader> ctrls,final Locatable loc,final Locatable mateLoc) {
		for(SamReader sr:ctrls) {
			final SAMFileHeader hdr=sr.getFileHeader();
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(hdr);
			final ContigNameConverter converter= ContigNameConverter.fromOneDictionary(dict);
			final String ctg = converter.apply(loc.getContig());
			if(StringUtils.isBlank(ctg)) continue;
			int tid = dict.getSequenceIndex(ctg);
			final QueryInterval[] intervals= new QueryInterval[] {
					new QueryInterval(tid,loc.getStart(),loc.getEnd()),
					new QueryInterval(tid,mateLoc.getStart(),mateLoc.getEnd())
				};
			try(CloseableIterator<SAMRecord> iter = sr.queryOverlapping(QueryInterval.optimizeIntervals(intervals))) {
				while(iter.hasNext()) {
					final SAMRecord rec=iter.next();
					if(!is_FF_or_RR_read(rec)) continue;
					if(!CoordMath.overlaps(
							loc.getStart(), loc.getEnd(),
							mateLoc.getStart(), mateLoc.getEnd())) continue;

					if(!CoordMath.overlaps(
						loc.getStart(), loc.getEnd(),
						mateLoc.getStart(), mateLoc.getEnd())) continue;
					}
				}
			}
		return false;
		}


	/** open iterator for main SamReader */
	private CloseableIterator<SAMRecord> makeIterator(final SamReader sr) {
		final SAMFileHeader header=sr.getFileHeader();
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);

		if(this.bedFile==null) {
			return sr.iterator();
			}
		else
			{
			try(BedLineReader blr = BedLineReader.openPath(this.bedFile)) {
				blr.setContigNameConverter(dict);
				return sr.queryOverlapping(blr.optimizeIntervals(dict));
				}
			}
		}
	

	/** fetch reads in a given region */
	private List<SAMRecord> getReads(final SamReader sr,final Locatable loc,Predicate<SAMRecord> accept) {
		try(CloseableIterator<SAMRecord> r= sr.queryOverlapping(loc.getContig(), loc.getStart(), loc.getEnd())) {
			final List<SAMRecord> L= new ArrayList<SAMRecord>(loc.getLengthOnReference()*40);
			while(r.hasNext()) {
				final SAMRecord rec=r.next();
				if(accept!=null && !accept.test(rec)) continue;
				L.add(rec);
				}
			return L;
			}
		}
	
	
	
	private Event challenge(final LookForwardIterator srcIter,final SamReader sr,SAMRecord rec1,Predicate<Locatable> accept,List<SamReader> controls) {
		int idx=0;
		final int min_clip_size = 10;

		/* for reads for breakend 1 */
		List<SAMRecord> reads = new ArrayList<SAMRecord>(this.max_distance*50);
		reads.add(rec1);
		for(;;) {
			final SAMRecord reci = srcIter.at(idx);
			if(reci==null || !rec1.contigsMatch(reci)) break;
			if(!reci.withinDistanceOf(rec1, this.max_distance)) break;
			reads.add(reci);
			idx++;
			}
		
		List<SAMRecord> suppl_reads = new ArrayList<SAMRecord>();
		List<SAMRecord> split_reads = new ArrayList<SAMRecord>();
		
		for(SAMRecord rec: reads) {
			if(!SAMRecordDefaultFilter.accept(rec)) continue;
			if(rec.getReadNegativeStrandFlag()) continue;
			if(!rec.getCigar().isRightClipped()) // expect a clip
			if(rec.getCigar().getLastCigarElement().getLength()< min_clip_size) continue;
			if(!CoordMath.overlaps(
				rec.getStart(), rec.getEnd(),
				rec1.getStart(), rec1.getMateAlignmentStart())) continue;
			boolean has_suppl_read=false;
			for(SAMRecord suppl: SAMUtils.getOtherCanonicalAlignments(rec)) {
				if(!suppl.contigsMatch(rec1)) continue;
				if(!suppl.getReadNegativeStrandFlag()) continue;
				if(suppl.getStart() < rec1.getMateAlignmentStart()) continue;//we expect a 2nd break-end AFTER the mate
				suppl_reads.add(suppl);
				has_suppl_read=true;
				}
			if(has_suppl_read) {
				split_reads.add(rec1);
				}
			}
		
		final int min_reads = 5;
		if(split_reads.size() < min_reads) {
			return null;
			}
		
		if(suppl_reads.size() < min_reads) {
			return null;
			}
		
		final Event evt=new Event();
		evt.start = (int) new Median().evaluate(
				split_reads.stream().mapToDouble(R->R.getEnd()).sorted().toArray()
				);
		
		
		
		final Locatable mateInterval = new SimpleInterval(
				rec1.getContig(),
				suppl_reads.stream().mapToInt(R->R.getStart()).min().getAsInt(),
				suppl_reads.stream().mapToInt(R->R.getEnd()).max().getAsInt()
				);
		if(mateInterval.getLengthOnReference()> this.max_distance) {
			return null;
			}
		
		
		reads = getReads(sr, mateInterval, null);
		
		suppl_reads.clear();
		split_reads.clear();
		
		for(SAMRecord rec: reads) {
			if(!SAMRecordDefaultFilter.accept(rec)) continue;
			if(!rec.getReadNegativeStrandFlag()) continue;
			if(!rec.getCigar().isLeftClipped()) // expect a clip
			if(rec.getCigar().getFirstCigarElement().getLength()< min_clip_size) continue;
			if(!CoordMath.overlaps(
				rec.getStart(), rec.getEnd(),
				rec1.getStart(), BAMUtils.getMateEnd(rec1))) continue;
			boolean has_suppl_read=false;
			for(SAMRecord suppl: SAMUtils.getOtherCanonicalAlignments(rec)) {
				if(!suppl.contigsMatch(rec1)) continue;
				if(suppl.getReadNegativeStrandFlag()) continue;
				if(suppl.getStart() > rec1.getEnd()) continue;//we expect a 2nd break-end AFTER the mate
				suppl_reads.add(suppl);
				has_suppl_read=true;
				}
			if(has_suppl_read) {
				split_reads.add(rec1);
				}
			}
		
		if(split_reads.size() < min_reads) {
			return null;
			}
		
		if(suppl_reads.size() < min_reads) {
			return null;
			}
		evt.end = (int) new Median().evaluate(
				split_reads.stream().mapToDouble(R->R.getStart()).sorted().toArray()
				);
		
		return evt;
		}
	
	private static boolean is_FF_or_RR_read(SAMRecord rec) {
		if(!SAMRecordDefaultFilter.accept(rec)) return false;
		if(!rec.getReadPairedFlag()) return false;
		if(rec.getMateUnmappedFlag()) return false;
		if(!rec.getReferenceIndex().equals(rec.getMateReferenceIndex()))return false;
		// we expect different strand for read and mate
		if(rec.getReadNegativeStrandFlag() != rec.getMateNegativeStrandFlag()) return false;
		return true;
		}
	
	@Override
	public int doWork(final List<String> args) {
		final List<SamReader> control_readers = new ArrayList<>()  ;

		if(this.max_distance<100) {
			LOG.error("max size inversion must be >=100");
			return -1;
			}
		try {
			final String input = oneAndOnlyOneFile(args);
			final Path inputBam = Paths.get(input);
			final Path referenceFaidx = findReferenceForBam(inputBam);
			final SamReaderFactory srf =  super.createSamReaderFactory().
					validationStringency(ValidationStringency.LENIENT).
					referenceSequence(referenceFaidx)
					;
			
			if(this.controlFiles!=null) {
				for(Path path: IOUtils.unrollPath(this.controlFiles)) {
					control_readers.add(srf
							.referenceSequence(findReferenceForBam(path))
							.open(path)
							);
					}
				}
			
			try(ReferenceSequenceFile referenceSequenceFile = new IndexedFastaSequenceFile(referenceFaidx)) {
				try(SamReader samReader1= srf.referenceSequence(referenceFaidx).open(inputBam)) {
					if(!samReader1.hasIndex()) {
						LOG.error("file is not indexed "+input);
						return -1;
						}
					final SAMFileHeader samHeader = samReader1.getFileHeader();
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samHeader);
					final String sampleName= BAMUtils.getSampleName(samHeader, inputBam);
					final Hyperlink hyperlink = Hyperlink.compile(dict);
					final String build =  SequenceDictionaryUtils.getBuildName(dict).orElse(null);

					
					final Predicate<Locatable> accept_reads;
					if(this.excludeBedFile!=null) {
						try(BedLineReader blr = BedLineReader.openPath(this.excludeBedFile)) {
							blr.setContigNameConverter(dict);
							final IntervalTreeMap<Boolean> xmap= blr.toBooleanIntervalTreeMap();
							accept_reads = REC->!xmap.containsOverlapping(REC);
							}
						}
					else
						{
						accept_reads = REC->true;
						}
					
					//open a second reader to fetch the mate
					try(SamReader samReader2= srf.referenceSequence(referenceFaidx).open(inputBam)) {
						try(JsonWriter js = new JsonWriter(super.openPathOrStdoutAsPrintWriter(this.outputFile))) {
							js.setIndent(" ");
							js.beginArray();
							
							/** loop over eall reads */
							try(LookForwardIterator iter1= new LookForwardIterator( makeIterator(samReader1))) {
								while(iter1.hasNext()) {
									final SAMRecord rec1 = iter1.next();
									if(!SAMRecordDefaultFilter.accept(rec1,this.min_mapq)) {
										continue;
										}
									if(rec1.getReadNegativeStrandFlag()) continue;//we expect forward when scanning from 5' to 3'
									if(!is_FF_or_RR_read(rec1)) continue;
									
									final Event event = challenge(iter1,samReader2,rec1,accept_reads,control_readers);
									
									if(event!=null ) {
										js.beginObject();
										js.name("bam");js.value(String.valueOf(inputBam));
										js.name("fasta");js.value(String.valueOf(referenceFaidx));
										js.name("build");js.value(build);
										js.name("sample");js.value(sampleName);
										js.name("chromosome");js.value(rec1.getContig());
										js.name("start");js.value(event.start);
										js.name("end");js.value(event.end);
										final Optional<String> href= hyperlink.apply(rec1.getContig(),event.start,event.end);
										if(href.isPresent()) {
											js.name("href");js.value(href.get());
											}
										js.beginObject();
										
										/* consumme remaining data */
										while(iter1.hasNext()) {
											final SAMRecord rec2 = iter1.peek();
											if(rec2==null) break;
											if(!rec2.contigsMatch(rec1)) break;
											if(rec2.getStart()> event.start) break;
											iter1.next();//consumme
											}
										}
									
									
									
									
									}
								}
							
							
							js.endArray();
							js.flush();
							}//end vcw
						} // end cram2
					} //end cram1
				} // end fasta
			return 0;
			} 
		catch (final Throwable e) {
			LOG.error(e);
			return -1;
			}
		finally {
			for(SamReader sr:control_readers) {
				try {sr.close(); } catch(IOException err2) {}
				}
			}
		}

	
	public static void main(final String[] args)
		{
		new SamInversions().instanceMainWithExit(args);
		}
	}
