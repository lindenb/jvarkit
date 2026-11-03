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
package com.github.lindenb.jvarkit.tools.samtransloc;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.log.ProgressDict;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.BAMUtils;
import com.github.lindenb.jvarkit.samtools.util.AbstractLocatable;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SequenceUtil;
/**

BEGIN_DOC

## input

input is a set of bam files or one file with the suffix '.list' containing the path to the bams.

## Example

```
$ java -jar dist/jvarkit.jar samtransloc src/test/resources/HG02260.transloc.chr9.14.bam  > out.json
```


END_DOC
*/
@Program(name="samtransloc",
	description="Explore balanced translocations between two chromosomes using discordant paired-end reads.",
	keywords={"sam","bam","sv","translocation"},
	jvarkit_hidden = false,
	creationDate="20171130",
	modificationDate="20261101",
	menu="BAM Manipulation"
	)
public class SamTranslocations extends Launcher {
	private static final Logger LOG = Logger.of(SamTranslocations.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=CRAM_INDEXED_REFENCE)
	private List<Path> faidx_list = new ArrayList<>();
	@Parameter(names={"-d","--distance"},
			description="Max distance between two read to test if they both end at the same ~ position. " +DistanceParser.OPT_DESCRIPTION,
			converter=DistanceParser.StringConverter.class,
			splitter=NoSplitter.class
			)
	private int fuzzy_distance = 150;
	@Parameter(names={"-B","--bed"},description="Optional BED file. SV should overlap this bed.")
	private Path bedFile = null;
	@Parameter(names={"--exclude-bed","--exclude"},description="skip reads pverlapping that bed file")
	private Path excludeBedFile = null;
	@Parameter(names={"--chrom-regex"},description="Only consider the chromosomes matching the following regular expression.")
	private String chromRegex="(chr)?[0-9XY]+";
	@Parameter(names={"--mapq"},description="min mapping quality.")
	private int min_mapq = 30;
	@Parameter(names={"--control","--controls"},description="BAM/CRAM file that can be used as a control. I discordand reads are found in any control, the candidate transloc is rejected. A file with the '.list' suffix is interpreted as a file containing the path to the bam files.")
	private Path controlFiles = null;

	
	
	/** open iterator for main SamReader */
	private CloseableIterator<SAMRecord> makeIterator(final SamReader sr) {
		final SAMFileHeader header=sr.getFileHeader();
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);

		if(this.bedFile==null && StringUtils.isBlank(this.chromRegex)) {
			return sr.iterator();
			}
		else if(this.bedFile==null) {
			return sr.queryOverlapping(
				QueryInterval.optimizeIntervals(
					dict.getSequences().stream()
						.filter(SSR->SSR.getContig().matches(this.chromRegex))
						.map(SSR->new QueryInterval(dict.getSequenceIndex(SSR.getContig()), 1, SSR.getEnd()))
						.toArray(X->new QueryInterval[X])
					)
				);
			}
		else
			{
			try(BedLineReader blr = BedLineReader.openPath(this.bedFile)) {
				blr.setContigNameConverter(dict);
				return sr.queryOverlapping(
						QueryInterval.optimizeIntervals(
							blr.stream()
								.filter(BED->BED.getContig().matches(this.chromRegex))
								.map(BED->new QueryInterval(dict.getSequenceIndex(BED.getContig()), BED.getStart(), BED.getEnd()))
								.toArray(X->new QueryInterval[X])
							)
						);
				}
			}
		}
	
	/* return true if read is paired and reads are mapped on discordant contigs */
	private static boolean isDiscordant(final SAMRecord rec) {
		return !rec.getReadUnmappedFlag()	&& 
			rec.getReadPairedFlag() && 
			!rec.getMateUnmappedFlag() && 
			!rec.getReferenceIndex().equals(rec.getMateReferenceIndex())
			;
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
	
	/** hold the count for a given strand in a breakend */
	private static class StrandCount {
		int count_discordant;// discordant with the expected chromosome
		int count_other_discordant; //not the same discordant chromosome
		int count_splitread; // number of split reads, whatever
		int count_splitread_discordant; // number of split reads mapping on the discordant chromosome
		int count_splitread_other_discordant; // number of split reads mapping on other discordant chromosome
		boolean isValid() {
			return count_discordant>4 && count_splitread>0 && count_discordant > count_other_discordant*2;
			}
		
		void json(final JsonWriter j) throws IOException{
			j.beginObject();
			j.name("valid");j.value(isValid());
			j.name("high_qual");j.value(isHighQual());
			j.name("discordant");j.value(count_discordant);
			j.name("other_discordant");j.value(count_other_discordant);
			j.name("splitread");j.value(count_splitread);
			j.name("splitread_discordant");j.value(count_splitread_discordant);
			j.name("splitread_other_discordant");j.value(count_splitread_other_discordant);
			j.endObject();
			}
		boolean isHighQual() {
			return count_discordant>=10 && count_other_discordant==0 && count_splitread_discordant>0 && count_splitread_other_discordant==0;
			}
		}
	
	/** class describind a Breakend (left or right) on a chromosome */
	private static class BreakEnd extends AbstractLocatable {
		String contig="";
		int start=0;
		int end= 0;
		double depth = 0.0;
		double mapq = 0.0;
		int count_insertions=0;
		
		StrandCount plus;
		StrandCount minus;
		
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
		
		void json(final JsonWriter j,final Hyperlink hyperlink) throws IOException {
			j.beginObject();
			j.name("chrom");j.value(getContig());
			j.name("start");j.value(getStart());
			j.name("end");j.value(getEnd());
			j.name("interval");j.value(getContig()+":"+getStart()+"-"+getEnd());
			j.name("length");j.value(getLengthOnReference());
			j.name("mean_depth");j.value(depth);
			j.name("mean_mapq");j.value(mapq);
			j.name("insertions");j.value(count_insertions);
			final String url = hyperlink.apply(this).orElse(null);
			if(url!=null) {
				j.name("href");j.value(url);
				}
			j.name("valid");j.value(isValid());
			j.name("high_qual");j.value(isHighQual());
			j.name("left");
			plus.json(j);
			j.name("right");
			minus.json(j);
			

			j.endObject();
			}
		
		boolean isValid() {
			final double factor=10.0;
			return this.plus.isValid() &&
					this.minus!=null &&
					this.minus.isValid() && 
					this.depth<=this.plus.count_discordant*factor && 
					this.depth<=this.minus.count_discordant*factor;
			}
		boolean isHighQual() {
			return isValid() &&
					plus.isHighQual() &&
					minus!=null &&
					minus.isHighQual()
					;
			}
		}
	
	/** a translocation event, with two breakend , one for each chromosome */
	private static class Event {
		/* where to advance the Iterator for next scan */
		int end_position = 0;
		BreakEnd bnd1; //chromosome1
		BreakEnd bnd2; //chromosome2
		boolean isValid() {
			return bnd1!=null && bnd1.isValid() && bnd2!=null && bnd2.isValid();
			}
		boolean isHighQual() {
			return isValid() && bnd1!=null && bnd1.isHighQual() && bnd2!=null && bnd2.isHighQual();
			}
		}
	
	/** count features on a given strand */
	private StrandCount makeStrandCount(
			final List<SAMRecord> reads_contig,
			final Locatable loc,
			final Locatable mate,
			boolean is_negative_strand
			) {
		final StrandCount evt = new StrandCount();
		final Predicate<SAMRecord> accept_strand = R->R.getReadNegativeStrandFlag()==is_negative_strand;
		final Predicate<SAMRecord> accept_clip = is_negative_strand?R->R.getCigar().isLeftClipped():R->R.getCigar().isRightClipped();

		evt.count_discordant = (int)reads_contig
				.stream()
				.filter(R->isDiscordant(R) && !R.isSecondaryOrSupplementary() && !R.getDuplicateReadFlag())
				.filter(R->R.getMateReferenceName().equals(mate.getContig()))
				.filter(accept_strand)
				.count();
		
		
		evt.count_other_discordant = (int)reads_contig
				.stream()
				.filter(R->isDiscordant(R) && !R.isSecondaryOrSupplementary() && !R.getDuplicateReadFlag())
				.filter(R->!R.getMateReferenceName().equals(mate.getContig()))
				.filter(accept_strand)
				.count();
		
		evt.count_splitread = (int)reads_contig
				.stream()
				.filter(R->!R.getReadUnmappedFlag() && !R.isSecondaryOrSupplementary() && !R.getDuplicateReadFlag())
				.filter(accept_clip)
				.filter(accept_strand)
				.count();
		
		evt.count_splitread_discordant = (int)reads_contig
				.stream()
				.filter(R->!R.getReadUnmappedFlag() && !R.isSecondaryOrSupplementary() && !R.getDuplicateReadFlag())
				.filter(accept_clip)
				.filter(accept_strand)
				.flatMap(R->SAMUtils.getOtherCanonicalAlignments(R).stream())
				.filter(suppl->suppl.withinDistanceOf(mate, fuzzy_distance))
				.count();

		evt.count_splitread_other_discordant = (int)reads_contig
				.stream()
				.filter(R->!R.getReadUnmappedFlag() && !R.isSecondaryOrSupplementary() && !R.getDuplicateReadFlag())
				.filter(accept_clip)
				.filter(accept_strand)
				.flatMap(R->SAMUtils.getOtherCanonicalAlignments(R).stream())
				.filter(suppl->!suppl.withinDistanceOf(mate, fuzzy_distance))
				.count();


		return evt;
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
	
	/** get avg depth in the interval */
	private int count_insertions(final List<SAMRecord> records,int start,int end) {
		final Set<Integer> positions=new HashSet<>();
		for(SAMRecord rec:records) {
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.isSecondaryOrSupplementary()) continue;
			if(rec.getReadFailsVendorQualityCheckFlag()) continue;
			if(rec.getDuplicateReadFlag()) continue;
			if(rec.getMappingQuality()<= this.min_mapq) continue;
			if(!CoordMath.overlaps(start, end, rec.getStart(), rec.getEnd())) continue;

			int ref=rec.getAlignmentStart();
			for(CigarElement ce:rec.getCigar()) {
				if(ref>end) break;
				if(ce.getOperator().equals(CigarOperator.I) && ref>=start) {
					positions.add(ref);
					}
				
				if(ce.getOperator().consumesReferenceBases()) {
					ref+=ce.getLength();
					}
				}
			}
		return positions.size();
		}
	
	private BreakEnd makeBreakEnd(
			final List<SAMRecord> reads_contig,
			final Locatable loc,
			final Locatable mate
			) {
		final BreakEnd evt = new BreakEnd();
		evt.contig = loc.getContig();
		evt.start= loc.getStart();
		evt.end= loc.getEnd();
		
		
		evt.start =   reads_contig
				.stream()
				.filter(R->isDiscordant(R) && !R.isSecondaryOrSupplementary() && !R.getDuplicateReadFlag())
				.filter(R->R.getMateReferenceName().equals(mate.getContig()))
				.filter(R->!R.getReadNegativeStrandFlag())
				.mapToInt(R->R.getStart())
				.max()
				.orElse(loc.getStart());
		
		evt.end =   reads_contig
				.stream()
				.filter(R->isDiscordant(R) && !R.isSecondaryOrSupplementary() && !R.getDuplicateReadFlag())
				.filter(R->R.getMateReferenceName().equals(mate.getContig()))
				.filter(R->R.getReadNegativeStrandFlag())
				.mapToInt(R->R.getStart())
				.min()
				.orElse(loc.getEnd());
		
		
		evt.plus = makeStrandCount(reads_contig,loc,mate,false);
		
		if(evt.plus.isValid()) {
			evt.minus = makeStrandCount(reads_contig,loc,mate,true);
			
			
			if(evt.start>evt.end) {
				int x= evt.start;
				evt.start=evt.end;
				evt.end=x;
				}
			
			evt.depth =meanDepth(reads_contig,evt.start,evt.end);
			evt.count_insertions = count_insertions(reads_contig,evt.start,evt.end);
			evt.mapq = meanMAPQ(reads_contig,evt.start,evt.end).orElse(0.0);
			}
		return evt;
		}
	
	private Event challenge(final SamReader sr,final SAMRecord rec,final List<SamReader> controls) {
		final Event evt=new Event();
		evt.end_position=rec.getAlignmentEnd();
		
		
		final Locatable loc = new SimpleInterval(
				rec.getContig(),
				rec.getStart(),
				rec.getEnd()+this.fuzzy_distance
				);
		final Locatable mate =  new SimpleInterval(
				rec.getMateReferenceName(),
				Math.max(1,rec.getMateAlignmentStart() - this.fuzzy_distance),
				rec.getMateAlignmentStart()+ rec.getReadLength()+ this.fuzzy_distance
				);
		
		if(discordantFoundInControl(controls, rec, mate)) {
			return evt;
			}
		
		List<SAMRecord> reads_contig = getReads(sr,loc,null);
		evt.bnd1 = makeBreakEnd(reads_contig , loc, mate);
		
		reads_contig = getReads(sr,mate,null);
		evt.bnd2 = makeBreakEnd(reads_contig , mate, loc);
		
		evt.end_position = Math.max( evt.bnd1.end,rec.getAlignmentEnd())+1;
		return evt;
		}
	
	private boolean discordantFoundInControl(final List<SamReader> ctrls,final Locatable loc,final Locatable mateLoc) {
		for(SamReader sr:ctrls) {
			final SAMFileHeader hdr=sr.getFileHeader();
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(hdr);
			final ContigNameConverter converter= ContigNameConverter.fromOneDictionary(dict);
			final String ctg = converter.apply(loc.getContig());
			if(StringUtils.isBlank(ctg)) continue;
			final String mateCtg2  = converter.apply(mateLoc.getContig());
			if(StringUtils.isBlank(mateCtg2)) continue;
			try(CloseableIterator<SAMRecord> iter = sr.queryOverlapping(ctg,loc.getStart(),loc.getEnd())) {
				while(iter.hasNext()) {
					final SAMRecord rec=iter.next();
					if(rec.getReadUnmappedFlag())  continue;
					if(rec.getDuplicateReadFlag()) continue;
					if(!rec.getReadPairedFlag()) continue;
					if(!rec.contigsMatch(loc)) continue;
					if(rec.getMateUnmappedFlag()) continue;
					if(rec.isSecondaryOrSupplementary())continue;
					if(!rec.getMateReferenceName().equals(mateCtg2)) continue;
					if(CoordMath.overlaps(loc.getStart(), loc.getEnd(),
						mateLoc.getStart(), mateLoc.getEnd())) return true;
					}
				}
			}
		return false;
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
	
	@Override
	public int doWork(final List<String> args) {
		final List<SamReader> control_readers = new ArrayList<>()  ;

		if(this.fuzzy_distance<=0) {
			LOG.error("fuzzy_distance <=0 ("+fuzzy_distance+")");
			return -1;
			}
		if(this.faidx_list.isEmpty()) {
			LOG.error("No reference list was set by user.");
			return -1;
			}
		
		try {
			final Path bamInput= Paths.get(oneAndOnlyOneFile(args));
			final Path faidx_reference = findReferenceForBam(bamInput);
			final SamReaderFactory srf=super.createSamReaderFactory()
					.referenceSequence(faidx_reference);
			
			
			if(this.controlFiles!=null) {
				for(Path path: IOUtils.unrollPath(this.controlFiles)) {
					control_readers.add(srf
							.referenceSequence(findReferenceForBam(path))
							.open(path)
							);
					}
				}
			
			try(SamReader sr1 = srf.referenceSequence(faidx_reference).open(bamInput))  {
				if(!sr1.isQueryable()) {
					LOG.error("cannot query "+bamInput);
					return -1;
					}
				final SAMFileHeader header = sr1.getFileHeader();
				final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(header);
				final ProgressDict<SAMRecord> progress = new ProgressDict<SAMRecord>(LOG,dict);
				progress.setEnabled(this.bedFile==null);
				final Hyperlink hyperlink = Hyperlink.compile(dict);
				final String build =  SequenceDictionaryUtils.getBuildName(dict).orElse(null);
				final String sampleName= BAMUtils.getSampleName(header,bamInput);
				
				final Set<Integer> enabled_tids =  StringUtils.isBlank(this.chromRegex)?
						dict.getSequences().stream().map(SR->SR.getSequenceIndex()).collect(Collectors.toSet()):
						dict.getSequences().stream().filter(SR->SR.getContig().matches(chromRegex)).map(SR->SR.getSequenceIndex()).collect(Collectors.toSet());
				if(enabled_tids.isEmpty()) {
					LOG.error("no contig name is enabled");
					return -1;
					}
				final boolean[] accepted_contig_array=new boolean[dict.size()];
				Arrays.fill(accepted_contig_array, false);
				for(int accept_tid: enabled_tids) {
					accepted_contig_array[accept_tid]=true;
					}
				
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
				
				try(SamReader sr2 = srf.referenceSequence(faidx_reference).open(bamInput))  {
					try(PeekableIterator<SAMRecord> iter1= new PeekableIterator<SAMRecord>(makeIterator(sr1))) {
						try(PrintWriter w =super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
							try(JsonWriter j = new JsonWriter(w)) {
								j.setIndent(" ");
								j.beginArray();
								while(iter1.hasNext()) {
									final SAMRecord rec = progress.apply( iter1.next());
									if(rec.getReadUnmappedFlag())  continue;
									if(rec.getDuplicateReadFlag()) continue;
									if(rec.getReadFailsVendorQualityCheckFlag()) continue;
									if(!rec.getReadPairedFlag()) continue;
									if(rec.getMateUnmappedFlag()) continue;
									if(rec.isSecondaryOrSupplementary())continue;
									if(rec.getReadNegativeStrandFlag()) continue; //scanning from 5' to 3' we expect first to find on plus strand
									if(rec.getMappingQuality() <= this.min_mapq) continue;
									// no, may be using a workflow .. if(rec.getMateReferenceIndex().compareTo(rec.getReferenceIndex())<0) continue;; // already done as reads are mapped on reference index
									if(!isDiscordant(rec)) continue;

									if(!accepted_contig_array[rec.getReferenceIndex()]) continue;
									if(!accepted_contig_array[rec.getMateReferenceIndex()]) continue;
									
									/* read in black region ?*/
									if(this.excludeBedFile!=null) {
										if(!accept_reads.test(rec)) continue;
										if(!accept_reads.test(new SimpleInterval(
												rec.getMateReferenceName(),
												rec.getMateAlignmentStart(),
												BAMUtils.getMateEnd(rec)))) continue;
										}
									
									final Event event = challenge(sr2,rec,control_readers);
									
									if(event.isValid()) {
										j.beginObject();
										j.name("sample"); j.value(sampleName);
										j.name("bam"); j.value(String.valueOf(bamInput));
										j.name("high_qual"); j.value(event.isHighQual());
										if(faidx_reference!=null) {
											j.name("reference");j.value(String.valueOf(faidx_reference));
											}
										if(build!=null) {
											j.name("build");j.value(build);
											}
										j.name("bnd1"); event.bnd1.json(j,hyperlink);
										j.name("bnd2"); event.bnd2.json(j,hyperlink);
										
										j.endObject();
										j.flush();
										
										
										while(iter1.hasNext()) {
											final SAMRecord rec2 = iter1.peek();
											if(rec2.getReadUnmappedFlag()) {
												iter1.next();
												continue;
												}
											
											if(!rec2.contigsMatch(rec)) break;
											
											if(rec2.getReadNegativeStrandFlag()) {
												iter1.next();
												continue;
												}
											
											if(rec2.getStart() > event.end_position) {
												break;
												}

											iter1.next();//consumme	
											}
										
										}
									}
								j.endArray();
								j.flush();
								}
							w.flush();
							}
						}
					}
				}
			
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			for(SamReader sr:control_readers) {
				try {sr.close(); } catch(IOException err2) {}
				}
			}
		}
	
	public static void main(final String[] args) {
		new SamTranslocations().instanceMainWithExit(args);

	}

}
