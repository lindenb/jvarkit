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
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
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
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.BAMUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

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
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC


END_DOC

**/

@Program(name="saminversions",
	description="Scan inversions in SAM using paired read on same strand.",
	keywords={"sam","bam","sv","inversion"},
	modificationDate="20201211",
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
	@Parameter(names={"-s","-supporting"},description="Don't print the variant if max FORMAT/DP < 's'")
	private int min_supporting_reads = 1;
	@Parameter(names={"-B","--bed"},description="Optional BED file. SV should overlap this bed.")
	private Path bedFile = null;
	@Parameter(names={"--exclude-bed","--exclude"},description="skip reads pverlapping that bed file")
	private Path excludeBedFile = null;
	@Parameter(names={"--mapq"},description="min mapping quality.")
	private int min_mapq = 30;
	@Parameter(names={"--control","--controls"},description="BAM/CRAM file that can be used as a control. I discordand reads are found in any control, the candidate inversion is rejected. A file with the '.list' suffix is interpreted as a file containing the path to the bam files.")
	private Path controlFiles = null;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();

	
	private static class BreakEnd {
		int count_same_strand_paired_reads=0;
		int count_split_reads=0;
		int start;
		int end;
		double mean_depth=0;
		double mean_qual = 0;
		boolean isValid() {
			return false;
			}
		}
	
	private static class Event {
		int end_position = 0;
		BreakEnd bnd1=null;
		BreakEnd bnd2=null;
		boolean found_in_controls;
		boolean isValid() {
			
			return !found_in_controls &&  bnd1!=null && bnd2!=null && bnd1.isValid() && bnd2.isValid();
		}
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
	
	
	
	private Event challenge(final SamReader sr,SAMRecord rec1,Predicate<Locatable> accept,List<SamReader> controls) {
		final Event evt=new Event();
		evt.end_position= rec1.getStart();
		
		final Locatable interval1 = new SimpleInterval(rec1.getContig(),rec1.getStart(),rec1.getStart()+this.max_distance);
		final Locatable mateInterval = new SimpleInterval(
				rec1.getContig(),
				Math.max(rec1.getStart(), rec1.getMateAlignmentStart()-this.max_distance),
				BAMUtils.getMateEnd(rec1)+this.max_distance
				);
		
		List<SAMRecord> reads = getReads(sr, interval1,null);
		
		List<SAMRecord> rr_ff_reads= reads.stream()
			.filter(R->is_FF_or_RR_read(R))
			.filter(R->!R.getReadNegativeStrandFlag())
			.filter(R->CoordMath.overlaps(R.getMateAlignmentStart(),BAMUtils.getMateEnd(R), mateInterval.getStart(), mateInterval.getEnd()))
			.collect(Collectors.toList())
			;
		
		evt.bnd1=new BreakEnd();
		evt.bnd1.count_same_strand_paired_reads =rr_ff_reads.size();
		evt.bnd1.start = (int) new Median().evaluate(
				rr_ff_reads.stream().mapToDouble(R->R.getStart()).sorted().toArray()
				);
		evt.bnd1.end = (int) new Median().evaluate(
				rr_ff_reads.stream().mapToDouble(R->R.getEnd()).sorted().toArray()
				);
		
		if(evt.bnd1.start > evt.bnd1.end) {
			int x= evt.bnd1.start;
			evt.bnd1.start = evt.bnd1.end;
			evt.bnd1.end = x;
			}
		
		List<SAMRecord> split_reads = reads.stream()
				.filter(R->SAMRecordDefaultFilter.accept(R))
				.filter(R->!R.getReadNegativeStrandFlag())
				.filter(R->R.getCigar().isRightClipped())
				.flatMap(R->SAMUtils.getOtherCanonicalAlignments(R).stream())
				.filter(suppl->suppl.contigsMatch(rec1))
				.filter(suppl->rec1.getReadNegativeStrandFlag()!= suppl.getReadNegativeStrandFlag())
				.filter(suppl->CoordMath.overlaps(suppl.getStart(), suppl.getEnd(), mateInterval.getStart(), mateInterval.getEnd()))
				.collect(Collectors.toList())
				;
		
		evt.bnd1.count_split_reads = split_reads.size();
		evt.bnd1.mean_depth =  meanDepth(reads, evt.bnd1.start,evt.bnd1.end);
		evt.bnd1.mean_qual =  meanMAPQ(reads, evt.bnd1.start,evt.bnd1.end).orElse(0.0);
		
		evt.bnd2=new BreakEnd();
		
		
		reads = getReads(sr, mateInterval,null);
		
		rr_ff_reads= reads.stream()
			.filter(R->is_FF_or_RR_read(R))
			.filter(R->R.getReadNegativeStrandFlag())
			.filter(R->CoordMath.overlaps(R.getMateAlignmentStart(),BAMUtils.getMateEnd(R), interval1.getStart(), interval1.getEnd()))
			.collect(Collectors.toList())
			;
		
		evt.bnd2.count_same_strand_paired_reads =rr_ff_reads.size();

		evt.bnd2.start = (int) new Median().evaluate(
				rr_ff_reads.stream().mapToDouble(R->R.getStart()).sorted().toArray()
				);
		evt.bnd2.end = (int) new Median().evaluate(
				rr_ff_reads.stream().mapToDouble(R->R.getEnd()).sorted().toArray()
				);
		
		if(evt.bnd2.start > evt.bnd2.end) {
			int x= evt.bnd2.start;
			evt.bnd2.start = evt.bnd2.end;
			evt.bnd2.end = x;
			}
		
		evt.bnd2.mean_depth =  meanDepth(reads, evt.bnd2.start,evt.bnd2.end);
		evt.bnd2.mean_qual =  meanMAPQ(reads, evt.bnd2.start,evt.bnd2.end).orElse(0.0);

		
		split_reads = reads.stream()
				.filter(R->SAMRecordDefaultFilter.accept(R))
				.filter(R->R.getReadNegativeStrandFlag())
				.filter(R->R.getCigar().isLeftClipped())
				.flatMap(R->SAMUtils.getOtherCanonicalAlignments(R).stream())
				.filter(suppl->suppl.contigsMatch(rec1))
				.filter(suppl->rec1.getReadNegativeStrandFlag()!= suppl.getReadNegativeStrandFlag())
				.filter(suppl->CoordMath.overlaps(suppl.getStart(), suppl.getEnd(), interval1.getStart(), interval1.getEnd()))
				.collect(Collectors.toList())
				;
		
		evt.bnd2.count_split_reads = split_reads.size();
		
		evt.found_in_controls = discordantFoundInControl(controls,interval1,mateInterval);
		
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
			final Allele REF = Allele.create("N", true);
			final Allele INV = Allele.create("<INV>", false);
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
						try(VariantContextWriter vcw = this.writingVariants.open(this.outputFile)) {

							final Set<VCFHeaderLine> meta=new HashSet<>();
							VCFStandardHeaderLines.addStandardFormatLines(meta,true,
									VCFConstants.GENOTYPE_KEY,
									VCFConstants.DEPTH_KEY
									
									);
							VCFStandardHeaderLines.addStandardInfoLines(meta,true,
									VCFConstants.DEPTH_KEY,
									VCFConstants.END_KEY
									);
							meta.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"SV length"));
							meta.add(new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String,"Structural variant type"));
							
							final VCFHeader header=new VCFHeader(meta,Collections.singleton(sampleName));
							JVarkitVersion.getInstance().addMetaData(this, header);

							header.setSequenceDictionary(dict);
							vcw.writeHeader(header);
	
							
							/** loop over eall reads */
							try(PeekableIterator<SAMRecord> iter1= new PeekableIterator<>( makeIterator(samReader1))) {
								while(iter1.hasNext()) {
									final SAMRecord rec1 = iter1.next();
									if(!SAMRecordDefaultFilter.accept(rec1,this.min_mapq)) {
										continue;
										}
									if(rec1.getReadNegativeStrandFlag()) continue;//we expect forward when scanning from 5' to 3'
									if(!is_FF_or_RR_read(rec1)) continue;
									
									final Event event = challenge(samReader2,rec1,accept_reads,control_readers);
									
									if(event.isValid()) {
										//TODO
									}
									
									
									while(iter1.hasNext()) {
										final SAMRecord rec2 = iter1.peek();
										if(!rec2.contigsMatch(rec1)) break;
										if(rec2.getStart() >= event.end_position) break;
										iter1.next();//consumme	
										}
									
									
									}
								}
							
							
							
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
