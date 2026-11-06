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

import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.function.Predicate;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.BAMUtils;
import com.github.lindenb.jvarkit.samtools.Read;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.SAMRecordRead;
import com.github.lindenb.jvarkit.samtools.util.AbstractLocatable;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.google.gson.stream.JsonWriter;

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
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

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
	private int min_mapq = 10;
	@Parameter(names={"--control","--controls"},description="BAM/CRAM file that can be used as a control. I discordand reads are found in any control, the candidate inversion is rejected. A file with the '.list' suffix is interpreted as a file containing the path to the bam files.")
	private Path controlFiles = null;
	private final int min_clip_size=1;
	private final int min_insert_size=1_000;
	
	private static class BreakEnd {
		int start;
		int end;
		boolean negativeStrand;
		}
	
	private static class JunctionCount {
		int count_split_read_proof=0;
		int count_overlap_properly_paired = 0;
		int count_overlap_inversion = 0;
		int count_overlap_other=0;
		boolean isValid() {
			return true;
			}
		void json(JsonWriter j) throws IOException {
			j.beginObject();
			j.name("valid");j.value(isValid());
			j.name("split_read");j.value(count_split_read_proof);
			j.name("overlap_properly_paired");j.value(count_overlap_properly_paired);
			j.name("overlap_inversion");j.value(count_overlap_inversion);
			j.name("overlap_other");j.value(count_overlap_other);
			j.endObject();
			}
		}
	
	private static class Inversion extends AbstractLocatable {
		String contig;
		int start;
		int end;
		final JunctionCount bnd1=new JunctionCount();
		final JunctionCount bnd2=new JunctionCount();
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
			return bnd1.isValid() && bnd2.isValid();
			}
		}
	
	
	private static class ControlBam implements Closeable {
		private final SamReader sr;
		private final SAMFileHeader hdr;
		private final SAMSequenceDictionary dict;
		private final ContigNameConverter converter;
		ControlBam(final SamReader sr) {
			this.sr = sr;
			this.hdr=sr.getFileHeader();
			this.dict = SequenceDictionaryUtils.extractRequired(hdr);
			this.converter= ContigNameConverter.fromOneDictionary(dict);
			}
		@Override
		public void close() throws IOException {
			sr.close();
			}
		}
	
	
	
	private Path findReferenceForBam(final Path bamPath) {
		final Optional<Path> optFasta= SequenceDictionaryUtils.findReferenceForBam(bamPath, faidx_list);
		return optFasta.orElseThrow(()-> new IllegalArgumentException("Cannot find REFerence for "+bamPath+"."));
		}
	
		
	private boolean discordantFoundInControl(final List<ControlBam> ctrls,final Inversion inv) {
		for(ControlBam control:ctrls) {
			final String ctg = control.converter.apply(inv.getContig());
			if(StringUtils.isBlank(ctg)) continue;
			for(int side=0;side<2;side++) {
				final int pos = (side==0?inv.start:inv.end);
				try(CloseableIterator<SAMRecord> iter = control.sr.queryOverlapping(ctg,Math.max(1, pos-min_insert_size),pos+min_insert_size)) {
					while(iter.hasNext()) {
						final Read rec= new SAMRecordRead(iter.next());
						
						}
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
	private void countOverlappingReads(final SamReader sr,final Inversion inv) {
		for(int side=0;side< 2;++side) {
			final int pos = (side==0?inv.start:inv.end);
			try(CloseableIterator<SAMRecord> r= sr.queryOverlapping(inv.getContig(), Math.max(1, pos-min_insert_size),pos+min_insert_size)) {
				while(r.hasNext()) {
					final Read rec= new SAMRecordRead(r.next());
					if(rec.isIgnorable()) continue;
					if(!rec.isPairedRead()) continue;
					if(rec.isMateUnmapped()) continue;
					if(!rec.isReadAndMateMappedOnSameChromosome()) continue;
					
					if(rec.getMateAlignmentStart()<=pos) continue;
					if(rec.getMateAlignmentStart()<=rec.getEnd()) continue;
					
					//read must be before pos
					if(rec.getEnd()>= pos) continue;
					// mate must be after pos
					if(rec.getMateAlignmentStart()<=pos) continue;

					final JunctionCount bnd = (side==0?inv.bnd1:inv.bnd2);
					
					if(rec.isPlusStrand() &&  rec.isMinusStrand()) {
						bnd.count_overlap_properly_paired++;
						}
					else if(side==0 && rec.isPlusStrand() && rec.isMatePlusStrand()) {
						bnd.count_overlap_inversion++;
						}
					else if(side==1 && rec.isMinusStrand() && rec.isMateMinusStrand()) {
						bnd.count_overlap_inversion++;
						}
					else
						{
						bnd.count_overlap_other++;
						}
					}
				}
			}
		}
	
	private boolean is_ff_rr_read(final Read rec1) {
		if(rec1.isIgnorable(this.min_mapq)) return false;
		if(!rec1.isPairedRead()) return false;
		if(rec1.isMateUnmapped()) return false;
		if(rec1.getMateReferenceIndex()!=rec1.getReferenceIndex()) return false;
		return true;
		}
	
	private void fillForwardCandidates(final Read rec1,final List<BreakEnd> candidates) {
		if(!is_ff_rr_read(rec1)) return ;
		if(rec1.isMinusStrand()) return;
	
		if(!(rec1.getEnd() <rec1.getMateAlignmentStart())) return;

	
		System.err.println("A0 ! "+rec1);
		
		if(!rec1.getCigar().isRightClipped()) return; // expect a clip
		if(rec1.getCigar().getLastCigarElement().getLength()< min_clip_size) return;
		for(Read suppl: rec1.getOtherCanonicalAlignments()) {
			
			if(!suppl.contigsMatch(rec1)) continue;
			System.err.println("A1 ! "+rec1+" "+suppl);
			if(suppl.isPlusStrand()) continue;
			System.err.println("A2 ! "+rec1+" "+suppl);
			if(!(rec1.getEnd() < suppl.getStart() &&  rec1.getMateAlignmentEnd()< suppl.getStart())) continue;
			
			final BreakEnd bnd=new BreakEnd();
			bnd.start = rec1.getEnd()+1;
			bnd.end = suppl.getStart()-1;
			bnd.negativeStrand=false;
			System.err.println("OK ! ");
			candidates.add(bnd);
			}
		}
	
	private void fillReverseCandidates(final Read rec1,final List<BreakEnd> candidates) {
		if(!is_ff_rr_read(rec1)) return ;
		if(rec1.isPlusStrand()) return;

		if(!rec1.getCigar().isLeftClipped()) return ; // expect a clip
		if(rec1.getCigar().getFirstCigarElement().getLength()< min_clip_size) return;
		if(!rec1.isMateMapped()) return ;
		if(rec1.getReferenceIndex()!=rec1.getMateReferenceIndex())return;
		if(!rec1.isMateMinusStrand()) return;
		if(!(rec1.getMateAlignmentEnd() < rec1.getStart())) return ;

		
		for(Read suppl: rec1.getOtherCanonicalAlignments()) {
			if(!suppl.contigsMatch(rec1)) continue;
			if(suppl.isMinusStrand()) continue;
			if(!(suppl.getEnd() < rec1.getMateAlignmentStart() && suppl.getEnd() < rec1.getStart())) continue;
			final BreakEnd bnd=new BreakEnd();
			bnd.start = suppl.getEnd()+1;
			bnd.end = rec1.getStart()-1;
			bnd.negativeStrand=true;
			
			candidates.add(bnd);
			}
		}

	private boolean are_close(int pos1,int pos2) {
		return Math.abs(pos1-pos2)<=this.max_distance;
		}
	
	@Override
	public int doWork(final List<String> args) {
		final List<ControlBam> control_readers = new ArrayList<>()  ;

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
					control_readers.add(new ControlBam(srf
							.referenceSequence(findReferenceForBam(path))
							.open(path)
							));
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
					final String buildName =  SequenceDictionaryUtils.getBuildName(dict).orElse(null);

					
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
							
							String prevContig = null;
							final List<BreakEnd> candidates = new ArrayList<>(100_000);
							/** loop over all reads */
							try(CloseableIterator<SAMRecord> iter1=  makeIterator(samReader1)) {
								for(;;) {
									final Read rec1 = iter1.hasNext()?new SAMRecordRead(iter1.next()):null;
									if(rec1!=null && rec1.isIgnorable(this.min_mapq)) continue;
									
									if(rec1==null || (prevContig!=null && !rec1.getContig().equals(prevContig))) {
										if(!candidates.isEmpty()) {
											Collections.sort(candidates,(A,B)->Integer.compare(A.start, B.start));
											
											int i=0;
											while(i< candidates.size()) {
												final BreakEnd bnd1 = candidates.remove(0);
												final Inversion inv = new Inversion();
												inv.contig = prevContig;
												
												final JunctionCount count_i = bnd1.negativeStrand?inv.bnd2:inv.bnd1;
												count_i.count_split_read_proof++;
												final List<Integer> starts = new ArrayList<>();
												starts.add(bnd1.start);
												final List<Integer> ends = new ArrayList<>();
												ends.add(bnd1.end);
												int j=0;
												while(j< candidates.size()) {
													final BreakEnd bnd2 = candidates.get(j);
													if(inv.start+this.max_distance < bnd2.start) break;
													if(are_close(bnd1.start, bnd2.start) && are_close(bnd1.start, bnd2.start) ) {
														starts.add(bnd2.start);
														ends.add(bnd2.end);
														final JunctionCount count_j = bnd2.negativeStrand?inv.bnd2:inv.bnd1;
														count_j.count_split_read_proof++;
														candidates.remove(j);
														}
													else
														{
														++j;
														}
													}
												inv.start =(int) new Median().evaluate(starts.stream().mapToDouble(P->P.doubleValue()).sorted().toArray());
												inv.end = (int)new Median().evaluate(ends.stream().mapToDouble(P->P.doubleValue()).sorted().toArray());
												
												countOverlappingReads(samReader2,inv);
												
												if(discordantFoundInControl(null,inv)) continue;
												if(!inv.isValid()) continue;
												js.beginObject();
												js.name("sample"); js.value(sampleName);
												js.name("bam"); js.value(String.valueOf(inputBam));
												js.name("reference");js.value(String.valueOf(referenceFaidx));
													
												if(buildName!=null) {
													js.name("build");js.value(buildName);
													}
												
												final String url = hyperlink.apply(inv).orElse(null);
												if(url!=null) {
													js.name("href");js.value(url);
													}
												
												js.name("left"); inv.bnd1.json(js);
												js.name("right"); inv.bnd2.json(js);
												
												js.endObject();
												js.flush();

												}
											
											candidates.clear();
											} //end candidates empty
										if(rec1==null) break;
										prevContig = rec1.getContig();
										} //end of contig==null
									
									
									if(rec1.isPlusStrand()) {
										fillForwardCandidates(rec1,candidates);
										}
									else if(rec1.isMinusStrand()) {
										fillReverseCandidates(rec1,candidates);
										}
									}//end for(;;
								
								} // end iterator
							
							
							js.endArray();
							js.flush();
					
							}//end js writer
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
			for(ControlBam sr:control_readers) {
				try {sr.close(); } catch(IOException err2) {}
				}
			}
		}

	
	public static void main(final String[] args)
		{
		new SamInversions().instanceMainWithExit(args);
		}
	}
