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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Predicate;

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
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.BAMUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.AbstractLocatable;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
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
	description="Scan inversions in SAM using split and FF/RR reads on same strand.",
	keywords={"sam","bam","sv","inversion"},
	modificationDate="20261108",
	creationDate="20200609",
	generate_doc=true
	)
public class SamInversions extends Launcher
	{
	private static final Logger LOG = Logger.of(SamInversions.class);
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--vcf-out","--out-vcf","--write-vcf"},description="save output in VCF format to that file")
	private Path outputVcfFile = null;
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
	@Parameter(names={"--no-clip"},description="Ignore the clipping status of the SAMRecord, just look at supplementary reads without looking at the clipping")
	private boolean ignore_clipping_status_flag = false;
	@Parameter(names={"--debug"},description="debug")
	private boolean debug_flag = false;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate= new WritingVariantsDelegate();
	
	
	private final int min_insert_size=1_000;
	
	private static enum EvidenceType {
		both,two5,two3,one5,one3
		}
	
	private static class Evidence {
		int start;
		int end;
		EvidenceType type;
		@Override
		public String toString() {
			return String.valueOf(start)+"-"+end+"("+type.name()+")";
			}
		}
	
	private static class JunctionCount {
		int count_split_read_proof=0;
		int count_overlap_properly_paired = 0;
		int count_overlap_inversion = 0;
		int count_overlap_other=0;
		boolean isHighQual() {
			return count_split_read_proof>0 && count_overlap_inversion>0;
		}
		boolean isValid() {
			return true;
			}
		void json(JsonWriter j) throws IOException {
			j.beginObject();
			j.name("valid");j.value(isValid());
			j.name("high_qual");j.value(isHighQual());
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
		final Counter<EvidenceType> evidence_count=new Counter<>();

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
		boolean isHighQual() {
			return bnd1.isHighQual() && bnd2.isHighQual() && evidence_count.getMaxCount().orElse(0L)>9;
		}
		boolean isValid() {
			return bnd1.isValid() && bnd2.isValid() && this.evidence_count.getTotal()>10L;
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
	
		
	private int countDiscordantFoundInControl(final List<ControlBam> ctrls,final  Inversion inv) {
		int n=0;
		for(ControlBam control:ctrls) {
			final String ctg = control.converter.apply(inv.getContig());
			boolean side5=false;
			boolean side3=false;
			if(StringUtils.isBlank(ctg)) continue;
			for(int side=0;side<2;++side) {
				final int pos = (side==0?inv.getStart():inv.getEnd());
				try(CloseableIterator<SAMRecord> iter = control.sr.queryOverlapping(ctg,Math.max(1, pos-this.max_distance),pos+max_distance)) {
					while(iter.hasNext()) {
						final SAMRecord rec= iter.next();
						if(!SAMRecordDefaultFilter.accept(rec)) continue;
						if(rec.getMateUnmappedFlag()) continue;
						if(!rec.getReferenceIndex().equals(rec.getMateReferenceIndex())) continue;
						if(side==0 && rec.getCigar().isRightClipped() && are_close(rec.getAlignmentEnd(),pos)) {
							side5=true;
							break;
							}
						if(side==0 && rec.getCigar().isLeftClipped() && are_close(rec.getAlignmentStart(),pos)) {
							side5=true;
							break;
							}
						if(side==1 && rec.getCigar().isRightClipped() && are_close(rec.getAlignmentEnd(),pos)) {
							side3=true;
							break;
							}
						if(side==1 && rec.getCigar().isLeftClipped() && are_close(rec.getAlignmentStart(),pos)) {
							side3=true;
							break;
							}
						}
					}
				}
			if(side5 && side3) n++;
			}
		return n;
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
	
	private String toString(final SAMRecord rec) {
		return rec.getContig()+":"+rec.getStart()+"-"+rec.getEnd()+"("+(rec.getReadNegativeStrandFlag()?"-":"+")+"):"+rec.getCigarString()+" <=> "+
					rec.getMateReferenceName()+":"+rec.getMateAlignmentStart()+"-"+BAMUtils.getMateEnd(rec)+"("+(rec.getMateNegativeStrandFlag()?"-":"+")+")"+
					(rec.hasAttribute(SAMTag.MC)? ":"+rec.getAttribute(SAMTag.MC):"")+
					(rec.hasAttribute(SAMTag.SA)? " SA:"+rec.getAttribute(SAMTag.SA):"")
					;
		}

	/** fetch reads in a given region */
	private void countOverlappingReads(final SamReader sr,final Inversion inv) {
		for(int side=0;side< 2;++side) {
			final int pos = (side==0?inv.start:inv.end);
			try(CloseableIterator<SAMRecord> r= sr.queryOverlapping(inv.getContig(), Math.max(1, pos-min_insert_size),pos+min_insert_size)) {
				while(r.hasNext()) {
					final SAMRecord rec= r.next();
					if(!SAMRecordDefaultFilter.accept(rec, this.min_mapq)) continue;
					if(!rec.getReadPairedFlag()) continue;
					if(rec.getMateUnmappedFlag()) continue;
					
					if(!(rec.getAlignmentEnd()<pos && pos < rec.getMateAlignmentStart())) continue;
					
					final JunctionCount bnd = (side==0?inv.bnd1:inv.bnd2);
					
					if(!rec.getReadNegativeStrandFlag() &&  rec.getMateNegativeStrandFlag()) {
						bnd.count_overlap_properly_paired++;
						}
					else if(side==0 && !rec.getReadNegativeStrandFlag() && !rec.getMateNegativeStrandFlag()) {
						bnd.count_overlap_inversion++;
						}
					else if(side==1 && rec.getReadNegativeStrandFlag() && rec.getMateNegativeStrandFlag()) {
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
	
	
	

	private boolean are_close(int pos1,int pos2) {
		return Math.abs(pos1-pos2)<=this.max_distance;
		}
	
	@Override
	public int doWork(final List<String> args) {
		final List<ControlBam> control_readers = new ArrayList<>()  ;
		VariantContextWriter vcw = null;
		
		if(this.max_distance<0) {
			LOG.error("max size inversion must be >0");
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
							
							if(this.outputVcfFile!=null) {
								vcw = this.writingVariantsDelegate.dictionary(dict)
										.open(this.outputVcfFile);
								final Set<VCFHeaderLine> metaData = new HashSet<>();
								metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY,true));
								metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY,true));
								metaData.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer, "SV length"));
								metaData.add(new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String, "SV type"));
								for(EvidenceType et: EvidenceType.values()) {
									metaData.add(new VCFInfoHeaderLine(et.name(),1, VCFHeaderLineType.Integer, "Evidence type count "+et.name()));
									}
								final VCFHeader vcfHeader=new VCFHeader(metaData, Collections.singletonList(sampleName));
								vcfHeader.setSequenceDictionary(dict);
								JVarkitVersion.getInstance().addMetaData(this, vcfHeader);
								vcw.writeHeader(vcfHeader);
								}
							
							
							
							String prevContig = null;
							final List<Evidence> candidates = new ArrayList<>(100_000);
							/** loop over all reads */
							try(CloseableIterator<SAMRecord> iter1=  makeIterator(samReader1)) {
								for(;;) {
									final SAMRecord rec1 = iter1.hasNext()?iter1.next():null;
									if(rec1!=null && (!SAMRecordDefaultFilter.accept(rec1,this.min_mapq) || !accept_reads.test(rec1))) {
										continue;
									}
									
									if(rec1==null || prevContig==null || !rec1.getContig().equals(prevContig)) {
										if(prevContig!=null && !candidates.isEmpty()) {
											Collections.sort(candidates,(A,B)->Integer.compare(A.start, B.start));
											
											int i=0;
											while(i< candidates.size()) {
												final Evidence bnd1 = candidates.remove(0);
												final Inversion inv = new Inversion();
												inv.contig = prevContig;
												inv.evidence_count.incr(bnd1.type);
												final List<Integer> starts = new ArrayList<>();
												starts.add(bnd1.start);
												final List<Integer> ends = new ArrayList<>();
												ends.add(bnd1.end);
												int j=0;
												while(j< candidates.size()) {
													final Evidence bnd2 = candidates.get(j);
													if(are_close(bnd1.start, bnd2.start) && are_close(bnd1.end, bnd2.end) ) {
														starts.add(bnd2.start);
														ends.add(bnd2.end);
														inv.evidence_count.incr(bnd2.type);
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
												
												//if(discordantFoundInControl(null,inv)) continue;
												if(!inv.isValid()) continue;
												js.beginObject();
												js.name("sample"); js.value(sampleName);
												js.name("bam"); js.value(String.valueOf(IOUtils.toRealPath(inputBam)));
												js.name("reference");js.value(String.valueOf(IOUtils.toRealPath(referenceFaidx)));
												js.name("chromosome");js.value(inv.getContig());
												if(buildName!=null) {
													js.name("build");js.value(buildName);
													}
												js.name("start");js.value(inv.getStart());
												js.name("end");js.value(inv.getEnd());
												js.name("length");js.value(inv.getLengthOnReference());
												js.name("interval");js.value(inv.getContig()+":"+inv.getStart()+"-"+inv.getEnd());
												
												for(EvidenceType et:EvidenceType.values()) {
													js.name("count_"+et.name());
													js.value(inv.evidence_count.count(et));
													}
												js.name("count_sum");
												js.value(inv.evidence_count.getTotal());

												
												final String url = hyperlink.apply(inv).orElse(null);
												if(url!=null) {
													js.name("href");js.value(url);
													}
												
												final int count_overlapping_in_controls = countDiscordantFoundInControl(control_readers, inv);
												js.name("n-controls");js.value(count_overlapping_in_controls);

												
												js.name("high_qual");js.value(inv.isHighQual() && count_overlapping_in_controls==0);
												js.name("left"); inv.bnd1.json(js);
												js.name("right"); inv.bnd2.json(js);
												
												
												
												js.endObject();
												js.flush();

												
												if(vcw!=null) {
													final List<Allele> alleles = Arrays.asList(Allele.REF_N,Allele.create("<INV>",false));
													
													final VariantContextBuilder vcb = new VariantContextBuilder(getProgramName(), inv.getContig(), inv.getStart(), inv.getEnd(),alleles)
																.attribute(VCFConstants.END_KEY,inv.end)
																.attribute("SVLEN",inv.getLengthOnReference())
																.attribute(VCFConstants.SVTYPE, "INV")
																.genotypes(Collections.singletonList(new GenotypeBuilder(sampleName, alleles).make()))
																;
													for(EvidenceType et: EvidenceType.values()) {
														vcb.attribute(et.name(),(int) inv.evidence_count.count(et));
														}
													
													vcw.add(vcb.make());
													}
												
												}
											
											candidates.clear();
											} //end candidates empty
										if(rec1==null) break;
										prevContig = rec1.getContig();
										} //end of contig==null
									
									
									if(!rec1.getReadPairedFlag()) continue;
									if(rec1.getMateUnmappedFlag()) continue;
									// mate mapped on same contrig
									if(!rec1.getReferenceIndex().equals(rec1.getMateReferenceIndex())) continue;

									//read before mate */
								 	
									/* read and mate overlap the break points  */
									//
									// ====||======>             <==||===============
									// ====||==>             <======||====
									if(!rec1.getReadNegativeStrandFlag() &&
										(this.ignore_clipping_status_flag ||  rec1.getCigar().isRightClipped()) &&
										rec1.getMateNegativeStrandFlag() &&
										rec1.getEnd() < rec1.getMateAlignmentStart() &&
										SAMUtils.getOtherCanonicalAlignments(rec1)
											.stream()
											.anyMatch(suppl->suppl.contigsMatch(rec1) && suppl.getReadNegativeStrandFlag() && rec1.getEnd() < suppl.getStart() && are_close(suppl.getEnd(), rec1.getMateAlignmentStart()))) {
										final SAMRecord mateRecord = samReader2.queryMate(rec1);
										if(mateRecord!=null &&
											(this.ignore_clipping_status_flag || mateRecord.getCigar().isLeftClipped()) &&
											SAMUtils.getOtherCanonicalAlignments(mateRecord)
												.stream()
												.anyMatch(suppl->suppl.contigsMatch(rec1) && !suppl.getReadNegativeStrandFlag() && suppl.getEnd() < mateRecord.getStart() && are_close(rec1.getAlignmentEnd(), suppl.getStart()))) {
													final Evidence proof=new Evidence();
													proof.start = rec1.getEnd()+1;
													proof.end = mateRecord.getStart()-1;
													proof.type = EvidenceType.both;
													if(debug_flag) {
														LOG.info("found proof double "+toString(rec1)+" "+proof);
														}
													candidates.add(proof);
													
												}
											
										continue;
										}
										
										// read overlap 5' breakpoint
										// *****
										// =====||===>    <==========         ||
										// =====||        ==========>     <===||
										if( !rec1.getReadNegativeStrandFlag() &&
											!rec1.getMateNegativeStrandFlag() &&
											(this.ignore_clipping_status_flag || rec1.getCigar().isRightClipped()) &&
											rec1.getEnd() < rec1.getMateAlignmentStart()
											) {
											for(SAMRecord suppl: SAMUtils.getOtherCanonicalAlignments(rec1)) {
												if(suppl.contigsMatch(rec1) && suppl.getReadNegativeStrandFlag()  && BAMUtils.getMateEnd(suppl) < suppl.getStart()) {
													final Evidence proof=new Evidence();
													proof.start = rec1.getEnd()+1;
													proof.end = suppl.getEnd();
													proof.type = EvidenceType.two5;
													if(debug_flag) {
														LOG.info("found proof 5' "+toString(rec1)+" "+proof);
													}
													candidates.add(proof);
													}
												}
											}
										
										
										// read overlap 3' breakpoint
										//                              *****
										//    ||  ========>       <===||=====
										//    ||===>       <========  ||=====
										if( rec1.getReadNegativeStrandFlag() &&
											rec1.getMateNegativeStrandFlag() &&
											(this.ignore_clipping_status_flag || rec1.getCigar().isLeftClipped()) &&
											BAMUtils.getMateEnd(rec1) < rec1.getStart()
											) {
											for(SAMRecord suppl: SAMUtils.getOtherCanonicalAlignments(rec1)) {
												if(suppl.contigsMatch(rec1) && !suppl.getReadNegativeStrandFlag() && suppl.getEnd() < rec1.getMateAlignmentStart()) {
													final Evidence proof=new Evidence();
													proof.start = suppl.getStart();
													proof.end = rec1.getStart()-1;
													proof.type = EvidenceType.two3;
													
													if(debug_flag) {
														LOG.info("found proof 3' "+toString(rec1)+" "+proof);
													}
													candidates.add(proof);
													}
												}
											}
										
										
										// read overlap 5' breakpoint
										//             *****
										//         <==||====              ||
										//            ||               ===||
										if( rec1.getReadNegativeStrandFlag() &&
												(this.ignore_clipping_status_flag || rec1.getCigar().isRightClipped()) &&
												BAMUtils.getMateEnd(rec1) < rec1.getStart() &&
												rec1.getProperPairFlag()
												) {
												for(SAMRecord suppl: SAMUtils.getOtherCanonicalAlignments(rec1)) {
													if(suppl.contigsMatch(rec1) && !suppl.getReadNegativeStrandFlag() && rec1.getEnd()<suppl.getStart()) {
														final Evidence proof=new Evidence();
														proof.start = rec1.getEnd()+1;
														proof.end = suppl.getEnd()-1;
														proof.type = EvidenceType.one5;
														
														if(debug_flag) {
															LOG.info("found proof2 5' "+toString(rec1)+" "+proof);
														}
														candidates.add(proof);
														}
													}
												}

										
										// read overlap 3' breakpoint
										//                              *****
										//    ||                   ===||=====> <====
										//    ||<===                  ||=====> <====
										if( !rec1.getReadNegativeStrandFlag() &&
												(this.ignore_clipping_status_flag || rec1.getCigar().isLeftClipped()) &&
												rec1.getEnd() < rec1.getMateAlignmentStart() &&
												rec1.getProperPairFlag()
												) {
												for(SAMRecord suppl: SAMUtils.getOtherCanonicalAlignments(rec1)) {
													if(suppl.contigsMatch(rec1) && suppl.getReadNegativeStrandFlag() && suppl.getEnd() < rec1.getStart()) {
														final Evidence proof=new Evidence();
														proof.start = suppl.getStart();
														proof.end = rec1.getStart()-1;
														proof.type = EvidenceType.one3;
														if(debug_flag) {
															LOG.info("found proof2 3' "+toString(rec1)+" "+proof);
														}
														candidates.add(proof);
														}
													}
												}
										
										
							
									}//end for(;;
						
								} // end iterator
							
							
							js.endArray();
							js.flush();
					
							}//end js writer
						} // end cram2
					} //end cram1
				
			return 0;
			} 
		catch (final Throwable e) {
			LOG.error(e);
			return -1;
			}
		finally {
			if(vcw!=null) {
				vcw.close();
				}
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
