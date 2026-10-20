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
package com.github.lindenb.jvarkit.tools.structvar;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.PeekIterator;
import htsjdk.samtools.util.StringUtil;
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
	@Parameter(names={"-R","--reference"},description=CRAM_INDEXED_REFENCE,required = true)
	private Path referenceFaidx = null;
	@Parameter(names={"--max-size"},description="max distance between read and mate.",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int max_distance = 1_000 ;
	@Parameter(names={"-B","--bed","-r","--rgn"},description=IntervalListProvider.OPT_DESC,splitter= NoSplitter.class,converter=IntervalListProvider.StringConverter.class)
	private IntervalListProvider intervallistProvider = null;
	@Parameter(names={"-s","-supporting"},description="Don't print the variant if max FORMAT/DP < 's'")
	private int min_supporting_reads = 1;
	@Parameter(names={"--mapq"},description="min MAPQ")
	private int mapq = 30;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();

	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.max_distance<100) {
			LOG.error("max size inversion must be >=100");
			return -1;
			}
		try {
			final Allele REF = Allele.create("N", true);
			final Allele INV = Allele.create("<INV>", false);
			final String input = oneAndOnlyOneFile(args);
			final Path inputPath = Paths.get(input);
			final SamReaderFactory srf =  super.createSamReaderFactory().
					validationStringency(ValidationStringency.LENIENT).
					referenceSequence(this.referenceFaidx)
					;
			try(ReferenceSequenceFile referenceSequenceFile = new IndexedFastaSequenceFile(this.referenceFaidx)) {
				try(SamReader samReader1= srf.open(inputPath)) {
					if(!samReader1.hasIndex()) {
						LOG.error("file is not indexed "+input);
						return -1;
						}
					final SAMFileHeader samHeader = samReader1.getFileHeader();
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samHeader);
					final String sampleName= samHeader.getReadGroups().stream().
							map(RG->RG.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(inputPath));

					//open a second reader to fetch the mate
					try(SamReader samReader2= srf.open(inputPath)) {
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
	
							final List<Interval> candidate_intervals = new ArrayList<>();
							
							String prev_contig = null;
							/** loop over eall reads */
							try(SAMRecordIterator iter1= samReader1.iterator()) {
								for(;;) {
									final SAMRecord rec1 = iter1.hasNext()?iter1.next():null;
									if(rec1!=null && !SAMRecordDefaultFilter.accept(rec1)) {
										continue;
										}
									if(rec1==null || !rec1.getContig().equals(prev_contig)) {
										for(Interval rgn : candidate_intervals) {
											final VariantContextBuilder vcb = new VariantContextBuilder();
											
											
											}
										
										if(rec1==null) break;
										candidate_intervals.clear();
										prev_contig = rec1.getContig();
										LOG.info("contig "+ rec1.getContig());
										}
									final Cigar cigar = rec1.getCigar();
									
									
									
									
									for(SAMRecord suppl:SAMUtils.getOtherCanonicalAlignments(rec1)) {
										if(!suppl.contigsMatch(rec1)) continue;
										if(rec1.getReadNegativeStrandFlag()== suppl.getReadNegativeStrandFlag()) continue;
										if(!rec1.withinDistanceOf(suppl, this.max_distance)) continue;
										
										
										// #########oooooooooooooooooooooo#######
										//      =======>
										//      <=======
										
										if(cigar.isRightClipped() && 
											rec1.getEnd()< suppl.getStart()
											) {
											final Interval bounds= new Interval(
													rec1.getContig(),
													rec1.getEnd()+1,
													suppl.getEnd(),
													rec1.getReadNegativeStrandFlag(),
													rec1.getReadName()
													);
											candidate_intervals.add(bounds);
											}
										// #########oooooooooooooooooooooo###########
										//                           =======>
										//                           <=======
										if(cigar.isLeftClipped() && 
												suppl.getEnd() < rec1.getStart()
												) {
												final Interval bounds= new Interval(
														rec1.getContig(),
														suppl.getStart(),
														rec1.getStart()-1,
														rec1.getReadNegativeStrandFlag(),
														rec1.getReadName()
														);
												candidate_intervals.add(bounds);
												}
										}
									// in inversion read and it's mate will be in same strand
									if(rec1.getReadPairedFlag() &&
											!rec1.getMateUnmappedFlag() &&
											rec1.getReferenceIndex().equals(rec1.getMateReferenceIndex()) &&
											rec1.getAlignmentStart() < rec1.getMateAlignmentStart() &&
											rec1.getReadNegativeStrandFlag()==rec1.getMateNegativeStrandFlag() &&
											Math.abs(rec1.getInferredInsertSize()) <= this.max_distance
											) {
										final SAMRecord mateRec = samReader2.queryMate(rec1);
										if(mateRec!=null) {
											final Interval interval = new Interval(rec1.getContig(),
													Math.min(rec1.getStart(), mateRec.getStart()),
													Math.max(rec1.getEnd(), mateRec.getEnd()),
													false,
													rec1.getReadName()
													);
											// TODO
											}
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
		
		}

	
	public static void main(final String[] args)
		{
		new SamInversions().instanceMainWithExit(args);
		}
	}
