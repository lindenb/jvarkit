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
package com.github.lindenb.jvarkit.tools.retrocopy;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.Paranoid;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.stream.HtsCollectors;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/*

STAR --genomeDir referenceDir --readFilesIn S4_L002_R1_001.fastq.gz S4_L002_R2_001.fastq.gz \
	--readFilesCommand gunzip -c  --outStd SAM --outReadsUnmapped None --outSAMattributes All |\
	
	

A00797:7:HJTYVDMXX:2:1101:14181:2378	99	9	123795763	255	17S68M26N66M	=	123795764	178	GTCCAAAAGGCAGCATTAAGTAGAGGGGATGACATGTGCAAAGGCACGGAGGTAGGAAAGCACTGGACATGCCAGACCATGGCTAGCACAAGATGAAGTTGGAGAAATGATAGCAGTAATAAAATAATAGGTTTAATTCTGCACTATGGAA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFFFF:	NH:i:1	HI:i:1	AS:i:267	nM:i:0	NM:i:0	MD:Z:134	jM:B:c,0	jI:B:i,123795831,123795856	MC:Z:67M26N84M
A00797:7:HJTYVDMXX:2:1101:14181:2378	147	9	123795764	255	67M26N84M	=	123795763	-178	AGTAGAGGGGATGACATGTGCAAAGGCACGGAGGTAGGAAAGCACTGGACATGCCAGACCATGGCTAGCACAAGATGAAGTTGGAGAAATGATAGCAGTAATAAAATAATAGGTTTAATTCTGCACTATGGAACATGTAGTATTTTCAAGA	:FFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFF,FFFF:FFFFFFFFFFFFFFF:FFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFF,FF:F:FFFFFF	NH:i:1	HI:i:1	AS:i:267	nM:i:0	NM:i:0	MD:Z:151	jM:B:c,0	jI:B:i,123795831,123795856	MC:Z:17S68M26N66M
A00797:7:HJTYVDMXX:2:1101:25120:2378	163	12	7312757	255	64M69N69M18S	=	7312981	368	TTTGCCTGGTCTGTGGCAAGAGCCTCCCACCTGCCTCCCACCTGCCTCCCACCTGTCCTCCCACCTCTCCTCCCACCTGCCTCCCACCTGTCCGCCCACCGGCCTCCCACCTGCCACCCACCAGGCCGCCACCGCTCAACCAGCCAGGCCC	FFFFFF,FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFF:FF,:,FF,:FFFFF,:FF,FFFFFF,F,F:,FFF,FF,:::,:::,:F,::F,,,F,FF,F,F,,FF,FF,:,:,,:F,,,F,:,F,:	NH:i:1	HI:i:1	AS:i:259	nM:i:8	NM:i:8	MD:Z:93T6T14T6T1C1T0C0G4	jM:B:c,2	jI:B:i,7312821,7312889	MC:Z:144M7S
A00797:7:HJTYVDMXX:2:1101:3974:2394	83	16	33884936	255	103M1107N42M6S	=	33884745	-1443	CAATCATTGAATGGAATGGAAAGGAATCGTCATCAAATTAAATAGAATAGAATCATCGAATGGAATCTAATGGAATCATCATCGAATGGAATCCAATGGAATCATCGAATGGAATCATCATCGAACGGAATCCAATGGAATATAATCCGGC	FF,FFF:FFF,FFFFFFFFFFF:FFFFFFFFFFFFF:F,FFFFFFFF:FFF,:F,FFFFFFFFF,FFFFF,FFF:FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NH:i:1	HI:i:1	AS:i:283	nM:i:1	NM:i:1	MD:Z:38G106	jM:B:c,0	jI:B:i,33885039,33886145	MC:Z:151M
A00797:7:HJTYVDMXX:2:1101:26214:2394	163	X	137457323	255	116M80N35M	=	137457628	436	GCCCATTCTTTGGGTCTAGTGCTGGAACTGTTATTAACCCACTCATTTTGCTTTGTGCCACTATTAGGCCATTCTCACTCTGCTAATAAAGACATACCCAAGACTGGGTAATTTACAAAGGCACGTCTTACATGGTGACATATACCCTGAG	FFFFFFFFFF:FFFFFF:F:FFFFFFFFFFFFFFFFFFFFFF:F:FF,FFFFFFFFFFFFFF:F,FFFFFFFFFFFFFFF:FFF:FFFFFFFFF:FFFFFFFF:FF,FF:FFFFFF:FFFFF,FFFFF:F,FFFFFFFFFFFFFFF,:,:F	NH:i:1	HI:i:1	AS:i:270	nM:i:1	NM:i:0	MD:Z:151	jM:B:c,0	jI:B:i,137457439,137457518	MC:Z:131M20S
A00797:7:HJTYVDMXX:2:1101:19153:2409	163	12	132640585	255	2S101M175642N48M	=	132816467	175969	CCACACACCACAATGCACACTACATACAAACCACTACACACACTACACACACACACCACTACACGCACACTACACACACCACTACATATACCACTACACACACCACACACACTTCCCGGGGCACTCCCTATACCCATTACACCCACACAAC	FFFFFF,FFFF,FFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFFFFF,FFF:FFFF,FF,:F:FFFFFFFFFFF,FFFFFF:F:F:F,:,,FFF,::FF:,FF	NH:i:1	HI:i:1	AS:i:210	nM:i:7	NM:i:6	MD:Z:6A0C2C1A107C19A8	jM:B:c,0	jI:B:i,132640686,132816327	MC:Z:87M64S 
 
 */
/**
BEGIN_DOC



END_DOC

*/
@Program(name="starretrocopy",
description="Scan retrocopies from the star-aligner/bwa output",
keywords={"sam","bam","cigar","clip","sv","retrocopy","star"},
creationDate="20190710",
modificationDate="20191008"
)
public class StarRetroCopy extends Launcher
	{
	private static final Logger LOG = Logger.build(StarRetroCopy.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-gtf","--gtf"},description=GtfReader.OPT_DESC,required=true)
	private Path gtfPath = null;
	@Parameter(names={"--mapq","-mapq"},description="Min mapping quality")
	private int min_read_mapq = 1;
	@Parameter(names={"--bam","--save-bam"},description="Optional: save matching reads in this bam file. Can be used to visualize the evidences in IGV.")
	private Path saveBamTo = null;
	@Parameter(names={"--min-depth","-D"},description="In a transcript one must found at least 'D' reads with a clip-length> 'min-cigar-size'.")
	private int min_depth=1;
	@Parameter(names={"--low-depth","-d"},description="Min number of reads to set FILTER=PASS.",hidden=true)
	private int low_depth_threshold=10;
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION+" Can be used for reading CRAM.")
	private Path faidx = null;
	@Parameter(names={"--no-ji","--no-star","--bwa"},description="Do not use the 'Ji' Attribute generated by STAR aligner to get the intron bounds but use the CIGAR string (e.g: output of bwa)")
	private boolean use_cigar_string = false;
	@Parameter(names={"--bai","-bai","--with-bai"},description="Use random access BAM using the bai and using the GTF data. May be slow at startup")
	private boolean use_bai = false;
	@Parameter(names={"-n","--min-cigar-size"},description="Minimal cigar element length for clipped sequence.")
	private int minCigarSize = 6;
	@Parameter(names={"--transcript-id"},description="Limit to those transcript ID (one id per line)")
	private Path transcriptIDPath=null;
	@Parameter(names={"--bedpe","-P","-J"},description="Optional. Save possible sites of insertion in this Bed-PE file.")
	private Path saveBedPeTo=null;

	
	private abstract static class AbstractLoc implements Locatable{
		int start;
		int end;
		@Override public int getStart() { return this.start;}
		@Override public int getEnd() { return this.end;}
		
		@Override
		public String toString() {
			return getContig()+":"+getStart()+"-"+getEnd();
		}
	}
	
	private static class Transcript extends AbstractLoc {
		String contig;
		String transcript_id;
		final Map<String,String> attributes = new HashMap<>();
		

		
		final List<Segment> segments = new ArrayList<>();
		@Override
		public String getContig() {
			return contig;
			}
	
		}
	
	private static class Segment extends AbstractLoc {
		Transcript transcript;
		int match=0;
		@Override
		public String getContig() {
			return this.transcript.getContig();
			}
		
		}
	
	private final IntervalTreeMap<List<Segment>> intronTreeMap = new IntervalTreeMap<>();
	private final List<Transcript> all_transcripts = new ArrayList<>();	
	private final Paranoid paranoid = Paranoid.createThrowingInstance();
	

	private static final String ATT_LOW_DEPTH_FILTER="LowQual";
	private final static String ATT_INTRONS_COUNT="COUNT_INTRONS";
	private final static String ATT_INTRONS_CANDIDATE_COUNT="MATCHING_INTRONS_COUNT";
	private final static String ATT_INTRONS_CANDIDATE_FRACTION="MATCHING_INTRONS_FRACTION";
	private final static String ATT_NOT_ALL_INTRONS="NOT_ALL_INTRONS";
	private static final String ENSEMBL_TRANSCRIPT_ATTS[]=new String[] {"gene_version","transcript_id","transcript_version","gene_name","gene_source","gene_biotype","transcript_name","transcript_source","transcript_biotype","tag","ccds_id","havana_transcript","havana_transcript_version","tag"};
	private final static String INTRON_START="BEG";
	private final static String INTRON_END="END";

	
	
	private void loadGTF(final SAMSequenceDictionary dict) throws IOException {
		final Set<String> retain_transcript_ids ;
		if(transcriptIDPath!=null)
			{
			retain_transcript_ids  = Files.
					lines(this.transcriptIDPath).
					filter(L->!L.startsWith("#")).
					filter(L->!StringUtils.isBlank(L)).
					collect(Collectors.toSet());
			}
		else
			{
			retain_transcript_ids = null;
			}
		
		final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(dict);
		final GTFCodec gtfCodec = new GTFCodec();
		final Map<String,Transcript> id2transcript = new HashMap<>(50_000);
		//Emmanuelle je t'aime contre le mur
		try(BufferedReader br=IOUtils.openPathForBufferedReading(this.gtfPath))
			{
			br.lines().
			filter(S->!S.startsWith("#")).
			filter(S->!StringUtils.isBlank(S)).
			map(S->gtfCodec.decode(S)).
			filter(S->converter.apply(S.getContig())!=null).
			filter(T->T.getType().equals("transcript")|| T.getType().equals("exon")).
			forEach(T->{
				final String transcript_id = T.getAttribute("transcript_id");
				if(StringUtils.isBlank(transcript_id)) throw new RuntimeIOException("no transcript_id in "+ T);
				
				if(retain_transcript_ids!=null && !retain_transcript_ids.contains(transcript_id)) {
					return;
					}
				
				Transcript transcript = id2transcript.get(transcript_id);
				
				if(transcript==null) {
					transcript = new Transcript();
					transcript.contig = converter.apply(T.getContig());
					transcript.transcript_id = transcript_id;
					id2transcript.put(transcript_id,transcript);
					}
				
				if(T.getType().equals("transcript")) {
					transcript.start = T.getStart();
					transcript.end = T.getEnd();
					for(final String att:ENSEMBL_TRANSCRIPT_ATTS) {
						final String v=T.getAttribute(att);
						if(StringUtils.isBlank(v)) continue;
						transcript.attributes.put(att, v);
						}
					}
				else
					{
					final Segment exon = new Segment();
					exon.transcript = transcript;
					exon.start = T.getStart();
					exon.end = T.getEnd();
					transcript.segments.add(exon);
					}
				});
			}
		id2transcript.values().forEach(T->{
			Collections.sort(T.segments,(A,B)->Integer.compare(A.start,B.start));
			final List<Segment> introns = new ArrayList<>();
			for(int i=0;i+1< T.segments.size();i++)
				{
				final Segment ex1 = T.segments.get(i);
				final Segment ex2 = T.segments.get(i+1);
				final Segment intron = new Segment();
				intron.transcript = T;
				intron.start = ex1.getEnd()+1;
				intron.end= ex2.getStart()-1;
				if(intron.start>intron.end) throw new IllegalStateException("overlapping exons "+ex1+" "+ex2+ " in "+T.transcript_id);
				introns.add(intron);
				final Interval r=new Interval(intron);
				List<Segment> L=this.intronTreeMap.get(r);

				if(L==null) {
					L=new ArrayList<>();
					this.intronTreeMap.put(r,L);
					}
					
				L.add(intron);
				}
			T.segments.clear();
			T.segments.addAll(introns);
			});
		
		this.all_transcripts.addAll(
				id2transcript.values().
					stream().
					filter(T->!T.segments.isEmpty()).
					collect(Collectors.toList())
				);
		
		
		}
	
	private boolean isWithinIntronBound(final Segment seg,final Interval interval,int side) {
		final int distance;
		if(side==0 ) //match cigar right of record
			{
			distance =  Math.abs(seg.getStart() - interval.getStart());
			}
		else // left of cigar match 3' of intron
			{
			distance =  Math.abs(interval.getEnd() - seg.getEnd());
			}
		//System.err.println(seg+" "+interval+" "+side+":"+distance);
		return distance <=10;
		}
	
	
	private void saveInsertion(final PrintWriter pw,final SAMRecord rec) {
		if(this.saveBedPeTo==null) return;
		
		final List<SimpleInterval> rec2list = new ArrayList<>();
		
		final Predicate<Locatable> shouldBeSaved = REC2->{
			if(!rec.getContig().equals(REC2.getContig()))
				{
				return true;
				}
			if(rec.overlaps(REC2)) return false;
			
			return this.intronTreeMap.getOverlapping(rec).
				stream().
				flatMap(L->L.stream()).
				map(L->L.transcript).
				noneMatch(P->P.overlaps(REC2));
			};
		
		if(rec.getReadPairedFlag() &&
			!rec.getMateUnmappedFlag() && 
			!rec.getProperPairFlag())
			{
			final int mateStart = rec.getMateAlignmentStart();
			final int mateEnd;
			if(SAMUtils.hasMateCigar(rec)) {
				mateEnd = SAMUtils.getMateAlignmentEnd(rec);
				}
			else
				{
				mateEnd = mateStart;
				}
			final SimpleInterval rec2 = new SimpleInterval(rec.getMateReferenceName(),mateStart,mateEnd);
			if(shouldBeSaved.test(rec2))
				{
				rec2list.add(rec2);
				}
			}
		/* add supplementary alignments that could be an evidence */
		SAMUtils.getOtherCanonicalAlignments(rec).
			stream().
			filter(shouldBeSaved).
			forEach(R-> rec2list.add(new SimpleInterval(R)));
						
		if(rec2list.isEmpty()) return;
		
		for(final SimpleInterval rec2: rec2list) {
			final String distanceStr;
			if(rec.overlaps(rec2)) {
				distanceStr = "0";
				}
			else  if(rec.contigsMatch(rec2)) {
				distanceStr = String.valueOf(Math.abs(Math.min(
						rec.getStart()-rec2.getEnd(),
						rec2.getStart()-rec.getEnd())));
				}
			else
				{
				distanceStr = ".";
				}
			
			// save insertions
			pw.print(rec.getContig());
			pw.print("\t");
			pw.print(rec.getStart()-1);
			pw.print("\t");
			pw.print(rec.getEnd());
			pw.print("\t");
			pw.print(rec2.getContig());
			pw.print("\t");
			pw.print(rec2.getStart()-1);
			pw.print("\t");
			pw.print(rec2.getEnd());
			pw.print("\t");
			pw.print(rec.getReadName()+(rec.getReadPairedFlag()?(rec.getFirstOfPairFlag()?"/1":"/2"):""));//name
			pw.print("\t");
			pw.print(".");//score
			pw.print("\t");
			pw.print(".");//strand 1
			pw.print("\t");
			pw.print(".");//strand 2
			pw.print("\t");

			// "Any number of additional, user-defined fields ..."
			pw.print(distanceStr);
			pw.println();
			}
		}

	
	@Override
	public int doWork(final List<String> args) {
		
		if(this.min_depth <1) {
			LOG.error("Bad min depth");
			return -1;
			}
		PrintWriter saveInsertionsPw = null;
		SamReader sr = null;
		VariantContextWriter vcw0=null;
		CloseableIterator<SAMRecord> iter = null;
		SAMFileWriter sfw = null;
		try {
			/* load the reference genome */
			/* create a contig name converter from the REF */
	
			
			// open the sam file
			final String input = oneFileOrNull(args);
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.faidx!=null) srf.referenceSequence(this.faidx);
			if(input!=null)
				{
				sr = srf.open(SamInputResource.of(input));
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(sr.getFileHeader());
				/* READ KNOWGENES FILES */
				loadGTF(dict);

				if(this.use_bai && !sr.hasIndex()) {
					LOG.warning("Cannot used bai because input is not indexed");
					}
				
				// if there is a bai, only query the bam in the regions of splicing
				if(this.use_bai && sr.hasIndex())
					{
					LOG.info("building intervals...");
					final QueryInterval intervals[] = this.intronTreeMap.keySet().
							stream().
							flatMap(intron->{
								// we need the reads overlapping the exon bounds
								final int tid=dict.getSequenceIndex(intron.getContig());
								final int extend = 1 +Math.max(0,this.minCigarSize);
								final QueryInterval q1=new QueryInterval(tid,Math.max(1, intron.getStart()-extend),intron.getStart()+extend);
								final QueryInterval q2=new QueryInterval(tid,Math.max(1, intron.getEnd()-extend),intron.getEnd()+extend);
								return Arrays.stream(new QueryInterval[]{q1,q2});
							}).
							sorted().
							collect(HtsCollectors.optimizedQueryIntervals());
					
					LOG.debug("Query bam using "+intervals.length+" random access intervals. Please wait...");
					iter = sr.queryOverlapping(intervals);
					}
				else
					{
					iter= sr.iterator();
					}
				}
			else
				{
				sr = srf.open(SamInputResource.of(stdin()));
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(sr.getFileHeader());
				/* READ GTF FILES */
				loadGTF(dict);
				iter= sr.iterator();
				}
			final SAMFileHeader samFileHeader = sr.getFileHeader();
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(samFileHeader);
			
				
			
			/* save gene writer */
			if(this.saveBedPeTo!=null) {
				saveInsertionsPw = super.openPathOrStdoutAsPrintWriter(this.saveBedPeTo);
				}
			else
				{
				saveInsertionsPw = NullOuputStream.newPrintWriter();
				}
			
			if(this.saveBamTo!=null) {
				sfw = new SAMFileWriterFactory().
						makeSAMOrBAMWriter(samFileHeader, true, this.saveBamTo);
				}
			
			final String sample = samFileHeader.getReadGroups().stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtils.isBlank(S)).
					findFirst().
					orElse("SAMPLE");
			
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.newInstance().dictionary(samFileHeader).logger(LOG).build();
			final String SAM_ATT_JI="jI";
			while(iter.hasNext()) {
				final SAMRecord rec = progress.apply(iter.next());
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.getMappingQuality() < this.min_read_mapq) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				if(rec.getReadFailsVendorQualityCheckFlag()) continue;
				
				boolean save_read_to_bam = false;

				/* save read if it is not properly mapped (problem with size) and he and his mate surround an intron */
				if(rec.getReadPairedFlag() && 
					!rec.getMateUnmappedFlag() &&
					rec.getReferenceIndex().equals(rec.getMateReferenceIndex()) &&
					!rec.getProperPairFlag() &&/* MUST NOT be proper pair */
					rec.getReadNegativeStrandFlag()!=rec.getMateNegativeStrandFlag() )
					{
					final SimpleInterval intronInterval;
					if(rec.getEnd()+1 < rec.getMateAlignmentStart()) {
						intronInterval = new SimpleInterval(
								rec.getContig(),
								rec.getEnd()+1,
								rec.getMateAlignmentStart()-1
								);
						}
					else if(SAMUtils.hasMateCigar(rec) && 
							SAMUtils.getMateAlignmentEnd(rec) +1 < rec.getAlignmentStart())
						{
						intronInterval = new SimpleInterval(
								rec.getContig(),
								SAMUtils.getMateAlignmentEnd(rec) +1,
								rec.getAlignmentStart()-1
								);
						}
					else
						{
						intronInterval = null;
						}
					if(intronInterval!=null) {
						if(this.intronTreeMap.getOverlapping(intronInterval).
							stream().
							flatMap(L->L.stream()).
							anyMatch(S->intronInterval.contains(S))) {
							save_read_to_bam = true;
							}
						}
					}
				
				
				/* WE use STAR DATA */
				if(!this.use_cigar_string)
					{
					if(!rec.hasAttribute(SAM_ATT_JI)) continue;
						
					final Object tagValue = rec.getAttribute(SAM_ATT_JI);
					paranoid.assertTrue((tagValue instanceof int[]));
					final int bounds[]= (int[])tagValue;
					// jI:B:i,-1
					if(bounds.length==1 && bounds[0]<0) continue;
					if(bounds.length%2!=0) {
						LOG.warn("bound.length%2!=0 with "+ rec.getSAMString());
						continue;
						}
	
					for(int i=0;i< bounds.length;i+=2) {
						int intron_start = bounds[i];
						int intron_end = bounds[i+1];
						
						final  Interval r = new Interval(rec.getContig(),intron_start,intron_end);
						final List<Segment> introns= this.intronTreeMap.get(r);// don't use overlapping : with STAR it is strictly the intron boundaries
						if(introns==null) continue;
						save_read_to_bam = true;
						for(final Segment intron:introns) {
							intron.match++;
							}
						}
					}
				else /* WE use other bam like bwa-mem */
					{
					final Cigar cigar = rec.getCigar();
					if(cigar==null || cigar.numCigarElements()<2) continue;
					
					int ref1=rec.getAlignmentStart();
					for(final CigarElement ce: cigar.getCigarElements())
						{
						final CigarOperator op = ce.getOperator();
						final  int ref_end = ref1 + (op.consumesReferenceBases()?ce.getLength():0);
						if(op.equals(CigarOperator.N) || op.equals(CigarOperator.D))
							{
							final  Interval r = new Interval(rec.getContig(),ref1,ref_end-1);
							final List<Segment> introns= this.intronTreeMap.get(r);
							if(introns==null) continue;
							save_read_to_bam = true;
							for(final Segment intron:introns) {
								intron.match++;
								}
							}
						ref1 = ref_end;
						}
					/** 2019-07-29. I tried using SA:Z:tag doesn't work well , so let's look a the clipping only */
					if(cigar.isClipped()) {
						for(int side=0;side < 2;side++) {
							final  Interval r;
							if(side==0 && 
								cigar.isRightClipped() &&
								cigar.getLastCigarElement().getLength() >= this.minCigarSize
								)
								{
								r = new Interval(
										rec.getContig(),
										rec.getEnd()+1,
										rec.getEnd()+1+this.minCigarSize
										);
								}
							else if(side==1 && 
									cigar.isLeftClipped() && 
									cigar.getFirstCigarElement().getLength() >= this.minCigarSize
									)
								{
								r = new Interval(rec.getContig(),
										Math.max(1,rec.getStart()-(1+this.minCigarSize)),
										Math.max(1,rec.getStart()-(1)))
										;
								}
							else
								{
								continue;
								}
							final int final_side = side;
							final List<Segment> introns= this.intronTreeMap.getOverlapping(r).
									stream().
									flatMap(L->L.stream()).
									filter(SEG->isWithinIntronBound(SEG,r,final_side)).
									collect(Collectors.toList());

							if(introns.isEmpty()) continue;
							//System.err.println("SA for "+r+" "+rec.getReadName()+" "+introns.size());

							save_read_to_bam = true;
							for(final Segment intron:introns) {
								intron.match++;
								}
							}
							
						}
					}
				
				if(save_read_to_bam) {
					saveInsertion(saveInsertionsPw,rec);
					if(sfw!=null) sfw.addAlignment(rec);
					}
				}
			final ContigDictComparator contigCmp = new ContigDictComparator(refDict);
			
			this.all_transcripts.removeIf(T->T.segments.stream().noneMatch(S->S.match>=min_depth));

			final int max_introns = this.all_transcripts.stream().mapToInt(K->K.segments.size()).max().orElse(1);
			final List<String> intron_names = IntStream.range(0, max_introns).
					mapToObj(IDX->String.format("%s_INTRON_%04d", sample,1+IDX)).
					collect(Collectors.toList())
					;
			/** build vcf header */
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY,true));
			metaData.add(new VCFInfoHeaderLine(VCFConstants.SVTYPE, 1, VCFHeaderLineType.String,"Variation type"));
			metaData.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"Variation Length"));
			
			for(final String att:ENSEMBL_TRANSCRIPT_ATTS)
				{
				metaData.add(new VCFInfoHeaderLine(att, 1, VCFHeaderLineType.String,"Value for the attribute '"+att+"' in the gtf"));
				}
			
			
			//metaData.add(new VCFFormatHeaderLine(ATT_COUNT_SUPPORTING_READS, 2,VCFHeaderLineType.Integer,"Count supporting reads [intron-left/intron-right]"));
			//metaData.add(new VCFInfoHeaderLine(ATT_RETRO_DESC, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
			//		"Retrocopy attributes: transcript-id|strand|exon-left|exon-left-bases|exon-right-bases|exon-right"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_COUNT, 1, VCFHeaderLineType.Integer,"Number of introns for the Transcript"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_CANDIDATE_COUNT, 1, VCFHeaderLineType.Integer,"Number of introns found retrocopied for the transcript"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_CANDIDATE_FRACTION, 1, VCFHeaderLineType.Float,"Fraction of introns found retrocopied for the transcript"));
			metaData.add(new VCFFormatHeaderLine(INTRON_START,1,VCFHeaderLineType.Integer, "Introns start"));
			metaData.add(new VCFFormatHeaderLine(INTRON_END,1,VCFHeaderLineType.Integer, "Introns end"));

			
			
			metaData.add(new VCFFilterHeaderLine(ATT_LOW_DEPTH_FILTER+this.low_depth_threshold,"Number of read is lower than :"+this.min_depth));
			metaData.add(new VCFFilterHeaderLine(ATT_NOT_ALL_INTRONS,"Not all introns were found retrocopied"));

			
			
			
			
			final VCFHeader header=new VCFHeader(metaData, intron_names); 
				
			JVarkitVersion.getInstance().addMetaData(this, header);
			header.setSequenceDictionary(refDict);
			
			/* open vcf for writing*/
			vcw0=super.openVariantContextWriter(this.outputFile);
			final VariantContextWriter vcw=vcw0;
			vcw.writeHeader(header);
			
			
			Collections.sort(this.all_transcripts,(A,B)->{
				int i= contigCmp.compare(A.getContig(), B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				return  Integer.compare(A.getEnd(), B.getEnd());
			});
			
			final Allele ref= Allele.create((byte)'N', true);
			final Allele alt= Allele.create("<RETROCOPY>", false);
			
			for(final Transcript kg: this.all_transcripts)
				{
			
				
				// ok good candidate
				final VariantContextBuilder vcb = new VariantContextBuilder();
				vcb.chr(kg.getContig());
				vcb.start(kg.getStart());
				vcb.stop(kg.getEnd());
				vcb.id(kg.transcript_id);
				final List<Allele> alleles = Arrays.asList(ref,alt);

				final int max_depth = kg.segments.stream().mapToInt(X->X.match).max().orElse(0);
				vcb.attribute(VCFConstants.DEPTH_KEY,max_depth);
				vcb.log10PError(max_depth/-10.0);

				boolean filter_set=false;
				if(max_depth < this.low_depth_threshold)
					{
					vcb.filter(ATT_LOW_DEPTH_FILTER+this.low_depth_threshold);
					filter_set = true;
					}
				
				vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,2);
				vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,1);
				vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,0.5);
				vcb.attribute(VCFConstants.SVTYPE,"DEL");
				vcb.attribute(VCFConstants.END_KEY,kg.getEnd());
				vcb.attribute("SVLEN",kg.getLengthOnReference());
				for(final String att: kg.attributes.keySet()) {
					vcb.attribute(att,VCFUtils.escapeInfoField(kg.attributes.get(att)));

				}
				
				vcb.alleles(alleles);
				
				// introns sequences
				vcb.attribute(ATT_INTRONS_CANDIDATE_COUNT,kg.segments.stream().filter(I->I.match>0).count());
				vcb.attribute(ATT_INTRONS_COUNT,kg.segments.size());
				vcb.attribute(ATT_INTRONS_CANDIDATE_FRACTION,kg.segments.stream().filter(I->I.match>0).count()/(float)kg.segments.size());
				if(kg.segments.stream().filter(I->I.match>0).count()!=kg.segments.size()) {
					vcb.filter(ATT_NOT_ALL_INTRONS);
					filter_set=true;
					}
				
				final List<Genotype> genotypes = new ArrayList<>(kg.segments.size());
				/* build genotypes */
				for(int i=0;i< kg.segments.size();i++)
					{
					final Segment intron = kg.segments.get(i);
					final GenotypeBuilder gb= new GenotypeBuilder(intron_names.get(i), Arrays.asList(ref,alt));
					gb.DP(intron.match);
					gb.attribute(INTRON_START, intron.start);
					gb.attribute(INTRON_END, intron.end);
					genotypes.add(gb.make());
					}
				
				
				vcb.genotypes(genotypes);
				
				
				if(!filter_set) {
					vcb.passFilters();
					}
				
				vcw.add(vcb.make());
				}
			progress.close();
			vcw.close();
			iter.close();
			iter=null;
			sr.close();
			sr=null;
			
			saveInsertionsPw.flush();
			saveInsertionsPw.close();
			saveInsertionsPw=null; 
			
			if(sfw!=null) {
				sfw.close();
				sfw=null;
				}
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sr);
			CloserUtil.close(vcw0);
			CloserUtil.close(sfw);
			CloserUtil.close(saveInsertionsPw);
			}
		}
	
	public static void main(final String[] args) {
		new StarRetroCopy().instanceMainWithExit(args);
	}
}
