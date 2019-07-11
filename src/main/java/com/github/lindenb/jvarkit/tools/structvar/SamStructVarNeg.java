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

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.OptionalInt;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.iterator.FilterIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**

BEGIN_DOC


```
$ java -jar dist/svneg.jar --bams jeter.list src/test/resources/HG02260.transloc.chr9.14.bam

##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SR,Number=4,Type=Integer,Description="Supporting reads: contig1-forward,contig1-reverse,contig2-forward,contig2-reverse">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=CHROM2,Number=1,Type=String,Description="other chromosome">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="other position">
##INFO=<ID=STDDEV_POS1,Number=1,Type=Integer,Description="std deviation to position 1">
##INFO=<ID=STDDEV_POS2,Number=1,Type=Integer,Description="std deviation to position 2">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variation type">
##svneg.meta=compilation:20190415102101 githash:4c6799f7 htsjdk:2.19.0 date:20190415102110 cmd:--bams jeter.list src/test/resources/HG02260.transloc.chr9.14.bam
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG02260
9	137230996	9:137230996:14:79839048	N	<TRANSLOC>	18	.	AC=1;AF=0.500;AN=2;CHROM2=14;DP=18;POS2=79839048;STDDEV_POS1=120;STDDEV_POS2=187;SVTYPE=BND	GT:DP:SR	0/1:18:5,13,13,5
14	79839119	14:79839119:9:137230968	N	<TRANSLOC>	18	.	AC=1;AF=0.500;AN=2;CHROM2=9;DP=18;POS2=137230968;STDDEV_POS1=154;STDDEV_POS2=145;SVTYPE=BND	GT:DP:SR	0/1:18:13,5,5,13


```

END_DOC
*/
@Program(name="svneg",
	description="Find Structural Variation by Negative Comparaison",
	keywords={"sam","bam","sv","translocation"},
	creationDate="20190413",
	modificationDate="20190701"
	)
public class SamStructVarNeg extends Launcher {
	private static final Logger LOG = Logger.build(SamStructVarNeg.class).make();
	
	
	
	private class StructuralVariant
		implements Locatable
		{
		final String contig;
		final int start;
		String contig2;
		StructuralVariantType svType;
		int start2;
		int count = 0;
		StructuralVariant(final SAMRecord vc) {
			this.contig = vc.getContig();
			this.start = vc.getStart() - vc.getStart()%bin_size;
			//
			this.contig2 = this.contig;
			this.start2 = this.start;
			}
		@Override
		public String getContig() {
			return this.contig;
			}
		@Override
		public int getStart() {
			return this.start;
			}
		@Override
		public int getEnd() {
			return this.start + bin_size;
			}
		}
	private final List<StructuralVariant> buffer = new ArrayList<>();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	
	
	@Parameter(names={"-bin","--bin"},
			description="bin size. Break the genome in parts of 'bin' size." +DistanceParser.OPT_DESCRIPTION,
			converter=DistanceParser.StringConverter.class,
			splitter=NoSplitter.class
			)
	private int bin_size = 1_000;
	@Parameter(names={"-m","--min"},description="Min number of events to validate the translocation")
	private int min_number_of_events=3;
	@Parameter(names={"--max-controls"},description="Maximum number controls bam matching the event")
	private int max_number_of_controls = 1;
	@Parameter(names={"-b","--controls"},description="Control Bams. One path per lines",required=true)
	private Path controlBamPaths = null;
	@Parameter(names={"-c","--cases"},description="Control Bams. One path per lines",required=true)
	private Path casesBamPaths = null;
	
	@Parameter(names={"--mapq"},description="min mapping quality.")
	private int min_mapq = 0;
	@Parameter(names={"-x","--exclude"},description="Optional BED file Excluding. SV shouldn't overlap this bed.")
	private Path excludeBedFile = null;
	
	@Parameter(names={"-D"},description="Presence of a discordant read in the control, whatever is the contig or the distance to theoritical mate is enought to invalidate the candiate. Make things quicker but less sensitive.")
	private boolean presence_of_discordant_in_ctrl_is_enough = false;
	
	
	private final IntervalTreeMap<Interval> excludeBedMap = new IntervalTreeMap<Interval>();
	
	private class BamResource implements Closeable{
		private final int index;
		private final Path bamPath;
		private SamReader samReader = null;
		private ContigNameConverter ctgNameConverter;
		private final Set<String> unseenContig = new HashSet<>();
		
		BamResource(final int index,final Path bamPath) {
			this.index  = index;
			this.bamPath = bamPath;
			
			}
		
		private void open() throws IOException {
			if( this.samReader != null ) return;
			LOG.info("opening ["+(index+1)+"/"+ controlBams.size()+"] "+bamPath);
			final SamReaderFactory srf = SamReaderFactory.
					makeDefault().
					validationStringency(ValidationStringency.SILENT);
			this.samReader = srf.open(this.bamPath);
			if(!samReader.hasIndex()) {
				this.samReader.close();
				throw new IOException("BAM is not indexed " + this.bamPath);
				}
			final SAMFileHeader header = samReader.getFileHeader();
			this.ctgNameConverter = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(header));
			if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
				this.samReader.close();
				throw new JvarkitException.BamBadSortOrder(SAMFileHeader.SortOrder.coordinate, header.getSortOrder());
				}		
			}
		
		@Override
		public void close() throws IOException {
			if(this.samReader!=null) CloserUtil.close(this.samReader);
			this.samReader = null;
			}
		}
	
	
	/* for the data I tested there was a bug in the sorting (reads not sorted on FLAG see htsjdk.samtools.SAMRecordCoordinateComparator)
	 but we just need fileOrderCompare */
	
	
	
	private final List<BamResource> controlBams = new ArrayList<>(); 	
	private final List<BamResource> casesBams = new ArrayList<>(); 	
		
	
	private class StructVarIterator 
		extends AbstractIterator<List<StructuralVariant>>
		implements CloseableIterator<List<StructuralVariant>>
		{
		double sum_insert_size = 0.0;
		long count_insert_size = 0L;
		
		final CloseableIterator<SAMRecord> iter;
		final List<StructuralVariant> buffer = new ArrayList<>();
		StructVarIterator(final CloseableIterator<SAMRecord> iter) {
			this.iter = iter;
			}
		
		private int convertBin(int pos) {
			return pos - pos%bin_size;
		}
		
		@Override
		protected List<StructuralVariant> advance() {
			final short SA = SAMTag.SA.getBinaryTag();
			for(;;) {
				
				final SAMRecord rec;
				rec= this.iter.hasNext()?this.iter.next():null;
				if(rec!=null && (rec.getReadUnmappedFlag() || rec.getMappingQuality() < min_mapq || excludeBedMap.containsOverlapping(rec))) continue;
				if(rec==null ||
					(!buffer.isEmpty() && !buffer.get(0).getContig().equals(rec.getContig())) ||
					(!buffer.isEmpty() && buffer.get(0).start!=convertBin(rec.getStart())))
					{
					if(!this.buffer.isEmpty())
						{
						final List<StructuralVariant> copy = new ArrayList<>(this.buffer);
						this.buffer.clear();
						return copy;
						}
					if(rec==null) {
						close();
						return null;
						}
					
					}
				if(rec.getReadPairedFlag())
					{
					// mate mapped
					if(!rec.getMateUnmappedFlag()) {
						if(rec.getReferenceIndex().equals(rec.getMateReferenceIndex())) {
							
							if(rec.getProperPairFlag() && 
								rec.getFirstOfPairFlag())
								{
								this.sum_insert_size += Math.abs(rec.getInferredInsertSize());
								this.count_insert_size++;
								}
							if(this.count_insert_size > 1000L && 
									!rec.getProperPairFlag() && 
									!rec.getReadNegativeStrandFlag() &&
									rec.getMateNegativeStrandFlag() &&
									rec.getInferredInsertSize() > 2.0*(this.sum_insert_size/this.count_insert_size))
								{
								
								}
							if(!rec.getReadNegativeStrandFlag() &&
									rec.getCigar()!=null &&
									rec.getCigar().isRightClipped())
								{
								
								}
							}
						else /* discordant reads not same chromosomes */
							{
							
							}
						}
					}
				}
			}
		@Override
		public void close() {
			this.iter.close();
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.bin_size <= 0) {
			LOG.error("fuzzy_distance <=0 (" + this.bin_size + ")");
			return -1;
		}
		if(this.controlBamPaths==null ) {
			LOG.error("Control BAM is undefined");
			return -1;
			}
		if(this.max_number_of_controls <1 ) {
			LOG.error("Bad number for max_number_of_controls");
			return -1;
			}
		
		CloseableIterator<SAMRecord> iter0 = null;
		StructVarIterator iter =null;
		VariantContextWriter out = null;
		
		if(!args.isEmpty()) {
			LOG.error("illegal number of arguments");
			return -1;
			}
		
		try {
			
			
			
			try(BufferedReader br = IOUtils.openPathForBufferedReading(this.controlBamPaths)) {
				br.lines().
					filter(L->!L.startsWith("#")).
					filter(L->!StringUtils.isBlank(L)).
					map(L->Paths.get(L)).
					forEach(P->{
						this.controlBams.add(new BamResource(this.controlBams.size(),P));
					});
				}
			
			try(BufferedReader br = IOUtils.openPathForBufferedReading(this.casesBamPaths)) {
				br.lines().
					filter(L->!L.startsWith("#")).
					filter(L->!StringUtils.isBlank(L)).
					map(L->Paths.get(L)).
					forEach(P->{
						this.casesBams.add(new BamResource(this.casesBams.size(),P));
					});
				}
			
			
			if(this.controlBams.isEmpty()) {
				LOG.error("No control bam was defined");
				return -1;
				}
			if(this.casesBams.isEmpty()) {
				LOG.error("No case bam was defined");
				return -1;
				}
			
			if(this.controlBams.size() < this.max_number_of_controls) {
				LOG.error("Number of bam  is lower than ");
				return -1;
				}
			
			final SamReader caseSamReader = this.casesBams.get(0).samReader;
			
			iter = new StructVarIterator(caseSamReader.iterator());
			
				
			
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(caseSamReader.getFileHeader());
			
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.
					newInstance().
					dictionary(refDict).
					logger(LOG).
					build();
			
			if(this.excludeBedFile!=null) {
				final BedLineCodec bedCodec = new BedLineCodec();
				try(BufferedReader br=com.github.lindenb.jvarkit.io.IOUtils.openPathForBufferedReading(this.excludeBedFile)) {
					br.lines().
					filter(line->!(line.startsWith("#") ||  com.github.lindenb.jvarkit.util.bio.bed.BedLine.isBedHeader(line) ||  line.isEmpty())).
					map(line->bedCodec.decode(line)).
					filter(B->B!=null).
					map(B->B.toInterval()).
					filter(L->L.getStart()<L.getEnd()).
					forEach(B->{
						this.excludeBedMap.put(B,B);							
						});	
					}
				}
		
			
			
			final Set<VCFHeaderLine> metaData=new HashSet<>();
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
			final VCFFormatHeaderLine supportingReadsFormat = new VCFFormatHeaderLine("SR",
					4,
					VCFHeaderLineType.Integer,
					"Supporting reads: contig1-forward,contig1-reverse,contig2-forward,contig2-reverse"
					);
			metaData.add(supportingReadsFormat);
			final VCFInfoHeaderLine stdDevContig1Info = new VCFInfoHeaderLine("STDDEV_POS1",
					1,
					VCFHeaderLineType.Integer,
					"std deviation to position 1"
					);
			metaData.add(stdDevContig1Info);
			final VCFInfoHeaderLine stdDevContig2Info = new VCFInfoHeaderLine("STDDEV_POS2",
					1,
					VCFHeaderLineType.Integer,
					"std deviation to position 2"
					);
			metaData.add(stdDevContig2Info);
			final VCFInfoHeaderLine chrom2Info = new VCFInfoHeaderLine("CHROM2",
					1,
					VCFHeaderLineType.String,
					"other chromosome"
					);
			metaData.add(chrom2Info);
			final VCFInfoHeaderLine pos2Info = new VCFInfoHeaderLine("POS2",
					1,
					VCFHeaderLineType.Integer,
					"other position"
					);
			metaData.add(pos2Info);
			

			final VCFHeader vcfHeader= new VCFHeader(metaData, Collections.singletonList(sampleName));
			vcfHeader.setSequenceDictionary(refDict);
			JVarkitVersion.getInstance().addMetaData(this, vcfHeader);
			
			out = VCFUtils.createVariantContextWriterToPath(this.outputFile);
			out.writeHeader(vcfHeader);
			
			
			final Allele REF=Allele.create("N", true);
			final Allele ALT=Allele.create("<TRANSLOC>", false);

			
			int prev_tid=-1;
			while(iter.hasNext())
				{
				final List<StructuralVariant> candidates = iter.next();
			
				int count_cases = 0;
				
				for(int i=1 /* start from 1 , 0 is the current case */;i< this.casesBams.size();++i)
					{
					for(final StructuralVariant sv:candidates) {
						final SAMRecordIterator sri2 = this.casesBams.get(i).samReader.query(sv.getContig(), sv.getStart(), sv.getEnd(), false);
						StructVarIterator iter2= new StructVarIterator(sri2);
						
						iter2.close();
						sri2.close();
						}
					}
				
				
			
				/* fill with current buffer until end */
				int buffer_index = 0;
				while(buffer_index < buffer.size())
					{
					final SAMRecord rec2 = buffer.get(buffer_index);
					
					if(!rec2.getReferenceIndex().equals(rec.getReferenceIndex())) {
						break;
						}
					if(rec2.getAlignmentStart() > end) {//not unclipped to avoid side effect
						break;
						}
					/* same mate contig for rec and rec2 */
					if(!rec2.getMateReferenceIndex().equals(rec.getMateReferenceIndex())) {
						buffer_index++;
						continue;
						}
					
					final int mate1 = rec.getMateNegativeStrandFlag()?
								rec.getMateAlignmentStart():
								getMateAlignmenEnd(rec)
								;
					final int mate2 = rec2.getMateNegativeStrandFlag()?
							rec2.getMateAlignmentStart():
							getMateAlignmenEnd(rec2)
							;

					if(Math.abs(mate1-mate2) > this.fuzzy_distance) {
						buffer_index++;
						continue;
						}
					buffer.remove(buffer_index);
					candidates.add(rec2);						
					}
				
					
				
				
				
				if(candidates.isEmpty() || candidates.size()<this.min_number_of_events) {
					continue;
					}
				

				final int pos_ctg1=(int)candidates.stream().mapToInt(
						SR->SR.getReadNegativeStrandFlag()?
								SR.getAlignmentStart():
								SR.getAlignmentEnd()).
						average().
						orElse(-1);
				
				if(pos_ctg1<1) {
					//putbackInBuffer(candidates, buffer);
					continue;
					}
				final int stddev_ctg1 = (int)candidates.stream().mapToInt(
						SR->SR.getReadNegativeStrandFlag()?SR.getAlignmentStart():SR.getAlignmentEnd()).
						map(X->Math.abs(X-pos_ctg1)).
						average().
						orElse(0.0);
				
				final int pos_ctg2=(int)candidates.stream().
					mapToInt(SR->SR.getMateNegativeStrandFlag()?
							SR.getMateAlignmentStart():
							getMateAlignmenEnd(SR)).
					average().orElse(-1);
				
				if(pos_ctg2<1) {
					//putbackInBuffer(candidates, buffer);
					continue;
					}
				
				final int stddev_ctg2 = (int)candidates.stream().
					mapToInt(SR->SR.getMateNegativeStrandFlag()?SR.getMateAlignmentStart():getMateAlignmenEnd(SR)).
					map(X->Math.abs(X-pos_ctg2)).
					average().
					orElse(0.0);
				
				
				int next_end = 0;
				int control_count = 0;
				for(final ControlBam controlBam : this.controlBams) {
					final OptionalInt optEnd = controlBam.overlap(
						new Interval(
							rec.getContig(),
							pos_ctg1 - stddev_ctg1 - this.fuzzy_distance ,
							pos_ctg1 + stddev_ctg1 + this.fuzzy_distance
							),
						new Interval(
								rec.getMateReferenceName(),
								Math.max(0, pos_ctg2 - stddev_ctg2 - this.fuzzy_distance) ,
								pos_ctg2 + stddev_ctg2 + this.fuzzy_distance
								)
						);
					/* control contains SV */
					if(optEnd.isPresent()) {
						next_end = Math.max(next_end, optEnd.getAsInt());
						control_count ++;
						if(control_count >= this.max_number_of_controls) break;
						}
					}
				if(control_count == 0 ) {
					LOG.info("all bams for "+rec.getContig()+":"+rec.getStart()+"-"+rec.getEnd());
					}
				
				if(control_count >= this.max_number_of_controls) {
					ignore_to_position = next_end;
					putbackInBuffer(candidates, buffer);
					continue;
					}
				
				//check in both sides
				final int count_plus =  (int)candidates.stream().filter(SR->!SR.getReadNegativeStrandFlag()).count();
				final int count_minus = (int)candidates.stream().filter(SR-> SR.getReadNegativeStrandFlag()).count();
				if(count_plus<this.min_number_of_events || count_minus<this.min_number_of_events) {
					putbackInBuffer(candidates, buffer);
					continue;
					}
								
				
				
				final VariantContextBuilder vcb=new VariantContextBuilder(null,
						rec.getContig(), 
						pos_ctg1,
						pos_ctg1,
						Arrays.asList(REF,ALT)
						);
					
				final GenotypeBuilder gb=new GenotypeBuilder(sampleName,Arrays.asList(REF,ALT));
				final int sn_contig1_count_plus =  (int)candidates.stream().filter(SR->!SR.getReadNegativeStrandFlag()).count();
				final int sn_contig1_count_minus = (int)candidates.stream().filter(SR-> SR.getReadNegativeStrandFlag()).count();
				final int sn_contig2_count_plus =  (int)candidates.stream().filter(SR->!SR.getMateNegativeStrandFlag()).count();
				final int sn_contig2_count_minus = (int)candidates.stream().filter(SR-> SR.getMateNegativeStrandFlag()).count();
				
						
				gb.DP(candidates.size());
				gb.attribute(supportingReadsFormat.getID(),
					new int[] {
						sn_contig1_count_plus,
						sn_contig1_count_minus,
						sn_contig2_count_plus,
						sn_contig2_count_minus
						}
					);
				
				vcb.genotypes(Collections.singletonList( gb.make()));
				
				vcb.id(rec.getReferenceName()+":"+pos_ctg1+":"+rec.getMateReferenceName()+":"+pos_ctg2);
				vcb.attribute(stdDevContig1Info.getID(), stddev_ctg1);
				vcb.attribute(stdDevContig2Info.getID(), stddev_ctg2);
				vcb.attribute(chrom2Info.getID(), rec.getMateReferenceName());
				vcb.attribute(pos2Info.getID(),pos_ctg2);
				
				
				
				vcb.attribute(VCFConstants.SVTYPE, StructuralVariantType.BND.name());
				vcb.attribute(VCFConstants.DEPTH_KEY,candidates.size());
				vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,1);
				vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,2);
				vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,0.5);
				
				vcb.log10PError(candidates.size()/-10.0);
				vcb.alleles(Arrays.asList(REF,ALT));
				
				out.add(vcb.make());
				
				buffer.removeIf(REC->REC.getReferenceIndex().equals(rec.getReferenceIndex()) && REC.getAlignmentEnd() <= rec.getEnd());
				}
			progress.close();
			iter.close();
			
			CloserUtil.close(iter); iter=null;
			CloserUtil.close(samReader);iter=null;
			out.close();
			out=null;
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			for(final ControlBam b:this.controlBams) CloserUtil.close(b);
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			CloserUtil.close(out);
			}
		}
	
	public static void main(final String[] args) {
		new SamStructVarNeg().instanceMainWithExit(args);
		}

}
