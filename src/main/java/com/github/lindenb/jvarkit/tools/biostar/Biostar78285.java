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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.iterator.FilterIterator;
import com.github.lindenb.jvarkit.util.iterator.MergingIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
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

/**

BEGIN_DOC

## Example

```bash
$ java -jar dist/biostar78285.jar -m 300   -R  ref.fa S*.bam 

##fileformat=VCFv4.2
##Biostar78285.SamFilter=record.getMappingQuality()<1 || record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag() || record.isSecondaryOrSupplementary()
##FILTER=<ID=DP_LT_300,Description="All  genotypes have DP< 300">
##FORMAT=<ID=DF,Number=1,Type=Integer,Description="Number of Reads on plus strand">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of Reads on minus strand">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AVG_DP,Number=1,Type=Float,Description="Mean depth">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=FRACT_DP_LT_300,Number=1,Type=Float,Description="Fraction of genotypes having DP< 300">
##INFO=<ID=GC_PERCENT,Number=1,Type=Integer,Description="GC% window_size:20">
##INFO=<ID=MAX_DP,Number=1,Type=Integer,Description="Max depth">
##INFO=<ID=MEDIAN_DP,Number=1,Type=Float,Description="Median depth">
##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description="Min depth">
##INFO=<ID=NUM_DP_LT_300,Number=1,Type=Integer,Description="Number of genotypes having DP< 300">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	1	.	G	.	.	DP_LT_300	AVG_DP=4.25;DP=17;FRACT_DP_LT_300=1.0;GC_PERCENT=38;MAX_DP=5;MEDIAN_DP=4.50;MIN_DP=3;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:5:5:0	./.:5:5:0	./.:3:3:0	./.:4:4:0
rotavirus	2	.	G	.	.	DP_LT_300	AVG_DP=9.50;DP=38;FRACT_DP_LT_300=1.0;GC_PERCENT=40;MAX_DP=14;MEDIAN_DP=8.50;MIN_DP=7;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:14:14:0	./.:9:9:0	./.:8:8:0	./.:7:7:0
rotavirus	3	.	C	.	.	DP_LT_300	AVG_DP=12.25;DP=49;FRACT_DP_LT_300=1.0;GC_PERCENT=39;MAX_DP=18;MEDIAN_DP=11.50;MIN_DP=8;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:18:18:0	./.:11:11:0	./.:12:12:0	./.:8:8:0
rotavirus	4	.	T	.	.	DP_LT_300	AVG_DP=16.25;DP=65;FRACT_DP_LT_300=1.0;GC_PERCENT=37;MAX_DP=22;MEDIAN_DP=17.00;MIN_DP=9;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:22:22:0	./.:16:16:0	./.:18:18:0	./.:9:9:0
rotavirus	5	.	T	.	.	DP_LT_300	AVG_DP=20.75;DP=83;FRACT_DP_LT_300=1.0;GC_PERCENT=40;MAX_DP=27;MEDIAN_DP=21.50;MIN_DP=13;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:27:27:0	./.:18:18:0	./.:25:25:0	./.:13:13:0
rotavirus	6	.	T	.	.	DP_LT_300	AVG_DP=24.25;DP=97;FRACT_DP_LT_300=1.0;GC_PERCENT=42;MAX_DP=33;MEDIAN_DP=25.50;MIN_DP=13;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:30:30:0	./.:21:21:0	./.:33:33:0	./.:13:13:0
rotavirus	7	.	T	.	.	DP_LT_300	AVG_DP=28.00;DP=112;FRACT_DP_LT_300=1.0;GC_PERCENT=40;MAX_DP=38;MEDIAN_DP=30.00;MIN_DP=14;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:38:38:0	./.:23:23:0	./.:37:37:0	./.:14:14:0
rotavirus	8	.	A	.	.	DP_LT_300	AVG_DP=30.50;DP=122;FRACT_DP_LT_300=1.0;GC_PERCENT=42;MAX_DP=41;MEDIAN_DP=33.00;MIN_DP=15;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:41:41:0	./.:26:26:0	./.:40:40:0	./.:15:15:0
rotavirus	9	.	A	.	.	DP_LT_300	AVG_DP=34.75;DP=139;FRACT_DP_LT_300=1.0;GC_PERCENT=44;MAX_DP=48;MEDIAN_DP=37.50;MIN_DP=16;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:46:46:0	./.:29:29:0	./.:48:48:0	./.:16:16:0
rotavirus	10	.	T	.	.	DP_LT_300	AVG_DP=40.75;DP=163;FRACT_DP_LT_300=1.0;GC_PERCENT=43;MAX_DP=56;MEDIAN_DP=43.00;MIN_DP=21;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:56:56:0	./.:35:35:0	./.:51:51:0	./.:21:21:0
rotavirus	11	.	G	.	.	DP_LT_300	AVG_DP=44.75;DP=179;FRACT_DP_LT_300=1.0;GC_PERCENT=45;MAX_DP=58;MEDIAN_DP=49.50;MIN_DP=22;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:58:58:0	./.:42:42:0	./.:57:57:0	./.:22:22:0
rotavirus	12	.	C	.	.	DP_LT_300	AVG_DP=48.75;DP=195;FRACT_DP_LT_300=1.0;GC_PERCENT=43;MAX_DP=66;MEDIAN_DP=53.00;MIN_DP=23;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:66:66:0	./.:46:46:0	./.:60:60:0	./.:23:23:0
rotavirus	13	.	T	.	.	DP_LT_300	AVG_DP=53.50;DP=214;FRACT_DP_LT_300=1.0;GC_PERCENT=42;MAX_DP=73;MEDIAN_DP=58.50;MIN_DP=24;NUM_DP_LT_300=4	GT:DF:DP:DR	./.:73:73:0	./.:51:51:0	./.:66:66:0	./.:24:24:0
```

## History

* 20180227: moved output to VCF, printing everything, adding optional BED

END_DOC

*/
@Program(name="biostar78285",
	biostars=78285,
	keywords={"sam","bam","depth","coverage"},
	description="Extract BAMs coverage as a VCF file.")
public class Biostar78285 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar78285.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	@Parameter(names={"-f","--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-B","--bed","--capture"},description="Limit analysis to this bed file")
	private File captureBed = null;
	@Parameter(names={"--partition"},description="When using display READ_GROUPS, how should we partition the ReadGroup ? "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition= SAMRecordPartition.sample;
	@Parameter(names={"-m","--min-depth"},description="Min depth tresholds.")
	private List<Integer> minDepthTresholds = new ArrayList<>();
	@Parameter(names={"-R","--reference"},description="Optional. "+ INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File refFile = null;
	@Parameter(names={"-gcw","--gc-percent-window","--gcw"},description="GC% window size. (if REF is defined)")
	private int gc_percent_window=20;

	private static class PosInfo
		{
		final String sample;
		int dp=0;
		int negative_strand =0;
		PosInfo(final String sample) {
			this.sample = sample;
			}
		}
	
	
    @Override
	public int doWork(final List<String> args) {
    	if(this.gc_percent_window<1) {
    		LOG.error("Bad GC% window size:"+this.gc_percent_window);
    		return -1;
    	}
    	
		final List<File> bamFiles = IOUtil.unrollFiles(args.stream().map(F->new File(F)).collect(Collectors.toCollection(HashSet::new)), ".bam");
		SAMSequenceDictionary dict=null;
		final List<SamReader> samReaders  = new ArrayList<>(); 
		final List<CloseableIterator<SAMRecord>> samIterators  = new ArrayList<>(); 
		final TreeSet<String> samples = new TreeSet<>(); 
		final String DEFAULT_PARTITION = "UNDEFINED_PARTITION";
		IndexedFastaSequenceFile indexedFastaSequenceFile = null;
    	VariantContextWriter out=null;
		try
			{
			
			final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			
			for(final File bamFile: bamFiles) {
				LOG.info("Opening "+bamFile);
				final SamReader samReader = samReaderFactory.open(bamFile);
				samReaders.add(samReader);
				final SAMFileHeader header=samReader.getFileHeader();
				if(header == null)
		    		{
		    		LOG.error("No header in "+bamFile);
		    		return -1;
		    		}
				JvarkitException.BamBadSortOrder.verify(SortOrder.coordinate, header);
		    	samples.addAll(header.getReadGroups().stream().map(RG->this.partition.apply(RG, DEFAULT_PARTITION)).collect(Collectors.toSet()));
		    	
		    	final SAMSequenceDictionary currDict =header.getSequenceDictionary();
		    	if(currDict==null)
		    		{
		    		LOG.error("SamFile doesn't contain a SAMSequenceDictionary : "+bamFile);
		    		return -1;
		    		}
		    	if(dict==null)
		    		{
		    		dict = currDict;
		    		}
		    	else if(!SequenceUtil.areSequenceDictionariesEqual(dict, currDict))
		    		{
		    		LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(dict,currDict));
		    		return -1;
		    		}
		    	
				}
			if(samReaders.isEmpty()) {
				LOG.error("no bam");
				return -1;
				}	
			if(dict==null) {
				LOG.error("no dictionary");
				return -1;
				}
			
			
			
			final QueryInterval intervals[];
			if(this.captureBed!=null)
				{
				LOG.info("Opening "+this.captureBed);
				ContigNameConverter.setDefaultAliases(dict);
				final List<QueryInterval> L = new ArrayList<>();
				final BedLineCodec codec= new BedLineCodec();
				final LineIterator li = IOUtils.openFileForLineIterator(this.captureBed);
				while(li.hasNext()) {
					final BedLine bed = codec.decode(li.next());
					if(bed==null) continue;
					final QueryInterval q= bed.toQueryInterval(dict);
					L.add(q);
					}
				CloserUtil.close(li);
				intervals = QueryInterval.optimizeIntervals(L.toArray(new QueryInterval[L.size()]));
				}
			else
				{
				intervals = null;
				}
			for(final SamReader samReader : samReaders)
				{
				LOG.info("querying "+samReader.getResourceDescription());
				final CloseableIterator<SAMRecord> iter;
				if(intervals==null)
					{
					iter = samReader.iterator();
					}
				else
					{
					iter = samReader.queryOverlapping(intervals);
					}
				
				samIterators.add(new FilterIterator<SAMRecord>(iter,R->!R.getReadUnmappedFlag() && !filter.filterOut(R)));
				}
			
			if(this.refFile!=null) {
				LOG.info("opening "+refFile);
				indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.refFile);
				final SAMSequenceDictionary refdict = indexedFastaSequenceFile.getSequenceDictionary();
				ContigNameConverter.setDefaultAliases(refdict);
				if(refdict==null) {
					throw new JvarkitException.FastaDictionaryMissing(this.refFile);
					}
				 if(!SequenceUtil.areSequenceDictionariesEqual(dict, refdict))
		    		{
		    		LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(dict,refdict));
		    		return -1;
		    		}
				}
			
			out = openVariantContextWriter(this.outputFile);
			
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true, 
					VCFConstants.DEPTH_KEY,
					VCFConstants.GENOTYPE_KEY
					);
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true, 
					VCFConstants.DEPTH_KEY
					);
			
			metaData.add(new VCFFormatHeaderLine("DF",1,VCFHeaderLineType.Integer,"Number of Reads on plus strand"));
			metaData.add(new VCFFormatHeaderLine("DR",1,VCFHeaderLineType.Integer,"Number of Reads on minus strand"));
			
			metaData.add(new VCFInfoHeaderLine("AVG_DP",1,VCFHeaderLineType.Float, "Mean depth"));
			metaData.add(new VCFInfoHeaderLine("MEDIAN_DP",1,VCFHeaderLineType.Float, "Median depth"));
			metaData.add(new VCFInfoHeaderLine("MIN_DP",1,VCFHeaderLineType.Integer, "Min depth"));
			metaData.add(new VCFInfoHeaderLine("MAX_DP",1,VCFHeaderLineType.Integer, "Max depth"));
			metaData.add(new VCFHeaderLine(Biostar78285.class.getSimpleName()+".SamFilter",this.filter.toString()));
			for(final Integer treshold: this.minDepthTresholds)
				{
				metaData.add(new VCFFilterHeaderLine("DP_LT_"+treshold, "All  genotypes have DP< "+treshold));
				metaData.add(new VCFInfoHeaderLine("NUM_DP_LT_"+treshold,1,VCFHeaderLineType.Integer, "Number of genotypes having DP< "+treshold));
				metaData.add(new VCFInfoHeaderLine("FRACT_DP_LT_"+treshold,1,VCFHeaderLineType.Float, "Fraction of genotypes having DP< "+treshold));
				}
			
			if(indexedFastaSequenceFile!=null)
				{
				metaData.add(new VCFInfoHeaderLine("GC_PERCENT",1,VCFHeaderLineType.Integer, "GC% window_size:"+this.gc_percent_window));
				}
			
			final List<Allele> refAlleles = Collections.singletonList(Allele.create("N", true));
			final List<Allele> NO_CALLS = Arrays.asList(Allele.NO_CALL,Allele.NO_CALL);
			final VCFHeader vcfHeader = new VCFHeader(metaData,samples);
			vcfHeader.setSequenceDictionary(dict);
			out.writeHeader(vcfHeader);
			
			
			final SAMRecordCoordinateComparator samRecordCoordinateComparator = new SAMRecordCoordinateComparator();
			final PeekableIterator<SAMRecord> peekIter = new PeekableIterator<>(
					new MergingIterator<>(
						(R1,R2)->samRecordCoordinateComparator.fileOrderCompare(R1, R2),
						samIterators
						)
					);
			
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(dict);
			for(final SAMSequenceRecord ssr: dict.getSequences()) {
				
				final IntervalTree<Boolean> capturePos;
				if(intervals!=null)
					{
					if(!Arrays.stream(intervals).anyMatch(I->I.referenceIndex == ssr.getSequenceIndex())) {
						continue;
						}
					capturePos = new IntervalTree<>();
					Arrays.stream(intervals).
							filter(I->I.referenceIndex == ssr.getSequenceIndex()).
							forEach(I->capturePos.put(I.start,I.end,true));
							;
					}
				else
					{
					capturePos = null;
					}
				
				final GenomicSequence genomicSequence;
				if(indexedFastaSequenceFile!=null && indexedFastaSequenceFile.getSequenceDictionary().getSequence(ssr.getSequenceName())!=null) {
					genomicSequence = new GenomicSequence(indexedFastaSequenceFile, ssr.getSequenceName());
				} else {
					genomicSequence = null;
				}
					
				
				final List<SAMRecord> buffer = new ArrayList<>(); 
				for( int ssr_pos=1;ssr_pos <= ssr.getSequenceLength();++ssr_pos)
					{
					if(capturePos!=null && !capturePos.overlappers(ssr_pos, ssr_pos).hasNext()) continue;
					
					progress.watch(ssr.getSequenceName(), ssr_pos);
					
					while(peekIter.hasNext())
						{
						final SAMRecord rec = peekIter.peek();
						if(rec.getReadUnmappedFlag())
							{
							peekIter.next();//consumme
							continue;
							}
						if(this.filter.filterOut(rec))
							{
							peekIter.next();//consumme
							continue;
							}
						if(rec.getReferenceIndex()< ssr.getSequenceIndex())
							{
							throw new IllegalStateException("should not happen");
							}
						if(rec.getReferenceIndex() > ssr.getSequenceIndex())
							{
							break;
							}
						if(rec.getAlignmentEnd() < ssr_pos)
							{
							throw new IllegalStateException("should not happen");
							}
						if(rec.getAlignmentStart() > ssr_pos)
							{
							break;
							}
						buffer.add(peekIter.next());
						}
					
					int x=0;
					while(x< buffer.size())
						{
						final SAMRecord R = buffer.get(x);
						if(R.getReferenceIndex()!=ssr.getSequenceIndex() || 
							R.getAlignmentEnd()<ssr_pos)
							{
							buffer.remove(x);
							}
						else
							{
							x++;
							}
						}
					
					final Map<String,PosInfo> count = 
							samples.stream().
							map(S->new PosInfo(S)).
							collect(Collectors.toMap(P->P.sample, Function.identity()));
						
					for(final SAMRecord rec:buffer)
						{
						if(rec.getReferenceIndex()!=ssr.getSequenceIndex()) throw new IllegalStateException("should not happen");
						if(rec.getAlignmentEnd() < ssr_pos) continue;
						if(rec.getAlignmentStart() > ssr_pos) continue;
						final Cigar cigar=rec.getCigar();
						if(cigar==null) continue;
			    		int refpos = rec.getAlignmentStart();
			    		final String sample = this.partition.getPartion(rec,DEFAULT_PARTITION); 
			    		for(final CigarElement ce:cigar.getCigarElements())
			    			{
			    			if(refpos > ssr_pos) break;
			    			final CigarOperator op=ce.getOperator();
			    			if(op.consumesReferenceBases())
			    				{	
			    				if(op.consumesReadBases())
			    					{
			    					if(refpos<=ssr_pos &&  ssr_pos <= refpos+ce.getLength())
			    		    			{
			    						final PosInfo posInfo = count.get(sample);
			    						if(posInfo!=null) {
				    						posInfo.dp++;
				    						if(rec.getReadNegativeStrandFlag()) {
				    							posInfo.negative_strand++;
				    							}
				    						}
			    						break;
		    		    				}
			    					}
			    				refpos += ce.getLength();
			    				}		    				
			    			}
						}
					final VariantContextBuilder vcb = new VariantContextBuilder();
					final Set<String> filters = new HashSet<>();
					
					
					
					
					vcb.chr(ssr.getSequenceName());
					vcb.start(ssr_pos);
					vcb.stop(ssr_pos);
					if(genomicSequence==null) {
						vcb.alleles(refAlleles);
						} else {
						vcb.alleles(Collections.singletonList(Allele.create((byte)genomicSequence.charAt(ssr_pos-1),true)));
						final GenomicSequence.GCPercent gcp = genomicSequence.getGCPercent(
							Math.max((ssr_pos-1)-this.gc_percent_window,0),
							Math.min(ssr_pos+this.gc_percent_window,ssr.getSequenceLength())
							);
						if(!gcp.isEmpty()) {
							vcb.attribute("GC_PERCENT",gcp.getGCPercentAsInteger());
							}
						}
					
					vcb.attribute(VCFConstants.DEPTH_KEY, count.values().stream().mapToInt(S->S.dp).sum());
					vcb.genotypes(count.values().stream().
							map(C->new GenotypeBuilder(C.sample, NO_CALLS).
									DP(C.dp).
									attribute("DR",C.negative_strand).
									attribute("DF",C.dp-C.negative_strand).
									make()).
							collect(Collectors.toList()));
					
					for(final Integer treshold: this.minDepthTresholds)
						{
						final int count_lt = (int) count.values().stream().
								filter(S->S.dp<treshold).
								count()
								;
						if(count_lt == samples.size())
							{
							filters.add("DP_LT_"+treshold);
							}
						vcb.attribute("NUM_DP_LT_"+treshold, count_lt);
						if(!samples.isEmpty())
							{
							vcb.attribute("FRACT_DP_LT_"+treshold, count_lt/(float)samples.size());
							}	
						}
					if(!samples.isEmpty())
						{
						final int array[] =  count.values().stream().
							mapToInt(S->S.dp).
							toArray();
						vcb.attribute("AVG_DP",Percentile.average().evaluate(array));
						vcb.attribute("MEDIAN_DP",Percentile.median().evaluate(array));
						vcb.attribute("MIN_DP",(int)Percentile.min().evaluate(array));
						vcb.attribute("MAX_DP",(int)Percentile.max().evaluate(array));
						}
					
					if(filters.isEmpty())
						{
						vcb.passFilters();
						}
					else
						{
						vcb.filters(filters);
						}
					out.add(vcb.make());
					}
				
				}
			progress.finish();
			peekIter.close();

			out.close();
			out = null;
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(samIterators);
			CloserUtil.close(samReaders);
			CloserUtil.close(indexedFastaSequenceFile);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar78285().instanceMainWithExit(args);

	}

}
