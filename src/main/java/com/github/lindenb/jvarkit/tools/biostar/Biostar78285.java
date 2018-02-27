/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.iterator.MergingIterator;
import com.github.lindenb.jvarkit.util.iterator.FilterIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
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
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**

BEGIN_DOC

## Example

```bash
$ java -jar dist/biostar78285.jar -m 5 -m 10 ~/src/gatk-ui/testdata/S*.bam 
##fileformat=VCFv4.2
##FILTER=<ID=DP_LT_10,Description="All  genotypes have DP< 10">
##FILTER=<ID=DP_LT_5,Description="All  genotypes have DP< 5">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=FRACT_DP_LT_10,Number=1,Type=Float,Description="Fraction of  genotypes having DP< 10">
##INFO=<ID=FRACT_DP_LT_5,Number=1,Type=Float,Description="Fraction of  genotypes having DP< 5">
##INFO=<ID=NUM_DP_LT_10,Number=1,Type=Integer,Description="Number of  genotypes having DP< 10">
##INFO=<ID=NUM_DP_LT_5,Number=1,Type=Integer,Description="Number of  genotypes having DP< 5">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	1	.	N	.	.	DP_LT_10	DP=17;FRACT_DP_LT_10=1.0;FRACT_DP_LT_5=0.5;NUM_DP_LT_10=4;NUM_DP_LT_5=2	GT:DP	./.:5	./.:5	./.:3	./.:4
rotavirus	2	.	N	.	.	DP_LT_10	DP=21;FRACT_DP_LT_10=1.0;FRACT_DP_LT_5=0.5;NUM_DP_LT_10=4;NUM_DP_LT_5=2	GT:DP	./.:9	./.:4	./.:5	./.:3
rotavirus	3	.	N	.	.	DP_LT_10;DP_LT_5	DP=11;FRACT_DP_LT_10=1.0;FRACT_DP_LT_5=1.0;NUM_DP_LT_10=4;NUM_DP_LT_5=4	GT:DP	./.:4	./.:2	./.:4	./.:1
rotavirus	4	.	N	.	.	DP_LT_10	DP=16;FRACT_DP_LT_10=1.0;FRACT_DP_LT_5=0.5;NUM_DP_LT_10=4;NUM_DP_LT_5=2	GT:DP	./.:4	./.:5	./.:6	./.:1
rotavirus	5	.	N	.	.	DP_LT_10	DP=18;FRACT_DP_LT_10=1.0;FRACT_DP_LT_5=0.5;NUM_DP_LT_10=4;NUM_DP_LT_5=2	GT:DP	./.:5	./.:2	./.:7	./.:4
rotavirus	6	.	N	.	.	DP_LT_10	DP=14;FRACT_DP_LT_10=1.0;FRACT_DP_LT_5=0.75;NUM_DP_LT_10=4;NUM_DP_LT_5=3	GT:DP	./.:3	./.:3	./.:8	./.:0
rotavirus	7	.	N	.	.	DP_LT_10	DP=15;FRACT_DP_LT_10=1.0;FRACT_DP_LT_5=0.75;NUM_DP_LT_10=4;NUM_DP_LT_5=3	GT:DP	./.:8	./.:2	./.:4	./.:1
rotavirus	8	.	N	.	.	DP_LT_10;DP_LT_5	DP=10;FRACT_DP_LT_10=1.0;FRACT_DP_LT_5=1.0;NUM_DP_LT_10=4;NUM_DP_LT_5=4	GT:DP	./.:3	./.:3	./.:3	./.:1
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

	

    @Override
	public int doWork(final List<String> args) {
		final List<File> bamFiles = IOUtil.unrollFiles(args.stream().map(F->new File(F)).collect(Collectors.toCollection(HashSet::new)), ".bam");
		SAMSequenceDictionary dict=null;
		final List<SamReader> samReaders  = new ArrayList<>(); 
		final List<CloseableIterator<SAMRecord>> samIterators  = new ArrayList<>(); 
		final TreeSet<String> samples = new TreeSet<>(); 
		final String DEFAULT_PARTITION = "UNDEFINED_PARTITION";
		
    	VariantContextWriter out=null;
		try
			{
			
			final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			
			for(final File bamFile: bamFiles) {
				final SamReader samReader = samReaderFactory.open(bamFile);
				samReaders.add(samReader);
				final SAMFileHeader header=samReader.getFileHeader();
				if(header == null)
		    		{
		    		LOG.error("No header in "+bamFile);
		    		return -1;
		    		}
		    	if(header.getSortOrder()!=SortOrder.coordinate)
		    		{
		    		LOG.error("Sam file "+bamFile+" is not sorted on coordinate :"+header.getSortOrder());
		    		return -1;
		    		}
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
			if(dict==null) {
				LOG.error("no bam or no dictionary");
				return -1;
				}
			
			final QueryInterval intervals[];
			if(captureBed!=null)
				{
				final List<QueryInterval> L = new ArrayList<>();
				final BedLineCodec codec= new BedLineCodec();
				final LineIterator li = IOUtils.openFileForLineIterator(this.captureBed);
				while(li.hasNext()) {
					final BedLine bed = codec.decode(li.next());
					if(bed==null) continue;
					final int tid;
					if((tid=dict.getSequenceIndex(bed.getContig()))==-1)
						{
						LOG.error("not in dictionary :"+bed);
						return -1;
						}
					final QueryInterval q= new QueryInterval(tid, bed.getStart(),  bed.getEnd());
					L.add(q);
					}
				CloserUtil.close(li);
				intervals = L.toArray(new QueryInterval[L.size()]);
				}
			else
				{
				intervals = null;
				}
			for(final SamReader samReader : samReaders)
				{
				final CloseableIterator<SAMRecord> iter;
				if(intervals==null)
					{
					iter = samReader.iterator();
					}
				else
					{
					iter = samReader.queryOverlapping(intervals);
					}
				samIterators.add(iter);
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
			
			for(final Integer treshold: this.minDepthTresholds)
				{
				metaData.add(new VCFFilterHeaderLine("DP_LT_"+treshold, "All  genotypes have DP< "+treshold));
				metaData.add(new VCFInfoHeaderLine("NUM_DP_LT_"+treshold,1,VCFHeaderLineType.Integer, "Number of  genotypes having DP< "+treshold));
				metaData.add(new VCFInfoHeaderLine("FRACT_DP_LT_"+treshold,1,VCFHeaderLineType.Float, "Fraction of  genotypes having DP< "+treshold));
				}
			
			final List<Allele> refAlleles = Collections.singletonList(Allele.create("N", true));
			final List<Allele> NO_CALLS = Arrays.asList(Allele.NO_CALL,Allele.NO_CALL);
			final VCFHeader vcfHeader = new VCFHeader(metaData,samples);
			vcfHeader.setSequenceDictionary(dict);
			out.writeHeader(vcfHeader);
			
			
			final SAMRecordCoordinateComparator samRecordCoordinateComparator = new SAMRecordCoordinateComparator();
			final PeekableIterator<SAMRecord> peekIter = new PeekableIterator<>(new FilterIterator<>(
					new MergingIterator<>(
						(R1,R2)->samRecordCoordinateComparator.fileOrderCompare(R1, R2),
						samIterators
						),R->!R.getReadUnmappedFlag() && !this.filter.filterOut(R))
					);
			
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
				final List<SAMRecord> buffer = new ArrayList<>(); 
				for(int pos=1;pos <= ssr.getSequenceLength();++pos)
					{
					if(capturePos!=null && !capturePos.overlappers(pos, pos).hasNext()) continue;
	
					while(peekIter.hasNext())
						{
						final SAMRecord rec = peekIter.peek();
						if(rec.getReadUnmappedFlag())
							{
							peekIter.next();//consumme
							continue;
							}
						if(rec.getReferenceIndex()< ssr.getSequenceIndex())
							{
							peekIter.next();//consumme
							continue;
							}
						if(rec.getReferenceIndex() > ssr.getSequenceIndex())
							{
							break;
							}
						if(rec.getAlignmentEnd() < pos)
							{
							peekIter.next();//consumme
							continue;
							}
						if(rec.getAlignmentStart() > pos)
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
							R.getAlignmentEnd()<pos)
							{
							buffer.remove(x);
							}
						else
							{
							x++;
							}
						}
					
					
					final Counter<String> count = new Counter<>();
					for(final SAMRecord rec:buffer)
						{
						if(rec.getReferenceIndex()!=ssr.getSequenceIndex()) continue;
						if(rec.getAlignmentEnd() < pos) continue;
						if(rec.getAlignmentStart() > pos) continue;
						final Cigar cigar=rec.getCigar();
						if(cigar==null) continue;
			    		int refpos = rec.getAlignmentStart();
			    		final String sample = this.partition.getPartion(rec,DEFAULT_PARTITION); 
			    		for(final CigarElement ce:cigar.getCigarElements())
			    			{
			    			if(refpos > pos) break;
			    			final CigarOperator op=ce.getOperator();
			    			if(op.consumesReferenceBases())
			    				{	
			    				if(op.consumesReadBases())
			    					{
			    					if(refpos>=pos &&  pos <= refpos+ce.getLength())
			    		    			{
			    						count.incr(sample);
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
					vcb.start(pos);
					vcb.stop(pos);
					vcb.alleles(refAlleles);
					vcb.attribute(VCFConstants.DEPTH_KEY, count.getTotal());
					vcb.genotypes(samples.stream().
							map(S->new GenotypeBuilder(S, NO_CALLS).
									DP((int)count.count(S)).
									make()).
							collect(Collectors.toList()));
					
					for(final Integer treshold: this.minDepthTresholds)
						{
						final int count_lt = (int) samples.stream().
								filter(S->count.count(S)<treshold).
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
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar78285().instanceMainWithExit(args);

	}

}
