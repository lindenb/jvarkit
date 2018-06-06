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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.AbstractUnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.DataSerializable;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.iterator.MergingIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
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


/**
BEGIN_DOC


END_DOC

 */
@Program(name="naivecnvdetector",
	description="experimental CNV detection for multiple samples. Doesn't work for now.",
	keywords= {"cnv","bam","sam"},
	generate_doc=false
	)
public class NaiveCnvDetector extends Launcher
	{
	private static final Logger LOG = Logger.build(NaiveCnvDetector.class).make();
	/** reference */
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;	
	
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File refFile=null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-B","--bed","--capture"},description="Bed capture")
	private File capturebedFile=null;
	/** size of a window */
	@Parameter(names={"-w"},description="window size")
	private int windowSize=1000;
	@Parameter(names={"-s"},description="window shift")
	private int windowShift=500;
	@Parameter(names={"-filter","--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION)
	private String filterStr  = SamRecordJEXLFilter.DEFAULT_FILTER;
	@Parameter(names={"-uvsd","--mean-sample-depth"},description="How to normalize the coverage between all samples in a window.")
	private UnivariateStatistic univariateAllSamplesDepthInRegion = UnivariateStatistic.median;
	@Parameter(names={"-t","--num-threads"},description="Number of threads.")
	private int numThreads=1;
	@Parameter(names={"-et","--del-treshold"},description="Deletion Treshold")
	private double delTreshold=0.2;
	@Parameter(names={"-ut","--dup-treshold"},description="Duplication Treshold")
	private double dupTreshold=2.0;
	@Parameter(names={"-gct","--gc-treshold"},description="Segment having GC<x or GC>(1.0-x) will be discarded.")
	private double gcTreshold=0.05;
	@Parameter(names={"-af","--affected"},description="A file containing the 'affected' sample. Will be used to produce a fisher test case/control.")
	private File affectedSamplesFile = null;
	@Parameter(names={"--ignoreContigPattern"},description="Ignore those Chromosomes")
	private String ignoreContigPattern = "^(chr)?(GL.*|M|MT|NC_.*)$";	
	@Parameter(names={"--chrom"},description="Limit to that chromosome")
	private String limitToChrom=null;
	@Parameter(names={"--weirdDepth"},description="Treat normalized depth greater than this value as 'weird' and discard the sliding windows at this place.")
	private int weirdDepth=500;
	@Parameter(names={"--univariateDepth"},description="How to calculate depth in a sliding window for ONE sample")
	private UnivariateStatistic univariateOneSampleDepth = UnivariateStatistic.mean;
	@Parameter(names={"--no-variant"},description="Show No-variant")
	private boolean showNonVariant=false;
	@Parameter(names={"-nafd","--no-affected-for-depth"},description="When calculating the average depth in one window, do not include the affected samples. use when the number of affected ~= non-affected.")
	private boolean non_affected_for_depth=false;
	@Parameter(names={"-md","--min-dp"},description="At least one sample must have a normalized-depth greater than this value.")
	private int min_depth = 20;
	
	
	
	enum UnivariateStatistic
		{
		mean{
			@Override
			 	AbstractUnivariateStatistic  create() {
				return new Mean();
				}
			
		},median{
			@Override
			AbstractUnivariateStatistic  create() {
				return new Median();
				}
		};
		
		abstract  AbstractUnivariateStatistic  create();
		public double evaluate(final double array[]) {
			return create().evaluate(array); 
			}
		public double evaluate(final double array[],int start,int len) {
			return create().evaluate(array,start,len); 
			}
		}
	
	
	enum UnivariateInerpolation
		{
		loess,neville,difference,linear,spline,identity
		}
	
	private class IterableOfBed implements Iterable<Interval>
		{
		final File bedFile;
		IterableOfBed(final File bedFile) {
			this.bedFile = bedFile;
			}
		@Override
		public Iterator<Interval> iterator() {
			return new MyIter();
			}
		private class MyIter extends AbstractIterator<Interval> {
			final BedLineCodec codec=new BedLineCodec();
			final BufferedReader br;
			MyIter() {
				this.br = IOUtil.openFileForBufferedReading(bedFile);
				}
			@Override
			protected Interval advance() {
				try {
					String line;
					while((line=this.br.readLine())!=null)
						{
						final BedLine bed = this.codec.decode(line);
						if(bed==null) continue;
						if(!StringUtil.isBlank(limitToChrom) && !bed.getContig().equals(limitToChrom))
							{
							continue;
							}
						if(bed.getContig().matches(ignoreContigPattern)) {
							continue;
							}
						return bed.toInterval();
						}
					this.br.close();
					return null;
				} catch (IOException e) {
					throw new RuntimeIOException(e);
					}
				}
			}
		}	
	
	private class IterableOfDict implements Iterable<Interval>
		{
		final SAMSequenceDictionary dict;
		IterableOfDict(final SAMSequenceDictionary dict) {
			this.dict = dict;
			}
		@Override
		public Iterator<Interval> iterator() {
			return new MyIter();
			}
		private class MyIter extends AbstractIterator<Interval> {
			int tid=0;
			@Override
			protected Interval advance() {
				while(tid<dict.size())
					{
					final SAMSequenceRecord ssr = dict.getSequence(tid);
					if(!StringUtil.isBlank(limitToChrom) && !ssr.getSequenceName().equals(limitToChrom))
						{
						tid++;
						continue;
						}
					if(ssr.getSequenceName().matches(ignoreContigPattern)) {
						tid++;
						continue;
						}
					final Interval rgn =new Interval(
							ssr.getSequenceName(),
							1, ssr.getSequenceLength()
							);
					tid++;
					return rgn;
					}
				return null;
				}
			}
		}
	
	
	private class CovFileIterator
		extends AbstractIterator<CaptureData>
		implements AutoCloseable
		{
		final File covFile;
		final FileInputStream fis;
		final DataInputStream dis;
		CovFileIterator(final File covFile) {
			this.covFile = covFile;
			try {
				fis = new FileInputStream(this.covFile);
				dis = new DataInputStream(fis);
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		protected CaptureData advance() {
			try {
				final CaptureData cd = new CaptureData();
				cd.readDataInputStream(dis);
				return cd;
				}
			catch(final EOFException err) {
				close();
				return null;
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		
		@Override
		public void close() {
			CloserUtil.close(this.dis);
			CloserUtil.close(this.fis);
			this.covFile.delete();
			}
		}
	
	private class SampleInfo
		implements Callable<File> 
		{
		final File bamFile;
		SAMSequenceDictionary dict;
		String sampleName;
		int idx;
		double mean_cov;
		
		boolean affected=false;
		public SampleInfo(final File bamFile) {
			this.bamFile = bamFile;
			}
		
		@Override
		public File call() throws Exception {
			SamReader samReader = null;
			File tmpFile = null;
			FileOutputStream fos = null;
			DataOutputStream daos = null;
			try
				{
				final SamRecordFilter samRecordFilter = SamRecordJEXLFilter.create(NaiveCnvDetector.this.filterStr);
				samReader = SamReaderFactory.makeDefault().
						validationStringency(ValidationStringency.LENIENT).
						open(this.bamFile);
				
				
				if(!samReader.hasIndex())
					{
					throw new JvarkitException.UserError("Bam file is not indexed "+this.bamFile);
					}
				
				final SAMFileHeader header = samReader.getFileHeader();
				this.dict = header.getSequenceDictionary();
				this.sampleName =  samReader.getFileHeader().getReadGroups().
						stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtil.isBlank(S)).
						findFirst().
						orElse(bamFile.getName());
				
				tmpFile = File.createTempFile("cnv_"+this.sampleName,".cov");
				fos = new FileOutputStream(tmpFile);
				daos = new DataOutputStream(fos);
				
				double depth_sum =0;
				long depth_count =0L;
				
				
				final Iterable<Interval> intervalSource ;
				if(NaiveCnvDetector.this.capturebedFile==null)
					{
					intervalSource =  new IterableOfDict(this.dict);
					}
				else
					{
					intervalSource = new IterableOfBed(NaiveCnvDetector.this.capturebedFile);	
					}
				
				for(final Interval bed:intervalSource) {
					if(bed.getLengthOnReference()==0) continue;
					final SAMSequenceRecord ssr = this.dict.getSequence(bed.getContig());
					if(ssr==null) 
						{
						throw new JvarkitException.ContigNotFoundInDictionary(bed.getContig(), this.dict);
						}
					final int tid = ssr.getSequenceIndex();

					final QueryInterval qi=new QueryInterval(tid, bed.getStart(),bed.getEnd());
					if(bed.getLengthOnReference()> 1_000_000) {
						LOG.info("warning allocating memory for sizeof(int)*"+bed.getLengthOnReference());
						System.gc();	
						}
					
					final double cov[]= new double[bed.getLengthOnReference()];
					Arrays.fill(cov, 0.0);
					final SAMRecordIterator iter = samReader.query(
							new QueryInterval[] {qi},false);
					while(iter.hasNext())
						{
						final SAMRecord rec = iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(samRecordFilter.filterOut(rec)) continue;
						final Cigar cigar = rec.getCigar();
						if(cigar==null || cigar.isEmpty()) continue;
						int refpos= rec.getAlignmentStart();
						for(final CigarElement ce:cigar) {
							final CigarOperator op = ce.getOperator();
							if(op.consumesReferenceBases())
								{
								if(op.consumesReadBases())
									{
									for(int x=0;x<ce.getLength();++x)
										{
										final int array_index = (refpos+x)-bed.getStart();
										if(array_index<0) continue;
										if(array_index>= cov.length) break;
										cov[array_index]++;
										}
									}
								refpos+=ce.getLength();
								}
							}
						
						}
					iter.close();
					
					for(final double dp:cov)
						{
						depth_sum+=dp;
						depth_count++;
						}
					int x=0;
					while(x <= cov.length)
						{
						final int len = Math.min(windowSize ,cov.length-x);
						if(len> windowSize/2) {
							final CaptureData captureData = new CaptureData();
							captureData.tid = tid;
							captureData.sample_idx = this.idx;
							captureData.start = bed.getStart()+x;
							captureData.len = len;
							captureData.depth = (float)univariateOneSampleDepth.evaluate(cov,x,len);
							captureData.writeDataOutputStream(daos);
							}
						
						x+= NaiveCnvDetector.this.windowShift;
						}
					}
				this.mean_cov = depth_sum / depth_count;
				fos.flush();
				daos.flush();
				}
			finally
				{
				CloserUtil.close(daos);
				CloserUtil.close(fos);
				CloserUtil.close(samReader);
				}
			return tmpFile;
			}
		}
	
	
	private static class CaptureData
		implements DataSerializable
		{
		int tid;
		int start;
		int len;
		int sample_idx;
		float depth;
		float norm_depth;
		
		int comparePosition(final CaptureData o) {
			int i = Integer.compare(this.tid, o.tid);
			if(i!=0) return i;
			i = Integer.compare(this.start, o.start);
			return i;
			}
		
		/* int compareGC(final CaptureData o) {
			return Float.compare(this.gc_percent, o.gc_percent);
			}*/
	
		@Override
		public void readDataInputStream(final DataInputStream dai) throws IOException {
			this.tid  = dai.readInt();
			this.start = dai.readInt();
			this.len = dai.readInt();
			this.sample_idx = dai.readInt();
			this.depth = dai.readFloat();
			// don't save this.norm_depth = dai.readFloat();
			//this.gc_percent =  dai.readFloat();
			}
		@Override
		public void writeDataOutputStream(final DataOutputStream dao) throws IOException {
			dao.writeInt(this.tid);
			dao.writeInt(this.start);
			dao.writeInt(this.len);
			dao.writeInt(this.sample_idx);
			dao.writeFloat(this.depth);
			//don't read dao.writeFloat(this.norm_depth);
			//dao.writeFloat(this.gc_percent);
			}
		}
	
	
	
	@Override
	public int doWork(final List<String> args) {		
		if(this.refFile==null)
			{
			LOG.error("Undefined REF file");
			return -1;
			}
		
		final List<File> bamFiles = IOUtil.unrollFiles(args.stream().map(S->new File(S)).collect(Collectors.toList()),
				".bam");
	
		if(bamFiles.isEmpty()) {
			LOG.info("No input bam defined");
			return -1;
		}
		if(bamFiles.size()<3) {
			LOG.info("Warning: low number of bam files : "+bamFiles.size());
			return -1;
		}
		if(this.windowSize<10) {
			LOG.error("low window size");
			return -1;
		}
		if(this.windowShift<10) {
			LOG.error("low window shift");
			return -1;
		}
		
		
		try
			{
			/* loading REF Reference */
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.refFile);
			final SAMSequenceDictionary dict = this.indexedFastaSequenceFile.getSequenceDictionary();
			if(dict==null || dict.isEmpty())
				{
				throw new JvarkitException.FastaDictionaryMissing(refFile);
				}
			final Set<String> affectedSamplesNames;
			if(this.affectedSamplesFile!=null)
				{
				affectedSamplesNames = Files.
						lines(this.affectedSamplesFile.toPath()).
						filter(S->!(S.startsWith("#") || StringUtil.isBlank(S))).
						collect(Collectors.toSet());
				}
			else
				{
				affectedSamplesNames = Collections.emptySet();
				}
			
			final List<SampleInfo> sampleToIdx = new ArrayList<>();
			for(final File bamFile: bamFiles) {
				final SampleInfo si = new SampleInfo(bamFile);
				si.idx = sampleToIdx.size();
				sampleToIdx.add(si);
				}
			
			final ExecutorService executor = Executors.newFixedThreadPool(this.numThreads);
			final List<Future<File>> covFiles = executor.invokeAll(sampleToIdx);
			executor.shutdown();//first https://stackoverflow.com/questions/18425026
			executor.awaitTermination(10000L, TimeUnit.DAYS);
			
			final List<File> covTmpFiles = new ArrayList<>(covFiles.size());
			for(final Future<File> tmp: covFiles)
				{
				if(!tmp.isDone()) {
					LOG.error("not done ??");
					return -1;
				}
				final File tmpFile;
				tmpFile = tmp.get();
				covTmpFiles.add(tmpFile);
				}
			
			
			for(final SampleInfo sampleInfo: sampleToIdx) {
				if(!SequenceUtil.areSequenceDictionariesEqual(dict, sampleInfo.dict))
					{
					LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(dict,  sampleInfo.dict));
					return -1;
					}
				if(sampleToIdx.stream().
					filter(S2->S2.idx!=sampleInfo.idx).	
					anyMatch(S->S.sampleName.equals(sampleInfo.sampleName))) {
					LOG.info("duplicate sample name "+sampleInfo.sampleName);
					return -1;
					}
				sampleInfo.affected = affectedSamplesNames.contains(sampleInfo.sampleName);
				}
			
			if(sampleToIdx.stream().anyMatch(R->R.affected))
				{
				if(!sampleToIdx.stream().noneMatch(R->!R.affected)) {
					LOG.error("All samples have status affected");
					return -1;
					}
				}
			
			final double median_all_sample_depth = 
					this.univariateAllSamplesDepthInRegion.evaluate(sampleToIdx.
					stream().
					mapToDouble(S->S.mean_cov).
					toArray())
					;
			final boolean contains_affected = sampleToIdx.stream().anyMatch(S->S.affected);
			
			final VariantContextWriter vcw=super.openVariantContextWriter(this.outputFile);
			final Set<VCFHeaderLine> meta= new HashSet<>();
			final VCFFormatHeaderLine normDepthFormat = new VCFFormatHeaderLine("ND",1,VCFHeaderLineType.Float,"Normalized Depth");
			meta.add(normDepthFormat);
			final VCFInfoHeaderLine fisherInfo = new VCFInfoHeaderLine("FISHER",1,VCFHeaderLineType.Float,"Fisher Test");
			if(contains_affected) meta.add(fisherInfo);
			final VCFInfoHeaderLine gcInfo = new VCFInfoHeaderLine("GC",1,VCFHeaderLineType.Float,"GC ratio");
			meta.add(gcInfo);
			final VCFInfoHeaderLine mediaDepthInfo = new VCFInfoHeaderLine(
					this.univariateAllSamplesDepthInRegion.name().toUpperCase()+"_DEPTH",1,
					VCFHeaderLineType.Float,
					this.univariateAllSamplesDepthInRegion.name()+" Depth"
					);
			meta.add(mediaDepthInfo);
			
			final VCFFilterHeaderLine gcFilterLine = new VCFFilterHeaderLine("BAD_GC","BAD GC ratio");
			meta.add(gcFilterLine);
			final VCFFilterHeaderLine weirdDepthFilterLine = new VCFFilterHeaderLine("BAD_DEPTH","BAD Depth");
			meta.add(weirdDepthFilterLine);
			final VCFFilterHeaderLine lowDepthFilterLine = new VCFFilterHeaderLine("LOW_DEPTH","LOW Depth <" +this.min_depth);
			meta.add(lowDepthFilterLine);
			
			final Allele alleleDel = Allele.create("<DEL>", false);
			final Allele alleleDup = Allele.create("<DUP>", false);

			VCFStandardHeaderLines.addStandardFormatLines(meta, true,VCFConstants.DEPTH_KEY,VCFConstants.GENOTYPE_KEY);
			VCFStandardHeaderLines.addStandardInfoLines(meta, true,VCFConstants.DEPTH_KEY,VCFConstants.END_KEY);
			
			final VCFHeader vcfHeader = new VCFHeader(
					meta,
					sampleToIdx.stream().
						map(S->S.sampleName).
						collect(Collectors.toList())
					);
			vcfHeader.setSequenceDictionary(dict);
			vcw.writeHeader(vcfHeader);
			
			final EqualRangeIterator<CaptureData> equal_range_iter = new EqualRangeIterator<>(
					new MergingIterator<>((A,B)->A.comparePosition(B), covTmpFiles.stream().map(F->new CovFileIterator(F)).collect(Collectors.toList())),
					(A,B)->A.comparePosition(B)
					);
			while(equal_range_iter.hasNext())
				{
				boolean filter_flag=false;
				final List<CaptureData> row = equal_range_iter.next();
				if(row.isEmpty()) continue;
				if(row.size()!=sampleToIdx.size())
					{
					LOG.error("unexpected number of covdata: " +row.size()+" expected "+sampleToIdx.size());
					return -1;
					}
				final CaptureData first = row.get(0);
				
				final String contig = dict.getSequence(first.tid).getSequenceName();
				
				int total_gc =0;
				int total_N =0;
				final byte array[]= this.indexedFastaSequenceFile.
						getSubsequenceAt(contig,
								first.start,
								first.start+ first.len -1
								).getBases();
				
				for(int i=0;i < array.length;++i) 
					{
					switch (array[i]) {
						case 'c': case 'C': case 'g': case 'G': case 's': case 'S': {
							total_gc++;
							break;
							}
						case 'n': case 'N':
							total_N++;
							break;
						default:
							break;
						}
					}
				final double gc_percent;
				if(total_N>0)
					{
					continue;
					}
				else
					{
					gc_percent = total_gc / (double) array.length;
					}
				
				
				final Allele refAllele = Allele.create(array[0],true);
				
				for(final CaptureData item:row)
					{
					double sample_mean_depth = sampleToIdx.get(item.sample_idx).mean_cov;
					double ratio = sample_mean_depth / median_all_sample_depth;
					item.norm_depth = (float)(item.depth*ratio);
					//System.err.println("all: "+ median_all_sample_depth+ " d="+item.depth+" n:"+item.norm_depth+" r="+ratio);
					}
				final VariantContextBuilder vcb=new VariantContextBuilder();
				vcb.chr(contig);
				vcb.start(first.start);
				vcb.stop(first.start+first.len-1);
				vcb.attribute(VCFConstants.END_KEY, (first.start+first.len-1));
				vcb.attribute(VCFConstants.DEPTH_KEY,row.stream().mapToInt(X->(int)X.depth).sum());
				
				
				
				if(gc_percent< gcTreshold || gc_percent>(1.0-gcTreshold) ) {
					filter_flag = true;
					vcb.filter(gcFilterLine.getID());
					}
				
				//final double min_row_depth = row.stream().mapToDouble(R->R.norm_depth).min().orElse(0.0);
				final double median_row_depth = 
						new Median().evaluate(row.stream().
								filter(R->!non_affected_for_depth ||  !sampleToIdx.get(R.sample_idx).affected).
								mapToDouble(R->R.norm_depth).toArray()
						);
				
				if(median_row_depth<=0) continue;
				
				if(this.weirdDepth>0 && median_row_depth> this.weirdDepth)
					{
					filter_flag = true;
					vcb.filter(weirdDepthFilterLine.getID());
					}
				if( row.stream().allMatch(R->R.norm_depth < min_depth)) {
					filter_flag = true;
					vcb.filter(lowDepthFilterLine.getID());
					}
				
				final Set<Allele> alleles=new HashSet<>();
				alleles.add(refAllele);
				
				final List<Allele> HOM_REF= Arrays.asList(refAllele,refAllele);
				
				int affected_snv =0;
				int affected_wild =0;
				int healthy_snv =0;
				int healthy_wild =0;
				final List<Genotype> genotypes = new ArrayList<>(sampleToIdx.size());
				for(final CaptureData item:row)
					{
					final SampleInfo sampleInfo = sampleToIdx.get(item.sample_idx);
					final GenotypeBuilder gb = new GenotypeBuilder(sampleInfo.sampleName);
					final double depth = (item.norm_depth)/median_row_depth;
					
					gb.attribute(normDepthFormat.getID(), depth);
					gb.DP((int)item.depth);
					if(depth>= this.dupTreshold ) {
						if(sampleInfo.affected)
							{
							affected_snv++;
							}
						else
							{
							healthy_snv++;
							}
						alleles.add(alleleDup);
						gb.alleles(Arrays.asList(alleleDup));
						}
					else if(depth<this.delTreshold ) {
						if(sampleInfo.affected)
							{
							affected_snv++;
							}
						else
							{
							healthy_snv++;
							}
						alleles.add(alleleDel);
						gb.alleles(Arrays.asList(alleleDel));
						}
					else
						{
						if(sampleInfo.affected)
							{
							affected_wild++;
							}
						else
							{
							healthy_wild++;
							}
						gb.alleles(HOM_REF);
						}
					
					genotypes.add(gb.make());
					}
				
				if(!this.showNonVariant && alleles.size()<2) continue;
				
				if(contains_affected) {
					final FisherExactTest ft=FisherExactTest.compute(
							affected_snv,
							affected_wild,
							healthy_snv,
							healthy_wild
							);
					vcb.attribute(fisherInfo.getID(), ft.getAsDouble());
					}
				vcb.attribute(gcInfo.getID(),gc_percent);
				vcb.attribute(mediaDepthInfo.getID(), median_row_depth);
				if(!filter_flag) vcb.passFilters();
				vcb.alleles(alleles);
				vcb.genotypes(genotypes);
				vcw.add(vcb.make());
				}
			equal_range_iter.close();
			
			vcw.close();
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			}	
		}
	
	public static void main(final String[] args) {
		new NaiveCnvDetector().instanceMainWithExit(args);
		}
	}
