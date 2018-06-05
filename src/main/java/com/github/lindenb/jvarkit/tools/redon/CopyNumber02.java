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
package com.github.lindenb.jvarkit.tools.redon;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Identity;
import org.apache.commons.math3.analysis.interpolation.DividedDifferenceInterpolator;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.interpolation.NevilleInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.descriptive.AbstractUnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.DataSerializable;
import com.github.lindenb.jvarkit.io.DataSerializableCodec;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.tools.misc.GcPercentAndDepth;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter.OnNotFound;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;


/**
BEGIN_DOC


END_DOC

 */
@SuppressWarnings("unused")
@Program(name="copynumber02",
	description="experimental CNV detection. Doesn't work for now.",
	keywords= {"cnv","bam","sam"}
	)
public class CopyNumber02 extends Launcher
	{
	private static final Logger LOG = Logger.build(CopyNumber02.class).make();
	/** reference */
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;	
	/** global sam dict */
	private SAMSequenceDictionary samDictionary=null;
	
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File refFile=null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-B","--bed","--capture"},description="Bed capture",required=true)
	private File capturebedFile=null;
	
	/** size of a window */
	@Parameter(names={"-w"},description="window size")
	private int windowSize=1000;
	@Parameter(names={"-s"},description="window shift")
	private int windowShift=500;
	@Parameter(names={"-filter","--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-uvsd","--mean-sample-depth"},description="How to calculate mean depth between all samples.")
	private UnivariateStatistic univariateSampleDepth = UnivariateStatistic.mean;
	
	
	@Parameter(names={"--univariateDepth"},description="How to calculate depth in a sliding window")
	private UnivariateStatistic univariateDepth = UnivariateStatistic.mean;
	@Parameter(names={"--univariateGC"},description="Loess needs only one GC value: we need to merge Depth with same GC%. How do we merge ?")
	private UnivariateStatistic univariateGCLoess = UnivariateStatistic.median;
	@Parameter(names={"--univariateMid"},description="Vertical Depth normalization")
	private UnivariateStatistic univariateMid = UnivariateStatistic.median;
	@Parameter(names={"--univariateSmooth"},description="How to smooth data")
	private UnivariateStatistic univariateSmooth = UnivariateStatistic.mean;
	@Parameter(names={"-et","--del-treshold"},description="Deletion Treshold")
	double delTreshold=1.5;
	@Parameter(names={"-ut","--dup-treshold"},description="Duplication Treshold")
	double dupTreshold=0.5;
	@Parameter(names={"-gct","--gc-treshold"},description="Segment having GC<x or GC>(1.0-x) will be discarded.")
	double gcTreshold=0.05;

	
	
	@Parameter(names={"--ignoreContigPattern"},description="Ignore those Chromosomes")
	private String ignoreContigPattern = "^(chr)?(GL.*|M|MT|NC_.*)$";
	
	@Parameter(names={"--chrom"},description="Limit to that chromosome")
	private String limitToChrom=null;

@Parameter(names={"--gcDepthInterpolation"},description="Method to interpolate GC% and depth. See https://commons.apache.org/proper/commons-math/javadocs/api-3.0/org/apache/commons/math3/analysis/interpolation/UnivariateInterpolator.html ")
	private UnivariateInerpolation gcDepthInterpolation=UnivariateInerpolation.loess;
	@Parameter(names={"--weirdDepth"},description="Treat depth greater than this value as 'weird' and discard the sliding windows at this place.")
	private int weirdDepth=500;

	
	
	private ContigNameConverter sam2faiContigNameConverter=null;
	
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
			int pos=1;
			@Override
			protected Interval advance() {
				while(tid<dict.size())
					{
					SAMSequenceRecord ssr = dict.getSequence(tid);
					if(ssr.getSequenceName().matches(ignoreContigPattern)) {
						tid++;
						pos=1;
						continue;
						}
					if(pos+windowSize > ssr.getSequenceLength()) {
						tid++;
						pos=1;
						continue;
						}
					Interval rgn =new Interval(ssr.getSequenceName(), pos, pos+windowSize);
					pos+=windowShift;
					return rgn;
					}
				return null;
				}
			}
		}
	
	
	private class DiskSorter<T>
		implements Iterable<T>,AutoCloseable
		{
		private final Class<T> clazz;
		private final SortingCollection<T> sortingCol;
		private final SortingCollection.Codec<T> codec;
		private final Comparator<T> comparator;
		private CloseableIterator<T> iter = null;
		
		DiskSorter(final Class<T> clazz,SortingCollection.Codec<T> codec,
				final Comparator<T> comparator) 
			{
			this.clazz = clazz;
			this.comparator = comparator;
			this.codec = codec;
			this.sortingCol = SortingCollection.newInstance(
				clazz,
				codec,
				comparator,
				50000
				);
			this.sortingCol.setDestructiveIteration(true);
			
			}
		
		void add(final T rec) {
			this.sortingCol.add(rec);
			}
		@Override
		public CloseableIterator<T> iterator() {
			if(iter==null)
				{
				sortingCol.doneAdding();
				iter = this.sortingCol.iterator();
				}
			return iter;
			}
		
		DiskSorter<T> sorter(final Comparator<T> cmp2) {
			return new DiskSorter<>(this.clazz, this.codec, cmp2);
			}
				
		@Override
		public void close() {
			CloserUtil.close(iter);
			iter=null;
			sortingCol.cleanup();
			}
		}
	
				
	
	
	
	
		
	
	private UnivariateInterpolator createInterpolator()
		{	
		switch(this.gcDepthInterpolation)
			{
			case loess: return new LoessInterpolator(0.5,4);
			case neville: return new NevilleInterpolator();
			case difference  : return new DividedDifferenceInterpolator();
			case linear: return new LinearInterpolator();
			case spline : return new SplineInterpolator();
			case identity : return new UnivariateInterpolator()
					{
					@Override
					public UnivariateFunction interpolate(double[] arg0, double[] arg1)
							throws MathIllegalArgumentException, DimensionMismatchException {
						return new Identity();
						}
					};
			default: throw new IllegalStateException("Not implemented");
			}
		}
	
	private class SampleInfo
		{
		String name;
		int idx;
		double mean_cov;
		boolean affected=false;
		}
	
	private class CaptureData
		implements DataSerializable
		{
		int tid;
		int start;
		int sample_idx;
		float depth;
		float norm_depth;
		float gc_percent;
		
		int comparePos(final CaptureData o) {
			int i = Integer.compare(this.tid, o.tid);
			if(i!=0) return i;
			i = Integer.compare(this.start, o.start);
			if(i!=0) return i;
			i = Integer.compare(this.sample_idx, o.sample_idx);
			return i;
			}
		int compareGC(final CaptureData o) {
			return Float.compare(this.gc_percent, o.gc_percent);
			}
		@Override
		public void readDataInputStream(DataInputStream dai) throws IOException {
			this.tid  = dai.readInt();
			this.start = dai.readInt();
			this.sample_idx = dai.readInt();
			this.depth = dai.readFloat();
			this.norm_depth = dai.readFloat();
			this.gc_percent =  dai.readFloat();
			}
		@Override
		public void writeDataOutputStream(DataOutputStream dao) throws IOException {
			dao.writeInt(this.tid);
			dao.writeInt(this.start);
			dao.writeInt(this.sample_idx);
			dao.writeFloat(this.depth);
			dao.writeFloat(this.norm_depth);
			dao.writeFloat(this.gc_percent);
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
		
		
		
		SamReader samReader = null;
		try
			{
			DiskSorter<CaptureData> sorting1 = new DiskSorter<>(
					CaptureData.class, 
					new DataSerializableCodec<>(()->new CaptureData()),
					(A,B)->A.comparePos(B));
			
			/* loading REF Reference */
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.refFile);
			final SAMSequenceDictionary dict = this.indexedFastaSequenceFile.getSequenceDictionary();
			if(dict==null || dict.isEmpty())
				{
				throw new JvarkitException.FastaDictionaryMissing(refFile);
				}
			
			
			final List<SampleInfo> sampleToIdx = new ArrayList<>();
			for(final File bamFile: bamFiles) {
				
				LOG.info("processing "+bamFile);
				samReader = super.createSamReaderFactory().open(bamFile);
				if(!samReader.hasIndex())
					{
					LOG.error("file is not indexed "+bamFile);
					return -1;
					}
				final SAMFileHeader header = samReader.getFileHeader();
				final SAMSequenceDictionary dict2 = header.getSequenceDictionary();
				if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict2))
					{
					LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(dict, dict2));
					return -1;
					}
				final SampleInfo sampleInfo  = new SampleInfo();
				sampleInfo.idx = sampleToIdx.size();
				sampleInfo.name =  samReader.getFileHeader().getReadGroups().
						stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtil.isBlank(S)).
						findFirst().
						orElse(bamFile.getName());
				
				if(sampleToIdx.stream().anyMatch(S->S.name.equals(sampleInfo.name))) {
					LOG.info("duplicate sample name "+sampleInfo.name);
					return -1;
				}
				
				
				sampleToIdx.add(sampleInfo);
				double depth_sum =0;
				long depth_count =0L;
				
				
				final Iterable<Interval> intervalSource ;
				if(this.capturebedFile==null)
					{
					intervalSource =  new IterableOfDict(dict);
					}
				else
					{
					intervalSource = new IterableOfBed(this.capturebedFile);	
					}
				
				for(final Interval bed:intervalSource) {
					int tid = dict.getSequenceIndex(bed.getContig());
					if(tid<0) 
						{
						LOG.error(JvarkitException.ContigNotFoundInDictionary.getMessage(bed.getContig(), dict));
						return -1;
						}
					final QueryInterval qi=new QueryInterval(tid, bed.getStart(),bed.getEnd());
					final double cov[]= new double[bed.getLengthOnReference()];
					Arrays.fill(cov, 0.0);
					final SAMRecordIterator iter = samReader.query(
							new QueryInterval[] {qi},false);
					while(iter.hasNext())
						{
						final SAMRecord rec = iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(filter.filterOut(rec)) continue;
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
					while(x+this.windowSize <= cov.length)
						{
						final CaptureData captureData = new CaptureData();
						captureData.gc_percent = -1f;
						captureData.tid = tid;
						captureData.start = bed.getStart()+x;
						captureData.depth = (float)univariateDepth.evaluate(cov,x,this.windowSize);
						captureData.sample_idx = sampleInfo.idx;
						
						sorting1.add(captureData);
						x+=this.windowShift;
						}
					
					}
				samReader.close();
				
				sampleInfo.mean_cov = depth_sum / depth_count;
				
				}
			final double median_all_sample_depth = 
					this.univariateSampleDepth.evaluate(sampleToIdx.
					stream().
					mapToDouble(S->S.mean_cov).
					toArray())
					;

			
			EqualRangeIterator<CaptureData> equal_range_iter = null;
			while(equal_range_iter.hasNext())
				{
				final List<CaptureData> row = equal_range_iter.next();
				
				if(row.isEmpty()) continue;
				CaptureData first = row.get(0);
				
				int total_gc =0;
				int total_N =0;
				byte array[]= this.indexedFastaSequenceFile.
						getSubsequenceAt(dict.getSequence(first.tid).getSequenceName(),
								first.start,
								first.start+ this.windowSize -1
								).getBases();
				
				for(int i=0;i < array.length;++i) 
					{
					switch (array[i]) {
						case 'c':
						case 'C':
						case 'g':
						case 'G':
						case 's':
						case 'S': {
							total_gc++;
							break;
							}
						case 'n':
						case 'N':
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
				
				if(gcTreshold< gc_percent || gc_percent>(1.0-gcTreshold) ) {
					continue;
					}

				
				for(final CaptureData item:row)
					{
					double sample_mean_depth = sampleToIdx.get(item.sample_idx).mean_cov;
					double ratio = sample_mean_depth / median_all_sample_depth;
					item.depth = (float)(item.depth*ratio);
					item.gc_percent = (float)gc_percent;
					}
				}
			
			final VariantContextWriter vcw=super.openVariantContextWriter(this.outputFile);
			final Set<VCFHeaderLine> meta= new HashSet<>();
			final VCFFormatHeaderLine normDepthFormat = new VCFFormatHeaderLine("CN",1,VCFHeaderLineType.Float,"Normalized Depth");
			meta.add(normDepthFormat);
			final VCFInfoHeaderLine fisherInfo = new VCFInfoHeaderLine("FISHER",1,VCFHeaderLineType.Float,"Fisher Test");
			meta.add(fisherInfo);

			
			VCFStandardHeaderLines.addStandardFormatLines(meta, true,VCFConstants.DEPTH_KEY,VCFConstants.GENOTYPE_KEY);
			VCFStandardHeaderLines.addStandardInfoLines(meta, true,VCFConstants.DEPTH_KEY,VCFConstants.END_KEY);
			
			
			final VCFHeader vcfHeader = new VCFHeader(
					meta,
					sampleToIdx.stream().map(S->S.name).collect(Collectors.toList())
					);
			vcw.writeHeader(vcfHeader);
			
			while(equal_range_iter.hasNext())
				{
				final List<CaptureData> row = equal_range_iter.next();
				if(row.isEmpty()) continue;
				final CaptureData first = row.get(0);
				final String contig = dict.getSequence(first.tid).getSequenceName();
				final VariantContextBuilder vcb=new VariantContextBuilder();
				vcb.chr(contig);
				vcb.start(first.start);
				vcb.stop(first.start+this.windowSize-1);
				vcb.attribute(VCFConstants.END_KEY, (first.start+this.windowSize-1));
				vcb.attribute(VCFConstants.DEPTH_KEY,row.stream().mapToInt(X->(int)X.depth).sum());
				
				final double min_depth = row.stream().mapToDouble(R->R.norm_depth).min().orElse(0.0);
				final double median_depth = 
						new Median().evaluate(row.stream().mapToDouble(R->R.norm_depth-min_depth).toArray()
						);
				
				Set<Allele> alleles=new HashSet<>();
				Allele refAllele = Allele.create(this.indexedFastaSequenceFile.getSubsequenceAt(contig, first.start,first.start).getBases()[0],true);
				alleles.add(refAllele);
				List<Allele> HOM_REF= Arrays.asList(refAllele,refAllele);
				Allele alleleDel = Allele.create("<DEL>", false);
				Allele alleleDup = Allele.create("<DUP>", false);
				
				int affected_snv =0;
				int affected_wild =0;
				int healthy_snv =0;
				int healthy_wild =0;
				final List<Genotype> genotypes = new ArrayList<>(sampleToIdx.size());
				for(final CaptureData item:row)
					{
					final SampleInfo sampleInfo = sampleToIdx.get(item.sample_idx);
					final GenotypeBuilder gb = new GenotypeBuilder(sampleInfo.name);
					final double depth = (item.norm_depth-min_depth)/median_depth;
					
					gb.attribute(normDepthFormat.getID(), depth);
					gb.DP((int)item.depth);
					if(depth>= dupTreshold * median_depth) {
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
					else if(depth<delTreshold*median_depth) {
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
				final FisherExactTest ft=FisherExactTest.compute(
						affected_snv,
						affected_wild,
						healthy_snv,
						healthy_wild
						);
				vcb.attribute(fisherInfo.getID(), ft.getAsDouble());
				vcb.alleles(alleles);
				vcb.genotypes(genotypes);
				}
			
			
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
			CloserUtil.close(samReader);
			CloserUtil.close(this.indexedFastaSequenceFile);
			}	
		}
	
	public static void main(String[] args) {
		new CopyNumber02().instanceMainWithExit(args);
		}
	}
