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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.redon;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
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

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
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
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.samtools.CoverageFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.stream.HtsCollectors;
import com.github.lindenb.jvarkit.tools.misc.GcPercentAndDepth;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence.GCPercent;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;


/**
BEGIN_DOC


END_DOC

 */
@SuppressWarnings("unused")
@Program(name="copynumber01",
	description="experimental CNV detection. Doesn't work for now.",
	keywords= {"cnv","bam","sam"},
	creationDate="20140201",
	modificationDate="20210315",
	generate_doc=false
	)
public class CopyNumber01 extends Launcher
	{
	private static final Logger LOG = Logger.build(CopyNumber01.class).make();
	
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refFile=null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"--bed","--capture"},description="Exome Capture as BED",required=false)
	private Path bedFile=null;

	
	/** size of a window */
	@Parameter(names={"-w"},description="window size. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int windowSize=1000;
	@Parameter(names={"-s"},description="window shift. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int windowShift=500;
	
	@Parameter(names={"--univariateDepth"},description="How to calculate depth in a sliding window")
	private UnivariateStatistic univariateDepth = UnivariateStatistic.mean;
	@Parameter(names={"--univariateGC"},description="Loess needs only one GC value: we need to merge Depth with same GC%. How do we merge ?")
	private UnivariateStatistic univariateGCLoess = UnivariateStatistic.median;
	@Parameter(names={"--univariateMid"},description="Vertical Depth normalization")
	private UnivariateStatistic univariateMid = UnivariateStatistic.median;
	@Parameter(names={"--univariateSmooth"},description="How to smooth data")
	private UnivariateStatistic univariateSmooth = UnivariateStatistic.mean;
	
	@Parameter(names={"--smooth"},description="Smooth window")
	private int SMOOTH_WINDOW=5;
	@Parameter(names={"--blackListedBed"},description="Black Listed Bed Regions")
	private File blackListedBedFile=null;
	@Parameter(names={"--gcDepthInterpolation"},description="Method to interpolate GC% and depth. See https://commons.apache.org/proper/commons-math/javadocs/api-3.0/org/apache/commons/math3/analysis/interpolation/UnivariateInterpolator.html ")
	private UnivariateInerpolation gcDepthInterpolation=UnivariateInerpolation.loess;
	@Parameter(names={"--weirdDepth"},description="Treat depth greater than this value as 'weird' and discard the sliding windows at this place.")
	private int weirdDepth=500;
	@Parameter(names={"--mapq"},description="Min mapping quality")
	private int mappingQuality = 1;
	@DynamicParameter(names = "-D", description = "style. Undocumented.",hidden=true)
	private Map<String, String> dynaParams = new HashMap<>();

	
	private enum UnivariateStatistic
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
		}
	
	
	enum UnivariateInerpolation
		{
		loess,neville,difference,linear,spline,identity
		}
	
	/* fact: Y=depth X=GC% */
	private static class GCAndDepth extends SimpleInterval
		{
		GCAndDepth(final String ctg,int start,int end) {
			super(ctg,start,end);
			}
		
		double raw_depth=0;
		double norm_depth=0;
		double gc= -1;
		
		/** GC % */
		double getX() 	{ return gc; }
		/** DEPTH */
		double getY() { return norm_depth; }
		
		@Override
		public String toString() {
			return ""+getStart()+"-"+getEnd()+" GC%="+gc+" depth:"+this.raw_depth+" norm-depth:"+this.norm_depth;
			}
		}

	
	private UnivariateInterpolator createInterpolator()
		{	
		switch(this.gcDepthInterpolation)
			{
			case loess: return new LoessInterpolator(
					Double.parseDouble(dynaParams.getOrDefault("loess.bandwith","0.75")),
					Integer.parseInt(dynaParams.getOrDefault("loess.robustness","4"))
					);
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
			default: throw new IllegalStateException("Not implemented "+ this.gcDepthInterpolation);
			}
		}
	

	
	
	@Override
	public int doWork(final List<String> args) {
		ReferenceSequenceFile indexedFastaSequenceFile = null;
		try
			{
			/* loading REF Reference */
			indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile);
			final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(indexedFastaSequenceFile);
			final Comparator<Locatable> locComparator = new ContigDictComparator(dict).createLocatableComparator();
			
			final List<Locatable> intervals = new ArrayList<>();
			if(this.bedFile==null) {
				for(final Locatable loc: dict.getSequences()) {
					intervals.add(loc);
					}
				}
			else
				{
				final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(dict);
				final BedLineCodec codec = new BedLineCodec();
				try(BufferedReader br =  IOUtils.openPathForBufferedReading(this.bedFile)) {
					br.lines().
						filter(L->!StringUtil.isBlank(L) || L.startsWith("#")).
						map(L->codec.decode(L)).
						filter(L->!StringUtil.isBlank(converter.apply(L.getContig()))).
						forEach(B->{
							final String ctg = converter.apply(B.getContig());
							intervals.add(new SimpleInterval(ctg,B.getStart(),B.getEnd()));
						});
					}
				}
			
			if(intervals.isEmpty()) {
				LOG.error("no interval defined.");
				return -1;
				}
			
			final Map<String,List<GCAndDepth>> contig2items = new HashMap<>(dict.size());
			
			//split intervals
			for(final Locatable loc:intervals) {
				List<GCAndDepth> L = contig2items.get(loc.getContig());
				if(L==null) {
					L  = new ArrayList<>();
					contig2items.put(loc.getContig(),L);
					}
				int pos = loc.getStart();
				while(pos <  loc.getEnd())
					{
					final int pos_end = Math.min(pos + this.windowSize,loc.getEnd());
					
					if( (pos_end - pos) < this.windowSize*0.8) {
						break;
						}
					
					final GCAndDepth dataRow = new GCAndDepth(loc.getContig(),pos,pos_end);
					L.add(dataRow);
					
					pos += this.windowShift;
					}
				}
			intervals.clear();//free memory
			
			//fill gc percent
			LOG.info("fill gc%");
			for(final String ctg:contig2items.keySet()) {
				final List<GCAndDepth> items = contig2items.get(ctg);
				/* sort interval */
				Collections.sort(items,locComparator);
				
				final GenomicSequence gseq = new GenomicSequence(indexedFastaSequenceFile, ctg);
				for(final GCAndDepth dataRow:items) {
					final GCPercent gc=gseq.getGCPercent(dataRow.getStart(),dataRow.getEnd());
					if(gc.isEmpty()) continue;
					dataRow.gc = gc.getGCPercent();
					}
				}
			//remove strange gc
			LOG.info("remove string gc%");
			for(List<GCAndDepth> items:contig2items.values()) {
				items.removeIf(B->B.gc<=0.0);
				}
			
			final CoverageFactory coverageFactory = new CoverageFactory().setMappingQuality(this.mappingQuality);
			
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				/* header */
				pw.println("#CHROM\tSTART\tEND\tSample\tIDX\tGC\tRAW-DEPTH\tNORM-DEPTH");

				for(final Path bamPath: IOUtils.unrollPaths(args)) {
					//open this samReader
					try(SamReader samReader = super.createSamReaderFactory().referenceSequence(this.refFile).open(bamPath)) {
						if(!samReader.hasIndex())
							{
							LOG.error("file is not indexed "+bamPath);
							return -1;
							}
						final SAMFileHeader header = samReader.getFileHeader();
	
						SequenceUtil.assertSequenceDictionariesEqual(dict, SequenceDictionaryUtils.extractRequired(header));
					
						final String sampleName = header.getReadGroups().
								stream().
								map(RG->RG.getSample()).
								filter(S->!StringUtil.isBlank(S)).
								collect(Collectors.toSet()).
								stream().
								collect(HtsCollectors.toSingleton());
						
						/* loop over contigs */
						for(final SAMSequenceRecord ssr:dict.getSequences()) {
							if(!contig2items.containsKey(ssr.getSequenceName())) continue;
							
							final Function<Integer, Long> genomicIndex = START ->
								{
								long n= START;
								for(int i=0;i< ssr.getSequenceIndex();++i) n+=  dict.getSequence(i).getSequenceLength();
								return n;
								};
						

							/* create a **COPY** of the intervals */
							final List<GCAndDepth> items = new ArrayList<>(contig2items.get(ssr.getContig()));
							if(items.isEmpty()) continue;
							
							// get coverage
							final CoverageFactory.SimpleCoverage coverage = coverageFactory.getSimpleCoverage(samReader,items, sampleName);
							// fill coverage
							for(final GCAndDepth gc:items) {
								final OptionalDouble optCov;
								switch(this.univariateDepth) {
									case median : optCov =  coverage.getMedian(gc); break;
									case mean : optCov =  coverage.getAverage(gc); break;
									default: throw new IllegalStateException();
									}
								
								gc.raw_depth = optCov.orElse(-1.0);
								gc.norm_depth = gc.raw_depth;
								}
							items.removeIf(V->V.raw_depth<0);
							items.removeIf(V->V.raw_depth>this.weirdDepth);
							if(items.isEmpty()) continue;
							
							Collections.sort(items,(a,b)->{
								final int i = Double.compare(a.getX(), b.getX());
								if(i!=0) return i;
								return Double.compare(a.getY(), b.getY());
								});
						
						
							double x[]= items.stream().mapToDouble(R->R.getX()).toArray();
							double y[]= items.stream().mapToDouble(R->R.getY()).toArray();
							
							// get min GC
							final double min_x=x[0];
							// get max GC
							final double max_x=x[x.length-1];
							
							/* merge adjacent x (GC%) having same values */
							int i=0;
							int k=0;
							while(i  < x.length)
								{
								int j = i+1;
								while(j< x.length && Double.compare(x[i],x[j])==0)
									{
									++j;
									}
								x[k] = x[i];
								y[k] = this.univariateGCLoess.create().evaluate(y, i, j-i);
								++k;
								i=j;
								}
					
							/* reduce size of x et y */
							if(k != x.length)
								{
								x = Arrays.copyOf(x, k);
								y = Arrays.copyOf(y, k);
								}
							
							
							final UnivariateInterpolator interpolator = createInterpolator();
							UnivariateFunction  spline = null;
							try {
								spline = interpolator.interpolate(x, y);
								}
							catch(final org.apache.commons.math3.exception.NumberIsTooSmallException err)
								{
								spline=null;
								LOG.error("Cannot use "+interpolator.getClass().getName()+":"+err.getMessage());
								}
							//min depth cal
							int points_removed=0;
							i=0;
							while(i<items.size())
								{
								final GCAndDepth r= items.get(i);
								if(spline==null)
									{
									++i;
									}
								else if(r.getX()< min_x || r.getX()> max_x)
									{
									items.remove(i);
									++points_removed;
									}
								else
									{
									final double norm = spline.value(r.getX());
									if(Double.isNaN(norm) || Double.isInfinite(norm)  )
										{
										LOG.info("NAN "+r);
										items.remove(i);
										++points_removed;
										continue;
										}
									r.norm_depth = norm; 
									++i;
									}
								}
							if(items.isEmpty()) continue;
							spline=null;
							final double min_norm_depth  = items.stream().mapToDouble(G->G.norm_depth).min().orElse(Double.NaN);
							if(Double.isNaN(min_norm_depth)) continue;
							
							//fit to min, fill new y for median calculation
							for(final GCAndDepth gc:items) {
								gc.norm_depth -= min_norm_depth;
								}
							
							
							//normalize on median
							y= items.stream().mapToDouble(G->G.getY()).toArray();

							final double median_depth =  this.univariateMid.create().evaluate(y, 0, y.length);
							for(i=0;median_depth>0 && i< items.size();++i)
								{
								final GCAndDepth gc = items.get(i);
								gc.norm_depth /= median_depth;
								}
														
							//restore genomic order
							Collections.sort(items,locComparator);
							
							//  smoothing values with neighbours 
							y= items.stream().mapToDouble(V->V.getY()).toArray();
							
							for(i=0;i< items.size();++i)
								{
								final GCAndDepth gc = items.get(i);
								final int left = Math.max(i,i-SMOOTH_WINDOW);
								final int right = Math.min(y.length-1,i+SMOOTH_WINDOW);
								gc.norm_depth = this.univariateSmooth.create().evaluate(y, left,(right-left)+1);
								}

							/* print data */
							for(final GCAndDepth r:items)
								{
								pw.print(ssr.getSequenceName());
								pw.print('\t');
								pw.print( r.getStart() - 1);
								pw.print('\t');
								pw.print( r.getEnd());
								pw.print('\t');
								pw.print( sampleName);
								pw.print('\t');
								pw.print(genomicIndex.apply(r.getStart() -1));
								pw.print('\t');
								pw.printf("%.3f",r.gc);
								pw.print('\t');
								pw.printf("%.3f",r.raw_depth);
								pw.print('\t');
								pw.printf("%.3f",r.norm_depth);
								pw.println();
								}
							pw.flush();
							}
						
						}// samReader
					}//end loop bamPath
				pw.flush();
				}//output
	
				
				
				
			
			
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(indexedFastaSequenceFile);
			}	
		}
	
	public static void main(String[] args) {
		new CopyNumber01().instanceMainWithExit(args);
		}
	}
