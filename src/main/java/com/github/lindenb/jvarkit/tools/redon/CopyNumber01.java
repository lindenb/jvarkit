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
package com.github.lindenb.jvarkit.tools.redon;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;

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
import org.apache.commons.math3.util.Precision;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.CoverageFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.stream.HtsCollectors;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence.GCPercent;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Example:

```
$ java -jar dist/copynumber01.jar  -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam
[INFO][CopyNumber01]sorting...
[INFO][CopyNumber01]fill gc%
[INFO][CopyNumber01]remove high/low gc%
[INFO][CopyNumber01]Getting coverage for RF01 N=6
[INFO][CopyNumber01]Getting coverage for RF02 N=4
[INFO][CopyNumber01]Getting coverage for RF03 N=4
[INFO][CopyNumber01]Getting coverage for RF04 N=4
[INFO][CopyNumber01]Getting coverage for RF05 N=2
[INFO][CopyNumber01]Getting coverage for RF06 N=2
[INFO][CopyNumber01]Getting coverage for RF07 N=1
[INFO][CopyNumber01]Getting coverage for RF08 N=1
[INFO][CopyNumber01]Getting coverage for RF09 N=1
[INFO][CopyNumber01]removed 0. now N=25
[INFO][CopyNumber01]median norm depth : 8.331950991034539
#CHROM	START	END	Sample	IDX	GC	RAW-DEPTH	NORM-DEPTH
RF01	0	1001	S1	0	0.321	6.410	1.015
RF01	500	1501	S1	500	0.349	8.446	1.015
RF01	1000	2001	S1	1000	0.371	9.479	1.015
RF01	1500	2501	S1	1500	0.374	9.445	1.015
RF01	2000	3001	S1	2000	0.354	7.921	1.015
RF01	2500	3302	S1	2500	0.331	5.766	1.015
RF02	0	1001	S1	3302	0.347	7.380	0.998
RF02	500	1501	S1	3802	0.348	9.189	0.998
RF02	1000	2001	S1	4302	0.341	8.672	0.998
RF02	1500	2501	S1	4802	0.344	7.977	0.998
RF03	0	1001	S1	5989	0.314	7.060	0.931
RF03	500	1501	S1	6489	0.332	9.967	0.931
RF03	1000	2001	S1	6989	0.319	9.193	0.931
RF03	1500	2501	S1	7489	0.315	7.012	0.931
RF04	0	1001	S1	8581	0.352	6.119	1.021
RF04	500	1501	S1	9081	0.344	9.554	1.021
RF04	1000	2001	S1	9581	0.374	10.000	1.021
RF04	1500	2362	S1	10081	0.359	7.396	1.021
RF05	0	1001	S1	10943	0.326	8.827	0.980
RF05	500	1501	S1	11443	0.327	8.996	0.980
RF06	0	1001	S1	12522	0.363	8.571	1.035
RF06	500	1356	S1	13022	0.421	8.384	1.035
RF07	0	1001	S1	13878	0.329	7.900	0.996
RF08	0	1001	S1	14952	0.358	7.876	1.008
RF09	0	1001	S1	16011	0.374	7.919	1.089
```

END_DOC

 */
@Program(name="copynumber01",
	description="experimental CNV detection.",
	keywords= {"cnv","bam","sam"},
	creationDate="20140201",
	modificationDate="20210318"
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
	@Parameter(names={"-w","--win-size"},description="window size. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int windowSize=1_000;
	@Parameter(names={"-s","--win-shift"},description="window shift. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int windowShift=500;
	@Parameter(names={"--win-min"},description="Discard window where length on reference is lower than 'x'. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int windowMin=100;
	
	@Parameter(names={"--univariateDepth"},description="How to calculate depth in a BAM interval.")
	private UnivariateStatistic univariateDepth = UnivariateStatistic.mean;
	@Parameter(names={"--univariateGC"},description="Loess needs only one GC value: we need to merge Depth with same GC%. How do we merge ?")
	private UnivariateStatistic univariateGCLoess = UnivariateStatistic.median;
	@Parameter(names={"--univariateMid"},description="Depth normalization. Used when we want to normalize the depths between 0.0 and 1.0")
	private UnivariateStatistic univariateMid = UnivariateStatistic.median;
	@Parameter(names={"--univariateSmooth"},description="How to smooth data with the --smooth option.")
	private UnivariateStatistic univariateSmooth = UnivariateStatistic.mean;
	
	@Parameter(names={"--smooth"},description="Smooth normalized depth window. smooth normalized depth with the 'n' neightbours")
	private int smooth_window = 5;
	@Parameter(names={"--smooth-distance"},description="When using --smooth. Only merge if windows are within that distance." + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int smoothDistance = 1_000;
	@Parameter(names={"--gcDepthInterpolation"},description="Method to interpolate GC% and depth. See https://commons.apache.org/proper/commons-math/javadocs/api-3.0/org/apache/commons/math3/analysis/interpolation/UnivariateInterpolator.html ")
	private UnivariateInerpolation gcDepthInterpolation=UnivariateInerpolation.loess;
	@Parameter(names={"--min-depth"},description="Treat depth lower than this value as 'weird' and discard the sliding windows at this place.")
	private int weirdMinDepth= 0;
	@Parameter(names={"--max-depth"},description="Treat depth greater than this value as 'weird' and discard the sliding windows at this place.")
	private int weirdMaxDepth= 500;
	@Parameter(names={"--min-gc"},description="Min GC%")
	private double minGC = 0.0;
	@Parameter(names={"--max-gc"},description="Max GC%")
	private double maxGC = 1.0;
	@Parameter(names={"--mapq"},description="Min mapping quality")
	private int mappingQuality = 1;
	@Parameter(names={"--sex"},description="Sexual contigs, comma or space separated")
	private String sexContigStr="chrX,chrY,X,Y";
	@DynamicParameter(names = "-D", description = "style. Undocumented.",hidden=true)
	private Map<String, String> dynaParams = new HashMap<>();

	private final Set<String> sexContigSet = new HashSet<>();
	
	
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
		}
	
	
	enum UnivariateInerpolation
		{
		loess,neville,difference,linear,spline,identity
		}
	private static interface XY {
		public double getX();
		public double getY();
	}
	
	private static class XYImpl implements XY {
		private final double x;
		private final double y;
		XYImpl(double x,double y) {
			this.x=x;this.y=y;
			}
		@Override public double getX() 	{ return x; }
		@Override public double getY() { return y; }
		@Override
		public String toString()
			{
			return "("+x+","+y+")";
			}
		}
	
	
	/* fact: Y=depth X=GC% */
	private static class GCAndDepth extends SimpleInterval implements XY 
		{
		GCAndDepth(final String ctg,int start,int end) {
			super(ctg,start,end);
			}
		
		double raw_depth=0;
		double norm_depth=0;
		double gc= -1;
		
		/** GC % */
		@Override public double getX() 	{ return gc; }
		/** DEPTH */
		@Override public double getY() { return norm_depth; }
		
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
					public UnivariateFunction interpolate(double[] x, double[] arg1)
							throws MathIllegalArgumentException, DimensionMismatchException {
						return new Identity();
						}
					};
			default: throw new IllegalStateException("Not implemented "+ this.gcDepthInterpolation);
			}
		}
	
	private boolean isSex(final String s) {
		return sexContigSet.contains(s);
		}
	
	private long getGenomicIndex(final SAMSequenceDictionary dict,final String ctg,int start) {
		long n= start;
		for(final SAMSequenceRecord ssr:dict.getSequences()) {
			if(ctg.equals(ssr.getSequenceName())) return n;
			n+=  ssr.getSequenceLength();
			}
		throw new IllegalArgumentException(ctg);
		}
	
	@Override
	public int doWork(final List<String> args) {
		ReferenceSequenceFile indexedFastaSequenceFile = null;
		try
			{
			this.sexContigSet.addAll(Arrays.stream(this.sexContigStr.split("[, \t]")).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toSet()));
			
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
				try(BedLineReader br =  new BedLineReader(this.bedFile)) {
					br.stream().
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
			
			LOG.info("intervals N=" + intervals.size()+" mean-size:" + intervals.stream().mapToInt(R->R.getLengthOnReference()).average().orElse(0.0));
				
			
			final List<GCAndDepth> user_items = new ArrayList<>();	
			//split intervals
			for(final Locatable loc:intervals) {
				int pos = loc.getStart();
				while(pos <  loc.getEnd())
					{
					final int pos_end = Math.min(pos + this.windowSize,loc.getEnd());
					
					final GCAndDepth dataRow = new GCAndDepth(loc.getContig(),pos,pos_end);
					if(dataRow.getLengthOnReference() < this.windowMin) {
						break;
						}
					user_items.add(dataRow);
					pos += this.windowShift;
					}
				}
			intervals.clear();//free memory
			LOG.info("sorting N=" + user_items.size());
			Collections.sort(user_items,locComparator);

			//fill gc percent
			LOG.info("fill gc% N=" + user_items.size());
			for(final String ctg: user_items.stream().map(T->T.getContig()).collect(Collectors.toSet())) {
				
				final GenomicSequence gseq = new GenomicSequence(indexedFastaSequenceFile, ctg);
				for(final GCAndDepth dataRow: user_items) {
					if(!dataRow.getContig().equals(ctg)) continue;
 					final GCPercent gc=gseq.getGCPercent(dataRow.getStart(),dataRow.getEnd());
					if(gc.isEmpty()) continue;
					dataRow.gc = gc.getGCPercent();
					}
				}
			//remove strange gc
			user_items.removeIf(B->B.gc < this.minGC);
			user_items.removeIf(B->B.gc > this.maxGC);
			LOG.info("remove high/low gc% N="+ user_items.size());
			
						
			if(user_items.stream().allMatch(P->isSex(P.getContig()))) {
				LOG.error("All chromosomes are defined as sexual. Cannot normalize");
				return -1;
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
						
						final List<GCAndDepth> bam_items = new ArrayList<>(user_items.size());

						
						/* loop over contigs */
						for(final SAMSequenceRecord ssr:dict.getSequences()) {
							/* create a **COPY** of the intervals */
							final List<GCAndDepth> ctgitems = user_items.
									stream().
									filter(T->T.contigsMatch(ssr)).
									collect(Collectors.toList());
							if(ctgitems.isEmpty()) continue;
							LOG.info("Getting coverage for "+ ssr.getSequenceName()+" N="+ctgitems.size());

							
							// get coverage
							final CoverageFactory.SimpleCoverage coverage = coverageFactory.getSimpleCoverage(samReader,ctgitems, sampleName);
							// fill coverage
							for(final GCAndDepth gc:ctgitems) {
								final OptionalDouble optCov;
								switch(this.univariateDepth) {
									case median : optCov =  coverage.getMedian(gc); break;
									case mean : optCov =  coverage.getAverage(gc); break;
									default: throw new IllegalStateException();
									}
								
								gc.raw_depth = optCov.orElse(-1.0);
								gc.norm_depth = gc.raw_depth;
								}
							ctgitems.removeIf(V->V.raw_depth < 0);
							ctgitems.removeIf(V->V.raw_depth > this.weirdMaxDepth);
							ctgitems.removeIf(V->V.raw_depth < this.weirdMinDepth);
							if(ctgitems.isEmpty()) continue;
							bam_items.addAll(ctgitems);
							}
						
						double y[] = bam_items.stream().filter(R->!isSex(R.getContig())).mapToDouble(R->R.raw_depth).toArray(); 
						LOG.info("median raw depth "+ new Median().evaluate(y,0,y.length));
						

						Collections.sort(bam_items,(a,b)->{
							final int i = Double.compare(a.getX(), b.getX());
							if(i!=0) return i;
							return Double.compare(a.getY(), b.getY());
							});
					
						
					
						double x[] = bam_items.stream().filter(R->!isSex(R.getContig())).mapToDouble(R->R.getX()).toArray();
						       y   = bam_items.stream().filter(R->!isSex(R.getContig())).mapToDouble(R->R.getY()).toArray();
						
						// get min GC
						final double min_x=x[0];
						// get max GC
						final double max_x=x[x.length-1];
						LOG.info("min/max gc "+ min_x+" "+ max_x);
						/* merge adjacent x (GC%) having same values */
						int i=0;
						int k=0;
						while(i  < x.length)
							{
							int j = i+1;
							while(j< x.length && Precision.equals(x[i],x[j]))
								{
								++j;
								}
							x[k] = x[i];
							y[k] = this.univariateGCLoess.create().evaluate(y, i, j-i);
							++k;
							i=j;
							}
						LOG.info("merged n="+(x.length-k)+" items.");
						/* reduce size of x et y */
						final List<XY> xyL = new ArrayList<>(k);
						for(int t=0; t<k; ++t) {
							xyL.add(new XYImpl(x[t], y[t]));
							}
						/* sort on Y */
						Collections.sort(xyL,(a,b)->{
							final int d = Double.compare(a.getX(), b.getX());
							if(d!=0) return d;
							return Double.compare(a.getY(), b.getY());
							});

						x  = xyL.stream().mapToDouble(R->R.getX()).toArray();
						y  = xyL.stream().mapToDouble(R->R.getY()).toArray();

						
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
						while(i<bam_items.size())
							{
							final GCAndDepth r= bam_items.get(i);
							if(spline==null)
								{
								++i;
								}
							else if(r.getX()< min_x || r.getX()> max_x)
								{
								bam_items.remove(i);
								++points_removed;
								}
							else
								{
								final double norm;
								if(this.gcDepthInterpolation.equals(UnivariateInerpolation.identity)) {
									norm = r.getY();
									}
								else
									{
									norm = spline.value(r.getX());
									}
								
								if(Double.isNaN(norm) || Double.isInfinite(norm)  )
									{
									LOG.info("NAN "+r);
									bam_items.remove(i);
									++points_removed;
									continue;
									}
								r.norm_depth = norm;
								++i;
								}
							}
						LOG.info("removed "+points_removed+". now N="+bam_items.size());
						if(bam_items.isEmpty()) continue;
						spline=null;
						// DO NOT NORMALIZE ON MINIMUM DEPTH, STUPID.
						
						
						//normalize on median
						y= bam_items.stream().mapToDouble(G->G.getY()).toArray();

						final double median_depth =  this.univariateMid.create().evaluate(y, 0, y.length);
						LOG.info("median norm depth : " + median_depth);

						for(i=0;median_depth>0 && i< bam_items.size();++i)
							{
							final GCAndDepth gc = bam_items.get(i);
							gc.norm_depth /= median_depth;
							}

						//restore genomic order
						Collections.sort(bam_items,locComparator);
						// smoothing values with neighbours
						y= bam_items.stream().mapToDouble(V->V.getY()).toArray();

						for(i=0;i< bam_items.size();++i)
							{
							final GCAndDepth gc = bam_items.get(i);
							int left=i;
							for(int j=Math.max(0,i-this.smooth_window);j<=i;++j) {
								final GCAndDepth gc2 = bam_items.get(j);
								if(!gc2.withinDistanceOf(gc,this.smoothDistance)) continue;
								left=j;
								break;
								}
							int right=i;
							for(int j=i;j<=i+this.smooth_window && j< bam_items.size();++j) {
								final GCAndDepth gc2 = bam_items.get(j);
								if(!gc2.withinDistanceOf(gc,this.smoothDistance)) break;
								right=j;
								// no break;
								}
							gc.norm_depth = this.univariateSmooth.create().evaluate(y, left,(right-left)+1);
							}

						/* print data */
						for(final GCAndDepth r:bam_items)
							{
							pw.print(r.getContig());
							pw.print('\t');
							pw.print( r.getStart() - 1);
							pw.print('\t');
							pw.print( r.getEnd());
							pw.print('\t');
							pw.print( sampleName);
							pw.print('\t');
							pw.print(getGenomicIndex(dict,r.getContig(),r.getStart()) -1);
							pw.print('\t');
							pw.printf("%.3f",r.gc);
							pw.print('\t');
							pw.printf("%.3f",r.raw_depth);
							pw.print('\t');
							pw.printf("%.3f",r.norm_depth);
							pw.println();
							}
						pw.flush();
						
						
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
