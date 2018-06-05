/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
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
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
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
	
	
	
	@Parameter(names={"--ignoreContigPattern"},description="Ignore those Chromosomes")
	private String ignoreContigPattern = "^(chr)?(GL.*|M|MT|NC_.*)$";
	@Parameter(names={"--smooth"},description="Smooth window")
	private int SMOOTH_WINDOW=5;
	@Parameter(names={"--blackListedBed"},description="Black Listed Bed Regions")
	private File blackListedBedFile=null;
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
	
	
	/* fact: Y=depth X=GC% */
	private class GCAndDepth
		{
		int start;
		int end;
		double raw_depth=0;
		double norm_depth=0;
		double gc=0;
		
		/** GC % */
		double getX()
			{
			return gc;
			}
		/** DEPTH */
		double getY()
			{
			return norm_depth;
			}
		
		
		
		
		@Override
		public String toString() {
			return ""+start+"-"+end+" GC%="+gc+" depth:"+this.raw_depth+" norm-depth:"+this.norm_depth;
			}
		}

			
			
	
	
	/** constructor */
	private CopyNumber02()
		{
		}
		
	
	private boolean ignoreChromosomeName(final String chrom)
		{
		return Pattern.compile(ignoreContigPattern).matcher(chrom).matches();
		}
	
	
	private class ContigProcessor
		{
		private final String sampleName;
		private final SamReader samReader;
		private final List<GCAndDepth> items ;
		private final SAMSequenceRecord ssr;
		ContigProcessor(final SamReader samReader,
				final SAMSequenceRecord ssr,
				final String sampleName
				)
			{
			this.samReader = samReader;
			this.sampleName = sampleName;
			this.ssr = ssr;
			this.items = new ArrayList<>(100000);
			
			final int contig_len = this.getSequenceLength();
			int pos = 0;
			while(pos <  contig_len)
				{
				final int pos_end = Math.min(pos + CopyNumber02.this.windowSize,contig_len);
				
				if( (pos_end - pos) < CopyNumber02.this.windowSize*0.8)
					{
					break;
					}
				
				final GCAndDepth dataRow = new GCAndDepth();
				dataRow.start = pos+1;
				dataRow.end = pos_end;
				this.items.add(dataRow);
				
				pos += CopyNumber02.this.windowShift;
				}
			LOG.debug("number of items = " + this.items.size());
			}
		
		public String getContig()
			{
			return this.ssr.getSequenceName();
			}
		
		public int getSequenceLength() {
			return this.ssr.getSequenceLength();
			}
		
		void removeBlackListed()
			{
			if(blackListedBedFile==null) return;
			final BedLineCodec codec= new BedLineCodec();
			BufferedReader r=null;
			int countBefore = this.items.size();
			try
				{
				r = IOUtils.openFileForBufferedReading(blackListedBedFile);
				final IntervalTree<Boolean> badRgns = new IntervalTree<>();
				r.lines().map(L->codec.decode(L)).
					filter(B->B!=null).
					filter(B->B.getContig().equals(getContig())).
					forEach(B->{
						badRgns.put(B.getStart(),B.getEnd(), true);
						});
				r.close();r=null;
				
				int i=0;
				while(i<this.items.size())
					{
					GCAndDepth item = this.items.get(i);
					boolean ok=true;
					Iterator<IntervalTree.Node<Boolean>> it = badRgns.overlappers(item.start,item.end);
					while(it.hasNext())
						{
						final IntervalTree.Node<Boolean> node=it.next();
						if(node.getStart() < item.start && item.start< node.getEnd() && node.getEnd()<=item.end)
							{
							item.start=node.getEnd()+1;
							}
						else if(item.start <=node.getStart()  && node.getStart() <=item.end && item.end<node.getEnd())
							{
							item.end=node.getStart()-1;
							}
						else
							{
							ok=false;
							}
						if(item.end-item.start< 0.8*windowSize)
							{
							ok=false;
							}
						}
					if(!ok)
						{
						this.items.remove(i);
						}
					else
						{
						i++;
						}
					}
				if(countBefore!=this.items.size())
					{
					LOG.info("removed "+(countBefore-this.items.size())+" blacklisted regions.");
					}
				}
			catch(IOException err)
				{
				throw new RuntimeIOException(err);
				}
			finally
				{
				CloserUtil.close(r);
				}	
			}
		
		void fillGC()
				{
				LOG.debug("Fill GC% for "+this.getContig());
				final String faixContig = CopyNumber02.this.sam2faiContigNameConverter.apply(this.getContig());
				final GenomicSequence genomic = new GenomicSequence(CopyNumber02.this.indexedFastaSequenceFile,
						faixContig);
				final int chrom_end = this.getSequenceLength();
				int i=0;
				
				while( i < this.items.size()) 
					{
					final GCAndDepth dataRow = this.items.get(i);
					double total_gc = 0.0;
					double total_bases = 0.0;
					double total_N = 0.0;
					boolean foundN = true;
	
					for (int pos = dataRow.start;pos <= dataRow.end && pos <= chrom_end; ++pos ) {
						final char c = genomic.charAt(pos - 1);
						switch (c) {
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
						++total_bases;
					}
				if(total_N>0)
					{
					this.items.remove(i);
					}
				else
					{
					dataRow.gc = total_gc / (double) total_bases;
					i++;
					}
				}
			}
		
		/** get whole contig coverage */
		private void scanCoverage()
			{
			final int contig_length =  this.getSequenceLength();
			LOG.debug("Fill Depth for "+this.getContig()+" allocating sizeof(short)*"+contig_length);

			// region were we have items
			final IntervalTree<Boolean> itemRgns = new IntervalTree<>();
			
			{
			int i=0;
			
			while(i < this.items.size())
				{
				final GCAndDepth item_i = this.items.get(i);
				int end= item_i.end;
				int j=i+1;
				//reduce footprint: merge adjascents item to build itemRgns
				while(j<this.items.size())
					{
					final GCAndDepth item_j = this.items.get(j);

					if(!(item_i.end > item_j.start && item_i.end < item_j.end ))
						{
						break;
						}
					end = item_j.end;
					++j;
					}
				itemRgns.put(item_i.start, end, true);
				i=j;
				}
			}
			
			short array[]= new short[ 1 + contig_length];
			Arrays.fill(array, (short)0);
			final SAMRecordIterator sri=this.samReader.query(
					getContig(),
					1,
					array.length,
					false);
			while(sri.hasNext())
				{
				final SAMRecord rec = sri.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(CopyNumber02.this.filter.filterOut(rec)) continue;
				
				if(itemRgns.minOverlapper(rec.getStart(), rec.getEnd())==null) continue;
				
				final Cigar c= rec.getCigar();
				int refStart= rec.getAlignmentStart();
				for(final CigarElement ce:c.getCigarElements())
					{
					if(!ce.getOperator().consumesReferenceBases()) continue;
					if(ce.getOperator().consumesReadBases())
						{
						for(int x=0;x< ce.getLength();++x)
							{
							final int pos = refStart + x;
							if(pos >= 1 && pos <= contig_length && array[pos]<Short.MAX_VALUE)
								{
								array[ pos  ]++;
								}
							}
						}
					refStart+=ce.getLength();
					}		
				}
			sri.close();
			
			
			for(final GCAndDepth row:this.items)
				{
				final double sum_array[]= new double[1 + row.end - row.start];
				Arrays.fill(sum_array, 0.0);
				for(int i=0;i< sum_array.length;++i)
					{
					sum_array[ i ] = array[ i + (row.start-1) ];
					}
				row.raw_depth = CopyNumber02.this.univariateDepth.evaluate(sum_array);
				row.norm_depth = row.raw_depth;
				}
			array=null;
			System.gc();
			
			final int oldCount=this.items.size();
			if(this.items.removeIf(T->T.raw_depth >= weirdDepth))
				{
				LOG.info(""+(oldCount-items.size())+" items had a depth greater than treshold ("+weirdDepth+")");
				}
			
			}
		
		private void normalizeCoverage()
			{
			Collections.sort(this.items,(a,b)->{
				final int i = Double.compare(a.getX(), b.getX());
				if(i!=0) return i;
				return Double.compare(a.getY(), b.getY());
				});
			
			
			double x[]=new double[this.items.size()];
			double y[]=new double[this.items.size()];
	
			int i=0;
			for(i=0;i< this.items.size();++i)
				{
				final GCAndDepth r=this.items.get(i);
				x[i] = r.getX();
				y[i] = r.getY();
				}
			
			final double min_x=x[0];
			final double max_x=x[x.length-1];
			
			/* merge adjacent x (GC%) having same values */
			i=0;
			int k=0;
			while(i  < x.length)
				{
				int j = i+1;
				
				while(j< x.length && Double.compare(x[i],x[j])==0)
					{
					++j;
					}
				x[k] = x[i];
				y[k] = CopyNumber02.this.univariateGCLoess.create().evaluate(y, i, j-i);
				++k;
				i=j;
				}
	
			/* reduce size of x et y */
			if(k != x.length)
				{
				LOG.info("Compacting X from "+x.length+" to "+k);
				x = Arrays.copyOf(x, k);
				y  =Arrays.copyOf(y, k);
				}
			
			//min depth cal
			double min_depth=Double.MAX_VALUE;
			
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
				
			int points_removed=0;
			i=0;
			while(i<this.items.size())
				{
				final GCAndDepth r= this.items.get(i);
				if(spline==null)
					{
					min_depth=Math.min(min_depth,r.norm_depth);
					++i;
					}
				else if(r.getX()< min_x || r.getX()> max_x)
					{
					this.items.remove(i);
					++points_removed;
					}
				else
					{
					double norm = spline.value(r.getX());
					if(Double.isNaN(norm) || Double.isInfinite(norm)  )
						{
						LOG.info("NAN "+r);
						this.items.remove(i);
						++points_removed;
						continue;
						}
					r.norm_depth -= norm; 
					min_depth=Math.min(min_depth,r.norm_depth);
					++i;
					}
				}
			LOG.info("Removed "+points_removed+" because GC% is too small (Sexual chrom)" );
			spline=null;
			
			
			//fit to min, fill new y for median calculation
	
			y= new double[this.items.size()];
			for(i=0;i< this.items.size();++i)
				{
				final GCAndDepth gc= this.items.get(i);
				gc.norm_depth -= min_depth;
				y[i] = gc.norm_depth;
				}
			
			//normalize on median
			final double median_depth =  CopyNumber02.this.univariateMid.create().evaluate(y, 0, y.length);
			for(i=0;median_depth>0 && i< this.items.size();++i)
				{
				final GCAndDepth gc = this.items.get(i);
				gc.norm_depth /= median_depth;
				}
			
			// ok, we now have a normalization between 0 and 1. Let's multiply by 2 because we have a diploid
			for(final GCAndDepth gc:this.items)
				{
				gc.norm_depth*=2.0;
				}
			
			//restore genomic order
			Collections.sort(this.items,(A,B)->A.start-B.start);
			
			
			//  smoothing values with neighbours 
			y= new double[this.items.size()];
			for(i=0;i< this.items.size();++i)
				{
				y[i] = this.items.get(i).getY();
				}
			for(i=0;i< this.items.size();++i)
				{
				final GCAndDepth gc = this.items.get(i);
				final int left = Math.max(i,i-SMOOTH_WINDOW);
				final int right = Math.min(y.length-1,i+SMOOTH_WINDOW);
				gc.norm_depth = CopyNumber02.this.univariateSmooth.create().evaluate(y, left,(right-left)+1);
				}
			}
	

	
		private void saveCoverage( final PrintWriter pw )
			{
			LOG.info("Dumping coverage for "+getContig());
			
			final Function<Integer, Long> genomicIndex = START ->
					{
					long n= START;
					for(int i=0;i< ssr.getSequenceIndex();++i) n+=  samDictionary.getSequence(i).getSequenceLength();
					return n;
					};
			
			/* header */
			pw.println("#CHROM\tSTART\tEND\tIDX\tGC\tRAW-DP\tNORM-DP");
			
			/* get data */
			for(final GCAndDepth r:this.items)
				{
				pw.print(ssr.getSequenceName());
				pw.print('\t');
				pw.print( r.start);
				pw.print('\t');
				pw.print( r.end);
				pw.print('\t');
				pw.print(genomicIndex.apply(r.start));
				pw.print('\t');
				pw.print(r.gc);
				pw.print('\t');
				pw.print(r.raw_depth);
				pw.print('\t');
				pw.print(r.norm_depth);
				pw.println();
				}
			pw.flush();
			}
		
		public void run()
			{
			removeBlackListed();
			fillGC();
			scanCoverage();
			normalizeCoverage();
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
		{
		int tid;
		int start;
		int sample_idx;
		float depth;
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
		void read(DataInputStream dai) throws IOException
			{
			this.tid  = dai.readInt();
			this.start = dai.readInt();
			this.sample_idx = dai.readInt();
			this.depth = dai.readFloat();
			this.gc_percent =  dai.readFloat();
			}
		void write(DataOutputStream dao) throws IOException
			{
			dao.writeInt(this.tid);
			dao.writeInt(this.start);
			dao.writeInt(this.sample_idx);
			dao.writeFloat(this.depth);
			dao.writeFloat(this.gc_percent);
			}
		}
	
	private class CaptureDataCodec extends AbstractDataCodec<CaptureData>
		{
		@Override
		public CaptureData decode(DataInputStream dis) throws IOException {
			try
				{
				final CaptureData cd = new CaptureData();
				cd.read(dis);
				return cd;
				}
			catch(EOFException err)
				{
				return null;
				}
			}
		@Override
		public void encode(DataOutputStream dos, CaptureData object) throws IOException {
			object.write(dos);
			}
		@Override
		public AbstractDataCodec<CaptureData> clone() {
			return new CaptureDataCodec();
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
					new CaptureDataCodec(),
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
				
				
				final BedLineCodec bedLineCodec = new BedLineCodec();
				final BufferedReader br = IOUtils.openFileForBufferedReading(this.capturebedFile);
				String line;
				while((line=br.readLine())!=null) {
					final BedLine bed= bedLineCodec.decode(line);
					if(bed==null) continue;
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
				
				final VariantContextBuilder vcb=new VariantContextBuilder();
				vcb.chr(dict.getSequence(first.tid).getSequenceName());
				vcb.start(first.start);
				vcb.stop(first.start+this.windowSize-1);
				vcb.attribute("END", (first.start+this.windowSize-1));
				
				final double min_depth = row.stream().mapToDouble(R->R.depth).min().orElse(0.0);
				final double median_depth = 
						new Median().evaluate(row.stream().mapToDouble(R->R.depth-min_depth).toArray()
						);
				
				
				List<Genotype> genotypes = new ArrayList<>();
				for(final CaptureData item:row)
					{
					GenotypeBuilder gb = new GenotypeBuilder(sampleToIdx.get(item.sample_idx).name);
					double depth = (item.depth-min_depth)/median_depth;
					
					gb.attribute(normDepthFormat.getID(), depth);
					if(depth>1.5*median_depth) {
						
						}
					else if(depth<1.5*median_depth) {
						
						}
					else
						{
						
						}
					
					genotypes.add(gb.make());
					}
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
