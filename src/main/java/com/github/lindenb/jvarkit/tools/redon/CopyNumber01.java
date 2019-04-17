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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.redon;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
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
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.tools.misc.GcPercentAndDepth;
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


/**
BEGIN_DOC


END_DOC

 */
@SuppressWarnings("unused")
@Program(name="copynumber01",
	description="experimental CNV detection. Doesn't work for now.",
	keywords= {"cnv","bam","sam"},
	modificationDate="20190417"
	)
public class CopyNumber01 extends Launcher
	{
	private static final Logger LOG = Logger.build(CopyNumber01.class).make();
	/** reference */
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;	
	/** global sam dict */
	private SAMSequenceDictionary samDictionary=null;
	
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File refFile=null;
	@Parameter(names={"-o","--out"},description="output base name",required=true)
	private File archiveFile=null;
	/** size of a window */
	@Parameter(names={"-w"},description="window size")
	private int windowSize=1000;
	@Parameter(names={"-s"},description="window shift")
	private int windowShift=500;
	@Parameter(names={"-filter","--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	private SamRecordFilter filter  = SamRecordFilterFactory.getDefault();
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
	@Parameter(names={"--prefix"},description="File prefix")
	private String prefix="";
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
		}
	
	
	enum UnivariateInerpolation
		{
		loess,neville,difference,linear,spline,identity
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
	private CopyNumber01()
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
				final int pos_end = Math.min(pos + CopyNumber01.this.windowSize,contig_len);
				
				if( (pos_end - pos) < CopyNumber01.this.windowSize*0.8)
					{
					break;
					}
				
				final GCAndDepth dataRow = new GCAndDepth();
				dataRow.start = pos+1;
				dataRow.end = pos_end;
				this.items.add(dataRow);
				
				pos += CopyNumber01.this.windowShift;
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
				final String faixContig = CopyNumber01.this.sam2faiContigNameConverter.apply(this.getContig());
				final GenomicSequence genomic = new GenomicSequence(CopyNumber01.this.indexedFastaSequenceFile,
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
				if(CopyNumber01.this.filter.filterOut(rec)) continue;
				
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
				row.raw_depth = CopyNumber01.this.univariateDepth.evaluate(sum_array);
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
				y[k] = CopyNumber01.this.univariateGCLoess.create().evaluate(y, i, j-i);
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
			final double median_depth =  CopyNumber01.this.univariateMid.create().evaluate(y, 0, y.length);
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
				gc.norm_depth = CopyNumber01.this.univariateSmooth.create().evaluate(y, left,(right-left)+1);
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
	

	
	
	@Override
	public int doWork(final List<String> args) {		
		if(this.refFile==null)
			{
			LOG.error("Undefined REF file");
			return -1;
			}
		
		
		if(this.archiveFile==null)
			{
			LOG.error("Undefined output file.");
			return -1;
			}
	
		if(!StringUtil.isBlank(prefix))
			{
			if(!(prefix.endsWith(".") || prefix.endsWith("_")))
				{
				prefix = prefix +".";
				}
			}
		
		
		SamReader samReader = null;
		ArchiveFactory archive=null;
		try
			{
			final String input = oneAndOnlyOneFile(args);
			samReader = super.openSamReader(input);
			if(!samReader.hasIndex())
				{
				LOG.error("file is not indexed "+input);
				return -1;
				}
			
			final String sampleName =  samReader.getFileHeader().getReadGroups().
					stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtil.isBlank(S)).
					findFirst().
					orElse("SAMPLE");
			
			final SAMFileHeader header = samReader.getFileHeader();
			
			
			this.samDictionary= header.getSequenceDictionary();
			if(this.samDictionary==null || this.samDictionary.isEmpty()) {
				throw new JvarkitException.DictionaryMissing(input);
				}
			
		
			
			/* loading REF Reference */
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(refFile);
			final SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
			if(dict==null)
				{
				throw new JvarkitException.DictionaryMissing(refFile.getPath());
				}
			
			this.sam2faiContigNameConverter = ContigNameConverter.fromDictionaries(this.samDictionary, dict);
			
			archive = ArchiveFactory.open(this.archiveFile);

			
			PrintWriter mkw = archive.openWriter(this.prefix+"cnv01.mk");
			mkw.println("## "+this.getProgramCommandLine());
			mkw.println(".PHONY=all all2");
			mkw.println("TARGETS=");
			mkw.println("SCREEN_WIDTH?=2600");
			mkw.println("SCREEN_HEIGHT?=1000");
			mkw.println("GNUPLOT?=gnuplot");
			mkw.println("all : all2");
			
			final List<String> all_chrom_files= new ArrayList<>();
			double maxDepth=0.0;
			
			for(final SAMSequenceRecord ssr: this.samDictionary.getSequences()) {
				if(!StringUtil.isBlank(this.limitToChrom) && !this.limitToChrom.equals(ssr.getSequenceName()))
					{
					LOG.info("Skipping "+ssr.getSequenceName()+"....");
					continue;
					}
				
				if(this.sam2faiContigNameConverter.apply(ssr.getSequenceName())==null)
					{
					LOG.info("Ignoring "+ssr.getSequenceName()+" because it's not in REF");
					continue;
					}
				if(ignoreChromosomeName(ssr.getSequenceName()))
					{
					LOG.info("Ignoring "+ssr.getSequenceName());
					continue;
					}
				LOG.info("processing chromosome "+ssr.getSequenceName());
				ContigProcessor proc = new ContigProcessor(samReader, ssr,sampleName);
				
				
				proc.run();
				
				final String tsvFilename=this.prefix + proc.getContig()+"."+proc.sampleName+".bed";
				final String pngFilename="$(addsuffix .png,$(basename "+tsvFilename+"))";
				final double depth = proc.items.stream().mapToDouble(I->I.norm_depth).max().orElse(4);
				
				maxDepth  = Math.max(maxDepth, depth);
				
				all_chrom_files.add(tsvFilename);
				mkw.println("TARGETS+="+pngFilename);
				mkw.println(pngFilename+":"+tsvFilename);
				mkw.println("\techo 'set terminal png truecolor size ${SCREEN_WIDTH},${SCREEN_HEIGHT};"
						+ "set key autotitle columnhead;"
						+ "set datafile separator \"\t\";"
						+ "set title \""+proc.sampleName+" "+ssr.getSequenceName()+"\";"
						+ "set ylabel \"Normalized Depth - Number of copies\";"
						+ "set xlabel \"Position on "+ssr.getSequenceName()+"\";"
						//+ "set style circle radius 0.02;" 
						+ "set nokey;"
						+ "set yrange [0:"+Math.min(4,Math.ceil(depth))+"];"
						+ "set xrange [1:"+ssr.getSequenceLength()+"];"
						+ "set xtic rotate by 90 right;"
						+ "set output \"$@\";"
						//+ "plot \"$<\" using 1:2:1 w points linecolor variable "
						+ "plot \"$<\" using 2:7 w lines "
						+ "' | ${GNUPLOT}");
				
				
				PrintWriter pw= archive.openWriter(tsvFilename);
				proc.saveCoverage(pw);
				pw.flush();pw.close();
				proc=null;
				}
			
			
			samReader.close();
			
			
			final String pngFilename=this.prefix+"wholeGenome.png";
			mkw.println("TARGETS+="+pngFilename);
			mkw.println(pngFilename+":"+String.join(" ", all_chrom_files));
			mkw.println("\tgrep -v '^#' $^ | cut -f 1,4,7 | sed "
					+ samDictionary.getSequences().stream().map(SSR->" -e 's/^"+SSR.getSequenceName()+"\t/"+SSR.getSequenceIndex()+"\t/' ").collect(Collectors.joining())
					+ "| sort -t '\t' -k2,2n > $(addsuffix .tmp.tsv,$@)");
			mkw.println("\techo 'set terminal png truecolor size ${SCREEN_WIDTH},${SCREEN_HEIGHT};"
					+ "set datafile separator \"\t\";"
					+ "set title \""+ sampleName+" (Whole Genome)\";"
					+ "set ylabel \"Normalized Depth - Number of copies\";"
					+ "set xlabel \"Genomic Position\";"
					+ "set yrange [0:"+Math.min(4,Math.ceil(maxDepth))+"];"
					+ "set xrange [1:"+ this.samDictionary.getReferenceLength()+".0];"
					+ "set xtic rotate by 90 right;"
					+ "set nokey;"
					+ "set output \"$@\";"
					+ "plot \"$(addsuffix .tmp.tsv,$@)\" using 2:3:1 w points linecolor variable "
					+ "' | ${GNUPLOT}");
			mkw.println("\trm -f $(addsuffix .tmp.tsv,$@)");
			
			mkw.println("all2: ${TARGETS}");
			mkw.flush();
			mkw.close();
			mkw=null;
			
			archive.close();
			archive=null;
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
			CloserUtil.close(archive);
			CloserUtil.close(this.indexedFastaSequenceFile);
			}	
		}
	
	public static void main(String[] args) {
		new CopyNumber01().instanceMainWithExit(args);
		}
	}
