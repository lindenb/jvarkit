/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.coveragegrid;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.canvas.Canvas;
import com.github.lindenb.jvarkit.canvas.CanvasFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.Average;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.math.MinMaxInteger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
/**
BEGIN_DOC
## input

input a set of bam/cram files or one file with the suffix '.list' containing the path to the bams

## output

output is a HTML+SVG file

## example:

```
find dir -type f -name "*bam" > in.list 
java -jar dist/jvarkit.jar coveragegrid -R src/test/resources/rotavirus_rf.fa --region "RF01:100-200" in.list
```


END_DOC 
 */
@Program(
	name="coveragegrid",
	description="Display an image of depth to display any anomaly an intervals+bams as a grid image",
	keywords={"cnv","bam","depth","coverage","svg","postscript"},
	creationDate="20241009",
	modificationDate="20241009",
	jvarkit_amalgamion =  true,
	menu="CNV/SV"
	)
public class CoverageGrid extends Launcher {
	private enum PlotType {COVERAGE,MEDIAN_COVERAGE,PILEUP,PILEUP_PAIR,GRID};
	private static final Logger LOG = Logger.build( CoverageGrid.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refPath = null;
	@Parameter(names={"--regions","--region","--interval"},description = "Interval region",required=true)
	private String intervalStr=null;
	@Parameter(names={"--mapq"},description = "min mapping quality")
	private int min_mapq=1;
	@Parameter(names={"--dimension","--dim"},description = "Image Dimension. " + DimensionConverter.OPT_DESC, converter=DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,300);
	@Parameter(names={"--extend","-x"},description = "extends the interval x times")
	private double extendFactor= 3.0;
	@Parameter(names= {"--max-y"},description="Max normalized Y, when using --type MEDIAN_COVERAGE")
	private double max_normalized_y = 3.0;
	@Parameter(names={"--threads"},description="number of threads")
	private int nThreads = 1;
	@Parameter(names= {"--format"},description=CanvasFactory.OPT_FORMAT_DESC)
	private CanvasFactory.Format canvas_format=CanvasFactory.Format.SVG;
	@Parameter(names= {"--title"},description="Set Title for graphics")
	private String title="";
	@Parameter(names= {"--type"},description="Plot type. "
			+ "COVERAGE: DEPTH coverage, "
			+ "MEDIAN_COVERAGE: coverage normalize by external boundaries, "
			+ "PILEUP: show reads, "
			+ "PILEUP_PAIR: show paired fragments, "
			+ "GRID: matrix of read name co-occurence. Slow and memory consuming."
			)
	private PlotType plotType = PlotType.MEDIAN_COVERAGE;
	@DynamicParameter(names = "-D", description = "other parameters '-Dkey=value'. "+
			"-Dhide_insertions=true "+
			"-Dinversion_mode=true "+
			"-Doverlap_boundaries=true "
			)
	private Map<String, String> dynaParams = new HashMap<>();

	
	private final Color EXTERNAL_AREA_COLOR =new Color(230,230,230);

	
	private abstract class BamTask implements Runnable {
		String sampleName="";
		Path bamPath;
		Path referencePath;
		Locatable userInterval;
		Locatable extendedInterval;
		Rectangle boundaries;
		int status = -1;
		final boolean inversion_mode;
		final boolean overlap_boundaries;
		
		protected BamTask() {
			inversion_mode= CoverageGrid.this.dynaParams.getOrDefault("inversion_mode", "false").equals("true");
			overlap_boundaries= CoverageGrid.this.dynaParams.getOrDefault("overlap_boundaries", "false").equals("true");
			}
		
		protected SamReader openSamReader() throws IOException {
			final SamReaderFactory samReaderFactory = SamReaderFactory.
					makeDefault().
					disable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES).
					referenceSequence(this.referencePath).
					validationStringency(ValidationStringency.LENIENT);
			final SamReader sr= samReaderFactory.open(this.bamPath);
			if(!sr.hasIndex()) {
				this.status=-1;
				throw new IOException("index missng for "+bamPath);
				}
			final SAMFileHeader hdr = sr.getFileHeader();
			this.sampleName = hdr.getReadGroups().
					stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtils.isBlank(S)).
					findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(this.bamPath))
					;
			return sr;
			}
		protected double trimX(double coordX) {
			return Math.min(this.boundaries.getMaxX(), Math.max(this.boundaries.getX(),coordX));
			}
		
		protected double position2pixel(int pos) {
			return this.boundaries.getX() +
					((pos-this.extendedInterval.getStart())/(double)this.extendedInterval.getLengthOnReference())*this.boundaries.getWidth();
			}
		
		protected CloseableIterator<SAMRecord> query(final SamReader sr) {
			final int extend = 100;//give a chance to get clipping
			return sr.query(
					this.extendedInterval.getContig(),
					Math.max(1,this.extendedInterval.getStart() - extend),
					this.extendedInterval.getEnd() + extend, 
					false
					);
			}
		
		protected boolean acceptRead(final SAMRecord rec) {
			if(!SAMRecordDefaultFilter.accept(rec, CoverageGrid.this.min_mapq)) return false;
			if(inversion_mode) {
				boolean clip = rec.getCigar().isClipped();
				boolean bad_pair = rec.getReadPairedFlag() && !rec.getProperPairFlag();
				if(!(clip || bad_pair)) return false;
				}
			if(overlap_boundaries) {
				int start = rec.getUnclippedStart();
				int end = rec.getUnclippedEnd();
				if(rec.getReadPairedFlag() && 
				  !rec.getMateUnmappedFlag() && 
				  rec.getReferenceName().equals(rec.getMateReferenceName())) {
					start = Math.min(rec.getMateAlignmentStart(), start);
					end = Math.max(getMateEnd(rec), end);
					}
				if(!(
						(start <= userInterval.getStart() && userInterval.getStart()<=end ) ||
						(start <= userInterval.getEnd() && userInterval.getEnd()<=end ) 
					))
					{
					return false;
					}
				}
			return true;
			}
		
		protected void plotVerticalGuides(Canvas canvas) {
			// draw original region vertical bounds 
			if(!this.userInterval.equals(this.extendedInterval)) {
				for( int i=0;i<2;i++) {
					final int p = (i==0?userInterval.getStart():userInterval.getEnd());
					final double x = position2pixel(p);
					canvas.line(
							x,
							boundaries.getY(),
							x,
							this.boundaries.getMaxY(),
							FunctionalMap.of(
								Canvas.KEY_STROKE,Color.GREEN,
								Canvas.KEY_TITLE, StringUtils.niceInt(p),
								Canvas.KEY_STROKE_WIDTH,0.5
								)
							);
					}
				}
			}
		
		protected void plotFrame(Canvas canvas) {
			// draw frame
			canvas.rect(
				boundaries.getX(),
				boundaries.getY(),
				boundaries.getWidth(),
				boundaries.getHeight(),
				FunctionalMap.of(
						Canvas.KEY_FILL,null,
						Canvas.KEY_STROKE,Color.DARK_GRAY,
						Canvas.KEY_STROKE_WIDTH,0.5
						)
				);
			}
		protected void plotBackgroundFrame(Canvas canvas) {
			
			for(int i=0;i< 2;i++) {
				double x1 = trimX(position2pixel(i==0?this.extendedInterval.getStart():this.userInterval.getEnd()));
				double x2 = trimX(position2pixel(i==0?this.userInterval.getStart():this.extendedInterval.getEnd()));
				if(x1>=x2) continue;
				// draw frame
				canvas.rect(
					x1,
					boundaries.getY(),
					(x2-x1),
					boundaries.getHeight(),
					FunctionalMap.of(
						Canvas.KEY_FILL,EXTERNAL_AREA_COLOR,
						Canvas.KEY_STROKE,null
						)
					);
				}
			}
		
		abstract void plot(Canvas canvas);
		}	

	
	private class PileupTask extends BamTask {
		final boolean by_pair;
		PileupTask(boolean by_pair) {
			this.by_pair  = by_pair;
		}
		
		/** simplified read */
		private  class MiniRead implements Locatable {
			final int start;
			final int end;
			final Cigar cigar;
			final int flags;
			MiniRead(int start,int end,Cigar cigar,int flags) {
				this.start= start;
				this.end = end;
				this.cigar= cigar;
				this.flags = flags;
				}
			MiniRead(SAMRecord rec) {
				this(rec.getAlignmentStart(),rec.getAlignmentEnd(),rec.getCigar(),rec.getFlags());
				}
			public String getContig() {
				return PileupTask.this.extendedInterval.getContig();
				}
			public int getStart() {
				return start;
				}
			public int getEnd() {
				return end;
				}
			public Cigar getCigar() {
				return cigar;
				}
			public int getUnclippedStart() {
				return SAMUtils.getUnclippedStart(getStart(), cigar);
				}
			public int getUnclippedEnd() {
				return SAMUtils.getUnclippedEnd(getEnd(), cigar);
				}
			private boolean testFlag(SAMFlag flg) {
				return (this.flags & flg.intValue()) != 0;
				}
			public boolean isReadPaired() {
		        return testFlag(SAMFlag.READ_PAIRED);
				}
			}
		
		private class MiniReadPair implements Locatable,Comparable<MiniReadPair>,Iterable<MiniRead> {
			private final MiniRead R1 ;
			private MiniRead R2 = null;
			
			public MiniReadPair(final MiniRead rec)
				{
				this.R1 = rec;
				}
			
			
			@Override
			public int compareTo(final MiniReadPair o) {
				int i= Integer.compare(this.getStart(),o.getStart());
				if(i!=0) return i;
				i= Integer.compare(this.getEnd(),o.getEnd());
				if(i!=0) return i;
				return 0;
				}
			
			@Override
			public String getContig() {
				return R1.getContig();
				}
			@Override
			public int getStart() {
				return R2==null?
					R1.getUnclippedStart():
					Math.min(R1.getUnclippedStart(), R2.getUnclippedStart())
					;
				}
			@Override
			public int getEnd() {
				return R2==null?
						R1.getUnclippedEnd():
						Math.max(R1.getUnclippedEnd(), R2.getUnclippedEnd())
						;
				}
			@Override
			public Iterator<MiniRead> iterator() {
				return (R2==null?Arrays.asList(R1):Arrays.asList(R1,R2)).iterator();
				}
			boolean isStrangelyPaired() {
				if(!R1.isReadPaired()) return false;
				if(!R1.testFlag(SAMFlag.PROPER_PAIR)) return true;
				if(R1.testFlag(SAMFlag.READ_REVERSE_STRAND)==R1.testFlag(SAMFlag.MATE_REVERSE_STRAND)) return true;
				return false;
				}
			}

		
		final List<List<MiniReadPair>> rows = new ArrayList<>();
		
		void extractPileup() {
			final Map<String,MiniReadPair> readName2pairs = new HashMap<>();
			try(SamReader sr=openSamReader()) {
				 try(CloseableIterator<SAMRecord> iter=this.query(sr)) {
					 while(iter.hasNext()) {
						 final SAMRecord rec0=iter.next();
						 if(!acceptRead(rec0)) continue;
						 final Cigar cigar = rec0.getCigar();
						 if(cigar==null || cigar.isEmpty()) continue;
						 
						 final MiniRead rec = new MiniRead(rec0);
						 final String pair_name = this.by_pair ? rec0.getReadName(): rec0.getPairedReadName();
						 MiniReadPair pair=readName2pairs.get(pair_name);
						 if(pair==null) {
							pair =new MiniReadPair(rec);
							if( this.by_pair &&
								rec0.getReadPairedFlag() &&
								!rec0.getMateUnmappedFlag() && 
								rec0.getReferenceName().equals(rec0.getMateReferenceName())) {
								final int mate_start = rec0.getMateAlignmentStart();
								final int mate_end = getMateEnd(rec0);
								final int mate_flags = SAMFlag.READ_PAIRED.intValue()+
										(rec0.getFirstOfPairFlag()?SAMFlag.SECOND_OF_PAIR.intValue():SAMFlag.FIRST_OF_PAIR.intValue())+
										(rec0.getReadNegativeStrandFlag()?0:SAMFlag.READ_REVERSE_STRAND.intValue()) +
										(rec0.getReadNegativeStrandFlag()?SAMFlag.MATE_REVERSE_STRAND.intValue():0)  +
										(rec0.getProperPairFlag()?SAMFlag.PROPER_PAIR.intValue():0) 
										;
								final MiniRead mate = new MiniRead(
										mate_start,
										mate_end,
										new Cigar(Arrays.asList(new CigarElement(1+mate_end-mate_start, CigarOperator.M))),
										mate_flags
										);
								pair.R2 = mate;
								}
							readName2pairs.put(pair_name,pair);
						 	}
						 else {
							if(!rec0.getReadPairedFlag()) LOG.error("reads not paired with same name "+pair_name);
							pair.R2 = rec;
						 	}
					 	 }
					 }//end iterator
				 
				 final Pileup<MiniReadPair> pileup = new Pileup<>((L,R)->position2pixel(L.getEnd()+1) +1  < position2pixel(R.getStart()));
				 for(MiniReadPair pair: readName2pairs.values().stream().sorted().collect(Collectors.toList())) {
					 if(!pair.overlaps(extendedInterval)) continue;
					 pileup.add(pair);
				 	}
				
				 this.rows.addAll(pileup.getRows());
				 this.status=0;
				}//end samreader
			catch(IOException err) {
				 this.status=-1;
				LOG.error(err);
				}
			}
		@Override
		void plot(Canvas canvas) {
			 if(this.status!=0) return;
		     final double featureHeight= Math.min(20,(this.boundaries.getHeight()/Math.max(1.0,(double)this.rows.size())));
			 final boolean hide_insertions=  CoverageGrid.this.dynaParams.getOrDefault("hide_insertions", "false").equals("true");

		     double y= boundaries.getMaxY()-featureHeight;
		     
		     
		     plotBackgroundFrame(canvas);
			 plotVerticalGuides(canvas);

		     for(final List<MiniReadPair> row:this.rows) {
		    	final double h2= Math.max(featureHeight*0.9,featureHeight-2);
		    	for(final MiniReadPair pair: row) {		    
		    		final List<Integer> insertions = new ArrayList<>();
		    		/* draw rec itself */
		    		final double midy=y+h2/2.0;
		    		canvas.line(
		    				trimX(position2pixel(pair.getStart())),
		    				midy,
		    				trimX(position2pixel(pair.getEnd())),
		    				midy,
		    				FunctionalMap.of(
	    						Canvas.KEY_STROKE, pair.isStrangelyPaired()? Color.RED:Color.GRAY,
	    						Canvas.KEY_STROKE_WIDTH, 0.1
	    						)
		    				);
		    		
			    	for(MiniRead rec:pair) {
			    		final int unclipped_start = rec.getUnclippedStart();
			    		int ref1 = unclipped_start;
			    		for(final CigarElement ce: rec.getCigar().getCigarElements()) {
			    			if(ref1> this.extendedInterval.getEnd()) break;
			    			final CigarOperator op=ce.getOperator();
			    			Shape shape = null;
			    			Color fill=null;
			    			Color stroke=Color.DARK_GRAY;
			    			switch(op) {
			    				case P: break;
			    				case M://through
			    				case X://through
			    				case EQ://through
			    				case S: //through
			    				case H: 
			    						final double x1= trimX(position2pixel(ref1));
			    						final double x2 =trimX(position2pixel(ref1+ce.getLength()+1));
			    						if(x1 < x2) {
				    						shape = new Rectangle2D.Double(
					    						x1,
					    						y,
					    						(x2-x1),
					    						h2
					    						);
				    						}
			    						
			    						ref1+=ce.getLength();
			    						switch(op) {
			    							case H: case S: fill=Color.YELLOW;break;
			    							case X: fill=Color.RED;break;
			    							case EQ: case M: fill=(pair.isStrangelyPaired()?Color.PINK:Color.LIGHT_GRAY);break;
			    							default:break;
			    							}
			    						break;
			    				case N://through
			    				case D: shape=null;fill=null;stroke=null;ref1+=ce.getLength();break;
			    				case I: shape=null;fill=null;stroke=null;if(!hide_insertions) insertions.add(ref1);break;
			    				default: throw new IllegalStateException(""+op);
			    				}
			    			if(ref1 < this.extendedInterval.getStart()) continue;
			    			
			    			if(shape!=null && fill!=null) {
			    				canvas.shape(shape,
		    						FunctionalMap.of(
	    	    						Canvas.KEY_FILL,fill,
	    	    						Canvas.KEY_STROKE,(h2>4?stroke:null),
	    	    						Canvas.KEY_STROKE_WIDTH, 1
	    	    						)
		    						);
			    				}
			    			} // end loop cigar
			    		
			    		
			    		
			    		
			    		for(int px:insertions) {
			    			canvas.line(
				    				position2pixel(px),
				    				y-0.5,
				    				position2pixel(px),
				    				y+h2+0.5,
				    				FunctionalMap.of(
			    						Canvas.KEY_STROKE, Color.RED,
			    						Canvas.KEY_STROKE_WIDTH, 0.5
			    						)
				    				);
			    			}
			    		}
		    		}
		    	y-=featureHeight;
		     	}
			final int fontSize=Math.min(12,(int)(this.boundaries.getHeight()/10.0));
		 	canvas.text(
					this.sampleName,
					this.boundaries.getX()+1,
					this.boundaries.getY()+fontSize+3,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_STROKE, null,
						Canvas.KEY_FILL, Color.DARK_GRAY,
						Canvas.KEY_TITLE,this.sampleName+" "+this.bamPath
						)
					); 
		 	plotFrame(canvas);
			}
		
		@Override
		public void run() {
			this.extractPileup();
		}
	}
	
	
	private class CoverageTask extends BamTask {
		double[] coverage = null;
		final Average mean_depth= new Average();
		final boolean normalize_median;
		CoverageTask(boolean normalize_median) {
			this.normalize_median=normalize_median;
		}
		
		
		void extractCoverage() {
			
			final float[] array = new float[this.extendedInterval.getLengthOnReference()];
			Arrays.fill(array, 0);
			try(SamReader sr= openSamReader()) {
				try(CloseableIterator<SAMRecord> iter= this.query(sr)) {
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						if(!acceptRead(rec)) continue;
						for(AlignmentBlock ab: rec.getAlignmentBlocks()) {
							for(int x=0;x< ab.getLength();++x) {
								final int refpos1 = ab.getReferenceStart()+x;
								final int array_index = refpos1 - this.extendedInterval.getStart();
								if(array_index<0 || array_index>=array.length) continue;
								array[array_index]++;
								}
							}
						}
					}
				this.coverage = new double[boundaries.width];
				
				if(normalize_median) {
					final DiscreteMedian<Integer> discreteMedian= new DiscreteMedian<>();
					for(int i=0;i< array.length;i++) {
						this.mean_depth.accept(array[i]);
						final int refpos1 = this.extendedInterval.getStart() + i;
						if(CoordMath.overlaps(
								refpos1, refpos1,
								this.userInterval.getStart(),this.userInterval.getEnd()
								)) {
							continue;
							}
						discreteMedian.add((int)array[i]);
						}					
					final double median = Math.max(1.0,discreteMedian.getMedian().orElse(0.0));
					for(int i=0;i< array.length;i++) {
						array[i]=(float)(array[i]/median);
						}
					}
				else
					{
					for(int i=0;i< array.length;i++) {
						this.mean_depth.accept(array[i]);
						}
					}
				for(int i=0;i< this.coverage.length;i++) {
					      int x1 = (int)(((i+0)/(double) this.coverage.length)*array.length);
					final int x2 = (int)(((i+1)/(double) this.coverage.length)*array.length);
					final Average avg=new Average();
					while(x1<x2 && x1 < array.length) {
						avg.accept(array[x1]);
						x1++;
						}
					this.coverage[i]=avg.getAsDouble();
					}
				this.status=0;
				}
			catch(Throwable err) {
				LOG.error(err);
				this.status=-1;
				this.coverage=null;
				}
			}
		@Override
		public void run() {
			this.extractCoverage();
			}
		@Override
		void plot(Canvas canvas) {
			if(status!=0) return;
			final double bottom_y = this.boundaries.getMaxY();
			final double max_y_value = 
					this.normalize_median ?
						CoverageGrid.this.max_normalized_y:
						Math.max(1, this.mean_depth.get().orElse(1.0)*2.0)
						;
			final ToDoubleFunction<Double> cov2y = COV-> bottom_y - (COV/max_y_value)*this.boundaries.getHeight();
			final int fontSize=Math.min(12,(int)(this.boundaries.getHeight()/10.0));

			canvas.comment(this.sampleName);
			plotBackgroundFrame(canvas);
			
			canvas.text(
					this.sampleName+" avgDP:"+(int)this.mean_depth.getAsDouble(),
					this.boundaries.getX()+1,
					this.boundaries.getY()+fontSize+3,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_STROKE, null,
						Canvas.KEY_FILL, Color.DARK_GRAY,
						Canvas.KEY_TITLE,this.sampleName+" "+this.bamPath
						)
					);
			// draw horizonal lines
			for(double v=0.5;v<=1.5 && this.normalize_median ;v+=0.5) {
				
				double y = cov2y.applyAsDouble(v);
				if(y< this.boundaries.getY()) continue;
				canvas.line(
						this.boundaries.getX(),
						y,
						this.boundaries.getMaxX(),
						y,
						FunctionalMap.of(
							Canvas.KEY_STROKE,(v==1.0?Color.BLUE:Color.ORANGE),
							Canvas.KEY_TITLE, String.valueOf(v),
							Canvas.KEY_STROKE_WIDTH,0.5
							)
						);
				}

			

			final List<Point2D> points = new ArrayList<>(2+this.coverage.length);
			points.add(new Point2D.Double(this.boundaries.getX(), bottom_y));
			for(int i=0;i< this.coverage.length;i++) {
				double y = Math.max(this.boundaries.getY(), cov2y.applyAsDouble(this.coverage[i]));
				points.add(new Point2D.Double(this.boundaries.getX()+i, y));
				}
			points.add(new Point2D.Double(this.boundaries.getMaxX(), bottom_y));
			int i=1;
			while(i+2 < points.size())
				{
				final Point2D p1 = points.get(i+0);
				final Point2D p2 = points.get(i+1);
				final Point2D p3 = points.get(i+2);
				if(p1.getY()==p2.getY() && p2.getY()==p3.getY() && p1.getX() < p3.getX()) {
					points.remove(i+1);
					}
				else
					{
					i++;
					}
				}
			canvas.polygon(
					points,
					FunctionalMap.of(
							Canvas.KEY_STROKE, null,
							Canvas.KEY_FILL,Color.DARK_GRAY
							))
					;

			plotVerticalGuides(canvas);
			plotFrame(canvas);
			
			}
		}

	
	/**************************************************************************
	 * 
	 * GridTask : Co-Occurence of read names in memory
	 *
	 */
	private static class ReadName {
		final String name;
		final byte side;
		ReadName(final SAMRecord rec) {
			this.name= rec.getReadName();
			if(!rec.getReadPairedFlag()) {
				side=0;
				}
			else if(rec.getFirstOfPairFlag())
				{
				side=1;
				}
			else if(rec.getSecondOfPairFlag())
				{
				side=2;
				}
			else
				{
				throw new IllegalArgumentException();
				}
			}
		@Override
		public int hashCode() {
			return Objects.hash(name,side);
			}
		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null || !(obj instanceof ReadName)) return false;
			final ReadName other = (ReadName) obj;
			return side == other.side && this.name.equals(other.name);
			}
		
		}
	private class GridTask extends BamTask {
		private int square_size=1;
		private int[] counts;
		
		
		private int px2pos(int m) {
			return this.extendedInterval.getStart()+(int)((m/(double)this.square_size)*this.extendedInterval.getLengthOnReference());
			}
		private double pos2px(int p) {
			return ((p-this.extendedInterval.getStart())/(double)this.extendedInterval.getLengthOnReference())*this.square_size;
			}
		private int getCount(final Collection<Set<ReadName>> colx,final Collection<Set<ReadName>> coly) {
			int n=0;
			Iterator<ReadName> ix = colx.stream().flatMap(L->L.stream()).iterator();
			while(ix.hasNext()) {
				final ReadName rdName = ix.next();
				if(coly.stream().flatMap(L->L.stream()).anyMatch(NY->NY!=rdName /* NOT equals() , same read in memory */ && NY.name.equals(rdName.name))) {
					n++;
					}
				}
			return n;
			}
		@Override
		public void run() {
			try {
				final IntervalTreeMap<Set<ReadName>> treemap=new IntervalTreeMap<>();
				try(SamReader sr= openSamReader()) {
					try(CloseableIterator<SAMRecord> iter= this.query(sr)) {
						while(iter.hasNext()) {
							final SAMRecord rec = iter.next();
							 if(!acceptRead(rec)) continue;
							 final Interval r= new Interval(rec.getContig(), rec.getUnclippedStart(), rec.getUnclippedEnd(),rec.getReadNegativeStrandFlag(),rec.getReadName());
							 if(!this.extendedInterval.overlaps(r)) continue;
							
							 Set<ReadName> L = treemap.get(r);
							 if(L==null) {
								L=new HashSet<>();
								treemap.put(r, L); 
							 	}
							 L.add(new ReadName(rec));
						     }
						}
					}
				this.square_size = Math.max(Math.min(this.boundaries.width, this.boundaries.height),1);
				this.counts=new int[this.square_size*this.square_size];
				for(int x=0;x<this.square_size;x++) {
					final Interval rx = new Interval(this.extendedInterval.getContig(),px2pos(x),px2pos(x+1));
					final Collection<Set<ReadName>> colx = treemap.getOverlapping(rx);
					if(colx.isEmpty()) continue;
					for(int y=0;y<this.square_size;y++) {
						final Interval ry = new Interval(this.extendedInterval.getContig(),px2pos(y),px2pos(y+1));
						final Collection<Set<ReadName>> coly = treemap.getOverlapping(ry);
						if(coly.isEmpty()) continue;
						this.counts[y*this.square_size+x] = getCount(colx,coly);
					}
				}
				this.status=0;
				}
			catch(Throwable err) {
				this.status=-1;
				LOG.error(err);
			}
			
		}

		@Override
		void plot(final Canvas canvas) {
			if(status!=0) return;
			final int margin_top = this.boundaries.y +  (this.boundaries.height - this.square_size)/2;
			final int margin_left =  this.boundaries.x +   (this.boundaries.width - this.square_size)/2;

			final int fontSize=Math.min(12,(int)(this.boundaries.getHeight()/10.0));

			canvas.comment(this.sampleName);
			
			canvas.text(
					this.sampleName,
					this.boundaries.getX()+1,
					this.boundaries.getY()+fontSize+3,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_STROKE, null,
						Canvas.KEY_FILL, Color.DARK_GRAY,
						Canvas.KEY_TITLE,this.sampleName+" "+this.bamPath
						)
					);
			// area frames
			for(int side=0;side<2;++side) {
				final int p1 = side==0?this.extendedInterval.getStart():this.userInterval.getEnd();
				final int p2 = side==0?this.userInterval.getStart():this.extendedInterval.getEnd();
				final double v1 = pos2px(p1);
				final double v2 = pos2px(p2);
				if(v1>=v2) continue;
				canvas.rect(
						margin_left+v1,
						margin_top, 
						(v2-v1),
						this.square_size,
						FunctionalMap.of(
							Canvas.KEY_STROKE,null,
							Canvas.KEY_FILL,EXTERNAL_AREA_COLOR
							)
						);
				canvas.rect(
						margin_left,
						margin_top+v1, 
						this.square_size,
						(v2-v1),
						FunctionalMap.of(
							Canvas.KEY_STROKE,null,
							Canvas.KEY_FILL,EXTERNAL_AREA_COLOR
							)
						);
				}
						
			final MinMaxInteger mM = MinMaxInteger.of(Arrays.stream(this.counts));
			if(!mM.isEmpty()) {
				for(int y=0;y<this.square_size;y++) {
					int x=0;
					while(x<this.square_size) {
						if(this.counts[y*this.square_size+x]==0) {
							++x;
							continue;
							}
						final int x0=x;
						while(x<this.square_size && 
							counts[y*this.square_size+x]==counts[y*this.square_size+(x0)] ) {
							++x;
							}
						
						final float f = (float) mM.normalize(this.counts[y*this.square_size+x]);
						final Color c =  Color.getHSBColor(f, 0.85f, 1.0f);
						canvas.rect(
								margin_left+x0,
								margin_top+(this.square_size-(1.0/square_size))-y, /* inverse y coordinate so origin is at bottom*/
								(x-x0),
								1,
								FunctionalMap.of(
									Canvas.KEY_FILL, c,
									Canvas.KEY_STROKE,null
									)
								);
						}	
					}
				}
			
			// frame
			canvas.rect(
					margin_left,
					margin_top, 
					this.square_size,
					this.square_size,
					FunctionalMap.of(
						Canvas.KEY_FILL, null,
						Canvas.KEY_STROKE,Color.DARK_GRAY
						)
					);
			}
		}
	/********************************************************************************************************************/
	
/** get mate end position or mate-start */
private static int getMateEnd(final SAMRecord rec0) {
	return rec0.hasAttribute(SAMTag.MC)?SAMUtils.getMateAlignmentEnd(rec0):rec0.getMateAlignmentStart();
	}
	
@Override
public int doWork(final List<String> args) {
	try
		{
		if(this.extendFactor < 1.5) {
			System.err.println("extendFactor should be >=1.5");
			return -1;
			}
		if(this.max_normalized_y < 2.0) {
			System.err.println("max Y should be >= 2.0");
			return -1;
			}
		
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refPath);
		
		
		
		
		final List<Path> inputBams =  IOUtils.unrollPaths(args);
		
		if(inputBams.isEmpty()) {
			LOG.error("input bam file missing.");
			return -1;
			}
		
		final SimpleInterval user_interval = IntervalParserFactory.newInstance(dict).
				make().
				apply(this.intervalStr).
				orElseThrow(()->new IllegalArgumentException("Cannot parse interval:\""+this.intervalStr+"\""));
		
		final int mid_position = user_interval.getStart() + user_interval.getLengthOnReference()/2;
		int new_len = (int)(user_interval.getLengthOnReference()* this.extendFactor);
		int new_start = Math.max(1, mid_position-new_len/2);
		int new_end = Math.min(dict.getSequence(user_interval.getContig()).getEnd(), mid_position+new_len/2);
		final SimpleInterval extended_interval = new SimpleInterval(user_interval.getContig(),new_start,new_end);

		final List<BamTask> tasks = new ArrayList<>(inputBams.size());
		
		int marginTop=50;

		final int ncols = (int)Math.max(1,Math.floor(Math.sqrt(inputBams.size())));
		final int nrows = (int)Math.max(1, Math.ceil(inputBams.size()/(double)ncols));
		final double wi  = this.dimension.getWidth()/ncols;
		final double hi  = (this.dimension.getHeight()-marginTop)/nrows;
		
		LOG.info("tile:"+ncols+"x"+nrows+" w:"+(int)wi+" h:"+(int)hi);
		
		if(wi < 10 || hi < 10) {
			LOG.error("dimension is too small");
			return -1;
			}
		
		int nr=0;
		int nc=0;
		
		for(final Path bamPath: inputBams) {
			final BamTask task ;
				switch(this.plotType) {
				case GRID: task= new GridTask(); break;
				case PILEUP: task= new PileupTask(false); break;
				case PILEUP_PAIR: task= new PileupTask(true); break;
				case COVERAGE:  task = new CoverageTask(false); break;
				default: task = new CoverageTask(true); break;
				}
			
			task.bamPath = bamPath;
			task.referencePath = this.refPath;
			task.userInterval = user_interval;
			task.extendedInterval = extended_interval;
			task.boundaries = new Rectangle(
					(int)(nc*wi), 
					(int)(marginTop+nr*hi),
					(int)wi-1,
					(int)hi-1
					);
			tasks.add(task);
			nc++;
			if(nc==ncols) {
				nc=0;
				nr++;
				}
			}
	    final ExecutorService executorService = Executors.newFixedThreadPool(this.nThreads);
		for(int i=0; i < tasks.size();i+=this.nThreads) {
			final List<BamTask> batch = tasks.subList(i, Math.min(i + this.nThreads, tasks.size()));
			batch.forEach(executorService::execute);
			}
		executorService.shutdown();
		executorService.awaitTermination(365, TimeUnit.DAYS);
		
		final String title= SequenceDictionaryUtils.getBuildName(dict).orElse("") + " "+ 
				user_interval.toNiceString()+
			" length:"+ StringUtils.niceInt(user_interval.getLengthOnReference())+
			" extended to "+extended_interval.toNiceString()+
			" length:"+StringUtils.niceInt(extended_interval.getLengthOnReference())+
			" model:"+this.plotType.name()+" "+
			this.title
			;
		try(Canvas canvas = new CanvasFactory().
				setWidth(this.dimension.width).
				setHeight(this.dimension.height).
				setFormat(this.canvas_format).
				open(this.outputFile,FunctionalMap.of(
					Canvas.KEY_TITLE,title,
					Canvas.KEY_FILL, Color.BLACK,
					Canvas.KEY_STROKE,null
					) )
				)
			{
			final int fontSize=Math.min(14,marginTop);
			final String url =Hyperlink.compile(dict).apply(extended_interval).orElse(null);
			
			
			canvas.text(
					title,
					10,
					fontSize,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_FILL, Color.DARK_GRAY,
						Canvas.KEY_HREF, url,
						Canvas.KEY_TITLE, title
						)
					);
		
			
			for(BamTask task: tasks) {
				canvas.begin(FunctionalMap.make());
				task.plot(canvas);
				canvas.end();
				}
			
			}
			
	
		return 0;
		}
	catch(final Throwable err)
		{
		err.printStackTrace();
		LOG.error(err);
		return -1;
		}

	}

public static void main(final String[] args) {
	new CoverageGrid().instanceMainWithExit(args);
	}

}
