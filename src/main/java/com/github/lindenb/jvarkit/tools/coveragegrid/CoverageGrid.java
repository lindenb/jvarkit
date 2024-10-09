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
import java.awt.geom.Point2D;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.ToDoubleFunction;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.canvas.Canvas;
import com.github.lindenb.jvarkit.canvas.CanvasFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.Average;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
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
	keywords={"cnv","bam","depth","coverage","svg","postscrip"},
	creationDate="20241009",
	modificationDate="20241009",
	jvarkit_amalgamion =  true,
	menu="CNV/SV"
	)
public class CoverageGrid extends Launcher {
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
	@Parameter(names= {"--max-y"},description="Max normalized Y")
	private double max_normalized_y = 3.0;
	@Parameter(names={"--threads"},description="number of threads")
	private int nThreads = 1;
	@Parameter(names= {"--format"},description=CanvasFactory.OPT_FORMAT_DESC)
	private CanvasFactory.Format canvas_format=CanvasFactory.Format.SVG;


	

	

	private class CoverageTask implements Runnable {
		String sampleName="";
		Path bamPath;
		Path referencePath;
		Locatable userInterval;
		Locatable extendedInterval;
		Rectangle boundaries;
		double[] coverage = null;
		int status = -1;
		Average mean_depth= new Average();
		
		void extractCoverage() {
			final SamReaderFactory samReaderFactory = SamReaderFactory.
					makeDefault().
					disable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES).
					referenceSequence(this.referencePath).
					validationStringency(ValidationStringency.LENIENT);
			final float[] array = new float[this.extendedInterval.getLengthOnReference()];
			Arrays.fill(array, 0);
			try(SamReader sr=samReaderFactory.open(this.bamPath)) {
				if(!sr.hasIndex()) {
					LOG.error("index missng for "+bamPath);
					this.status=-1;
					return;
					}
				final SAMFileHeader hdr = sr.getFileHeader();
				this.sampleName = hdr.getReadGroups().
						stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtils.isBlank(S)).
						findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(this.bamPath))
						;
				try(CloseableIterator<SAMRecord> iter=sr.query(this.extendedInterval.getContig(), this.extendedInterval.getStart(), this.extendedInterval.getEnd(), false)) {
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						if(!SAMRecordDefaultFilter.accept(rec, CoverageGrid.this.min_mapq)) continue;
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
				
				final DiscreteMedian<Integer> discreteMedian= new DiscreteMedian<>();
				for(int i=0;i< array.length;i++) {
					this.mean_depth.accept(array[i]);
					final int refpos1 = this.extendedInterval.getStart()+i;
					if(this.userInterval.getStart()>=refpos1 && refpos1<=this.userInterval.getEnd()) {
						continue;
						}
					discreteMedian.add((int)array[i]);
					}
				
				
				final double median = Math.max(1.0,discreteMedian.getMedian().orElse(0.0));
				for(int i=0;i< array.length;i++) {
					array[i]=(float)(array[i]/median);
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
		void plot(Canvas canvas) {
			if(status!=0) return;
			final double bottom_y = this.boundaries.getMaxY();
			ToDoubleFunction<Double> cov2y = COV-> bottom_y - (COV/max_normalized_y)*this.boundaries.getHeight();
			ToDoubleFunction<Integer> pos2pixel = POS-> this.boundaries.getX() + ((POS-this.extendedInterval.getStart())/(double)this.extendedInterval.getLengthOnReference())*this.boundaries.getWidth();
			final int fontSize=Math.min(12,(int)(this.boundaries.getHeight()/10.0));

			canvas.comment(this.sampleName);
			
			
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
			for(double v=0.5;v<=1.5;v+=0.5) {
				
				double y = cov2y.applyAsDouble(v);
				if(y< this.boundaries.getY()) continue;
				canvas.line(
						this.boundaries.getX(),
						y,
						this.boundaries.getMaxX(),
						y,
						FunctionalMap.of(
							Canvas.KEY_STROKE,(v==1.0?Color.BLUE:Color.ORANGE),
							Canvas.KEY_TITLE, String.valueOf(v)
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
				Point2D p1 = points.get(i+0);
				Point2D p2 = points.get(i+1);
				Point2D p3 = points.get(i+2);
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
							Canvas.KEY_STROKE, Color.ORANGE,
							Canvas.KEY_FILL,Color.DARK_GRAY
							))
					;

			// draw original region vertical bounds 
			if(!this.userInterval.equals(this.extendedInterval)) {
				for( i=0;i<2;i++) {
					final int p = (i==0?userInterval.getStart():userInterval.getEnd());
					final double x = pos2pixel.applyAsDouble(p);
					canvas.line(
							x,
							boundaries.getY(),
							x,
							bottom_y,
							FunctionalMap.of(
								Canvas.KEY_STROKE,Color.GREEN,
								Canvas.KEY_TITLE, StringUtils.niceInt(p)
								)
							);
					}
				}
			// draw frame
			canvas.rect(
					boundaries.getX(),
					boundaries.getY(),
					boundaries.getWidth(),
					boundaries.getHeight(),
					FunctionalMap.of(
							Canvas.KEY_FILL,null,
							Canvas.KEY_STROKE,Color.DARK_GRAY)
					);
			}
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
		
		int mid_position = user_interval.getStart() + user_interval.getLengthOnReference()/2;
		int new_len = (int)(user_interval.getLengthOnReference()* this.extendFactor);
		int new_start = Math.max(1, mid_position-new_len/2);
		int new_end = Math.min(dict.getSequence(user_interval.getContig()).getEnd(), mid_position+new_len/2);
		final SimpleInterval extended_interval = new SimpleInterval(user_interval.getContig(),new_start,new_end);

		final List<CoverageTask> tasks = new ArrayList<>(inputBams.size());
		
		int marginTop=50;

		final int ncols = (int)Math.max(1,Math.floor(Math.sqrt(inputBams.size())));
		final int nrows = (int)Math.max(1, Math.ceil(inputBams.size()/(double)ncols));
		final double wi  = this.dimension.getWidth()/ncols;
		final double hi  = (this.dimension.getHeight()-marginTop)/nrows;
		int nr=0;
		int nc=0;
		
		for(final Path bamPath: inputBams) {
			final CoverageTask task = new CoverageTask();
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
			final List<CoverageTask> batch = tasks.subList(i, Math.min(i + this.nThreads, tasks.size()));
			batch.forEach(executorService::execute);
			}
		executorService.shutdown();
		executorService.awaitTermination(365, TimeUnit.DAYS);
		
		final String title= SequenceDictionaryUtils.getBuildName(dict).orElse("") + " "+ 
		user_interval.toNiceString()+
			" length:"+ StringUtils.niceInt(user_interval.getLengthOnReference())+
			" extended to "+extended_interval.toNiceString()+
			" length:"+StringUtils.niceInt(extended_interval.getLengthOnReference());
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
			System.err.println(canvas.getClass());
			final int fontSize=Math.min(14,marginTop);

			canvas.text(
					title,
					10,
					fontSize,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_FILL, Color.DARK_GRAY
						)
					);
		
			
			for(CoverageTask task: tasks) {
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
