package com.github.lindenb.jvarkit.tools.minicov;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.Closeable;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.function.ToDoubleFunction;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPInputStream;


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.canvas.Canvas;
import com.github.lindenb.jvarkit.canvas.CanvasFactory;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.TabixReader;

public class MiniCovPlotter extends Launcher {
	@Parameter(names = { "-o", "--out" }, description = "Ouput prefix",required = true)
	private String outputPrefix = null;
	@Parameter(names = { "-r", "--region","--interval" }, description ="CHR:START-END interval or use default.")
	private String intervalStr = null;
	@Parameter(names = { "-t", "--threads" }, description ="number of threads")
	private int nThreads=1;
	@Parameter(names = { "--from" }, description ="start from 0-based index in the data. inclusive")
	private int from_index=0;
	@Parameter(names = { "--to" }, description ="end at 0-based index in the data. inclusive. -1: any")
	private int to_index=-1;
	@Parameter(names = { "--batch-size","-B" }, description ="group coverage into batch of 'x' samples. -1 all batch")
	private int batch_size=-1;
	@Parameter(names={"--dimension","--dim"},description = "Image Dimension. " + DimensionConverter.OPT_DESC, converter=DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,300);
	@Parameter(names={"--format","-f"},description = CanvasFactory.OPT_FORMAT_DESC)
	private CanvasFactory.Format outputFormat =CanvasFactory.Format.PNG;
	@Parameter(names= {"--gtf"},description="indexed GTF file to show genes")
	private String gtfPath=null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	public MiniCovPlotter() {
		}

	private static class Record {
		String sample;
		Locatable loc;
		int[] array;
		}

	private static class RecordIterator extends AbstractIterator<Record> implements CloseableIterator<Record> {
		final BinaryCodec codec;
		boolean skip_array=false;
		volatile boolean cancel_flag=false;
		RecordIterator(final Path path) throws IOException {
			this.codec=new BinaryCodec(new GZIPInputStream(Files.newInputStream(path)));
			}
		@Override
		protected Record advance() {
			try {
				if(codec.knownAtEof() || this.cancel_flag) return null;
				final Record record=new Record();
				
				try {
					record.sample = codec.readNullTerminatedString();
					System.err.println(record.sample);
					}
				catch(RuntimeEOFException err) {
					return null;
					}
				final String chrom=codec.readNullTerminatedString();
				int start = codec.readInt();
				int end = codec.readInt();
				record.loc  = new SimpleInterval(chrom, start, end);
				System.err.println(record.loc);
				final int distance=end-start+1;
				int sizeof = (int)codec.readByte();
				if(sizeof!=1) throw new IllegalArgumentException(""+sizeof);
				if(skip_array) {
					record.array=null;
					codec.getInputStream().skip(distance);
					}
				else
					{
					record.array = new int[distance];
					for(int i=0;i< distance;i++) {
						record.array[i]= (int)codec.readUByte();
						}
					}
				if(codec.readByte()!=0) throw new IOException("expected 0");
				if(codec.readByte()!=0) throw new IOException("expected 0");
				return record;
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() {
			this.codec.close();
			}
		}
	
	private static class PosDepth implements Comparable<PosDepth> {
		int x;
		float y;
		int sample_idx;
		@Override
		public int compareTo(PosDepth o) {
			return Integer.compare(this.x, o.x);
			}
		}
	
	private static class PosDepthCodec extends AbstractDataCodec<PosDepth> {
		@Override
		public PosDepth decode(DataInputStream dis) throws IOException {
			PosDepth pd=new PosDepth();
			try {
				pd.x=dis.readInt();
				}
			catch(EOFException err) {
				return null;
				}
			pd.y=dis.readFloat();
			pd.sample_idx=dis.readInt();
			return pd;
			}
		@Override
		public void encode(DataOutputStream dos, PosDepth o) throws IOException {
			dos.writeInt(o.x);
			dos.writeFloat(o.y);
			dos.writeInt(o.sample_idx);
			}
		@Override
		public AbstractDataCodec<PosDepth> clone() {
			return new PosDepthCodec();
			}
		}
	
	private static class BoxPlot {
		int x;
		double mean;
		double median;
		double upperQuartile;
		double lowerQuartile;
		double upperWhisker;
		double lowerWhisker;
		List<Double> outliers;
		
		double calculateQuantile(final double[] values, final double p) {

	        double quantile = 0;

	        /* here are small symbols better to read than long ones */
	        final int s = values.length;

	        if (s > 1) {

	            /* Epsilon */
	            final double e = 0.000001;
	            final int k = (int) (p * (s));

	            if ((((s * p) % 1) < e) && (((s * p) % 1) > (-e))) {

	                /* s * p is element of N */
	                quantile = (((double) values[k - 1]) + ((double) values[k])) / 2;
	            } else {

	                /* s * p is not an element of N */
	                quantile = values[k];
	            }
	        } else if (s == 1) {
	            quantile = values[0];
	        }
	        return quantile;
	    }
		
		

	    private double calculateUpperWhisker(final double[] values) {

	        final double interQuartileRange = Math.abs(upperQuartile - lowerQuartile);

	        /* searching for the upperWhisker, starting at the upperQuartile and searching upwards */
	        double upperWhiskerCandidate = upperQuartile;
	        double upperWhisker = upperQuartile;

	        int i = (int) (values.length * 0.75);

	        while ((upperWhiskerCandidate <= (upperQuartile + (1.5d * interQuartileRange))) && (i <= values.length)) {
	            upperWhisker = upperWhiskerCandidate;
	            ++i;
	            if (i < values.length) {
	                upperWhiskerCandidate = values[i];
	            }
	        }
	        return upperWhisker;
	    }

		
	    private double calculateLowerWhisker(final double[] values) {
	        double lowerWhisker = 0d;

	        if (values.length > 0) {
	            final double interQuartileRange = Math.abs(upperQuartile - lowerQuartile);

	            /* searching for the lowerWhisker, starting at the lowerQuartile and searching downwards */
	            double lowerWhiskerCandidate = lowerQuartile;
	            lowerWhisker = lowerQuartile;

	            int i = (int) (values.length * 0.25);

	            while ((lowerWhiskerCandidate >= (lowerQuartile - (1.5d * interQuartileRange))) && (i >= 0)) {
	                lowerWhisker = lowerWhiskerCandidate;
	                --i;
	                if (i >= 0) {
	                    lowerWhiskerCandidate = values[i];
	                }
	            }
	        }
	        return lowerWhisker;
	    }
	    
	    private List<Double> calculateOutlier(final double[] values) {

	        final List<Double> outlier = new ArrayList<Double>();

	        if (values.length > 1) {

	            /* calculate upper outlier */
	            final int lastElement = values.length - 1;

	            if (upperWhisker < values[lastElement]) {

	                /* there are outlier */
	                int i = lastElement;

	                while (upperWhisker < values[i]) {
	                    outlier.add(values[i]);
	                    --i;
	                }
	            }

	            /* calculate lower outlier */
	            final double lowerWhisker = this.calculateLowerWhisker(values);
	            final int firstElement = 0;

	            if (lowerWhisker > values[firstElement]) {

	                /* there are outlier */
	                int i = firstElement;

	                while (lowerWhisker > values[i]) {
	                    outlier.add(values[i]);
	                    ++i;
	                }
	            }
	        }
	        return outlier;
	    }

	    
		void addAll(List<PosDepth> pds) {
			double[] values = pds.stream().mapToDouble(X->X.y).toArray();
			Arrays.sort(values);
			this.mean = Arrays.stream(values).average().orElse(0.0);
	        this.median = this.calculateQuantile(values,0.5);
	        this.upperQuartile =  this.calculateQuantile(values,0.75);
	        this.lowerQuartile = this.calculateQuantile(values,0.25);
	        this.upperWhisker = this.calculateUpperWhisker(values);
	        this.lowerWhisker = this.calculateLowerWhisker(values);
	        this.outliers = this.calculateOutlier(values);

		}
	}
	
	private  class Batch implements Closeable {
		int painted_count=0;
		final Canvas canvas;
		final Locatable userLoc;
		final double max_value=2.0;
		final int index_start;
		SortingCollection<PosDepth> sorter;
		Batch(final Locatable userLoc,int index) throws IOException {
			this.userLoc=userLoc;
			this.index_start=index;
			
			final String path= MiniCovPlotter.this.outputPrefix+ userLoc.getContig()+"_"+userLoc.getStart()+"_"+userLoc.getEnd()+"."+this.index_start+"."+ MiniCovPlotter.this.outputFormat.getSuffix();
			Path p = Paths.get(path);
			CanvasFactory cf = new CanvasFactory().setDimension(dimension.width+1, dimension.height+1);
			canvas = cf.open(p, FunctionalMap.make());
	
			
			
			//draw frame
			canvas.rect(0, 0,dimension.width+1,dimension.height+1,
					FunctionalMap.of(
						Canvas.KEY_FILL,Color.WHITE,
						Canvas.KEY_STROKE,Color.DARK_GRAY
					));
			
			this.sorter = SortingCollection.newInstance(
					PosDepth.class,
					new PosDepthCodec(),
					(A,B)->Integer.compare(A.x,B.x),
					writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths()
					);
			this.sorter.setDestructiveIteration(true);
			}
		
		
		protected void plotGenes() {
			if(MiniCovPlotter.this.gtfPath==null) return;
			final float midy= dimension.height-10;
			final GTFCodec codec=new GTFCodec();
			ToDoubleFunction<Integer> position2pixel= POS->{
				double p=(POS-this.userLoc.getStart())/(double)this.userLoc.getLengthOnReference();
				p=p*dimension.width;
				if(p<0) p=0;
				if(p>dimension.width) p=dimension.width;
				return p;
				};
			
			try(TabixReader tabix = new TabixReader(MiniCovPlotter.this.gtfPath)) {
				for(int step=0;step< 3;++step) {
					String expect;
					float featHeight;
					Color color;
					switch(step) {
						case 0: expect="gene";featHeight=1f;color=Color.MAGENTA;break;
						case 1: expect="exon";featHeight=4f;color=Color.MAGENTA;break;
						case 2: expect="CDS";featHeight=6f;color=Color.MAGENTA;break;
						default: expect=null;featHeight=0f;color=Color.MAGENTA;break;
						}
			
					final List<Rectangle2D> rects = new ArrayList<>();
					TabixReader.Iterator iter= tabix.query(
							this.userLoc.getContig(),
							this.userLoc.getStart(),
							this.userLoc.getEnd());
					for(;;) {
						final String line = iter.next();
						if(line==null) break;
						final GTFLine feat=codec.decode(line);
						if(feat==null) continue;
						if(!feat.getType().equals(expect)) continue;
						final double x1= position2pixel.applyAsDouble(feat.getStart());
						final double x2= position2pixel.applyAsDouble(feat.getEnd()+1);
						if(x1>=x2) continue;
						if(step>0 /* not a gene */ && (x2-x1)<2 /* too small to be displayed */) continue;
						rects.add(new Rectangle2D.Float((float)x1, midy-featHeight/2f, (float)(x2-x1), featHeight));
						}
					Collections.sort(rects,(A,B)->Double.compare(A.getX(),B.getX()));
					int k=0;
					while(k+1< rects.size()) {
						final Rectangle2D r1 = rects.get(k+0);
						final Rectangle2D r2 = rects.get(k+1);
						if(r2.getX() - r1.getMaxX() < 0.1) {
							final float x1 = (float)Math.min(r1.getX(), r2.getX());
							final float x2 = (float)Math.max(r1.getMaxX(), r2.getMaxX());
							rects.set(k,new Rectangle2D.Float(
								x1,
								(float)r1.getY(),
								(x2-x1),
								(float)r1.getHeight()
								));
							rects.remove(k+1);
							}
						else
							{
							k++;
							}
						}
					for(Rectangle2D rect:rects) {
						// draw feature
						canvas.shape(
							rect,
							FunctionalMap.of(
									Canvas.KEY_FILL, color,
									Canvas.KEY_STROKE,null
									)
							);
						}
					}
				}
			catch(Exception err) {
				err.printStackTrace();
				}
			}

		
		private List<Point2D> simplifyPolygon(final List<Point2D> points) {
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
			return points;
			}

	
		
		void add(final Record rec) {
			if(!this.userLoc.overlaps(rec.loc)) return;
			final Locatable loc = LocatableUtils.sharedInterval(this.userLoc, rec.loc);
			
			
			
			
			final DiscreteMedian<Integer> discreteMedian=new DiscreteMedian<>();
			Arrays.stream(rec.array).forEach(V->discreteMedian.add(V));
			final double median = Math.max(0, discreteMedian.getMedian().orElse(1.0));
			final int[] coverage = Arrays.copyOfRange(
					rec.array,
					loc.getStart()-rec.loc.getStart(),
					1+loc.getEnd()-rec.loc.getStart()
					);
			
			double[] y_array=new double[dimension.width];
			for(int i=0;i< y_array.length;++i) {
				int idx1 = (int)Math.floor(((i+0)/(double)y_array.length)*coverage.length);
				int idx2 = (int)Math.ceil(((i+1)/(double)y_array.length)*coverage.length);
				if(idx1==idx2) System.err.println("boum");
				y_array[i] = Arrays.stream(coverage,idx1, idx2).average().orElse(0);
				}
			
			final List<Point2D> points = new ArrayList<>(dimension.width+4);
			points.add(new Point2D.Double(0, dimension.height));
			for(int x=0;x < y_array.length;++x) {
				double v = y_array[x];
				if(median>0) v/=median;
				PosDepth pd=new PosDepth();
				pd.x=x;
				pd.y=(float)(dimension.height - ((float)v/this.max_value)*dimension.height);
				
				sorter.add(pd);
				
				
				points.add(new Point2D.Double(x,Math.min(dimension.height,pd.y)));
				}
			points.add(new Point2D.Double(y_array.length-1, dimension.height));
			
			float f_color= MiniCovPlotter.this.batch_size<1 ?(painted_count%255)/255f:painted_count/(float)MiniCovPlotter.this.batch_size;
			final Color c= Color.getHSBColor(f_color, 0.85f, 1.0f);
			float f_opacity =MiniCovPlotter.this.batch_size<1?0.1f:1f/MiniCovPlotter.this.batch_size;
			f_opacity=1f;
			//draw frame
			canvas.polygon(simplifyPolygon(points),
					FunctionalMap.of(
						Canvas.KEY_FILL,null,
						Canvas.KEY_STROKE,c,
						Canvas.KEY_STROKE_OPACITY,f_opacity,
						Canvas.KEY_STROKE_WIDTH,0.5
					));
			
			painted_count++;
			}
		
		@Override
		public void close() throws IOException {
			List<BoxPlot> boxPlots=new ArrayList<>(dimension.width);
			for(int i=0;i< dimension.width;i++) {
				BoxPlot bp=new  BoxPlot();
				bp.x=i;
				boxPlots.add(bp);
			}
			sorter.doneAdding();
			try(CloseableIterator<PosDepth> iter0  = sorter.iterator()) {
				
				try(EqualRangeIterator<PosDepth> iter1 = new EqualRangeIterator<>(iter0,(A,B)->Integer.compare(A.x, B.x))) {
					while(iter1.hasNext()) {
						List<PosDepth> row = iter1.next();
						PosDepth pd=row.get(0);
						BoxPlot boxplot= boxPlots.get(pd.x);
						boxplot.addAll(row);
						}
					}
				}
			
			
			
			
			
			
			
			
			for(int x=0;x<boxPlots.size();x++) {
				BoxPlot bp=boxPlots.get(x);
				for(double y:bp.outliers) {
					canvas.circle(
							x,y,2,
							FunctionalMap.of(
								Canvas.KEY_FILL,Color.RED,
								Canvas.KEY_STROKE,Color.BLACK
							));
					}
			}
			
			
			int fontSize=12;
			String title=this.userLoc.toString();
			canvas.text(
					title,
					1,
					fontSize+3,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_STROKE, null,
						Canvas.KEY_FILL, Color.DARK_GRAY,
						Canvas.KEY_TITLE,title
						)
					); 

			
			plotGenes();
			// draw horizonal lines
			for(double v=0.5;v<=1.5 ;v+=0.5) {
				double y = dimension.height - (v/max_value)*dimension.height;
				canvas.line(
						0,
						y,
						dimension.width,
						y,
						FunctionalMap.of(
							Canvas.KEY_STROKE,(v==1.0?Color.BLUE:Color.ORANGE),
							Canvas.KEY_TITLE, String.valueOf(v),
							Canvas.KEY_STROKE_WIDTH,0.5
							)
						);
				}

			canvas.close();
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			final Locatable userLoc;
			if(StringUtils.isBlank(this.intervalStr)) {
				userLoc =null;
				}
			else
				{
				userLoc = LocatableUtils.parse(this.intervalStr);
				}
			
			Batch batch=null;
			Path input= Paths.get(oneAndOnlyOneFile(args));
			IOUtil.assertFileIsReadable(input);
			int idx=0;
			try(RecordIterator iter=new RecordIterator(input)) {
				while(iter.hasNext()) {
					final Record record= iter.next();
					if(idx >= this.from_index) {
						if(batch!=null && this.batch_size>0 && batch.painted_count>=this.batch_size) {
							batch.close();
							batch=null;
							}
						if(batch==null) batch=new Batch(userLoc==null?record.loc:userLoc,idx);
						batch.add(record);
						}
					idx++;
					if(this.to_index>=0 && idx> this.to_index ) break;
					}
				}
			if(batch!=null) {
				batch.close();
				batch=null;
				}
			return 0;
			}
		catch(Throwable err) {
			err.printStackTrace();
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new MiniCovPlotter().instanceMainWithExit(args);
	}

}
