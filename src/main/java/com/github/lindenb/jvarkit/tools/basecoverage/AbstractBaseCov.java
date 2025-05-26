/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.basecoverage;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;
import java.util.function.Function;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.canvas.Canvas;
import com.github.lindenb.jvarkit.canvas.CanvasFactory;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;

public abstract class AbstractBaseCov extends Launcher {
	protected static final String FORMAT_NORM_DEPH= "MD";
	
	@Parameter(names={"--mapq"},description=" min mapping quality.")
	private int mapping_quality=1;
	@Parameter(names={"--runmed"},description=" moving median size")
	private int moving_median_length = 31;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@ParametersDelegate
	protected WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	
	protected double[] getCoverage(final SamReader sr,final Locatable queryInterval) {
		System.gc();
		final double[] coverage_d = new double[queryInterval.getLengthOnReference()];
		Arrays.fill(coverage_d,0);
		try(CloseableIterator<SAMRecord> it= sr.queryOverlapping(queryInterval.getContig(), queryInterval.getStart(), queryInterval.getEnd())) {
			while(it.hasNext()) {
				final SAMRecord rec =it.next();
				if(!SAMRecordDefaultFilter.accept(rec, this.mapping_quality)) continue;
				for(AlignmentBlock ab : rec.getAlignmentBlocks()) {
					for(int x=0;x < ab.getLength();++x) {
						final int array_index = (ab.getReferenceStart()+x) - queryInterval.getStart();
						if(array_index<0 ) continue;
						if(array_index>=coverage_d.length) break;
						coverage_d[array_index]++;
						}
					}
				}
			}
		return applyMovingMedian(coverage_d);
		}
	
	private double[] applyMovingMedian(final double[] source) {
		if(moving_median_length<=0) return source;
		final Median med = new Median();
		int half_size = this.moving_median_length/2;
		double[] coverage_f = new double[source.length];	
		for(int i=0;i< source.length;++i) {
			if(i-half_size>=0 && i-half_size+ this.moving_median_length <= source.length ) {
				coverage_f[i]=  med.evaluate(source, i-half_size,  this.moving_median_length);
			} else {
				coverage_f[i]=  source[i];
			}
			}
		return coverage_f;
		}
	
	
	protected OptionalDouble getMedian(final double[] coverage_int) {
		if(coverage_int.length==0) return OptionalDouble.empty();
		final Median med = new Median();
		return OptionalDouble.of(med.evaluate(coverage_int));
		}
	
	protected double[] normalizeOnMedian(final double[] coverage_d,double median_cov) {
		double[] coverage_norm = new double[coverage_d.length];	
		for(int array_index=0;array_index< coverage_d.length;++array_index) {
			coverage_norm[array_index]= (coverage_d[array_index]/median_cov);
			}
		return coverage_norm;
		}
	
	
	protected class Plot implements Closeable {
		final Canvas canvas;
		private final double max_y=2.5;
		private final int array_size;
		Plot(Path p,int length_on_ref) throws IOException {
			final CanvasFactory cf = new CanvasFactory();
			cf.setDimension(1000, 300);
			canvas = cf.open(p,FunctionalMap.make());
			canvas.rect(0, 0, canvas.getWidth()-1, canvas.getHeight()-1,FunctionalMap.of(Canvas.KEY_STROKE,Color.BLACK,Canvas.KEY_FILL,Color.WHITE));
			this.array_size=length_on_ref;
			}
		double coordX(double x) {
			return (x/(double)array_size)*canvas.getWidth();
			}
		double coordY(double v) {
			return canvas.getHeight() - (v/max_y)*canvas.getHeight();
			}
		Point2D toPoint(int i,double v) {
			return new Point2D.Double(coordX(i),coordY(v));
			}
		void simplify(List<Point2D> points, Function<List<Double>, Double> fun) {
			int i=0;
			while(i< points.size()) {
				List<Double> ys=new ArrayList<>();
				ys.add(points.get(i).getY());
				while(i+1 < points.size() && (int)points.get(i).getX()==(int)points.get(i+1).getX()) {
					ys.add(points.get(i+1).getY());
					points.remove(i+1);
					}
				points.set(i, new Point2D.Double(points.get(i).getX(),fun.apply(ys)));
				i++;
				}
			}
		@Override
		public void close() throws IOException {
			canvas.line(0, coordY(0.5), array_size, coordY(0.5), FunctionalMap.of(Canvas.KEY_STROKE,Color.ORANGE));
			canvas.line(0, coordY(1.5), array_size, coordY(1.5), FunctionalMap.of(Canvas.KEY_STROKE,Color.ORANGE));
			canvas.line(0, coordY(1.0), array_size, coordY(1.0), FunctionalMap.of(Canvas.KEY_STROKE,Color.BLUE));
			canvas.rect(0, 0, canvas.getWidth()-1, canvas.getHeight()-1,FunctionalMap.of(Canvas.KEY_STROKE,Color.BLACK,Canvas.KEY_FILL,null));
			canvas.close();
			}
		}
	
	protected AbstractBaseCov() {
		}

	}
