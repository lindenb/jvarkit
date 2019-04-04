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
*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.DoubleConsumer;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.regex.Pattern;
import java.util.stream.Collectors;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.chart.Axis;
import com.github.lindenb.jvarkit.chart.BarChart;
import com.github.lindenb.jvarkit.chart.BubbleChart;
import com.github.lindenb.jvarkit.chart.CategoryAxis;
import com.github.lindenb.jvarkit.chart.Chart;
import com.github.lindenb.jvarkit.chart.HeatMapChart;
import com.github.lindenb.jvarkit.chart.NumberAxis;
import com.github.lindenb.jvarkit.chart.PieChart;
import com.github.lindenb.jvarkit.chart.RExporter;
import com.github.lindenb.jvarkit.chart.ScatterChart;
import com.github.lindenb.jvarkit.chart.StackedBarChart;
import com.github.lindenb.jvarkit.chart.XYChart;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC 

## Plot types:

Year  X  Y
2018  1  2
2019  3  4


## Examples

```
$ gunzip -c in.vcf.gz | grep -v "^#" | cut -f 1 | sort | uniq -c | java -jar dist/simpleplot.jar -t SIMPLE_HISTOGRAM  -su 
$ gunzip -c in.vcf.gz | grep -v "^#" | cut -f 1 | sort | uniq -c | java -jar dist/simpleplot.jar -t PIE  -su
```

### HISTOGRAM

```
Year  X  Y
2018  1  2
2019  3  4
```

```
$ echo -e "Year\tX\tY\n2018\t1\t2\n2019\t3\t4" | java -jar dist/simpleplot.jar -t HISTOGRAM
```

produces a histogram with two series in the legend (2018 and 2019). On the X-axis: 4 items X-2018, X-2019, Y-2018, Y-2019



### STACKED_HISTOGRAM

```
Year  X  Y
2018  1  2
2019  3  4
```

```
$ echo -e "Year\tX\tY\n2018\t1\t2\n2019\t3\t4" | java -jar dist/simpleplot.jar -t STACKED_HISTOGRAM
``` 

produces a histogram with two series in the legend (2018 and 2019). On the X-axis: 2 items X( 2018 under-2019), Y (2018 under 2019)


### STACKED_HISTOGRAM_PIVOTED

```
Year  X  Y
2018  1  2
2019  3  4
```

```
$  echo -e "2018\t1\t2\n2019\t3\t4\n2020\t1\t10" | java -jar dist/simpleplot.jar -t STACKED_HISTOGRAM_PIVOTED -nh
``` 

produces a histogram with two series in the legend ($1 and $2). On the X-axis: 3 items 2018, 2019, 2020. Each with the stacked $1 and $2



### Example

plot base=function(position) in a fastq:

```
gunzip -c src/test/resources/S1.R1.fq.gz | \
	awk '(NR%4==2) {L=length($0);for(i=1;i<=L;i++) printf("%d\t%s\n",i,substr($0,i,1));}' |\
	sort | uniq -c |\
	java -jar dist/simpleplot.jar -su -t STACKED_XYV --xlabel "Position"
```


### Example

HeatMap

```
echo -e "A\tA\t1\nA\tB\t2\nB\tA\t3\nB\tB\t10" | java -jar dist/simpleplot.jar  -t HEATMAP
```


## History

  * 2019: removed jfx as openjdk doesn't support it... output is now R
 
 
END_DOC
 */
@Program(name="simpleplot",
description="simple figure plotter output is a R script",
keywords={"char","figure"})

public class SimplePlot extends Launcher {
	private static final Logger LOG = Logger.build(SimplePlot.class).make();
	
	
	
	private static class MinMax
		implements DoubleConsumer
		{
		Double min = null;
		Double max = null;
		@Override
		public void accept(double value) {
			if(min==null) {
				min=max=value;
				}
			else
				{
				min=Math.min(min, value);
				max=Math.max(max, value);
				}
			}
		boolean isEmpty() { return this.min==null;}
		double distance() { return max-min;}
		public NumberAxis toAxis() {
			return new NumberAxis(this.min, this.max, distance()/10.0);
		}
		}
	
	
	/** base generator for a Chart Maker */
	private abstract class ChartSupplier
		implements Supplier<Chart>
		{
		private SAMSequenceDictionary _dict = null;
		protected Pattern delimiter = Pattern.compile("[\t]");
		private ContigNameConverter _ctgConverter = null;
		
		protected <V> Map<String,V> createMapWithStringKey()
			{
			if(!input_is_sort_uniq)
				{
				return new TreeMap<>(new SmartComparator().caseSensitive());
				}
			else
				{
				return new LinkedHashMap<>();
				}
			}
		protected Object[] nextTriple() {
			for(;;)
				{
				final String line = lineSupplier.get();
				if(line==null) return null;
				final String x;
				final String y;
				final Double val;
				if(input_is_sort_uniq)
					{
					final Object array[] = parseSortUniq(line);
					if(array==null) continue;
					val = Double.class.cast(array[0]);
					String remain[] =this.delimiter.split(String.class.cast(array[1]));
					if(remain.length<2) continue;
					x=remain[0];
					y=remain[1];
					}
				else
					{
					final String tokens[] = this.delimiter.split(line);
					if(tokens.length<3) {
						LOG.warn("got "+ tokens.length +" tokens expected at least 3. Skeeping "+line);						
						continue;
						}
					x=tokens[0];
					y=tokens[1];
					val = parseDouble.apply(tokens[2]);
					}
				if(val==null || val.doubleValue()<0) continue;
				return new Object[] {x,y,val};
				}
			}
		
		protected <T> void updateAxisX(final Axis<T> axis) {
			if(!StringUtil.isBlank(xlabel))
			{
				axis.setLabel(xlabel);
			}
		}
		
		protected <T> void updateAxisY(final Axis<T> axis) {
			if(!StringUtil.isBlank(ylabel))
			{
				axis.setLabel(ylabel);
			}
		}
		
		/** 1st = count as Doube, 2nd=remains */
		protected Object[] parseSortUniq(final String line) {
			int i=0;
			//skip ws
			while(i<line.length() && Character.isWhitespace(line.charAt(i))) i++;
			
			//read count
			final StringBuilder sb = new StringBuilder();
			while(i<line.length() && Character.isDigit(line.charAt(i)))
				{
				sb.append(line.charAt(i));
				i++;
				}
			if(sb.length()==0) return null;
			Double val = this.parseDouble.apply(sb.toString());
			if(val==null || val.doubleValue()<0.0) return null;
			i++; //skip one whitespace
			if(i>=line.length()) return null;
			return new Object[] {val,line.substring(i)};
			}
		
		
		protected final Function<String,Double> parseDouble = (S)->{
			try {
				if(StringUtil.isBlank(S)) return null;
				if(S.equals(".") || S.equals("NA")) return null;
				final Double dbl = Double.parseDouble(S);
				return dbl;
			} catch(final NumberFormatException err) {
				return null;
			}
		};
		
		
		
		protected ContigNameConverter getContigNameConverter() {
			if(this._ctgConverter==null)  getDictionary();
			return this._ctgConverter;
		}
		
		protected SAMSequenceDictionary getDictionary()
			{
			if(this._dict==null) {
				this._dict = new SAMSequenceDictionary(
						SequenceDictionaryUtils.extractRequired(SimplePlot.this.faidx).
							getSequences().
							stream().
							filter(C->min_contig_size<0|| C.getSequenceLength()>= SimplePlot.this.min_contig_size).
							collect(Collectors.toList())
							);
				if(this._dict==null || this._dict.isEmpty()) {
					throw new JvarkitException.DictionaryMissing("empty dictionary extracted from "+SimplePlot.this.faidx);
					}
				this._ctgConverter =  ContigNameConverter.fromOneDictionary(this._dict);
				}
			return this._dict;
			}
		}
	
	private class BedGraphSupplier
		extends ChartSupplier
		{
		@Override
		public Chart get() {
				final SAMSequenceDictionary dict = this.getDictionary();
				final List<XYChart.Series<Number,Number>> chroms = dict.getSequences().
					stream().
					map(SSR->{
						final XYChart.Series<Number,Number> c= new XYChart.Series<>();
						c.setName(SSR.getSequenceName());
						return c;
						}).
					collect(Collectors.toList());
				final MinMax minmaxY = new MinMax();
				for(;;) {
					final String line = lineSupplier.get();
					if(line==null) break;
					if(StringUtil.isBlank(line) || BedLine.isBedHeader(line)  ) continue;
					final String tokens[] = super.delimiter.split(line);
					
					final String chrom = getContigNameConverter().apply(tokens[0]);
					if(chrom==null) continue;
					SAMSequenceRecord ssr = dict.getSequence(chrom);
					if(ssr==null) continue;
					
					final int start,end;
					final Double yVal;
					if(SimplePlot.this.input_is_chrom_position) {
						start = end = Integer.parseInt(tokens[1]);
						yVal =  this.parseDouble.apply(tokens[2]);
					} else
					{
						start = Integer.parseInt(tokens[1]);
						end = Integer.parseInt(tokens[3]);
						yVal =  this.parseDouble.apply(tokens[3]);
					}
					if(yVal==null || start>end || start<0 || end > ssr.getSequenceLength()) continue;
					int midX = (start + end)/2;
					
					
					long x = midX +
							dict.getSequences().
							subList(0, ssr.getSequenceIndex()).
							stream().
							mapToLong(C->C.getSequenceLength()).
							sum()
							;
										
					final XYChart.Data<Number,Number> data  = new XYChart.Data<>(x, yVal);
					chroms.get(ssr.getSequenceIndex()).getData().add(data);
					minmaxY.accept(yVal);
					
					}
				if(minmaxY.isEmpty() || chroms.isEmpty()) {
					LOG.error("no valid data");
					return null;
				}
		        final NumberAxis xAxis = new NumberAxis(1,dict.getReferenceLength(),dict.getReferenceLength()/10.0);
		        xAxis.setLabel("Genomic Index");
		        final NumberAxis yAxis = minmaxY.toAxis();
		       
			    final ScatterChart<Number,Number> sc = new ScatterChart<Number,Number>(xAxis,yAxis);
			    sc.getData().addAll(chroms);
			    
			    return sc;			
			    }
			}
	/** abstract supplier for Map<String,Double> */
	private abstract class AbstractKeyValueSupplier
		extends ChartSupplier
		{
		protected Map<String,Double> getKeyValueMap() {
			final Map<String,Double> key2count = super.createMapWithStringKey();
			for(;;) {
				String line = lineSupplier.get();
				if(line==null) break;
				if(StringUtil.isBlank(line) ||line.startsWith("#") ) continue;
				final String key;
				final Double val;
				if( input_is_sort_uniq) {
					final Object array[] = parseSortUniq(line);
					if(array==null) continue;
					val = Double.class.cast(array[0]);
					key=String.class.cast(array[1]);
					}
				else
					{
					final String tokens[]= super.delimiter.split(line,2);
					if(tokens.length<2 || key2count.containsKey(tokens[0])) continue;
					val = this.parseDouble.apply(tokens[1]);
					if(val==null || val.doubleValue()<0.0) continue;
					key=tokens[0];
					}
				
				key2count.put(key, val);
				}
			return key2count;
			}
		}
	
	
	/** supplier for pie char */
	private class PieChartSupplier extends AbstractKeyValueSupplier {
		@Override
		public Chart get() {
			final Map<String,Double> key2count = getKeyValueMap();			
			final List<PieChart.Data> pieChartData =
						key2count.entrySet().stream().
						map(KV->new PieChart.Data(KV.getKey(),KV.getValue())).
						collect(Collectors.toList())
						;
			final  PieChart chart = new PieChart(pieChartData);
			return chart;
			}
		}
	
	/** supplier for simple bar chart */
	private class SimpleHistogramSupplier extends AbstractKeyValueSupplier {
		@Override
		public Chart get() {
			final Map<String,Double> key2count = getKeyValueMap();			
			final List<StackedBarChart.Data<String,Number>> data =
				key2count.entrySet().stream().
				map(KV->new StackedBarChart.Data<String,Number>(KV.getKey(),KV.getValue())).
				collect(Collectors.toList())
				;
			final XYChart.Series<String, Number> series1 =
			            new XYChart.Series<String, Number>(data);
			final CategoryAxis xAxis = new CategoryAxis();
		    final NumberAxis yAxis = new NumberAxis();
			final StackedBarChart<String, Number> chart =
			            new StackedBarChart<String, Number>(xAxis, yAxis);
			chart.getData().add(series1);
			chart.setLegendVisible(false);
			return chart;
			}
		}
	
	private abstract class AbstractHistogramSupplier extends ChartSupplier
		{
		protected abstract <X,Y> XYChart<X,Y> create(Axis<X> xAxis,Axis<Y> yAxis);
		
		@Override
		public Chart get() {
			String header[]=null;
			
			if(!there_is_no_table_header) {
				String firstLine = lineSupplier.get();
				if(firstLine==null) {
					return null;
					}
				header = this.delimiter.split(firstLine);
				if(header.length<=1) return null;
				}
			final List<XYChart.Series<String, Number>> series = new ArrayList<>();

			for(;;)
				{
				String line = lineSupplier.get();
				if(line==null) break;
				if(StringUtil.isBlank(line)) continue;
				final String tokens[]=this.delimiter.split(line);
				if(header==null) // there_is_no_table_header
					{
					/* create the header */
					header= new String[tokens.length];
					for(int x=0;x<header.length;++x) {
						header[x]="$"+(x);
						}
					}
				
				
				final XYChart.Series<String, Number> serie = new XYChart.Series<>();
				serie.setName(tokens[0]);
				for(int x=1;x< header.length && x < header.length;++x)
				 	{
					Double yVal = parseDouble.apply(tokens[x]);
					if(yVal==null || yVal.doubleValue()<=0) yVal=0.0; 
					serie.getData().add(new XYChart.Data<String,Number>(header[x],yVal));
				 	}
				series.add(serie);
				}
			if(header==null) return null;
		    final CategoryAxis xAxis = new CategoryAxis();
		    final NumberAxis yAxis = new NumberAxis();
			final XYChart<String, Number> sbc =create(xAxis, yAxis);
			sbc.getData().addAll(series);
			sbc.setLegendVisible(series.size()>1);
			updateAxisX(xAxis);
			updateAxisY(yAxis);
			return sbc;
			}
		}
	
	
	private class HistogramSupplier extends AbstractHistogramSupplier {
		@Override
		protected <X, Y> XYChart<X, Y> create(Axis<X> xAxis, Axis<Y> yAxis) {
			return new BarChart<>(xAxis, yAxis);
			}
		}
	private class StackedHistogramSupplier extends AbstractHistogramSupplier {
		@Override
		protected <X, Y> XYChart<X, Y> create(Axis<X> xAxis, Axis<Y> yAxis) {
			return new StackedBarChart<>(xAxis, yAxis);
			}
		}
	
	private class StackedHistogramPivoted extends ChartSupplier
		{
		@Override
		public Chart get() {
			String header[]=null;
			
			if(!there_is_no_table_header) {
				String firstLine = lineSupplier.get();
				if(firstLine==null) {
					return null;
					}
				header = this.delimiter.split(firstLine);
				if(header.length<=1) return null;
				}
			List<XYChart.Series<String, Number>> series = null;
	
			for(;;)
				{
				String line = lineSupplier.get();
				if(line==null) break;
				if(StringUtil.isBlank(line)) continue;
				final String tokens[]=this.delimiter.split(line);
				if(header==null) // there_is_no_table_header
					{
					/* create the header */
					header= new String[tokens.length];
					for(int x=0;x<header.length;++x) {
						header[x]="$"+(x);
						}
					}
				if(series == null) //first time we see the serie/header
					{
					series= new ArrayList<>();
					for(int x=1;x<header.length;++x) {
						final XYChart.Series<String, Number> serie = new XYChart.Series<>();
						serie.setName(header[x]);
						series.add(serie);
						}
					}
				
				for(int x=1;x< header.length && x < header.length;++x)
				 	{
					Double yVal = parseDouble.apply(tokens[x]);
					if(yVal==null || yVal.doubleValue()<=0) yVal=0.0; 
					series.get(x-1).getData().add(new XYChart.Data<String,Number>(tokens[0],yVal));
				 	}
				}
			if(header==null || series==null || series.isEmpty()) return null;
		    final CategoryAxis xAxis = new CategoryAxis();
		    final NumberAxis yAxis = new NumberAxis();
			final XYChart<String, Number> sbc =new StackedBarChart<>(xAxis, yAxis);
			sbc.getData().addAll(series);
			updateAxisX(xAxis);
			updateAxisY(yAxis);
			return sbc;
			}
		}

	
	
	private class XYVHistogramSupplier extends ChartSupplier
		{
		boolean stacked=false;
		
		public XYVHistogramSupplier setStacked(boolean stacked) {
			this.stacked = stacked;
			return this;
			}
		private <X,Y> XYChart<X,Y> create(Axis<X> xAxis,Axis<Y> yAxis)
			{	
			return this.stacked?
				new StackedBarChart<>(xAxis, yAxis):
				new BarChart<>(xAxis, yAxis)
				;
			}

		
		
		@Override
		public Chart get() {
			Map<String,Map<String,Double>> xyv = createMapWithStringKey();
			final Set<String> allylabels = new LinkedHashSet<>();
			for(;;)
				{
				final Object[] triple=nextTriple();
				if(triple==null) break;
				
				final String keyx = String.class.cast(triple[0]);
				Map<String,Double> hashx = xyv.get(keyx);
				if(hashx==null) {
					hashx= createMapWithStringKey();;
					xyv.put(keyx,hashx);
					}
				
				final String keyy = String.class.cast(triple[1]);
				allylabels.add(keyy);
				if(hashx.containsKey(keyy)) {
					LOG.warn("Duplicate value for "+keyx+"/"+keyy);
					continue;
				}
				
				final Double yVal = Double.class.cast(triple[2]);
				hashx.put(keyy, yVal);
				}
			if(xyv.isEmpty()) return null;
		    final CategoryAxis xAxis = new CategoryAxis();
		    final NumberAxis yAxis = new NumberAxis();
		    final XYChart<String, Number> sbc =create(xAxis, yAxis);
		    for(final String keyy: allylabels)
				{
			    final XYChart.Series<String, Number> serie = new XYChart.Series<>();
				serie.setName(keyy);
				for(final String keyx:xyv.keySet())
					{
					Double yVal = xyv.get(keyx).getOrDefault(keyy, 0.0);
					
					serie.getData().add(new XYChart.Data<String,Number>(keyx,yVal));
					}
				sbc.getData().add(serie);
				}
			updateAxisX(xAxis);
			updateAxisY(yAxis);
			return sbc;
			}
		}

	@SuppressWarnings("unused")
	private class BubbleChartSupplier extends ChartSupplier {
		int num_cols=2;
		
		BubbleChartSupplier numCols(int n) {
			this.num_cols = n;
			return this;
		}
		
		@Override
		public Chart get() {
			final MinMax minmaxX = new MinMax();
			final MinMax minmaxY = new MinMax();
			final XYChart.Series<Number, Number> serie = new XYChart.Series<>();
			for(;;)
				{
				String line = lineSupplier.get();
				if(line==null) break;
				if(StringUtil.isBlank(line)) continue;
				String tokens[] = this.delimiter.split(line);
				if(tokens.length< this.num_cols) continue;
				Double vals[] = new Double[this.num_cols];
				int x=0;
				for(x=0;x<vals.length;++x)
					{
					vals[x] = parseDouble.apply(tokens[x]);
					if(vals[x]==null) break;
					}
				if(x!=vals.length) continue;
				minmaxX.accept(vals[0]);
				minmaxY.accept(vals[1]);
				
				serie.getData().add(new XYChart.Data<Number,Number>(vals[0],vals[1]));
				}
			if(minmaxX.isEmpty()) return null;
			if(minmaxY.isEmpty()) return null;
	        final NumberAxis xAxis = minmaxX.toAxis();
	        final NumberAxis yAxis = minmaxY.toAxis();
			final BubbleChart<Number,Number> blc = new
		            BubbleChart<Number,Number>(xAxis,yAxis);
			blc.getData().add(serie);
			updateAxisX(xAxis);
			updateAxisY(yAxis);
			return blc;
		}
	}
	

	private class HeatmapSupplier extends ChartSupplier {
		
		
		@Override
		public Chart get() {
			final HeatMapChart<Double> map=new HeatMapChart<>();
			for(;;) {
				final Object[] triple= nextTriple();
				if(triple==null) break;
				final String x=String.class.cast(triple[0]);
				final String y=String.class.cast(triple[1]);
				if(map.contains(x,y)) {
					LOG.warn("duplicate x/y "+x+"/"+y);
					}
				map.put(x,y,Double.class.cast(triple[2]));
				}
			updateAxisX(map.getXAxis());
			updateAxisY(map.getYAxis());
			return map;
		}
	}
	
	private enum PlotType {
		UNDEFINED,
		BEDGRAPH,
		PIE,
		SIMPLE_HISTOGRAM,
		HISTOGRAM,
		STACKED_HISTOGRAM,
		STACKED_HISTOGRAM_PIVOTED,
		XYV,
		STACKED_XYV,
		HEATMAP
	};
	

	@Parameter(names= {"--type","-t"},description = "type")
	private PlotType chartType = PlotType.UNDEFINED;
	@Parameter(names= {"-R","--reference"},description = Launcher.INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File faidx = null;
	@Parameter(names= {"--min-reference-size"},description ="When using a *.dict file, discard the contigs having a length < 'size'. Useful to discard the small unused contigs like 'chrM'. -1 : ignore.")
	private int min_contig_size = -1;
	@Parameter(names= {"-o","--out"},description = "Output file. If defined, save the picture in this file extension:png, jpg or R (experimental) and close the application.")
	private File outputFile = null;
	@Parameter(names= {"-chrompos","--chrom-position"},description = "When reading a genomic file. Input is not a BED file but a CHROM(tab)POS file. e.g: output of `samtools depth`")
	private boolean input_is_chrom_position=false;
	//@Parameter(names= {"--title-side"},description = "Title side")
	//private Side titleSide=Side.TOP;
	//@Parameter(names= {"--legend-side"},description = "Legend side")
	//private Side legendSide=Side.RIGHT;
	@Parameter(names= {"--title"},description = "Chart Title")
	private String chartTitle = "";
	@Parameter(names= {"--hide-legend"},description = "Hide Legend")
	private boolean hide_legend;
	@Parameter(names= {"-su","--sort-unique"},description = "For PIE or SIMPLE_HISTOGRAM the input is the output of a `sort | uniq -c` pipeline")
	private boolean input_is_sort_uniq=false;
	@Parameter(names= {"-xlab","-xlabel","--xlabel"},description = "X axis label.")
	private String xlabel=null;
	@Parameter(names= {"-ylab","-ylabel","--ylabel"},description = "Y axis label.")
	private String ylabel=null;
	@Parameter(names= {"-nh","--no-header"},description = "There is no header. Tested for Histograms")
	private boolean there_is_no_table_header;

	
	
	private Supplier<String> lineSupplier = ()->null;
	
	
	
	
	
	@Override
	public int doWork(final List<String> args) {
		Chart chartNode =null;
		try
			{
			
			
			final BufferedReader br= super.openBufferedReader(oneFileOrNull(args));
			
			
			this.lineSupplier = ()-> {
				try {for(;;) {
					final String L=br.readLine();
					if(L==null) return null;
					return L;
					}}
				catch(final IOException err) {
				throw new RuntimeIOException(err);
				}};
			
			switch(this.chartType) {
				case BEDGRAPH:
					{
						chartNode = new BedGraphSupplier().get();
					break;
					}
				case PIE: chartNode = new PieChartSupplier().get();break;
				case SIMPLE_HISTOGRAM : chartNode = new SimpleHistogramSupplier().get();break;
				case HISTOGRAM: chartNode = new HistogramSupplier().get();break;
				case STACKED_HISTOGRAM: chartNode = new StackedHistogramSupplier().get();break;
				case STACKED_HISTOGRAM_PIVOTED : chartNode = new StackedHistogramPivoted().get();break;
				case XYV:
				case STACKED_XYV:
					chartNode = new XYVHistogramSupplier().
						setStacked(this.chartType==PlotType.STACKED_XYV).
						get();
					break;
				case HEATMAP:
					chartNode = new HeatmapSupplier().get();
					break;					
				default:
					{
					LOG.error("Bad chart type : " + this.chartType);
					return -1;
					}
				}
			CloserUtil.close(br);
			
			if(chartNode==null) {
				LOG.error("No chart was generated");
				return -1;
			}
			
			if(StringUtil.isBlank(this.chartTitle)) {
				chartNode.setLegendVisible(false);
				}
			else
				{
				//chart.setTitleSide(this.titleSide);
				chartNode.setTitle(this.chartTitle);
				}
			//chart.setLegendSide(this.legendSide);
			chartNode.setLegendVisible(!this.hide_legend);
			

	  		try(PrintWriter pw = super.openFileOrStdoutAsPrintWriter(this.outputFile)) {
				final RExporter exporter=new RExporter();
				exporter.exportToR(pw,chartNode);
				pw.flush();
				}

			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			}
		
		
			
			
		
        return 0;
		}
	
	
	public static void main(final String[] args) {
		new SimplePlot().instanceMainWithExit(args);
		}
	}
