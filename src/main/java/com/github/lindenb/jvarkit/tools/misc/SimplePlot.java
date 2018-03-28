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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
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
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.embed.swing.SwingFXUtils;
import javafx.geometry.Side;
import javafx.scene.Scene;
import javafx.scene.SnapshotParameters;
import javafx.scene.chart.Axis;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.PieChart;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;
import javafx.scene.image.WritableImage;
import javafx.stage.Screen;
import javafx.stage.Stage;
/**
BEGIN_DOC 

## Examples

```
$ gunzip -c in.vcf.gz | grep -v "^#" | cut -f 1 | sort | uniq -c | java -jar dist/simpleplot.jar -t SIMPLE_HISTOGRAM  -su 
$ gunzip -c in.vcf.gz | grep -v "^#" | cut -f 1 | sort | uniq -c | java -jar dist/simpleplot.jar -t PIE  -su

$ echo -e "Year\tX\tY\n2018\t1\t2\n2019\t3\t4" | java -jar dist/simpleplot.jar -t HISTOGRAM
$ echo -e "Year\tX\tY\n2018\t1\t2\n2019\t3\t4" | java -jar dist/simpleplot.jar -t STACKED_HISTOGRAM


``` 
 
END_DOC
 */
@Program(name="simpleplot",
description="simple figure plotter using java JFX. You'd better use gnuplot or R.",
keywords={"char","figure","jfx"})

public class SimplePlot extends Application {
	private static final Logger LOG = Logger.build(SimplePlot.class).make();
	
	/** base generator for a Chart Maker */
	private abstract class ChartSupplier
		implements Supplier<Chart>
		{
		private SAMSequenceDictionary _dict = null;
		protected Pattern delimiter = Pattern.compile("[\t]");
		private ContigNameConverter _ctgConverter = null;
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
				if(SimplePlot.this.faidx == null) {
					throw new JvarkitException.DictionaryMissing("undefined reference file.");
					}
				
				this._dict = new SAMSequenceDictionary(
						SAMSequenceDictionaryExtractor.extractDictionary(SimplePlot.this.faidx).
							getSequences().
							stream().
							filter(C->min_contig_size<0|| C.getSequenceLength()>= SimplePlot.this.min_contig_size).
							collect(Collectors.toList())
							);
				if(this._dict==null || this._dict.isEmpty()) {
					throw new JvarkitException.DictionaryMissing("empty dictionary extracted from "+SimplePlot.this.faidx);
					}
				this._ctgConverter =  ContigNameConverter.fromOneDictionary(this._dict);
				this._ctgConverter.setOnNotFound(ContigNameConverter.OnNotFound.SKIP);
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
				Double minY=null;
				Double maxY=null;
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
					if(minY==null )
						{
						minY=maxY=yVal;
						}
					else
						{
						if(yVal<minY) minY=yVal;
						if(yVal>maxY) maxY=yVal;
						}
					}
				if(minY==null || maxY==null || chroms.isEmpty()) {
					LOG.error("no valid data");
					return null;
				}
		        final NumberAxis xAxis = new NumberAxis(1,dict.getReferenceLength(),dict.getReferenceLength()/10.0);
		        xAxis.setLabel("Genomic Index");
		        final NumberAxis yAxis = new NumberAxis(minY,maxY,(maxY-minY)/10.0); 
		       
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
			final Map<String,Double> key2count = new LinkedHashMap<>();
			for(;;) {
				String line = lineSupplier.get();
				if(line==null) break;
				if(StringUtil.isBlank(line) ||line.startsWith("#") ) continue;
				final String key;
				final Double val;
				if( input_is_sort_uniq) {
					int i=0;
					//skip ws
					while(i<line.length() && Character.isWhitespace(line.charAt(i))) i++;
					
					//read count
					StringBuilder sb = new StringBuilder();
					while(i<line.length() && Character.isDigit(line.charAt(i)))
						{
						sb.append(line.charAt(i));
						i++;
						}
					if(sb.length()==0) continue;
					val = this.parseDouble.apply(sb.toString());
					if(val==null || val.doubleValue()<0.0) continue;
					i++; //skip one whitespace
					if(i>=line.length()) continue;
					key=line.substring(i);
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
			final ObservableList<PieChart.Data> pieChartData =
					FXCollections.observableList(
							key2count.entrySet().stream().
							map(KV->new PieChart.Data(KV.getKey(),KV.getValue())).
							collect(Collectors.toList())
							);
			final  PieChart chart = new PieChart(pieChartData);
			return chart;
			}
		}
	
	/** supplier for simple bar chart */
	private class SimpleHistogramSupplier extends AbstractKeyValueSupplier {
		@Override
		public Chart get() {
			final Map<String,Double> key2count = getKeyValueMap();			
			final ObservableList<StackedBarChart.Data<String,Number>> data =
					FXCollections.observableList(
							key2count.entrySet().stream().
							map(KV->new StackedBarChart.Data<String,Number>(KV.getKey(),KV.getValue())).
							collect(Collectors.toList())
							);
			final XYChart.Series<String, Number> series1 =
			            new XYChart.Series<String, Number>(data);
			final CategoryAxis xAxis = new CategoryAxis();
		    final NumberAxis yAxis = new NumberAxis();
			final StackedBarChart<String, Number> chart =
			            new StackedBarChart<String, Number>(xAxis, yAxis);
			chart.getData().add(series1);
			
			return chart;
			}
		}
	
	private abstract class AbstractHistogramSupplier extends ChartSupplier
		{
		protected abstract <X,Y> XYChart<X,Y> create(Axis<X> xAxis,Axis<Y> yAxis);
		@Override
		public Chart get() {
			String firstLine = lineSupplier.get();
			if(firstLine==null) {
				return null;
				}
			String header[]=this.delimiter.split(firstLine);
			if(header.length<=1) return null;
			
			final List<XYChart.Series<String, Number>> series = FXCollections.observableArrayList();

			for(;;)
				{
				String line = lineSupplier.get();
				if(line==null) break;
				if(StringUtil.isBlank(line)) continue;
				final String tokens[]=this.delimiter.split(line);
				
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
		    final CategoryAxis xAxis = new CategoryAxis();
		    final NumberAxis yAxis = new NumberAxis();
			final XYChart<String, Number> sbc =create(xAxis, yAxis);
			sbc.getData().addAll(series);
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
	
	private enum PlotType {
		UNDEFINED,
		BEDGRAPH,
		PIE,
		SIMPLE_HISTOGRAM,
		HISTOGRAM,
		STACKED_HISTOGRAM
	};
	
	@ParametersDelegate
	private Launcher.UsageBuider usageBuider = null;
	@Parameter(description = "Files")
	private List<String> args = new ArrayList<>();
	@Parameter(names="--testng",description = "testng",hidden=true)
	private boolean this_is_a_unit_test = false;
	@Parameter(names= {"--type","-t"},description = "type")
	private PlotType chartType = PlotType.UNDEFINED;
	@Parameter(names= {"-R","--reference"},description = Launcher.INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File faidx = null;
	@Parameter(names= {"--min-reference-size"},description ="When using a *.dict file, discard the contigs having a length < 'size'. Useful to discard the small unused contigs like 'chrM'. -1 : ignore.")
	private int min_contig_size = -1;
	@Parameter(names= {"-o","--out"},description = "Output file. If defined, save the picture in this file and close the application.")
	private File outputFile = null;
	@Parameter(names= {"-chrompos","--chrom-position"},description = "When reading a genomic file. Input is not a BED file but a CHROM(tab)POS file. e.g: output of `samtools depth`")
	private boolean input_is_chrom_position=false;
	@Parameter(names= {"--title-side"},description = "Title side")
	private Side titleSide=Side.TOP;
	@Parameter(names= {"--legend-side"},description = "Legend side")
	private Side legendSide=Side.RIGHT;
	@Parameter(names= {"--title"},description = "Chart Title")
	private String chartTitle = "";
	@Parameter(names= {"--hide-legend"},description = "Hide Legend")
	private boolean hide_legend;
	@Parameter(names= {"-su","--sort-unique"},description = "For PIE or SIMPLE_HISTOGRAM the input is the output of a `sort | uniq -c` pipeline")
	private boolean input_is_sort_uniq=false;

	
	private Supplier<String> lineSupplier = ()->null;
	
	
	
	
	
	@Override
	public void start(final Stage primaryStage) throws Exception {
		Chart chart=null;
		try
			{
			
			final List<String> params = this.getParameters().getUnnamed();
			this.usageBuider = new Launcher.UsageBuider(SimplePlot.class);
			final JCommander jCommander = new JCommander(this);
			jCommander.parse(params.toArray(new String[params.size()]));
			if(this.usageBuider.shouldPrintUsage())
				{
				this.usageBuider.usage(jCommander);
				Platform.exit();
				return;
				}
			final BufferedReader br;
			if(this.args.isEmpty())
				{
				// open stdin
				br = IOUtils.openStdinForBufferedReader();
				}
			else if(this.args.size()==1)
				{
				br = IOUtils.openURIForBufferedReading(this.args.get(0));
				}
			else
				{
				LOG.error("Illegal Number of arguments: " + this.args);
				Platform.exit();
				return;
				}
			
			this.lineSupplier = ()-> {
				try {for(;;) {
					String L=br.readLine();
					if(L==null) return null;
					return L;
					}}
				catch(final IOException err) {
				throw new RuntimeIOException(err);
				}};
			
			switch(this.chartType) {
				case BEDGRAPH:
					{
					chart = new BedGraphSupplier().get();
					break;
					}
				case PIE: chart = new PieChartSupplier().get();break;
				case SIMPLE_HISTOGRAM : chart = new SimpleHistogramSupplier().get();break;
				case HISTOGRAM: chart = new HistogramSupplier().get();break;
				case STACKED_HISTOGRAM: chart = new StackedHistogramSupplier().get();break;
				default:
					{
					LOG.error("Bad chart type : " + this.chartType);
					Platform.exit();
					return;
					}
				}
			CloserUtil.close(br);
			}
		catch(final Exception err) {
			LOG.error(err);
			Platform.exit();
			return;
			}
		finally
			{
			
			}
		if(chart==null) {
			LOG.error("No chart was generated");
			Platform.exit();
			return;
		}
		
		if(StringUtil.isBlank(this.chartTitle)) {
			chart.setLegendVisible(false);
			}
		else
			{
			chart.setTitleSide(this.titleSide);
			chart.setTitle(this.chartTitle);
			}
		chart.setLegendSide(this.legendSide);
		chart.setLegendVisible(!this.hide_legend);
		
		
		
		if(this.outputFile!=null) 
			{
			chart.setAnimated(false);
      		LOG.info("saving as "+this.outputFile+" and exiting.");
			final Chart theChart=chart;
			primaryStage.setOnShown(WE->{
	       		 try
	       		 	{
	       			saveImageAs(theChart,this.outputFile);
	       		 	}
	       		 catch(final IOException err)
	       		 	{
	       			LOG.error(err);
	       			System.exit(-1);
	       		 	}
	       		 Platform.exit();
				});
			}
		final Screen scr = Screen.getPrimary();
		Scene scene  = new Scene(chart,
				scr.getBounds().getWidth()-100,
				scr.getBounds().getHeight()-100
				);
		primaryStage.setScene(scene);
        primaryStage.show();
		}
	
	/** save char in file */
	private void saveImageAs(
			final Chart   chart,
			final File file)
		 	throws IOException
	 	{
		final WritableImage image = chart.snapshot(new SnapshotParameters(), null);
		final String format=(file.getName().toLowerCase().endsWith(".png")?"png":"jpg");
        ImageIO.write(SwingFXUtils.fromFXImage(image, null), format, file);
	 	}

	
	public static void main(String[] args) {
		Application.launch(args);
		}
	}
