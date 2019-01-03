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
package com.github.lindenb.jvarkit.tools.bamstats01;

import java.io.Closeable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.Vector;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jfx.JFXChartExporter;
import com.github.lindenb.jvarkit.math.RangeOfIntegers;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.JfxLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.concurrent.Task;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.SnapshotParameters;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Alert;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Label;
import javafx.scene.image.WritableImage;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Screen;
import javafx.stage.Stage;
/*
BEGIN_DOC

## Examples

todo

## Screenshot

todo

END_DOC
*/
@Program(name="bamstatsjfx",
description="GUI: BAM statistics",
keywords={"bam","stats","jfx"},
generate_doc=false
)
public class BamStatsJfx extends JfxLauncher {
	private static final Logger LOG=Logger.build(BamStatsJfx.class).make();
	private final DecimalFormat niceIntFormat = new DecimalFormat("###,###");
	private final ReentrantLock lock = new ReentrantLock();
	private ProgressBar progressBar;

	@Parameter(names={"-s","--seconds"},description="Refresh screen every 's' seconds")
	private int refreshEverySeconds = 15;
	@Parameter(names={"-o","--output"},description="output images in directory or zip file (filename ends with '.zip') or 'R' file (filename ends with '.R' **UNDER CONSTRUCTION**) . If defined, the application will exit automatically")
	private File outputFile = null;
	@Parameter(names={"--stdin"},description="if there is no file argument. Read vcf from stdin instead of opening a FileOpen dialog")
	private boolean vcf_stdin = false;
	@Parameter(names={"-hr","--hr"},description="Ignore HOM_REF in genotype type.")
	private boolean ignore_HOM_REF = false;
	@Parameter(names={"-fgt","--fgt"},description="Ignore filtered **GENOTYPES**")
	private boolean ignore_filtered_genotypes = false;
	@Parameter(names={"--prefix"},description="Title Prefix")
	private String titlePrefix="";
	@Parameter(names={"--read-lengths"},description="tranches for Read lengths. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers readLengthTranches = new RangeOfIntegers(0,10,20,30,40,50,100,150,200,300,400,500,1000);
	@Parameter(names={"--max-read-length"},description="max read length for base composition. 0: ignore; <0 not limit")
	private int max_read_length = 500;
	@Parameter(names={"--clip-lengths"},description="tranches for Clip lengths. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers clipLengthTranches = new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10,20,30,40,50);

	 
	private volatile boolean stop = false;
	private abstract class ChartGenerator
		{
		Tab tab;
		boolean enabled = true;
		long nRecords = 0;
		boolean isEnabled() {
			return enabled;
			}
		String getChartTitle() {
			return "untitled";
		}
			
		String getTabTitle() {
			return getChartTitle();
		}
		Tab makeTab() {
			this.tab = new Tab(getTabTitle());
			this.tab.setClosable(true);
			this.tab.setOnClosed(AE->{this.enabled=false;});
			return this.tab;
			}
		
		void title(final Chart c,final String s) {
			c.setTitle(
					StringUtil.isBlank(titlePrefix)?
					s:
					titlePrefix+ " : " + s
					);
			}
		
		abstract void visit(final SAMRecord ctx) ;
		
		Chart makeChart() {
			return null;
			}
		
		void refresh() {
			if(!isEnabled()) return;
			final Chart chart= makeChart();
			if(chart==null) return;
			if(outputFile!=null) chart.setAnimated(false);
			this.tab.setContent(chart);
			}
		public String getFilename() {
			return getTabTitle().replaceAll("[^A-Za-z_0-9]+","")+".png";
			}
		
		protected void saveR(final JFXChartExporter exporter) {
		 	if(!isEnabled()) return;
		 	final Node content = this.tab.getContent();
		 	if(content==null || !(content instanceof Chart)) return;
		 	final Chart chart = Chart.class.cast(content);
		 	exporter.exportToR(chart);
			}
		protected void saveImageAs(final File dir)
			 	throws IOException
			 	{
			 	if(!isEnabled()) return;
			 	final Node content = this.tab.getContent();
			 	if(content==null) return;
	    		WritableImage image = content.snapshot(new SnapshotParameters(), null);
	    		final File file = new File(getFilename());
	            ImageIO.write(SwingFXUtils.fromFXImage(image, null), "png", file);
			 	}
		protected void saveImageAs(final ZipOutputStream zout)
			 	throws IOException
			 	{
			 	if(!isEnabled()) return;
			 	final Node content = this.tab.getContent();
			 	if(content==null) return;
	    		WritableImage image = content.snapshot(new SnapshotParameters(), null);
	            ImageIO.write(SwingFXUtils.fromFXImage(image, null), "png", zout);
			 	}
		}
	
	private class ReadLengthGenerator extends ChartGenerator
		{
		private final Counter<RangeOfIntegers.Range> len2count = new Counter<>();
		@Override
		String getChartTitle() {
			return "Read-Length";
			}
		
		
		@Override
		Chart makeChart() {

			if(this.len2count.isEmpty()) return null;
			
			final List<String> lengthCategories = this.len2count.keySet().
					stream().
					sorted().
					map(x->x.toString()).
					collect(Collectors.toList());
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(
					FXCollections.observableArrayList(lengthCategories)
					);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			
			final BarChart<String,Number> bc = 
		            new BarChart<>(xAxis,yAxis);
		
			
			for(final RangeOfIntegers.Range x:new TreeSet<>(this.len2count.keySet()))
				{
				
				series1.getData().add(new XYChart.Data<String,Number>(
						x.toString(),
						this.len2count.count(x)
					));
						
				}
			bc.getData().add(series1);
			bc.setCategoryGap(1);
			title(bc,
	        		this.getChartTitle()+
	        		" N. Primary SAMRecord "+niceIntFormat.format(nRecords)
	        		)
	        		;
	        xAxis.setLabel("Length");
	        yAxis.setLabel("Number of Reads");
	        xAxis.setTickLabelRotation(90);
	        bc.setLegendVisible(false);
	        return bc;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			this.nRecords ++ ;
			final RangeOfIntegers.Range r = readLengthTranches.getRange(rec.getReadLength());
			this.len2count.incr(r);
			}
		
		}
	private class ClippedReadsGenerator extends ChartGenerator
		{
		private final Counter<RangeOfIntegers.Range> clip5 = new Counter<>();
		private final Counter<RangeOfIntegers.Range> clip3 = new Counter<>();
		
		@Override
		String getChartTitle() {
			return "Clip";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(FXCollections.observableArrayList(clipLengthTranches.getRanges().stream().map(O->O.toString()).collect(Collectors.toList())));
			
			final StackedBarChart<String,Number> bc = 
		            new StackedBarChart<String,Number>(xAxis,yAxis);
			
	
			for(int side=0;side<2;++side)
				{
				final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
				series1.setName(side==0?"5'":"3'");
				for(final RangeOfIntegers.Range r: clipLengthTranches.getRanges()) {
					series1.getData().add(new XYChart.Data<String,Number>(
							r.toString(),
							(side==0?this.clip5:this.clip3).count(r)
							));
					}
				bc.getData().add(series1);
				}
			bc.setCategoryGap(1);
			bc.setLegendVisible(true);;
			title(bc,this.getChartTitle()+" N. Clipped Records: "+niceIntFormat.format(nRecords));
	        xAxis.setLabel("Clip-Length");
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			if(rec.getReadUnmappedFlag()) return ;
			final Cigar cigar = rec.getCigar();
			if(cigar==null || cigar.isEmpty()) return;
			if(!cigar.isClipped()) return;
			this.nRecords++;
			CigarElement ce =cigar.getFirstCigarElement();
			if(ce.getOperator().isClipping()) {
				this.clip5.incr(clipLengthTranches.getRange(ce.getLength()));
			}
			
			ce =cigar.getLastCigarElement();
			if(ce.getOperator().isClipping()) {
					this.clip3.incr(clipLengthTranches.getRange(ce.getLength()));
				}
			}
		}
	
	private class BaseCompositionGenerator extends ChartGenerator
		{
		private final List<Counter<Character>> pos2count = new ArrayList<>();
		BaseCompositionGenerator() {
			}
		
		@Override
		String getChartTitle() {
			return "Base Composition";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(
					FXCollections.observableArrayList(
						IntStream.
						range(0, pos2count.size()).
						mapToObj(x->niceIntFormat.format(x+1)).
						collect(Collectors.toList())
						)
					);
			final TreeSet<Character> bases = pos2count.stream().
					flatMap(F->F.keySet().stream()).
					collect(Collectors.toCollection(TreeSet::new));
			
			final StackedBarChart<String,Number> bc = 
		            new StackedBarChart<String,Number>(xAxis,yAxis);
		

			for(final Character base:bases)
				{
				final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
				series1.setName(base.toString());
				for(int x=0;x< this.pos2count.size();++x) {
					series1.getData().add(new XYChart.Data<String,Number>(
							niceIntFormat.format(x+1),
							this.pos2count.get(x).count(base)
							));
					}
				bc.getData().add(series1);
				}
			bc.setCategoryGap(1);
			title(bc,this.getChartTitle()+" N. Records: "+niceIntFormat.format(nRecords));
	        xAxis.setLabel("Position");
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Bases");
	        return bc;
			}
		
		private void visit(char c,int pos) {
			this.pos2count.get(pos).incr(Character.toUpperCase(c));
		}
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			this.nRecords++;
			final byte[] readBases = rec.getReadBases();
			
			while(this.pos2count.size()<=readBases.length && (max_read_length<0 || this.pos2count.size()<max_read_length)) {
				this.pos2count.add(new Counter<>());
				}
			
			if(rec.getReadNegativeStrandFlag())
				{
				for(int i=0;i< readBases.length && (max_read_length<0 || i<max_read_length);i++)
					{
					visit(AcidNucleics.complement((char)readBases[(readBases.length-1)-i]),i);
					}
				}
			else
				{
				for(int i=0;i< readBases.length && (max_read_length<0 || i<max_read_length);i++)
					{
					visit((char)readBases[i],i);
					}
				}
			}
		}

	private class QualGenerator extends ChartGenerator
		{
		private class Average {
			double sum=0.0;
			long n=0L;
			double get() {
				return n==0?0.0:sum/n;
			}
		}
		private final Map<Integer,Average> pos2qual = new HashMap<>();
		
		@Override
		String getChartTitle() {
			return "Quality";
			}
		
		@Override
		Chart makeChart() {

			if(this.pos2qual.isEmpty()) return null;
			
			final int max_length = pos2qual.keySet().stream().mapToInt(X->X.intValue()).max().getAsInt();
			
			final NumberAxis yAxis = new NumberAxis();
			final NumberAxis xAxis = new NumberAxis();
			
			final LineChart<Number,Number> bc = 
		            new LineChart<Number,Number>(xAxis,yAxis);
			final XYChart.Series<Number,Number> series1 = new XYChart.Series<>();
			for(int x=0;x<= max_length;++x)
				{
				series1.getData().add(new XYChart.Data<Number,Number>(
					x+1,
					pos2qual.get(x).get()
					));
				}
			bc.getData().add(series1);
						
			bc.setLegendVisible(false);
			title(bc,this.getChartTitle());
	        xAxis.setLabel("Position");
			xAxis.setLowerBound(0.0);
			xAxis.setUpperBound(max_length);
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Avg Qual.");
	        yAxis.setLowerBound(0);
	        return bc;
			}
		
		private void visit(int x,byte q) {
			Average a = this.pos2qual.get(x);
			if(a==null) {
				a = new Average();
				this.pos2qual.put(x, a);
				}
			a.n++;
			a.sum+=q;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			final byte[] readQuals = rec.getBaseQualities();
			if(readQuals.length==0) return;
			this.nRecords++;
			for(int i=0;i< readQuals.length ;i++)
				{
				int j = rec.getReadNegativeStrandFlag()? (readQuals.length-1-i):i;
				visit(j,readQuals[j]);	
				}
			}
		}
	
	private class GCPercentGenerator extends ChartGenerator
		{
		private final Counter<Integer> mappedgc2count = new Counter<>();
		private final Counter<Integer> unmappedgc2count = new Counter<>();
		GCPercentGenerator() {
			}
		
		@Override
		String getChartTitle() {
			return "GC Percent";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final NumberAxis xAxis = new NumberAxis();
			
			final LineChart<Number,Number> bc = 
		            new LineChart<Number,Number>(xAxis,yAxis);
	        for(int side=0;side< 2;++side) {
				final XYChart.Series<Number,Number> series1 = new XYChart.Series<>();
				series1.setName(side==0?"Mapped":"Unmapped");
				final Counter<Integer> rc =(side==0?this.mappedgc2count:this.unmappedgc2count);
				for(final Integer percent: new TreeSet<>(rc.keySet()))
					{
					series1.getData().add(new XYChart.Data<Number,Number>(
						percent,
						rc.count(percent)
						));
					}
				bc.getData().add(series1);
				}			
			bc.setLegendVisible(true);
			title(bc,this.getChartTitle());
	        xAxis.setLabel("GC%");
			xAxis.setLowerBound(0.0);
			xAxis.setUpperBound(100.0);
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Count");
	        yAxis.setLowerBound(0);
	        return bc;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			if(rec.isSecondaryOrSupplementary()) return;
			this.nRecords++;
			final byte[] readBases = rec.getReadBases();
			if(readBases.length==0) return;
			double ngc=0;
			for(int i=0;i< readBases.length;++i) {
				switch(readBases[i])
					{
					case 's':case 'S':
					case 'g':case 'G':
					case 'c':case 'C': ngc++;break;
					}
				}
			final int ngi = (int)((ngc/(double)readBases.length)*100.0);
			if(rec.getReadUnmappedFlag())
				{
				this.unmappedgc2count.incr(ngi);
				}
			else
				{
				this.mappedgc2count.incr(ngi);
				}
			}
		}
	
	
	private class ContigUsageGenerator extends ChartGenerator
		{
		private final Counter<String> count = new Counter<>();
		private final SAMSequenceDictionary dict;
		ContigUsageGenerator(final SAMSequenceDictionary dict) {
			this.dict = dict;
			}
		
		@Override
		String getChartTitle() {
			return "Contigs";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final ObservableList<String> contigs = 
					FXCollections.observableArrayList(this.dict.getSequences().stream().
						map(S->S.getSequenceName()).
						filter(S->count.count(S)>0).
						collect(Collectors.toList()))
						;
			if(this.count.count(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)>0L)
				{
				contigs.add(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
				}
			
			final CategoryAxis xAxis = new CategoryAxis(contigs);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final SAMSequenceRecord ssr: this.dict.getSequences())
				{
				long n = this.count.count(ssr.getSequenceName());
				if(n==0L) continue;
				series1.getData().add(new XYChart.Data<String,Number>(ssr.getSequenceName(),n));
				}
			long n=this.count.count(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
			if(n>0L)
				{
				series1.getData().add(new XYChart.Data<String,Number>(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME,n));
				}
			bc.getData().add(series1);
			title(bc,this.getChartTitle()+ " Num. Reads. "+niceIntFormat.format(this.nRecords));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("Contig");       
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			this.nRecords++;
			if(rec.getReadUnmappedFlag())
				{
				this.count.incr(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
				}
			else
				{
				this.count.incr(rec.getContig());
				}
			}
		}
	
	private class SamFlagUsagGenerator extends ChartGenerator
		{
		private final Counter<SAMFlag> count = new Counter<>();
		
		@Override
		String getChartTitle() {
			return "SamFlags";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final ObservableList<String> contigs = 
					FXCollections.observableArrayList(
						Arrays.asList(SAMFlag.values()).stream().map(S->S.name()).collect(Collectors.toList())
						)
						;
			
			final CategoryAxis xAxis = new CategoryAxis(contigs);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final SAMFlag f: SAMFlag.values())
				{
				series1.getData().add(new XYChart.Data<String,Number>(f.name(),count.count(f)));
				}
			bc.getData().add(series1);
			title(bc,this.getChartTitle()+ " Num. Reads. "+niceIntFormat.format(this.nRecords));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("Flag");
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final SAMRecord rec) {
			this.nRecords++;
			for(final SAMFlag f: rec.getSAMFlags())
				{
				count.incr(f);
				}
			}
		}
	
	private final List<ChartGenerator> chartGenerators = new Vector<>();
	private class VariantContextRunner
		extends Task<Void>
		implements Runnable,Closeable
		
		{
		final SamReader samReader;
		final SAMRecordIterator iter;
		
		VariantContextRunner(final SamReader samReader,final SAMRecordIterator iter)
			{
			this.samReader = samReader;
			this.iter = iter;
			}
		@Override
		protected Void call() throws Exception {
			refreshCharts();
			long last = -1L;
			
			Platform.runLater(()->{
				progressBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
				});
			while(!BamStatsJfx.this.stop && this.iter.hasNext()) {
				final SAMRecord rec = this.iter.next();
				try 
					{
					lock.lock();
					for(final ChartGenerator cg: chartGenerators) {
						if(!cg.isEnabled()) continue;
						cg.visit(rec);
						}
					}
				finally
					{
					lock.unlock();
					}
				final long now = System.currentTimeMillis();
				if(last<0L || now - last > refreshEverySeconds * 1000) {
					last=now;
					refreshCharts();
					}
				}
			
			refreshCharts();
			
			Platform.runLater(()->{
				progressBar.setProgress(1.0);
				});
			
			if(outputFile!=null && !stop)
	        	{
				LOG.info("saving as "+outputFile+ " in "+refreshEverySeconds+ " secs.");;
        		new java.util.Timer().schedule( 
    		        new java.util.TimerTask() {
    		            @Override
    		            public void run() {
    		            Platform.runLater(()->{
		        		 try
		        		 	{
		        			if(!stop) saveImagesAs(outputFile);
		        		 	}
		        		 catch(final IOException err)
		        		 	{
		        			LOG.error(err);
		        			System.exit(-1);
		        		 	}
		        		 LOG.info("Platform.exit called");
		        		 Platform.exit();
    		            	});}
    		            },refreshEverySeconds*1000);
	        	 	
	        	}
			return null;
			}
		@Override
		public void close() throws IOException {
			CloserUtil.close(this.samReader);
			CloserUtil.close(this.iter);
			}
		}
	
	 void refreshCharts() {
		Platform.runLater(()->{
			try {
				lock.lock();
				for(final ChartGenerator cg:chartGenerators) {
					if(!cg.isEnabled()) continue;
					cg.refresh();
					}
				}
			finally
				{
				lock.unlock();
				}
			});
		}
	
	@Override
	protected int doWork(final Stage primaryStage, final List<String> args) {
		
		try {
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
							validationStringency(ValidationStringency.LENIENT)
							;
			final SamReader samReader;
			final SAMRecordIterator samIterator;
			if(args.size()==1) {
				samReader = srf. open(new File(args.get(0)));
				samIterator = samReader.iterator();
				}
			else if(args.isEmpty() && this.vcf_stdin)
				{
				samReader = srf. open(SamInputResource.of(System.in));
				samIterator = samReader.iterator();
				}
			else if(args.isEmpty())
				{
				final FileChooser fc = new FileChooser();
				final File f = fc.showOpenDialog(null);
				if(f==null) return -1;
				samReader = srf. open(f);
				samIterator = samReader.iterator();
				}
			else
				{
				LOG.error("illegal number of arguments.");
				return -1;
				}
			final SAMFileHeader header =  samReader.getFileHeader();

			
			
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict!=null && !dict.isEmpty()) {
				chartGenerators.add(new ContigUsageGenerator(dict));
				}
			this.chartGenerators.add(new ReadLengthGenerator());
			if(max_read_length!=0) { /* negative : not limit */
				this.chartGenerators.add(new BaseCompositionGenerator());
				}
			this.chartGenerators.add(new GCPercentGenerator());
			this.chartGenerators.add(new QualGenerator());
			this.chartGenerators.add(new SamFlagUsagGenerator());
			this.chartGenerators.add(new ClippedReadsGenerator());
			
			
			
			final VariantContextRunner runner = new VariantContextRunner(samReader,samIterator);
			
			primaryStage.setOnShowing(AE->{
				new Thread(runner).start();
			});
			primaryStage.setOnCloseRequest(AE->{
				stop = true;
				CloserUtil.close(runner);
			});
			final BorderPane contentPane=new BorderPane();
			final TabPane tabPane = new TabPane();
			final MenuBar menuBar = new MenuBar();
		    Menu menu = new Menu("File");
		    MenuItem item;
		    item=new MenuItem("About...");
		    item.setOnAction(AE->doMenuAbout(AE));
		    menu.getItems().add(item);
		    menu.getItems().add(new SeparatorMenuItem());
		    item=new MenuItem("Save Current image as... (.R,.png,.jpg)");
		    item.setOnAction(AE->{doMenuSaveCurrentImage(tabPane.getSelectionModel().getSelectedItem());});
		    menu.getItems().add(item);
		    item=new MenuItem("Save All images in directory... ");
		    item.setOnAction(AE->doMenuSaveAllImages());
		    menu.getItems().add(item);
		    item=new MenuItem("Save All images as ... (.R,.zip)");
		    item.setOnAction(AE->doMenuSaveAllImagesInFile());
		    menu.getItems().add(item);
		    menu.getItems().add(new SeparatorMenuItem());
		    item=new MenuItem("Quit");
		    item.setOnAction(AE->{stop=true;Platform.exit();});
		    menu.getItems().add(item);
		    menuBar.getMenus().add(menu);
		    contentPane.setTop(menuBar);
		    
		    progressBar = new ProgressBar();
			
			
			contentPane.setCenter(tabPane);
			contentPane.setBottom(new HBox(new Label("Progress:"),this.progressBar));
			
			
			chartGenerators.removeIf(G->!G.isEnabled());
			for(final ChartGenerator g:chartGenerators) {
				tabPane.getTabs().add(g.makeTab());
				}
			final Screen scr = Screen.getPrimary();
			final Scene scene  = new Scene(
					contentPane,
					scr.getBounds().getWidth()-100,
					scr.getBounds().getHeight()-100
					);
			primaryStage.setScene(scene);
			primaryStage.setTitle(getClass().getSimpleName());
			
	        primaryStage.show();
			}
		catch(final Exception err) {
			err.printStackTrace();
			return -1;
			}
		finally
			{
			
			}
		return 0;
		}
	
	private void doMenuSaveAllImagesInFile()
		{
		final FileChooser fc = new FileChooser();
		final File f = fc.showSaveDialog(null);
		if(f==null) return;
		if(!(f.getName().endsWith(".R") || f.getName().endsWith(".zip"))){
			  final Alert alert=new Alert(AlertType.ERROR,
					  "Filename must end with .R or .zip:"+f);
	          alert.showAndWait();
	          return;
			}
		try {
			saveImagesAs(f);
			}
    	catch (final IOException e) {
            super.displayAlert(e);
        	}
		}
	
	private void doMenuSaveAllImages() {
		final DirectoryChooser fc = new DirectoryChooser();
		final File dir = fc.showDialog(null);
		if(dir==null) return;
		if(!dir.exists() || !dir.isDirectory()) {
			  final Alert alert=new Alert(AlertType.ERROR,
					  "Bad directory :"+dir);
	          alert.showAndWait();
	          return;
			}
		try {
			saveImagesAs(dir);
			}
    	catch (final IOException e) {
    		 super.displayAlert(e);
        	}
		}
	
	private void doMenuSaveCurrentImage(final Tab tab) {
		if(tab==null) return;
		final Node content =  tab.getContent();
	 	if(content==null) return;
		final FileChooser fc = new FileChooser();
		final File file = fc.showSaveDialog(null);
    	if(file==null) return;
    	
    	
    	PrintWriter pw= null;
    	try {
    		if(file.getName().endsWith(".R")) {
        		if(content instanceof Chart) {
            		pw = new PrintWriter(file);
            		JFXChartExporter exporter = new JFXChartExporter(pw);
            		exporter.exportToR(Chart.class.cast(content));
	        		pw.flush();
	        		pw.close();
        			}
        		}
    		else
	    		{
	    		final WritableImage image = content.snapshot(new SnapshotParameters(), null);
				
	    		final String format = file.getName().toLowerCase().endsWith("png")?"png":"jpg";
				ImageIO.write(SwingFXUtils.fromFXImage(image, null), format, file);
	    		}
    		}
    	catch (final IOException e) {
    		super.displayAlert(e);
            LOG.error(e);
        	}
    	finally
    		{
    		CloserUtil.close(pw);
    		}
		}
	
	private void saveImagesAs(final File out) throws IOException {
		PrintWriter pw =null;
		FileOutputStream fout = null;
		try {
			if(out.getName().endsWith(".R")) {
				pw = new PrintWriter(out);
				final JFXChartExporter chartExporter = new JFXChartExporter(pw);
				for(final ChartGenerator  cg:this.chartGenerators) {
					if(!cg.isEnabled()) continue;
					LOG.info("saving "+cg.getFilename());
					cg.saveR(chartExporter);
					}
				pw.flush();
				pw.close();
				}
			else if(out.getName().endsWith(".zip")) {
				fout = new FileOutputStream(out);
				final ZipOutputStream zout = new ZipOutputStream(fout);
				for(final ChartGenerator  cg:this.chartGenerators) {
					if(!cg.isEnabled()) continue;
					LOG.info("saving "+cg.getFilename());
					final ZipEntry zipEntry = new ZipEntry(cg.getFilename());
					zout.putNextEntry(zipEntry);
					cg.saveImageAs(zout);
					zout.closeEntry();
					}
				zout.finish();
				zout.close();
				fout.close();
				}
			else
				{
				IOUtil.assertDirectoryIsWritable(out);
				for(final ChartGenerator  cg:this.chartGenerators) {
					if(!cg.isEnabled()) continue;
					LOG.info("saving "+cg.getFilename());
					cg.saveImageAs(out);
					}
				}
			} 
	catch(final Exception err) {
		super.displayAlert(err);
		}
	finally
		{
		CloserUtil.close(pw);
		CloserUtil.close(fout);
		}
	}

	public static void main(final String[] args) {
		Application.launch(args);
	}

}
