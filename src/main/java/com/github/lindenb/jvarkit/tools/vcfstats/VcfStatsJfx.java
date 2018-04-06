package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.JfxLauncher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.concurrent.Task;
import javafx.scene.Scene;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedBarChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.stage.Screen;
import javafx.stage.Stage;

public class VcfStatsJfx extends JfxLauncher {
	private static final Logger LOG=Logger.build(VcfStatsJfx.class).make();
	private final ReentrantLock lock = new ReentrantLock();
	private int refreshEverySeconds = 2;
	
	private volatile boolean stop = false;
	private class ChartGenerator
		{
		Tab tab;
		boolean enabled = true;
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
			return this.tab;
			}
		void visit(final VariantContext ctx) {
			
			}
		Chart makeChart() {
			return null;
			}
		
		void refresh() {
			LOG.info("refresh()");
			final Chart chart= makeChart();
			if(chart==null) return;
			this.tab.setContent(chart);
			}
		}
	
	private class VariantTypeGenerator extends ChartGenerator
		{
		private final Counter<VariantContext.Type> count = new Counter<>();
		
		@Override
		String getChartTitle() {
			return "Variant Type";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();

			final CategoryAxis xAxis = new CategoryAxis(
					FXCollections.observableArrayList(
						Arrays.asList(VariantContext.Type.values()).
							stream().map(I->I.name()).
							collect(Collectors.toList()))
					);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final VariantContext.Type t:VariantContext.Type.values())
				{
				series1.getData().add(new XYChart.Data<String,Number>(t.name(),this.count.count(t)));
				}
			bc.getData().add(series1);
	        bc.setTitle(this.getChartTitle() +" N="+this.count.getTotal());
	        xAxis.setLabel("Variant Type");       
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			this.count.incr(ctx.getType());
			}
		}
	
	private class FilterUsageGenerator extends ChartGenerator
		{
		private final Counter<String> count = new Counter<>();
		FilterUsageGenerator(Collection<VCFFilterHeaderLine> filters) {
			count.initializeIfNotExists(VCFConstants.PASSES_FILTERS_v4);
			filters.forEach(F->count.initializeIfNotExists(F.getID()));
			}
		
		@Override
		String getChartTitle() {
			return "Filters";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
	
			final CategoryAxis xAxis = new CategoryAxis(
					FXCollections.observableArrayList(count.keySet())
					);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final String t:this.count.keySet())
				{
				series1.getData().add(new XYChart.Data<String,Number>(t,this.count.count(t)));
				}
			bc.getData().add(series1);
	        bc.setTitle(this.getChartTitle());
	        xAxis.setLabel("Filter");       
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			if(ctx.isFiltered())
				{
				for(String f:ctx.getFilters()) {
					this.count.incr(f);
					}
				}
			else {
				this.count.incr(VCFConstants.PASSES_FILTERS_v4);
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
			final CategoryAxis xAxis = new CategoryAxis(
					FXCollections.observableArrayList(this.dict.getSequences().stream().
							map(S->S.getSequenceName()).
							filter(S->count.count(S)>0).
							collect(Collectors.toList()))
					);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final SAMSequenceRecord ssr: this.dict.getSequences())
				{
				long n = this.count.count(ssr.getSequenceName());
				if(n==0L) continue;
				series1.getData().add(new XYChart.Data<String,Number>(ssr.getSequenceName(),n));
				}
			bc.getData().add(series1);
	        bc.setTitle(this.getChartTitle());
	        xAxis.setLabel("Contig");       
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			this.count.incr(ctx.getContig());
			}
		}

	private class GenotypeTypeGenerator extends ChartGenerator
		{
		
		private final Map<String,Counter<GenotypeType>> sample2count;
		private final List<String> samples;
		GenotypeTypeGenerator(final List<String> samples) {
			this.samples = samples;
			this.sample2count = new HashMap<>(samples.size());
			for(final String s:samples)
				{
				this.sample2count.put(s, new Counter<>());
				}
			}
		
		@Override
		String getChartTitle() {
			return "Sample/Genotype-Type";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(
					FXCollections.observableArrayList(this.samples)
					);
			final StackedBarChart<String,Number> bc = 
		            new StackedBarChart<String,Number>(xAxis,yAxis);
		

			for(final GenotypeType gt:GenotypeType.values())
				{
				final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
				series1.setName(gt.name());
				for(final String sn:this.samples) {
					series1.getData().add(new XYChart.Data<String,Number>(
							sn,this.sample2count.get(sn).count(gt)
							));
					}
				bc.getData().add(series1);
				}
			
	        bc.setTitle(this.getChartTitle());
	        xAxis.setLabel("Sample");       
	        yAxis.setLabel("GT");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			for(final String sn:this.samples)
				{
				this.sample2count.get(sn).incr(ctx.getGenotype(sn).getType());
				}
			}
		}

	
	private final List<ChartGenerator> chartGenerators = new Vector<>();
	private class VariantContextRunner
		extends Task<Void>
		implements Runnable,Closeable
		
		{
		final CloseableIterator<VariantContext> iter;
		final VCFHeader header;
		long last_refresh;
		final VCFFileReader vcfFileReader;
		
		VariantContextRunner(VCFFileReader vcfFileReader)
			{
			this.vcfFileReader = vcfFileReader;
			this.header = vcfFileReader.getFileHeader();
			this.iter = vcfFileReader.iterator();
			}
		@Override
		protected Void call() throws Exception {
			refreshCharts();
			long last = System.currentTimeMillis();
			System.err.println("Call called");
			while(!VcfStatsJfx.this.stop && this.iter.hasNext()) {
				final VariantContext ctx = this.iter.next();
				try 
					{
					lock.lock();
					for(final ChartGenerator cg: chartGenerators) {
						if(!cg.isEnabled()) continue;
						cg.visit(ctx);
						}
					}
				finally
					{
					lock.unlock();
					}
				final long now = System.currentTimeMillis();
				if(now - last > refreshEverySeconds * 1000) {
					LOG.info("now-last");
					last=now;
					refreshCharts();
					}
				}
			System.err.println("Done");
			refreshCharts();
			return null;
			}
		@Override
		public void close() throws IOException {
			CloserUtil.close(this.iter);
			CloserUtil.close(this.vcfFileReader);
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
			if(args.size()!=1) {
				System.err.println("VCF missing");
				return -1;
				}
			final VCFFileReader vcfFileReader = new VCFFileReader(new File(args.get(0)));
			final VCFHeader header =  vcfFileReader.getFileHeader();
			chartGenerators.add(new VariantTypeGenerator());
			if(!header.getFilterLines().isEmpty())
				{
				chartGenerators.add(new FilterUsageGenerator(header.getFilterLines()));
				}
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict!=null && !dict.isEmpty()) {
				chartGenerators.add(new ContigUsageGenerator(dict));
				}
			
			if(header.hasGenotypingData())
				{
				chartGenerators.add(new GenotypeTypeGenerator(header.getGenotypeSamples()));
				for(final String sn:header.getSampleNamesInOrder()) {
					
					}
				}
			
			final VariantContextRunner runner = new VariantContextRunner(vcfFileReader);
			
			primaryStage.setOnShowing(AE->{
				LOG.info("start");
				new Thread(runner).start();
			});
			primaryStage.setOnCloseRequest(AE->{
				stop = true;
				CloserUtil.close(runner);
				LOG.info("end");
			});
			
			final TabPane tabPane = new TabPane();
			chartGenerators.removeIf(G->!G.isEnabled());
			for(final ChartGenerator g:chartGenerators) {
				tabPane.getTabs().add(g.makeTab());
				}
			final Screen scr = Screen.getPrimary();
			Scene scene  = new Scene(tabPane,
					scr.getBounds().getWidth()-100,
					scr.getBounds().getHeight()-100
					);
			primaryStage.setScene(scene);
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

	public static void main(String[] args) {
		Application.launch(args);
	}

}
