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
package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.Closeable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.math.RangeOfIntegers;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.JfxLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.concurrent.Task;
import javafx.embed.swing.SwingFXUtils;
import javafx.geometry.Side;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.SnapshotParameters;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.Chart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.PieChart;
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

```
 java -jar dist/vcfstatsjfx.jar input.vcf.gz
```

```
grep -E '(^#|missense)' input.vcf | java -jar dist/vcfstatsjfx.jar --stdin 
```

## Screenshot

  *  https://twitter.com/yokofakun/status/983280288238317568


![https://video.twimg.com/tweet_video/DaVQGvXXkAAMSBw.mp4](https://video.twimg.com/tweet_video/DaVQGvXXkAAMSBw.mp4 "animation")


END_DOC
*/
@Program(name="vcfstatsjfx",
description="GUI: VCF statistics",
biostars= {308310},
keywords={"vcf","stats","jfx"}
)
public class VcfStatsJfx extends JfxLauncher {
	private static final Logger LOG=Logger.build(VcfStatsJfx.class).make();
	private final DecimalFormat niceIntFormat = new DecimalFormat("###,###");
	private final ReentrantLock lock = new ReentrantLock();
	private ProgressBar progressBar;

	@Parameter(names={"-s","--seconds"},description="Refresh screen every 's' seconds")
	private int refreshEverySeconds = 15;
	@Parameter(names={"--max-concordance"},description="Max number of concordance to display. disable if <=0 ")
	private int max_condordance = 100;
	@Parameter(names={"-o","--output"},description="output Directory or zip file. If defined, the application will exit automatically")
	private File outputFile = null;
	@Parameter(names={"--trancheAffected"},description="tranches for the number of affected. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers affectedTranches = new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10,20,50,100,200,300,400,500,1000);
	@Parameter(names={"--stdin"},description="if there is no file argument. Read vcf from stdin instead of opening a FileOpen dialog")
	private boolean vcf_stdin = false;
	@Parameter(names={"-hr","--hr"},description="Ignore HOM_REF in genotype type.")
	private boolean ignore_HOM_REF = false;
	@Parameter(names={"-fgt","--fgt"},description="Ignore filtered **GENOTYPES**")
	private boolean ignore_filtered_genotypes = false;
	@Parameter(names={"-ncl","--norm-contig-length"},description="For the 'contig' Panel, normalize on contig length.")
	private boolean normalize_on_contig_length = false;
	
	private volatile boolean stop = false;
	private class ChartGenerator
		{
		Tab tab;
		boolean enabled = true;
		long nVariants = 0;
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
		void visit(final VariantContext ctx) {
			
			}
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
		
		protected boolean isCalled(final Genotype g) {
			if(g==null || !g.isCalled() || !g.isAvailable()) return false;
			if(ignore_filtered_genotypes && g.isFiltered()) return false;
			return g.getAlleles().stream().anyMatch(A->A.isCalled() && !A.isReference());
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
	        bc.setTitle(this.getChartTitle() +" N="+niceIntFormat.format(this.count.getTotal()));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("Variant Type");       
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			this.count.incr(ctx.getType());
			}
		}
	private class StructuralVariantTypeGenerator extends ChartGenerator
		{
		private final Counter<StructuralVariantType> count = new Counter<>();
		
		@Override
		String getChartTitle() {
			return "SV Type";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
	
			final CategoryAxis xAxis = new CategoryAxis(
					FXCollections.observableArrayList(
						Arrays.asList(StructuralVariantType.values()).
							stream().map(I->I.name()).
							collect(Collectors.toList()))
					);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final StructuralVariantType t:StructuralVariantType.values())
				{
				series1.getData().add(new XYChart.Data<String,Number>(t.name(),this.count.count(t)));
				}
			bc.getData().add(series1);
	        bc.setTitle(this.getChartTitle() +" N="+niceIntFormat.format(this.count.getTotal()));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("SV Type");       
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			final StructuralVariantType svt=ctx.getStructuralVariantType();
			if(svt==null) return;
			this.count.incr(svt);
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
			final XYChart.Series<String,Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final String t:this.count.keySet())
				{
				series1.getData().add(new XYChart.Data<String,Number>(t,this.count.count(t)));
				}
			bc.getData().add(series1);
	        bc.setTitle(this.getChartTitle()+" N="+niceIntFormat.format(nVariants));
	        bc.setLegendVisible(false);
	        xAxis.setTickLabelRotation(90);
	        xAxis.setLabel("Filters");       
	        yAxis.setLabel("Variant Count N="+niceIntFormat.format(nVariants));
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			++nVariants;
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
			return "Contigs" +(normalize_on_contig_length?" Normalized":"");
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
				series1.getData().add(new XYChart.Data<String,Number>(ssr.getSequenceName(),
						normalize_on_contig_length ?
								1e6 * n/(double)ssr.getSequenceLength() : 
								n
						));
				}
			bc.getData().add(series1);
	        bc.setTitle(this.getChartTitle()+ " N="+niceIntFormat.format(nVariants));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("Contig");       
	        yAxis.setLabel(normalize_on_contig_length? "Variant(s) per Mbp":"Count");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			this.nVariants++;
			this.count.incr(ctx.getContig());
			}
		}
	private class NumAltsGenerator extends ChartGenerator
		{
		private final Counter<Integer> count = new Counter<>();
		NumAltsGenerator() {
			}
		
		@Override
		String getChartTitle() {
			return "Num. ALT";
			}
		
		@Override
		Chart makeChart() {
			if(this.count.getCountCategories()==0) return null;
			final NumberAxis yAxis = new NumberAxis();
			final int max_count  = count.getMostFrequent();
			final List<String> L = new ArrayList<>(max_count+1);
			for(int x=0;x<=max_count;++x) {
				L.add(niceIntFormat.format(x));
			}
			
			
			final CategoryAxis xAxis = new CategoryAxis(
					FXCollections.observableArrayList(L)
					);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(int x=0;x<=max_count;++x)
				{
				long n = this.count.count(x);
				series1.getData().add(new XYChart.Data<String,Number>(
						niceIntFormat.format(x),
						n)
						);
				}
			bc.getData().add(series1);
	        bc.setTitle(this.getChartTitle()+ " N="+niceIntFormat.format(nVariants));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("N. ALT");       
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			this.nVariants++;
			this.count.incr(ctx.getAlternateAlleles().size());
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
				if(ignore_HOM_REF && gt.equals(GenotypeType.HOM_REF)) continue;
				final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
				series1.setName(gt.name());
				for(final String sn:this.samples) {
					series1.getData().add(new XYChart.Data<String,Number>(
							sn,this.sample2count.get(sn).count(gt)
							));
					}
				bc.getData().add(series1);
				}
			bc.setCategoryGap(1);
	        bc.setTitle(this.getChartTitle()+" N. Variants "+niceIntFormat.format(nVariants));
	        xAxis.setLabel("Sample (N="+niceIntFormat.format(this.samples.size())+")");
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("GT");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			nVariants++;
			for(final String sn:this.samples)
				{
				final Genotype gt= ctx.getGenotype(sn);
				if(gt!=null && ignore_filtered_genotypes && gt.isFiltered()) continue;
				this.sample2count.get(sn).incr(gt.getType());
				}
			}
		}
	
	private class GenotypesConcordanceGenerator extends ChartGenerator
		{
		private final Counter<String> count=new Counter<>();
		
		@Override
		String getChartTitle() {
			return "Genotype Condordance";
			}
		
		@Override
		Chart makeChart() {
			if(this.count.isEmpty() || max_condordance<1) return null;
			
			final List<String> bestPairs = this.count.keySetDecreasing().stream().limit(max_condordance).collect(Collectors.toList());
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(FXCollections.observableArrayList(bestPairs));
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			
			final BarChart<String,Number> bc = 
		            new BarChart<>(xAxis,yAxis);
		
			bc.getData().add(series1);
			
			for(final String sn:bestPairs)
				{
				series1.getData().add(new XYChart.Data<String,Number>(
					sn,100.0*(this.count.count(sn)/(double)nVariants)
					));
						
				}
			bc.setCategoryGap(1);
	        bc.setTitle(
	        		this.getChartTitle()+
	        		" N. Variants "+niceIntFormat.format(nVariants)+
	        		(this.count.getCountCategories()<max_condordance?"":" (Best "+niceIntFormat.format(max_condordance)+")")
	        		)
	        		;
	        xAxis.setLabel("Pair");
	        yAxis.setLabel("Percent of Variants (N.= "+niceIntFormat.format(nVariants)+")");
	        xAxis.setTickLabelRotation(90);
	        bc.setLegendVisible(false);
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			nVariants++;
			for(int i=0;i+1 < ctx.getNSamples();i++)
				{
				final Genotype g1 = ctx.getGenotype(i);
				if(!isCalled(g1)) continue;
				for(int j=i+1; j < ctx.getNSamples();j++)
					{
					final Genotype g2 = ctx.getGenotype(j);
					if(!isCalled(g2)) continue;
					if(!g1.sameGenotype(g2)) continue;
					
					this.count.incr(g1.getSampleName()+" ~ "+g2.getSampleName());
					}
				}
			}
		}

	
	private class NumberOfNoCallGenerator extends ChartGenerator
		{
		private final int nSamples ;
		private final Counter<Integer> count=new Counter<>();
		NumberOfNoCallGenerator(int nSamples) {
			this.nSamples = nSamples;
		}
		@Override
		String getChartTitle() {
			return "NoCalls";
			}
		
		@Override
		Chart makeChart() {
			if(this.count.isEmpty()) return null;
			int max_no_call = this.count.keySet().stream().mapToInt(G->G.intValue()).max().orElse(0);
			if(max_no_call==0) return null;
			final List<String> L = new ArrayList<>();
			for(int x=1;x<=max_no_call;++x) {
				if(this.count.count(x) ==0) continue;
				L.add(niceIntFormat.format(x));
			}
			if(L.isEmpty()) return null;
			
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(FXCollections.observableArrayList(L));
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			
			final BarChart<String,Number> bc = 
		            new BarChart<>(xAxis,yAxis);
		
			bc.getData().add(series1);
			
			for(int x=1;x<=max_no_call;++x)
				{
				final long N=this.count.count(x);
				if(N==0) continue;
				series1.getData().add(new XYChart.Data<String,Number>(
					niceIntFormat.format(x),
					N
					));
						
				}
			bc.setCategoryGap(1);
	        bc.setTitle(
	        		this.getChartTitle() +
	        		" N. Variants "+niceIntFormat.format(nVariants) +
	        		" N. Samples "+niceIntFormat.format(nSamples)
	        		)
	        		;
	        yAxis.setLabel("Number of Variants (num variants "+niceIntFormat.format(this.nVariants)+")");
	        xAxis.setLabel("Number of  NO_CALL per genotypes (num samples "+niceIntFormat.format(this.nSamples)+")");
	        xAxis.setTickLabelRotation(90);
	        bc.setLegendVisible(false);
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			nVariants++;
			final int n = (int)ctx.getGenotypes().stream().
					filter(G->ignore_filtered_genotypes?!G.isFiltered():true).
					filter(G->G.isNoCall()).
					count();
			if(n==0) return;
			this.count.incr(n);
			}
		}

	
	
	private class AffectedSamplesGenerator extends ChartGenerator
		{
		private final Counter<RangeOfIntegers.Range> count = new Counter<>();
		
		@Override
		String getChartTitle() {
			return "Sample Affected";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			
			
		
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final List<RangeOfIntegers.Range> L = new ArrayList<>( affectedTranches.getRanges());
			while(!L.isEmpty() && this.count.count(L.get(0))==0L)
				{
				L.remove(0);
				}
			while(!L.isEmpty() && this.count.count(L.get(L.size()-1))==0L)
				{
				L.remove(L.size()-1);
				}
			if(L.isEmpty()) return null;
			
			final CategoryAxis xAxis = new CategoryAxis(
					FXCollections.observableArrayList(
							L.stream().
							map(R->R.toString()).
							collect(Collectors.toList()))
					);
			
			for(final RangeOfIntegers.Range range : L)
				{
				series1.getData().add(new XYChart.Data<String,Number>(
						range.toString(),
						this.count.count(range)
						));
				}
			final BarChart<String,Number> bc = 
		            new BarChart<String,Number>(xAxis,yAxis);
			bc.getData().add(series1);
	        bc.setTitle(this.getChartTitle()+" N. variants ="+ niceIntFormat.format(nVariants));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("Num. Samples Affected N. variants ="+ niceIntFormat.format(nVariants));       
	        yAxis.setLabel("Sample Count");
	        xAxis.setTickLabelRotation(90);
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			++nVariants;
			this.count.incr(
				VcfStatsJfx.this.affectedTranches.getRange(
					(int)ctx.getGenotypes().stream().
						filter(G->isCalled(G)).
						count()	)
				);
			}
		}
	
	private enum PredType {
		start_lost,
		stop_lost,
		stop_gained,
		missense,
		frameshift,
		inframe_indel,
		protein_altering,
		synonymous_variant,
		coding_seqence,
		utr,
		coding_transcript,
		splicing_variant,
		exon_variant,
		intron_variant,
		transcript,
		gene,
		intergenic,
		regulatory_region,
		undefined
		};
	private final Map<SequenceOntologyTree.Term,PredType> soTerm2type = new HashMap<>();
	
	private class PredictionGenerator extends ChartGenerator
		{
		private final VcfTools tools;
		private final String sampleName;
		private PredType best=PredType.undefined;
		private final SequenceOntologyTree tree;
		
		private Counter<PredType> count = new  Counter<>();

		PredictionGenerator(final VcfTools tools,final String sampleName) {
			this.tools = tools;
			this.sampleName = sampleName;
			this.tree = tools.getSequenceOntologyTree();
			}
		PredictionGenerator(final VcfTools tools) {
			this(tools,null);
			}
		
		@Override
		String getChartTitle() {
			return "Predictions "+(this.sampleName==null?"":this.sampleName);
			}
		@Override
		Chart makeChart() {
			final ObservableList<PieChart.Data> pieChartData =
                FXCollections.observableArrayList(
                		this.count.keySet().stream().
                			map(T-> new PieChart.Data(
                					T.name()+" "+niceIntFormat.format(this.count.count(T)), 
                					this.count.count(T))).
                			collect(Collectors.toList())
                			);
		           
		    final PieChart chart = new PieChart(pieChartData);
		    chart.setTitle(getChartTitle()+ " N="+niceIntFormat.format(nVariants));
		    chart.setLabelLineLength(10);
		    chart.setLabelsVisible(true);
		    chart.setLegendSide(Side.LEFT);
		    return chart;
			}
		private PredType toPredType(final SequenceOntologyTree.Term term) {
			if(term.isChildrenOf(tree.getTermByAcn("SO:0002012")))
				{
				return PredType.start_lost;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001578")))
				{
				return PredType.stop_lost;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001587")))
				{
				return PredType.stop_gained;
				}
			
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001583")))
				{
				return PredType.missense;
				}
			
			
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001819")))
				{
				return PredType.synonymous_variant;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001589")))
				{
				return PredType.frameshift;
				}			
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001820")))
				{
				return PredType.inframe_indel;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001818")))
				{
				LOG.warn("general "+term);
				return PredType.protein_altering;
				}
			
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001580")))
				{
				LOG.warn("general "+term);
				return PredType.coding_seqence;
				}
			
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001622")))
				{
				return PredType.utr;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001568")))
				{
				return PredType.splicing_variant;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001791")))
				{
				return PredType.exon_variant;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001627")))
				{
				return PredType.intron_variant;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001968")))
				{
				return PredType.coding_transcript;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001576")))
				{
				LOG.warn("general "+term);
				return PredType.transcript;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001564")))
				{
				LOG.warn("general "+term);
				return PredType.gene;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001566")))
				{
				return PredType.regulatory_region;
				}
			if(term.isChildrenOf(tree.getTermByAcn("SO:0001628")) ||
				term.isChildrenOf(tree.getTermByAcn("SO:0000605")))
				{
				return PredType.intergenic;
				}
			
			
			LOG.warn("cannot decode "+term);
			return PredType.undefined;
			}
		
		private void best(SequenceOntologyTree.Term term) {
			final PredType curr  ;
			if(soTerm2type.containsKey(term))
				{
				curr = soTerm2type.get(term);
				}
			else
				{
				curr = toPredType(term);
				soTerm2type.put(term, curr);
				}
			
			if(curr.compareTo(best)<0) this.best=curr;
			}
		
		@Override
		void visit(final VariantContext ctx)
			{
			this.best=PredType.undefined;
			if(this.sampleName!=null && !isCalled(ctx.getGenotype(this.sampleName)))
				{
				return;
				}
				
			this.nVariants++;
			this.tools.getAnnPredictions(ctx).stream().flatMap(P->P.getSOTerms().stream()).forEach(T->best(T));
			this.tools.getVepPredictions(ctx).stream().flatMap(P->P.getSOTerms().stream()).forEach(T->best(T));
			this.count.incr(this.best);
			this.best=PredType.undefined;
			}
		}
	
	
	private final List<ChartGenerator> chartGenerators = new Vector<>();
	private class VariantContextRunner
		extends Task<Void>
		implements Runnable,Closeable
		
		{
		final VcfIterator iter;
		
		VariantContextRunner(final VcfIterator vcfFileReader)
			{
			this.iter = vcfFileReader;
			}
		@Override
		protected Void call() throws Exception {
			refreshCharts();
			long last = -1L;
			
			Platform.runLater(()->{
				progressBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);
				});
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
				if(last<0L || now - last > refreshEverySeconds * 1000) {
					last=now;
					refreshCharts();
					}
				}
			
			refreshCharts();
			
			Platform.runLater(()->{
				progressBar.setProgress(1.0);
				});
			
			if(outputFile!=null)
	        	{
	        
				LOG.info("saving as "+outputFile);
        		new java.util.Timer().schedule( 
    		        new java.util.TimerTask() {
    		            @Override
    		            public void run() {
    		            Platform.runLater(()->{
		        		 try
		        		 	{
		        			saveImagesAs(outputFile);
		        		 	}
		        		 catch(IOException err)
		        		 	{
		        			LOG.error(err);
		        			System.exit(-1);
		        		 	}
		        		 Platform.exit();
    		            	});}
    		            },refreshEverySeconds*1000);
	        	 	
	        	}
			return null;
			}
		@Override
		public void close() throws IOException {
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
			final VcfIterator vcfIterator;
			if(args.size()==1) {
				vcfIterator = VCFUtils.createVcfIteratorFromFile(new File(args.get(0)));
				}
			else if(args.isEmpty() && this.vcf_stdin)
				{
				vcfIterator = VCFUtils.createVcfIteratorFromInputStream(System.in);
				}
			else if(args.isEmpty())
				{
				final FileChooser fc = new FileChooser();
				final File f = fc.showOpenDialog(null);
				if(f==null) return -1;
				vcfIterator = VCFUtils.createVcfIteratorFromFile(f);
				}
			else
				{
				LOG.error("illegal number of arguments.");
				return -1;
				}
			final VCFHeader header =  vcfIterator.getHeader();
			chartGenerators.add(new VariantTypeGenerator());
			if(!header.getFilterLines().isEmpty())
				{
				chartGenerators.add(new FilterUsageGenerator(header.getFilterLines()));
				}
			chartGenerators.add(new NumAltsGenerator());
			
			if(header.getInfoHeaderLine(VCFConstants.SVTYPE)!=null) {
				chartGenerators.add(new StructuralVariantTypeGenerator());
				}
			
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict!=null && !dict.isEmpty()) {
				chartGenerators.add(new ContigUsageGenerator(dict));
				}
			final VcfTools tools = new VcfTools(header);
			boolean hasPred = tools.getVepPredictionParser().isValid() || 
							tools.getAnnPredictionParser().isValid();
			
			if(hasPred)
				{
				chartGenerators.add(new PredictionGenerator(tools));
				}
			
			if(header.hasGenotypingData())
				{
				if(max_condordance>0) {
					chartGenerators.add(new GenotypesConcordanceGenerator());
				}
				
				
				chartGenerators.add(new GenotypeTypeGenerator(header.getGenotypeSamples()));
				chartGenerators.add(new AffectedSamplesGenerator());
				chartGenerators.add(new NumberOfNoCallGenerator(header.getNGenotypeSamples()));
				for(final String sn:header.getSampleNamesInOrder()) {
					if(hasPred) chartGenerators.add(new PredictionGenerator(tools,sn));
					}
				}
			
			final VariantContextRunner runner = new VariantContextRunner(vcfIterator);
			
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
		    MenuItem item=new MenuItem("Save Current image as...");
		    item.setOnAction(AE->{doMenuSaveCurrentImage(tabPane.getSelectionModel().getSelectedItem());});
		    menu.getItems().add(item);
		    item=new MenuItem("Save All images as...");
		    item.setOnAction(AE->doMenuSaveAllImages());
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
            LOG.error(e);
            final Alert alert=new Alert(AlertType.ERROR,e.getMessage());
            alert.showAndWait();
        	}
		}
	
	private void doMenuSaveCurrentImage(final Tab tab) {
		if(tab==null) return;
		final Node content =  tab.getContent();
	 	if(content==null) return;
		final FileChooser fc = new FileChooser();
		final File file = fc.showSaveDialog(null);
    	if(file==null) return;
    	try {
    		final WritableImage image = content.snapshot(new SnapshotParameters(), null);
			
    		final String format = file.getName().toLowerCase().endsWith("png")?"png":"jpg";
			ImageIO.write(SwingFXUtils.fromFXImage(image, null), format, file);
        	}
    	catch (final IOException e) {
            LOG.error(e);
            final Alert alert=new Alert(AlertType.ERROR,e.getMessage());
            alert.showAndWait();
        	}
		}
	
	private void saveImagesAs(final File out) throws IOException {
	if(out.getName().endsWith(".zip")) {
		final FileOutputStream fout = new FileOutputStream(out);
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

	public static void main(String[] args) {
		Application.launch(args);
	}

}
