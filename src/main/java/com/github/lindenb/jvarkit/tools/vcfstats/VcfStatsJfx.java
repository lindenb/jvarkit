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
package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.chart.BarChart;
import com.github.lindenb.jvarkit.chart.CategoryAxis;
import com.github.lindenb.jvarkit.chart.Chart;
import com.github.lindenb.jvarkit.chart.NumberAxis;
import com.github.lindenb.jvarkit.chart.RExporter;
import com.github.lindenb.jvarkit.chart.StackedBarChart;
import com.github.lindenb.jvarkit.chart.XYChart;
import com.github.lindenb.jvarkit.math.RangeOfDoubles;
import com.github.lindenb.jvarkit.math.RangeOfIntegers;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

/*
BEGIN_DOC

## Examples

```
 java -jar dist/vcfstatsjfx.jar input.vcf.gz  | Rscript -
 
```


## Screenshot

  *  https://twitter.com/yokofakun/status/983280288238317568


![https://video.twimg.com/tweet_video/DaVQGvXXkAAMSBw.mp4](https://video.twimg.com/tweet_video/DaVQGvXXkAAMSBw.mp4 "animation")


## History

   *  removed JFX/gui because openjdk doesn't support jfx :-(

END_DOC
*/
@Program(name="vcfstatsjfx",
description="VCF statistics",
biostars= {308310,353051},
keywords={"vcf","stats"}
)
public class VcfStatsJfx extends Launcher {
	private static final Logger LOG=Logger.build(VcfStatsJfx.class).make();
	private final DecimalFormat niceIntFormat = new DecimalFormat("###,###");

	@Parameter(names={"-s","--seconds"},description="Save Rscript screen every 's' seconds, if output was defined.")
	private int refreshEverySeconds = 15;
	@Parameter(names={"--max-concordance"},description="Max number of concordance to display. disable if <=0 ")
	private int max_condordance = 100;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"--trancheIndelSize"},description="tranches for the Indel size "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers indelTranches =new RangeOfIntegers(-1000,-100,-50,-20,-15,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,15,20,50,100,1000);
	@Parameter(names={"--trancheAffected"},description="tranches for the number of affected. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers affectedTranches = new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10,20,50,100,200,300,400,500,1000);
	@Parameter(names={"--stdin"},description="if there is no file argument. Read vcf from stdin instead of opening a FileOpen dialog")
	private boolean vcf_stdin = false;
	@Parameter(names={"-fgt","--fgt"},description="Ignore filtered **GENOTYPES**")
	private boolean ignore_filtered_genotypes = false;
	@Parameter(names={"-ncl","--norm-contig-length"},description="For the 'contig' Panel, normalize on contig length.")
	private boolean normalize_on_contig_length = false;
	@Parameter(names={"--altering","--damaging"},description="For Prediction, just display children of SO:0001818 ( protein_altering_variant )")
	private boolean only_protein_altering_variants = false;
	@Parameter(names={"--predictions-per-sample","-pps"},description="Show Predictions per sample.")
	private boolean enable_predictions_per_sample = false;
	@Parameter(names={"--prefix"},description="Title Prefix")
	private String titlePrefix="";

	
	private abstract class ChartGenerator
		{
		boolean enabled = true;
		long nVariants = 0;
		boolean isEnabled() {
			return enabled;
			}
		abstract String getChartTitle();
		
		void title(final Chart c,final String s) {
			c.setTitle(
					StringUtil.isBlank(titlePrefix)?
					s:
					titlePrefix+ " : " + s
					);
			}
		
		void visit(final VariantContext ctx) {
			
			}
		abstract Chart makeChart();
		
		
		// public String getFilename() { return getTabTitle().replaceAll("[^A-Za-z_0-9]+","")+".png";}
		
		protected boolean isCalled(final Genotype g) {
			if(g==null || !g.isCalled() || !g.isAvailable()) return false;
			if(ignore_filtered_genotypes && g.isFiltered()) return false;
			return g.getAlleles().stream().anyMatch(A->A.isCalled() && !A.isReference());
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
					Arrays.asList(VariantContext.Type.values()).
							stream().map(I->I.name()).
							collect(Collectors.toList())
					);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final VariantContext.Type t:VariantContext.Type.values())
				{
				series1.getData().add(new XYChart.Data<String,Number>(t.name(),this.count.count(t)));
				}
			bc.getData().add(series1);
			title(bc,this.getChartTitle() +" N="+niceIntFormat.format(this.count.getTotal()));
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
						Arrays.asList(StructuralVariantType.values()).
							stream().map(I->I.name()).
							collect(Collectors.toList())
					);
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final StructuralVariantType t:StructuralVariantType.values())
				{
				series1.getData().add(new XYChart.Data<String,Number>(t.name(),this.count.count(t)));
				}
			bc.getData().add(series1);
			title(bc,this.getChartTitle() +" N="+niceIntFormat.format(this.count.getTotal()));
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
	
			final CategoryAxis xAxis = new CategoryAxis( count.keySet());
			final XYChart.Series<String,Number> series1 = new XYChart.Series<>();
			final BarChart<String,Number> bc = 
			            new BarChart<String,Number>(xAxis,yAxis);
			
			for(final String t:this.count.keySet())
				{
				series1.getData().add(new XYChart.Data<String,Number>(t,this.count.count(t)));
				}
			bc.getData().add(series1);
			title(bc,this.getChartTitle()+" N="+niceIntFormat.format(nVariants));
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
					this.dict.getSequences().stream().
							map(S->S.getSequenceName()).
							filter(S->count.count(S)>0).
							collect(Collectors.toList())
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
			title(bc,this.getChartTitle()+ " N="+niceIntFormat.format(nVariants));
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
			return "Num. ALTernate Alleles";
			}
		
		@Override
		Chart makeChart() {
			if(this.count.getCountCategories()==0) return null;
			final NumberAxis yAxis = new NumberAxis();
			final int max_count  = count.keySet().stream().mapToInt(I->I.intValue()).max().orElse(1);
			final List<String> L = new ArrayList<>(max_count+1);
			for(int x=0;x<=max_count;++x) {
				L.add(niceIntFormat.format(x));
			}
			
			
			final CategoryAxis xAxis = new CategoryAxis(L);
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
			title(bc,this.getChartTitle()+ " N="+niceIntFormat.format(nVariants));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("N. ALT");       
	        yAxis.setLabel("Count");
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			this.nVariants++;
			final int n_alts = ctx.getAlternateAlleles().size();
			this.count.incr(n_alts);			}
		}
	
	private abstract class AbstractAverageGTGenerator extends ChartGenerator
		{
		protected class Depth{long sum=0L;long count=0L;
		double avg() { return count==0L?0.0:sum/(double)count;}
		}
		private final Map<String,Depth> sample2count;
		protected AbstractAverageGTGenerator(final List<String> samples) {
			this.sample2count = new HashMap<>(samples.size());
			for(final String s:samples)
				{
				this.sample2count.put(s, new Depth());
				}
			}
		@Override
		abstract String getChartTitle();
		
		@Override
		Chart makeChart() {
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final List<String> categories = new ArrayList<>(this.sample2count.size());
			for(final String sn: this.sample2count.keySet().
					stream().
					sorted((A,B)->Double.compare(sample2count.get(A).avg(),sample2count.get(B).avg())).
					collect(Collectors.toList())
					)
				{
				final Depth dp = this.sample2count.get(sn);
				final String key = sn+" "+niceIntFormat.format(dp.count);
				categories.add(key);
				series1.getData().add(new XYChart.Data<String,Number>(
						key,
						dp.avg()
						));
				}
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(categories);
			final BarChart<String,Number> bc = new BarChart<>(xAxis,yAxis);
	
			bc.getData().add(series1);
			bc.setCategoryGap(1);
			title(bc,getChartTitle()+ " N="+niceIntFormat.format(nVariants));
		    bc.setLegendVisible(false);
		    bc.setVerticalGridLinesVisible(false);
	        xAxis.setLabel("Sample N="+niceIntFormat.format(this.sample2count.size()));
	        yAxis.setLabel(getYLabel());
	        xAxis.setTickLabelRotation(90);
		    return bc;
			}
		
		abstract String getYLabel();
		
		abstract void visitGenotype(final Genotype gt);
		
		@Override
		void visit(final VariantContext ctx) {
			nVariants++;
			for(final Genotype gt:ctx.getGenotypes())
				{
				if(gt.isFiltered() && ignore_filtered_genotypes) continue;
				visitGenotype(gt);
				}
			}
	
		}

	private class AverageDepthGenerator extends AbstractAverageGTGenerator
		{
		AverageDepthGenerator(final List<String> samples) {
			super(samples);
			}
		@Override
		String getChartTitle() {
			return "Depth/Sample";
			}
		@Override
		String getYLabel() {
			return "Average Depth";
			}
		@Override
		void visitGenotype(final Genotype gt) {
			if(!gt.hasDP()) return;
			final Depth dp=super.sample2count.get(gt.getSampleName());
			dp.count++;
			dp.sum+=gt.getDP();
			}
		}
	
	private class AverageGQGenerator extends AbstractAverageGTGenerator
		{
		AverageGQGenerator(final List<String> samples) {
			super(samples);
			}
		@Override
		String getChartTitle() {
			return "Genotype Quality/Sample";
			}
		@Override
		String getYLabel() {
			return "Average GQ";
			}
		@Override
		void visitGenotype(final Genotype gt) {
			if(!gt.hasGQ()) return;
			final Depth dp=super.sample2count.get(gt.getSampleName());
			dp.count++;
			dp.sum+=gt.getGQ();
			}
		}
	private class AveragePhasedGTGenerator extends AbstractAverageGTGenerator
		{
		AveragePhasedGTGenerator(final List<String> samples) {
			super(samples);
			}
		@Override
		String getChartTitle() {
			return "Fraction of Phased Genotype/Sample";
			}
		@Override
		String getYLabel() {
			return "Fraction of Phased Genotypes";
			}
		@Override
		void visitGenotype(final Genotype gt) {
			if(!gt.hasGQ()) return;
			if(gt.getPloidy()<2) return;
			final Depth dp=super.sample2count.get(gt.getSampleName());
			dp.count++;
			dp.sum+=gt.isPhased()?1:0;
			}
		}
	
	/** per sample/ genotype type */
	private abstract class AbstractGenotypeTypeGenerator extends ChartGenerator
		{
		private final Map<String,Counter<GenotypeType>> sample2count;
		private final List<String> samples;
		AbstractGenotypeTypeGenerator(final List<String> samples) {
			this.samples = samples;
			this.sample2count = new HashMap<>(samples.size());
			for(final String s:samples)
				{
				this.sample2count.put(s, new Counter<>());
				}
			}
		
		@Override
		String getChartTitle() {
			return "Sample/Genotype-Type"+(IsIgnoringHomRef()?" (HOM_REF Excluded)":"");
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(
					IsIgnoringHomRef()?
						this.samples.stream().sorted((S1,S2)->
							{
							return Long.compare(
									sample2count.get(S1).getTotal(),
									sample2count.get(S2).getTotal()
									);
							}).collect(Collectors.toList()):
						this.samples)
					;
			final StackedBarChart<String,Number> bc = 
		            new StackedBarChart<String,Number>(xAxis,yAxis);
		

			for(final GenotypeType gt:GenotypeType.values())
				{
				if(IsIgnoringHomRef() && gt.equals(GenotypeType.HOM_REF)) continue;
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
			title(bc,this.getChartTitle()+" N. Variants "+niceIntFormat.format(nVariants));
	        xAxis.setLabel("Sample (N="+niceIntFormat.format(this.samples.size())+")");
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("GT");
	        return bc;
			}
		
		abstract boolean IsIgnoringHomRef();
		
		@Override
		void visit(final VariantContext ctx) {
			nVariants++;
			for(final String sn:this.samples)
				{
				final Genotype gt= ctx.getGenotype(sn);
				if(gt!=null && ignore_filtered_genotypes && gt.isFiltered()) continue;
				if(gt!=null && IsIgnoringHomRef() && gt.isHomRef()) continue;
				this.sample2count.get(sn).incr(gt.getType());
				}
			}
		}
	
	
	/** per sample/ genotype type */
	private class DeNovoGenerator extends ChartGenerator
		{
		private final Counter<String> sample2hiConfDeNovo = new Counter<>();
		private final Counter<String> sample2loConfDeNovo = new Counter<>();
		
		
		@Override
		String getChartTitle() {
			return "GATK DeNovo/Sample";
			}
		
		@Override
		Chart makeChart() {
			final Set<String> deNovoSamples = new TreeSet<>(sample2hiConfDeNovo.keySet());
			deNovoSamples.addAll(sample2loConfDeNovo.keySet());
			
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(
					deNovoSamples
					);
			final StackedBarChart<String,Number> bc = 
		            new StackedBarChart<String,Number>(xAxis,yAxis);
		
			
			for(int y=0;y< 2;++y)
				{
				Counter<String> counter = (y==0?this.sample2hiConfDeNovo:this.sample2loConfDeNovo);
				final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
				series1.setName(y==0?"hiConfDeNovo":"loConfDeNovo");
				for(final String sn:deNovoSamples) {
					series1.getData().add(new XYChart.Data<String,Number>(
							sn,counter.count(sn)
							));
					}
				bc.getData().add(series1);
				}
			bc.setCategoryGap(1);
			title(bc,this.getChartTitle());
	        xAxis.setLabel("Sample");
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("DeNovo Variants");
	        return bc;
			}
		
		
		@Override
		void visit(final VariantContext ctx) {
			this.nVariants++;
			ctx.getAttributeAsList("hiConfDeNovo").stream().map(O->String.valueOf(O)).forEach(S->sample2hiConfDeNovo.incr(S));
			ctx.getAttributeAsList("loConfDeNovo").stream().map(O->String.valueOf(O)).forEach(S->sample2loConfDeNovo.incr(S));
			}
		}

	

	/** chrX/chrY */
	private class MaleFemaleGenerator extends ChartGenerator
		{
		private class CountXY
			{
			final String name;
			CountXY(final String name) {
				this.name = name;
				}
			
			long x = 0L;
			long y = 0L;
			
			double score() {
				double n=x+y;
				if(n==0.0) return 0;
				return x/n;
				}
			}
		private final Map<String,CountXY> sample2count = new HashMap<>();
		
		@Override
		String getChartTitle() {
			return "ChrX/ChrY";
			}
		
		@Override
		Chart makeChart() {
			final List<String> orderedNames = this.sample2count.values().
					stream().
					sorted((A,B)->Double.compare(A.score(), B.score())).
					map(A->A.name).
					collect(Collectors.toList());
			
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(
					orderedNames
					);
			final StackedBarChart<String,Number> bc = 
		            new StackedBarChart<String,Number>(xAxis,yAxis);
		
			
			for(int y=0;y< 2;++y)
				{
				final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
				series1.setName(y==0?"chrX":"chrY");
				for(final String sn:orderedNames) {
					final CountXY countxy = this.sample2count.get(sn);
					series1.getData().add(new XYChart.Data<String,Number>(
							sn,(y==0?countxy.x:countxy.y)
							));
					}
				bc.getData().add(series1);
				}
			bc.setCategoryGap(1);
			title(bc,this.getChartTitle());
	        xAxis.setLabel("Sample");
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("ChrX/ChrY");
	        return bc;
			}
		
		
		@Override
		void visit(final VariantContext ctx) {
			this.nVariants++;
			final int where ;
			if(ctx.getContig().equals("X") || ctx.getContig().equals("chrX"))
				{
				where = 23;
				}
			else if(ctx.getContig().equals("Y") || ctx.getContig().equals("chrY"))
				{
				where = 24;
				}
			else
				{
				return;
				}
			for(final Genotype gt:ctx.getGenotypes())
				{
				if(ignore_filtered_genotypes && gt.isFiltered()) continue;
				CountXY count = this.sample2count.get(gt.getSampleName());
				if(count==null) {
					count = new CountXY(gt.getSampleName());
					this.sample2count.put(gt.getSampleName(),count);
					}
				switch(where) {
					case 23: count.x++;break;
					case 24: count.y++;break;
					default: break;
					}
				}
			}
		}
	
	
	/** AD allele depth */
	private abstract class AlleleDepthGenerator extends ChartGenerator
		{
		private final RangeOfDoubles rangeOfDoubles = new RangeOfDoubles(new double[] {0.2,0.3,0.4,0.6,0.7,0.8});
		private final Map<String,Counter<RangeOfDoubles.Range>> sample2count = new TreeMap<>();
		
		
		@Override
		String getChartTitle() {
			return "AL/(REF+ALT) (SNP Diallelic "+getGenotypeType()+" Genotypes)";
			}
		
		
		@Override
		Chart makeChart() {
			
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(
						this.sample2count.keySet().stream().sorted((S1,S2)->
							{
							return Long.compare(
									sample2count.get(S1).getTotal(),
									sample2count.get(S2).getTotal()
									);
							}).collect(Collectors.toList())
					);
			final StackedBarChart<String,Number> bc = 
		            new StackedBarChart<String,Number>(xAxis,yAxis);
		
			
			for(final RangeOfDoubles.Range range:rangeOfDoubles.getRanges())
				{
				final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
				series1.setName(range.toString());
				for(final String sn:this.sample2count.keySet()) {
					final Counter<RangeOfDoubles.Range> counter = this.sample2count.get(sn);
					final long n  = counter==null?0L:counter.count(range);
					
					series1.getData().add(new XYChart.Data<String,Number>(
							sn,n
							));
					}
				bc.getData().add(series1);
				}
			bc.setCategoryGap(1);
			title(bc,this.getChartTitle());
	        xAxis.setLabel("Sample");
	        bc.setVerticalGridLinesVisible(false);
	        xAxis.setTickLabelRotation(90);
	        yAxis.setLabel("ALT/(ALT+REF)");
	        return bc;
			}
		
		protected abstract GenotypeType getGenotypeType();
		
		@Override
		void visit(final VariantContext ctx) {
			if(!ctx.isBiallelic()) return ;
			if(!ctx.isSNP()) return ;
			this.nVariants++;
			for(final Genotype gt: ctx.getGenotypes())
				{
				if(!gt.hasAD()) continue;
				if(!getGenotypeType().equals(gt.getType())) continue;
				final int ads[] = gt.getAD();
				if(ads==null || ads.length!=2 || ads[0]+ads[1]<=0.0) continue;
				final RangeOfDoubles.Range range = rangeOfDoubles.getRange((double)ads[1]/(double)(ads[0]+ads[1]));
				Counter<RangeOfDoubles.Range> counter = this.sample2count.get(gt.getSampleName());
				if(counter==null) {
					counter = new Counter<>();
					this.sample2count.put(gt.getSampleName(),counter);
					}
				counter.incr(range);
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
			final CategoryAxis xAxis = new CategoryAxis(bestPairs);
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
			title(bc,
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
			final CategoryAxis xAxis = new CategoryAxis(L);
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
			title(bc,
	        		this.getChartTitle() +
	        		" N. Variants "+niceIntFormat.format(nVariants) +
	        		" N. Samples "+niceIntFormat.format(nSamples)
	        		)
	        		;
	        yAxis.setLabel("Number of Variants (num variants "+niceIntFormat.format(this.nVariants)+")");
	        xAxis.setLabel("Number of  NO_CALL Genotypes per Variant (num samples "+niceIntFormat.format(this.nSamples)+")");
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

	
	/*****************************************************************************/
	private class AffectedSamplesGenerator extends ChartGenerator
		{
		private final Counter<RangeOfIntegers.Range> count = new Counter<>();
		private final int nSamples;
		AffectedSamplesGenerator(final int nSamples) {
			this.nSamples = nSamples;
			}
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
							L.stream().
							map(R->R.toString()).
							collect(Collectors.toList()))
					;
			
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
			title(bc,this.getChartTitle()+" N. variants ="+ niceIntFormat.format(nVariants));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("Num. Samples Affected N. samples ="+ niceIntFormat.format(nSamples));       
	        yAxis.setLabel("Variant Count N="+ niceIntFormat.format(nVariants));
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
	
	
	/*****************************************************************************/
	private class VariantSizesGenerator extends ChartGenerator
		{
		private final Counter<RangeOfIntegers.Range> count = new Counter<>();
		
		@Override
		String getChartTitle() {
			return "Variant lengths";
			}
		
		@Override
		Chart makeChart() {
			final NumberAxis yAxis = new NumberAxis();		
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			final List<RangeOfIntegers.Range> L = new ArrayList<>( indelTranches.getRanges());
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
							L.stream().
							map(R->R.toString()).
							collect(Collectors.toList()))
					;
			
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
			title(bc,this.getChartTitle()+" N. variants ="+ niceIntFormat.format(nVariants));
	        bc.setLegendVisible(false);
	        xAxis.setLabel("Variant Size.");       
	        yAxis.setLabel("Variant Count.");
	        xAxis.setTickLabelRotation(90);
	        return bc;
			}
		
		@Override
		void visit(final VariantContext ctx) {
			++nVariants;
			this.count.incr(
				VcfStatsJfx.this.indelTranches.getRange(1 + ctx.getEnd()-ctx.getStart())
				);
			}
		}


	
	private class PredictionGenerator extends ChartGenerator
		{
		private final VcfTools tools;
		//private Counter<String> count = new Counter<>();
		private final String sampleName;
		private SequenceOntologyTree.Term bestTerm = null;
		private final SequenceOntologyTree tree;
		private final SequenceOntologyTree.DamagingComparator damagingComparator;
		private final SequenceOntologyTree.Term protein_altering_variant_term;
		private Counter<SequenceOntologyTree.Term> countTerms = new  Counter<>();

		PredictionGenerator(final VcfTools tools,final String sampleName) {
			this.tools = tools;
			this.sampleName = sampleName;
			this.tree = tools.getSequenceOntologyTree();
			this.damagingComparator =new SequenceOntologyTree.DamagingComparator(this.tree);
			this.protein_altering_variant_term = this.tree.getTermByAcn("SO:0001818");
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
			if(this.countTerms.isEmpty()) return null;
			final XYChart.Series<String, Number> series1 = new XYChart.Series<>();
			
			
			for(final SequenceOntologyTree.Term t: this.countTerms.keySetDecreasing())
				{
				series1.getData().add(new XYChart.Data<String,Number>(
						t.getLabel(),
						this.countTerms.count(t)
						));
						
				}
			final NumberAxis yAxis = new NumberAxis();
			final CategoryAxis xAxis = new CategoryAxis(this.countTerms.keySet().stream().map(T->T.getLabel()).collect(Collectors.toList()));
			final BarChart<String,Number> bc = new BarChart<>(xAxis,yAxis);

			bc.getData().add(series1);
			bc.setCategoryGap(1);
			title(bc,getChartTitle()+ " N="+niceIntFormat.format(nVariants));
		    bc.setLegendVisible(false);
	        
	        xAxis.setLabel("SOTerm");
	        yAxis.setLabel("Number of Variants (N.= "+niceIntFormat.format(nVariants)+")");
	        xAxis.setTickLabelRotation(90);
		    return bc;
			}
		
		private void best(SequenceOntologyTree.Term term) {
			if(this.bestTerm==null || this.damagingComparator.compare(term, this.bestTerm)<0)
				{
				if(bestTerm!=null) {
					//System.err.println(term.getLabel()+" is better than "+bestTerm.getLabel());
				}
				this.bestTerm = term;
				}
			}
		

		private boolean accept(final SequenceOntologyTree.Term t) {
			if(!only_protein_altering_variants) return true;
			return t.isChildrenOf(protein_altering_variant_term);
		}
		
		@Override
		void visit(final VariantContext ctx)
			{
			this.bestTerm = null;
			if(this.sampleName!=null && !isCalled(ctx.getGenotype(this.sampleName)))
				{
				return;
				}
			this.nVariants++;
			this.tools.getAnnPredictions(ctx).stream().flatMap(P->P.getSOTerms().stream()).filter(T->accept(T)).forEach(T->best(T));
			this.tools.getVepPredictions(ctx).stream().flatMap(P->P.getSOTerms().stream()).filter(T->accept(T)).forEach(T->best(T));
			if(this.bestTerm==null) return;
			this.countTerms.incr(this.bestTerm);
			this.bestTerm=null;
			}
		}
	
	
	private final List<ChartGenerator> chartGenerators = new ArrayList<>();
	
	
	private void save() throws IOException {
		try (final PrintWriter pw = super.openFileOrStdoutAsPrintWriter(this.outputFile)) {
		final RExporter exporter = new RExporter();
		for(final ChartGenerator cg: this.chartGenerators) {
			LOG.info("saving "+cg.getChartTitle());
			if(!cg.isEnabled()) continue;
			final Chart chart = cg.makeChart();
			if(chart==null) continue;
			exporter.exportToR(pw, chart);
			}
		pw.flush();
		}
	}


	
	@Override
	public int doWork(final List<String> args) {
		VCFIterator iter = null;
		try {
			iter =  super.openVCFIterator(oneFileOrNull(args));
			
			final VCFHeader header =  iter.getHeader();
			this.chartGenerators.add(new VariantTypeGenerator());
			if(!header.getFilterLines().isEmpty())
				{
				this.chartGenerators.add(new FilterUsageGenerator(header.getFilterLines()));
				}
			this.chartGenerators.add(new NumAltsGenerator());
			this.chartGenerators.add(new VariantSizesGenerator());
			
			if(header.getInfoHeaderLine(VCFConstants.SVTYPE)!=null) {
				this.chartGenerators.add(new StructuralVariantTypeGenerator());
				}
			
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict!=null && !dict.isEmpty()) {
				this.chartGenerators.add(new ContigUsageGenerator(dict));
				if(SequenceDictionaryUtils.hasXY(dict))
					{
					chartGenerators.add(new MaleFemaleGenerator());
					}
				}
			final VcfTools tools = new VcfTools(header);
			boolean hasPred = tools.getVepPredictionParser().isValid() || 
							tools.getAnnPredictionParser().isValid();
			
			
			if(header.getInfoHeaderLine("hiConfDeNovo")!=null && 
				header.getInfoHeaderLine("loConfDeNovo")!=null)
				{
				chartGenerators.add(new DeNovoGenerator());
				}
			
			if(header.getNGenotypeSamples()>0  &&
				header.getFormatHeaderLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS)!=null)
				{
				if(header.getNGenotypeSamples()>1) chartGenerators.add(new AlleleDepthGenerator() {
					@Override
					protected GenotypeType getGenotypeType() {
						return GenotypeType.HOM_REF;
						}
					});
				chartGenerators.add(new AlleleDepthGenerator() {
					@Override
					protected GenotypeType getGenotypeType() {
						return GenotypeType.HET;
						}
					});
				}
			
			
			
			if(hasPred)
				{
				chartGenerators.add(new PredictionGenerator(tools));
				}
			
			if(header.hasGenotypingData())
				{
				if(max_condordance>0) {
					chartGenerators.add(new GenotypesConcordanceGenerator());
				}
				
				
				chartGenerators.add(new AbstractGenotypeTypeGenerator(header.getGenotypeSamples()){
					@Override
					boolean IsIgnoringHomRef() {
						return false;
						}
					});
				chartGenerators.add(new AbstractGenotypeTypeGenerator(header.getGenotypeSamples()){
					@Override
					boolean IsIgnoringHomRef() {
						return true;
						}
					});
				chartGenerators.add(new AffectedSamplesGenerator(header.getNGenotypeSamples()));
				chartGenerators.add(new NumberOfNoCallGenerator(header.getNGenotypeSamples()));
				
				if(header.getFormatHeaderLine(VCFConstants.DEPTH_KEY)!=null) {
					chartGenerators.add(new AverageDepthGenerator(header.getGenotypeSamples()));
				}
				
				if(header.getFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY)!=null) {
					chartGenerators.add(new AverageGQGenerator(header.getGenotypeSamples()));
					}
				chartGenerators.add(new AveragePhasedGTGenerator(header.getGenotypeSamples()));
				
				for(final String sn:header.getSampleNamesInOrder()) {
					if(hasPred && this.enable_predictions_per_sample ) {
						this.chartGenerators.add(new PredictionGenerator(tools,sn));
						}
					}
				}
			chartGenerators.removeIf(G->!G.isEnabled());
			
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			long last = -1L;
			while(iter.hasNext())
				{
				final VariantContext ctx = progress.apply(iter.next());
				
				for(final ChartGenerator cg: this.chartGenerators) {
					if(!cg.isEnabled()) continue;
					cg.visit(ctx);
					}
				if(this.outputFile!=null)
					{
					final long now = System.currentTimeMillis();
					if(last<0L || now - last > this.refreshEverySeconds * 1000) {
						last=now;
						save();
						}
					}
				}
			progress.close();
			
			save();
			
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	
	
	public static void main(final String[] args) {
		new VcfStatsJfx().instanceMainWithExit(args);
	}

}
