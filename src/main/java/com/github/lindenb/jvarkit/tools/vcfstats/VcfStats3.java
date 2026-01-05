/*
The MIT License (MIT)
Copyright (c) 2026 Pierre Lindenbaum

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

import java.awt.geom.Point2D;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.DoubleConsumer;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gatk.GATKConstants;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.multiqc.MultiqcCustomContent;
import com.github.lindenb.jvarkit.pedigree.SamplePopulation;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.google.gson.stream.JsonWriter;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

## Example

```
java -jar dist/jvarkit.jar vcfstats src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz |  R --no-save 
```

END_DOC
 */
@Program(name="vcfstats2",
	description="Produce VCF statitics",
	keywords={"vcf","stats","multiqc"},
	modificationDate = "20240726",
	creationDate = "20131212",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfStats3 extends Launcher {
	private static final Logger LOG = Logger.of(VcfStats3.class);
	private static int ID_GENERATOR=0;

	@Parameter(names={"-o","--output"},description="output directory",required=true)
	private Path outputDirectory = null;
	
	@Parameter(names={"--pipe"},description="write input VCF to stdout")
	private boolean enabled_pipe_to_stdout = false;

	@Parameter(names={"--prefix"},description="prefix for output files")
	private String prefix = "vcfstats.";

	
	@Parameter(names={"--title"},description="main section title")
	private String main_section_title = "";
	@Parameter(names={"--description"},description="main section description")
	private String main_section_description = "";
	
	@Parameter(names={"--sample2population","--sample2pop"},description=SamplePopulation.OPT_DESC)
	private Path sample2catPath = null;
	@Parameter(names={"--other-samples"},description="if sample doesn't belong to a population in sample2pop, "
			+ "create a new population named 'x' and insert those lonely samples in that population")
	private String otherPopulationName="other";

	
	@Parameter(names={"-exclude","--exclude"},description="name of modules to be excluded",hidden=true)
	private String moduleExcludeStr = "";
	@Parameter(names={"--list"},description="list available modules and exit",help = true)
	private boolean list_modules = false;	
	

	
	private final SamplePopulation samplePopulation=new SamplePopulation();
	
	/** manage x =  ((int)(x/factor))*factor   */
	private class DoubleRounder implements DoubleUnaryOperator {
		final DecimalFormat decimalFormat;
		DoubleRounder(int n) {
			String str="#.";
			while(n>0) {
				str+="#";
				--n;
				}
			this.decimalFormat = new DecimalFormat(str);
			this.decimalFormat.setRoundingMode(java.math.RoundingMode.CEILING);
			}
		@Override
		public double applyAsDouble(double v) {
			return Double.parseDouble(this.decimalFormat.format(v));
			}
		}
	
	private class BinRounder implements DoubleUnaryOperator {
		final double[] bins;
		BinRounder(double...values) {
			bins = Arrays.copyOf(values, values.length);
			Arrays.sort(bins);
			if(bins.length==0) throw new IllegalArgumentException();
			}
		@Override
		public double applyAsDouble(double v) {
			for(int i=0;i< this.bins.length;i++) {
				if(v <= this.bins[i]) return this.bins[i];
				}
			return this.bins[this.bins.length-1];
			}
		}
	
	private static class DataPoint implements DoubleConsumer {
		private long count=0L;
		private double sum = 0.0;
		
		@Override
		public void accept(double value) {
			count++;
			sum+=value;
			}
		
		public long getCount() {
			return this.count;
			}
	
		public  OptionalDouble getAverage() {
			return isEmpty()?OptionalDouble.empty():OptionalDouble.of(this.sum/this.count);
			}
		public boolean isEmpty() {
			return this.getCount()==0;
			}
		}

	
	/** Analyzer **/
	private interface Analyzer {
		void init(VCFHeader h);
		String getName();
		String getDescription();
		void visit(final VariantContext ctx);
		void saveAsMultiqcJson() throws IOException;
		public boolean isEnabled();
		}
	
	/** AbstractAnalyzer **/
	private abstract class AbstractAnalyzer implements Analyzer {
		protected boolean enabled=true;
		protected Map<String,String> attributes=new HashMap<>();

		@Override
		public String getName() {
			return attributes.getOrDefault("name", getClass().getSimpleName());
			}
		
		@Override
		public String getDescription() {
			return attributes.getOrDefault("description", getName());
			}
		
		@Override
		public boolean isEnabled() {
			return enabled;
			}
		
		public AbstractAnalyzer attribute(String key,String value) {
			this.attributes.put(key, value);
			return this;
			}
		boolean acceptVariant(final VariantContext ctx) {
			return true;
			}
		boolean acceptGenotype(Genotype g) {
			return true;
			}
		
		protected abstract boolean hasAnyDataToSave();
		
		@Override
		public void saveAsMultiqcJson() throws IOException {
			if(!hasAnyDataToSave()) {
				LOG.warn("no data to save for "+getName());
				return;
			}
			final String fname= VcfStats3.this.prefix + getName().replaceAll("[^A-Za-z0-9_]", "_").replaceAll("[_]+","_");
			final Path jsonPath = outputDirectory.resolve(fname+"_mqc.json");
			LOG.info("saving "+jsonPath);
			try(JsonWriter w = new JsonWriter(Files.newBufferedWriter(jsonPath))) {
				saveAsMultiqcJson(w);
				w.flush();
				}
			}
		
		protected MultiqcCustomContent makeMultiqcCustomContent() {
			return new MultiqcCustomContent().
					attribute(MultiqcCustomContent.KEY_PARENT_ID, StringUtils.md5(StringUtils.ifBlank(VcfStats3.this.main_section_title,VcfStats3.class.getName()))).
					attribute(MultiqcCustomContent.KEY_PARENT_NAME, StringUtils.ifBlank(VcfStats3.this.main_section_title,VcfStats3.class.getName())).
					attribute(MultiqcCustomContent.KEY_PARENT_DESC,VcfStats3.this.main_section_description).
					attribute(MultiqcCustomContent.KEY_SECTION_ID,String.valueOf(ID_GENERATOR++)).
					attribute(MultiqcCustomContent.KEY_SECTION_NAME,getName()).
					attribute(MultiqcCustomContent.KEY_SECTION_DESC,getDescription())
					;
			}
		
		protected void boxplot(final JsonWriter w,final Map<String,List<Number>> content) throws IOException {
			makeMultiqcCustomContent().
				writeBoxplot(w, content);
			}
		
		protected <T extends Number> List<Point2D> toPoints(Counter<T> counter) {
			return counter.entrySet().
					stream().
					map(KV->new Point2D.Double(KV.getKey().doubleValue(), KV.getValue())).
					sorted((A,B)->Double.compare(A.getX(), B.getX())).
					collect(Collectors.toList());
			}
		
		protected abstract void saveAsMultiqcJson(JsonWriter w) throws IOException;
		
		protected Optional<Genotype> getSingleton(VariantContext ctx) {
			Genotype single = null;
			for(Genotype g:ctx.getGenotypes()) {
				if(!acceptGenotype(g)) continue;
				if(!g.hasAltAllele()) continue;
				if(single!=null) return Optional.empty();
				single=g;
				}
			return Optional.ofNullable(single);
			}
		}
	
	/** Abstract AbstractDensityPlot ******************************************************************************************/
	private abstract class AbstractDensityPlot extends AbstractAnalyzer {
		private final Map<String,Counter<Double>> pop2counter=new HashMap<>();
		private final String all_names= "ALL";
		@Override
		public void init(VCFHeader h) {
			if(!h.hasGenotypingData()) {
				pop2counter.put(all_names, new Counter<>());
				}
			else
				{
				for(SamplePopulation.Population pop: VcfStats3.this.samplePopulation.getPopulations()) {
					pop2counter.put(pop.getNameAndCount(), new Counter<>());
					}
				}
			}
		@Override
		protected boolean hasAnyDataToSave() {
			return pop2counter.entrySet().stream().anyMatch(KV->!KV.getValue().isEmpty());
			}
		protected abstract List<Double> toDouble(VariantContext ctx);
		
		protected boolean acceptValue(double v) {
			return true;
			}
		
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			final List<Double> values = toDouble(ctx);
			if(ctx.hasGenotypes()) {
				for(SamplePopulation.Population pop: VcfStats3.this.samplePopulation.getPopulations()) {
					if(pop.getSampleNames().stream().map(SN->ctx.getGenotype(SN)).filter(G->acceptGenotype(G)).anyMatch(G->G.hasAltAllele())) {
						for(double v  : values) {
							if(!acceptValue(v)) continue;
							pop2counter.get(pop.getNameAndCount()).incr(v);
							}
						}
					}
				}
			else
				{
				for(double v  : values) {
					if(!acceptValue(v)) continue;
					pop2counter.get(all_names).incr(v);
					}
				}
			}
		
		@Override
		protected void saveAsMultiqcJson(JsonWriter w) throws IOException {
			final Map<String,List<Point2D>> hash=new HashMap<>(this.pop2counter.size());
			for(String popName:this.pop2counter.keySet()) {
				final List<Point2D> points = toPoints(this.pop2counter.get(popName));
				if(points.isEmpty()) continue;
 				hash.put(popName, points);
				}
			makeMultiqcCustomContent().writeLineGraph(w,hash);
			}
		}
	
	private class QUALDensityPlot extends AbstractDensityPlot {
		private final DoubleUnaryOperator converter = new BinRounder(0,1,10,20,30,40,50,60,70,80,90,100,200,300,400,500,1000,2000,3000,4000,5000,10_000);
		@Override
		public void init(VCFHeader h) {
			super.init(h);
			}
		@Override
		protected List<Double> toDouble(VariantContext ctx) {			
			return Arrays.asList(converter.applyAsDouble(ctx.getPhredScaledQual()));
			}
		@Override
		boolean acceptVariant(final VariantContext ctx) {
			if(!ctx.hasLog10PError()) return false;
			return super.acceptVariant(ctx);
			}
		
		}
	private abstract class AbstractFORMATPerSampleBoxPlot extends AbstractAnalyzer {
		final Map<String,DataPoint> sn2data =new HashMap<>();
		@Override
		public void init(final VCFHeader h) {
			if(this.enabled) {
				this.enabled=h.hasGenotypingData() && h.getFormatHeaderLine(getFormatKey())!=null;
				}
			if(this.enabled) {
				for(final String sn: h.getSampleNamesInOrder()) {
					sn2data.put(sn, new DataPoint());
					}
				}
			}
		protected abstract String getFormatKey();
		protected abstract OptionalDouble extractValueFromGenotype(Genotype gt);
		@Override
		protected boolean hasAnyDataToSave() {
			return sn2data.values().stream().anyMatch(data->data.getCount()>0);
			}
		
		protected OptionalDouble extractDataPoint(DataPoint data) {
			return data.getAverage();
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			for(Genotype gt:ctx.getGenotypes()) {
				if(!acceptGenotype(gt)) continue;
				OptionalDouble optV = extractValueFromGenotype(gt);
				if(!optV.isPresent()) continue;
				sn2data.get(gt.getSampleName()).accept(optV.getAsDouble());
				}
			}
		@Override
		protected void saveAsMultiqcJson(JsonWriter w) throws IOException {
			final AutoMap<String,Number,List<Number>> content= AutoMap.makeList();
			for(String sample: sn2data.keySet()) {
				for(SamplePopulation.Population pop : samplePopulation.getSampleByName(sample).getPopulations()) {
					final OptionalDouble v= extractDataPoint(sn2data.get(sample));
					if(!v.isPresent()) continue;
					content.insert(pop.getName()+ "(N="+pop.size()+")",v.orElse(-1));
					}
				}
			makeMultiqcCustomContent().writeBoxplot(w, content);
			}
		}
	
	private class DepthPerSampleBoxPlot extends AbstractFORMATPerSampleBoxPlot {
		@Override
		protected String getFormatKey() {
			return VCFConstants.DEPTH_KEY;
			}
		@Override
		protected OptionalDouble extractValueFromGenotype(Genotype gt) {
			return gt.hasDP()?OptionalDouble.of(gt.getDP()):OptionalDouble.empty();
			}
		}
	
	private class GQPerSampleBoxPlot extends AbstractFORMATPerSampleBoxPlot {
		@Override
		protected String getFormatKey() {
			return VCFConstants.GENOTYPE_QUALITY_KEY;
			}
		@Override
		protected OptionalDouble extractValueFromGenotype(Genotype gt) {
			return gt.hasGQ()?OptionalDouble.of(gt.getGQ()):OptionalDouble.empty();
			}
		}
	private class FILTERedPerSampleBoxPlot extends AbstractFORMATPerSampleBoxPlot {
		@Override
		protected String getFormatKey() {
			return VCFConstants.GENOTYPE_FILTER_KEY;
			}
		@Override
		protected OptionalDouble extractValueFromGenotype(Genotype gt) {
			return OptionalDouble.of(gt.isFiltered()?1:0);
			}
		}
	
	private class HetAFPerSampleBoxPlot extends AbstractFORMATPerSampleBoxPlot {
		@Override
		protected String getFormatKey() {
			return VCFConstants.ALLELE_FREQUENCY_KEY;
			}
		@Override
		boolean acceptGenotype(Genotype g) {
			return g.isHet() && g.hasExtendedAttribute(getFormatKey());
			}
		@Override
		protected OptionalDouble extractValueFromGenotype(Genotype gt) {
			Object o=gt.getExtendedAttribute(getFormatKey());
			if(o==null) return OptionalDouble.empty();
			if(o instanceof Number)  {
				Number num = Number.class.cast(o);
				return OptionalDouble.of(num.doubleValue());
				}
			if(o instanceof List)  {
				List<?> L = List.class.cast(o);
				if(!L.isEmpty()) return OptionalDouble.of(Number.class.cast(L.get(0)).doubleValue());
				}
			return OptionalDouble.empty();
			}
		}
	
	private class INFODensityPlot extends AbstractDensityPlot {
		private String infoKey="";
		private final DoubleUnaryOperator normalizer;
		protected INFODensityPlot(final String infoKey, final DoubleUnaryOperator normalizer) {
			super();
			this.infoKey = infoKey;
			this.normalizer=normalizer;
			}
		@Override
		public void init(final VCFHeader h) {
			super.init(h);
			if(this.enabled) {
				this.enabled = h.getFormatHeaderLine(this.infoKey)!=null;
				}
			}
		@Override
		boolean acceptVariant(final VariantContext ctx) {
			if(!ctx.hasAttribute(this.infoKey)) return false;
			return super.acceptVariant(ctx);
			}
		@Override
		protected List<Double> toDouble(VariantContext ctx) {
			return ctx.getAttributeAsDoubleList(this.infoKey, 0).stream().
					mapToDouble(V->normalizer.applyAsDouble(V)).
					boxed().
					collect(Collectors.toList());
			}
		}
	
	
	

	/** Abstract BoxPlot ******************************************************************************************/

	
	private abstract class AbstractFractionOfVariants extends AbstractAnalyzer {
		protected VCFHeader vcfHeader= null;
		protected final Counter<String> counter = new Counter<>();
		protected long n_variants = 0L;
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData();
			this.vcfHeader = h;
			}
		@Override
		protected boolean hasAnyDataToSave() {
			return n_variants>0L;
			}
		
		protected Map<String,List<Number>> groupByPopulation() {
			final AutoMap<String,Number,List<Number>> content= AutoMap.makeList();
			for(String sample: this.vcfHeader.getSampleNamesInOrder()) {
				for(SamplePopulation.Population pop : samplePopulation.getSampleByName(sample).getPopulations()) {
					content.insert(pop.getName()+ "(N="+pop.size()+")",this.counter.count(sample)/(double)n_variants);
					}
				}
			return content;
			}
		
		@Override
		protected void saveAsMultiqcJson(JsonWriter w) throws IOException {
			boxplot(w,groupByPopulation());
			}
		}
	
	private class Singleton extends AbstractFractionOfVariants {
		@Override
		public void init(final VCFHeader h) {
			super.init(h);
			if(this.enabled && h.getNGenotypeSamples()<=1) {
				this.enabled=false;
				}
			}
		
		
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			super.n_variants++;
			final Genotype single = getSingleton(ctx).orElse(null);
			if(single!=null) super.counter.incr(single.getSampleName());
			}
		}
	
	
	/** export number of variants per category */
	private class FiltersPerPop extends AbstractAnalyzer {
		private final Map<SamplePopulation.Population,Counter<String>> pop2data = new HashMap<>();
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData() && !h.getFilterLines().isEmpty();
			for(SamplePopulation.Population pop: VcfStats3.this.samplePopulation.getPopulations()) {
				pop2data.put(pop, new Counter<>());
				}
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			
			final Set<String> filters;
			if(ctx.isFiltered()) {
				filters = ctx.getFilters();
				}
			else
				{
				filters = Collections.singleton(VCFConstants.PASSES_FILTERS_v4);
				}
			
			for(SamplePopulation.Population pop: VcfStats3.this.samplePopulation.getPopulations()) {
				if(pop.getSampleNames().stream().map(SN->ctx.getGenotype(SN)).filter(G->acceptGenotype(G)).anyMatch(G->G.hasAltAllele())) {
					for(final String ft:filters) {
						pop2data.get(pop).incr(ft);
						}
					}
				}
			}
		@Override
		protected void saveAsMultiqcJson(final JsonWriter w) throws IOException {
			final Map<String,Map<String,Number>> hash=new HashMap<>();
			for(final SamplePopulation.Population pop:this.pop2data.keySet()) {
				hash.put(
					pop+" (N="+pop.size()+")",
					this.pop2data.get(pop).entrySet().stream().
						collect(
								Collectors.toMap(
										KV->KV.getKey(),
										KV->(long)KV.getValue()
										)
						));
				}
			
			makeMultiqcCustomContent().writeMultiBarGraph(w,hash);
			}
		@Override
		protected boolean hasAnyDataToSave() {
			return pop2data.entrySet().stream().anyMatch(KV->!KV.getValue().isEmpty());
			}
		}
	
	/** export number of variants per category */
	private class VariantTypePerPop extends AbstractAnalyzer {
		private final Map<String,Counter<VariantContext.Type>> pop2data = new HashMap<>();
		private final String allPopName="ALL";
		@Override
		public void init(final VCFHeader h) {
			if(h.hasGenotypingData()) {
				if(h.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY)==null) this.enabled=false;
				for(SamplePopulation.Population pop: VcfStats3.this.samplePopulation.getPopulations()) {
					pop2data.put(pop.getNameAndCount(), new Counter<>());
					}
				}
			else
				{
				pop2data.put(allPopName, new Counter<>());
				}
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			if(ctx.hasGenotypes()) {
				for(SamplePopulation.Population pop: VcfStats3.this.samplePopulation.getPopulations()) {
					if(pop.getSampleNames().stream().map(SN->ctx.getGenotype(SN)).filter(G->acceptGenotype(G)).anyMatch(G->G.hasAltAllele())) {
						pop2data.get(pop.getNameAndCount()).incr(ctx.getType());
						}
					}
				}
			else
				{
				pop2data.get(allPopName).incr(ctx.getType());
				}
			}
		@Override
		protected void saveAsMultiqcJson(final JsonWriter w) throws IOException {
			final Map<String,Map<String,Number>> hash=new HashMap<>();
			for(final String popName:this.pop2data.keySet()) {
				if(this.pop2data.get(popName).isEmpty()) continue;
				hash.put(
					popName,
					this.pop2data.get(popName).entrySet().stream().
						collect(
								Collectors.toMap(
										KV->KV.getKey().name(),
										KV->(long)KV.getValue()
										)
						));
				}
			
			makeMultiqcCustomContent().writeMultiBarGraph(w,hash);
			}
		@Override
		protected boolean hasAnyDataToSave() {
			return pop2data.entrySet().stream().anyMatch(KV->!KV.getValue().isEmpty());
			}
		}
	
	/***************************************************************************/
	private class AFHistogram extends AbstractAnalyzer {
		private final Map<SamplePopulation.Population,Counter<Double>> pop2data = new HashMap<>();
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData() && h.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY)!=null;
			for(SamplePopulation.Population pop: VcfStats3.this.samplePopulation.getPopulations()) {
				pop2data.put(pop, new Counter<>());
				}
			}
		OptionalDouble getAF(final Allele allele,final SamplePopulation.Population pop,VariantContext ctx) {
			double total=0;
			double sum=0;
			for(String sn: pop.getSampleNames()) {
				final Genotype g= ctx.getGenotype(sn);
				if(!acceptGenotype(g)) continue;
				total+=g.getPloidy();
				sum+= g.getAlleles().stream().filter(A->A.equals(allele)).count();
				}
			if(total==0) return OptionalDouble.empty();
			return OptionalDouble.of(sum/total);
			}
		
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			List<Allele> alts = ctx.getAlternateAlleles();
			for(SamplePopulation.Population pop : VcfStats3.this.samplePopulation.getPopulations()) {
				for(Allele alt: alts) {
					final OptionalDouble opt = getAF(alt,pop,ctx);
					if(!opt.isPresent()) continue;
					pop2data.get(pop).incr(opt.getAsDouble());
					}
				}
			}
		@Override
		protected boolean hasAnyDataToSave() {
			return pop2data.entrySet().stream().anyMatch(KV->!KV.getValue().isEmpty());
			}
		@Override
		protected void saveAsMultiqcJson(JsonWriter w) throws IOException {
			final Map<String,List<Point2D>> hash=new HashMap<>();
			for(final SamplePopulation.Population pop:this.pop2data.keySet()) {
				final List<Point2D> L = toPoints(this.pop2data.get(pop));
				if(L.isEmpty()) continue;
				hash.put(pop+" (N="+pop.size()+")", L);
				}
			makeMultiqcCustomContent().writeLineGraph(w, hash);
			}
		}
	
	/*************************************************************************************************/


	private void loadPhenotypes(final VCFHeader header) throws IOException {
		if(!header.hasGenotypingData()) return;
		final Set<String> sampleset = new HashSet<>(header.getGenotypeSamples());
		if(sample2catPath!=null) {
			samplePopulation.load(sample2catPath);
			samplePopulation.retain(header);
			}
		if(!StringUtils.isBlank(otherPopulationName)) {
			if(this.samplePopulation.hasCollection(this.otherPopulationName)) {
				throw new IllegalArgumentException("collection \""+this.otherPopulationName+"\" already defined in "+this.sample2catPath);
				}
			for(String sn:sampleset) {
				if(samplePopulation.hasSample(sn)) continue;
				samplePopulation.insert(sn, this.otherPopulationName);
				}
			}
		}
	
	
	
	@Override
	public int doWork(final List<String> args) {
		final List<Analyzer> modules =new ArrayList<>();
		modules.add(new Singleton());
		modules.add(new QUALDensityPlot());
		modules.add(new INFODensityPlot(GATKConstants.FS_KEY,new DoubleRounder(1)));
		modules.add(new INFODensityPlot(GATKConstants.MQ_KEY,new DoubleRounder(1)));
		modules.add(new INFODensityPlot(GATKConstants.MQRankSum_KEY,new DoubleRounder(10)));
		modules.add(new INFODensityPlot(GATKConstants.QD_KEY,new DoubleRounder(1)));
		modules.add(new INFODensityPlot(GATKConstants.ReadPosRankSum_KEY,new DoubleRounder(1)));
		modules.add(new INFODensityPlot(GATKConstants.SOR_KEY,new DoubleRounder(10)));
		
		modules.add(new DepthPerSampleBoxPlot());
		modules.add(new GQPerSampleBoxPlot());
		modules.add(new FILTERedPerSampleBoxPlot());
		modules.add(new FiltersPerPop());
		modules.add(new VariantTypePerPop());
		modules.add(new HetAFPerSampleBoxPlot());
		modules.add(new AFHistogram());
		
		final String input = oneFileOrNull(args);
		
		// remove modules
		modules.removeIf(M->Arrays.stream(this.moduleExcludeStr.split("[,; \t:]")).anyMatch(S->S.equalsIgnoreCase(M.getName())));
		try {
			if(this.list_modules) {
				try(PrintWriter w = super.openPathOrStdoutAsPrintWriter(null)) {
					for(Analyzer analyzer:modules) {
						w.print(analyzer.getName());
						w.print("\t");
						w.println(analyzer.getDescription());
						}	
					w.flush();
					}
				return 0;
				}
			final VariantContextWriter variantCtxWriter =  enabled_pipe_to_stdout?
					VCFUtils.createVariantContextWriterToStdout():
					null;
			
			try(VCFIterator iter= super.openVCFIterator(input)) {	
				final VCFHeader header=iter.getHeader();
				loadPhenotypes(header);
				if(variantCtxWriter!=null) {
					variantCtxWriter.writeHeader(header);
					}
				
				for(Analyzer analyzer:modules) {
					analyzer.init(header);
					if(!analyzer.isEnabled()) {
						LOG.warn("module "+analyzer.getName()+" will be disabled. ["+analyzer.getClass().getSimpleName()+"]");
						}
					}
				modules.removeIf(M->!M.isEnabled());
				if(modules.isEmpty()) {
					LOG.warn("no module was enabled");
					}
				while(iter.hasNext()) {
					final VariantContext ctx = iter.next();
					if(variantCtxWriter!=null) variantCtxWriter.add(ctx);
					for(Analyzer analyzer:modules) {
						analyzer.visit(ctx);
						}
					}
				if(variantCtxWriter!=null) {
					variantCtxWriter.close();
					}
				
				for(Analyzer analyzer:modules) {
					analyzer.saveAsMultiqcJson();
					}
					
				}
			return 0;
		} catch (final Throwable e) {
			e.printStackTrace();
			LOG.error(e);
			return -1;
			}
		
		}
	
			
	public static void main(final String[] args)
		{
		new VcfStats3().instanceMainWithExit(args);
		}
	}
