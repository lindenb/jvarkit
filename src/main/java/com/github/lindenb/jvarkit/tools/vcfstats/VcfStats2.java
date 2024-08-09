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
package com.github.lindenb.jvarkit.tools.vcfstats;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.DoubleConsumer;
import java.util.function.DoublePredicate;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.multiqc.MultiqcCustomContent;
import com.github.lindenb.jvarkit.pedigree.SamplePopulation;
import com.github.lindenb.jvarkit.tools.multiqc.MultiqcPostProcessor;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
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
public class VcfStats2 extends Launcher {
	private static final Logger LOG = Logger.build(VcfStats2.class).make();
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
		public  OptionalDouble getSum() {
			return isEmpty()?OptionalDouble.empty():OptionalDouble.of(this.sum);
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
		protected Predicate<Genotype> acceptGT = V->true;

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
		boolean acceptVariant(VariantContext ctx) {
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
			final String fname= VcfStats2.this.prefix + getName().replaceAll("[^A-Za-z0-9_]", "_").replaceAll("[_]+","_");
			final Path jsonPath = outputDirectory.resolve(fname+"_mqc.json");
			LOG.info("saving "+jsonPath);
			try(JsonWriter w = new JsonWriter(Files.newBufferedWriter(jsonPath))) {
				saveAsMultiqcJson(w);
				w.flush();
				}
			}
		
		protected MultiqcCustomContent makeMultiqcCustomContent() {
			return new MultiqcCustomContent().
					attribute(MultiqcCustomContent.KEY_PARENT_ID, StringUtils.md5(StringUtils.ifBlank(VcfStats2.this.main_section_title,VcfStats2.class.getName()))).
					attribute(MultiqcCustomContent.KEY_PARENT_NAME, StringUtils.ifBlank(VcfStats2.this.main_section_title,VcfStats2.class.getName())).
					attribute(MultiqcCustomContent.KEY_PARENT_DESC,VcfStats2.this.main_section_description).
					attribute(MultiqcCustomContent.KEY_SECTION_ID,String.valueOf(ID_GENERATOR++)).
					attribute(MultiqcCustomContent.KEY_SECTION_NAME,getName()).
					attribute(MultiqcCustomContent.KEY_SECTION_DESC,getDescription())
					;
			}
		
		protected void boxplot(final JsonWriter w,final Map<String,List<Number>> content) throws IOException {
			makeMultiqcCustomContent().
				writeBoxplot(w, content);
			}
		
		protected abstract void saveAsMultiqcJson(JsonWriter w) throws IOException;
		
		Optional<Genotype> getSingleton(VariantContext ctx) {
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
		private final Counter<Double> counter=new Counter<>();
		
		
		@Override
		protected boolean hasAnyDataToSave() {
			return !counter.isEmpty();
			}
		protected abstract List<Double> toDouble(VariantContext ctx);
		
		protected abstract double toDiscreteValue(double v);

		protected boolean acceptValue(double v) {
			return true;
			}
		
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			for(Double value:toDouble(ctx)) {
				if(!acceptValue(value)) continue;
				counter.incr(toDiscreteValue(value));
				}
			}
		
		@Override
		protected void saveAsMultiqcJson(JsonWriter w) throws IOException {
			makeMultiqcCustomContent().writeXY(w,
					this.counter.entrySet().stream().
						map(KV->new Point2D.Double(KV.getKey(),KV.getValue())).
						collect(Collectors.toList())
						);
			}
		}
	
	private class QUALDensityPlot extends AbstractDensityPlot {
		@Override
		public void init(VCFHeader h) {
			
			}
		
		@Override
		protected List<Double> toDouble(VariantContext ctx) {			
			return Arrays.asList(toDiscreteValue(ctx.getPhredScaledQual()));
			}
		@Override
		boolean acceptVariant(VariantContext ctx) {
			if(!ctx.hasLog10PError()) return false;
			return super.acceptVariant(ctx);
			}
		@Override
		protected double toDiscreteValue(double v) {
			return v;
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
	
	private abstract class AbstractINFODensityPlot extends AbstractDensityPlot {
		private String infoKey="";
		@Override
		public void init(final VCFHeader h) {
			this.infoKey = super.attributes.getOrDefault("info.key","");
			if(StringUtils.isBlank(infoKey)) throw new IllegalStateException("undefined info.key");
			this.enabled = h.hasGenotypingData() && h.getFormatHeaderLine(this.infoKey)!=null;
			}
		@Override
		boolean acceptVariant(VariantContext ctx) {
			if(!ctx.hasAttribute(this.infoKey)) return false;
			return super.acceptVariant(ctx);
			}
		@Override
		protected List<Double> toDouble(VariantContext ctx) {
			return ctx.getAttributeAsDoubleList(this.infoKey, 0);
			}
		}
	private abstract class DensityPlot extends AbstractINFODensityPlot {
		 DoubleUnaryOperator converter=DoubleUnaryOperator.identity();
		final DoublePredicate acceptValue;
		DensityPlot(final String key,final DoubleUnaryOperator converter, final DoublePredicate acceptValue) {
			attribute("info.key",key);
			attribute("name",key);
			attribute("description","Density of INFO/"+key);
			this.converter=converter;
			this.acceptValue=acceptValue;
			}
		@Override
		protected boolean acceptValue(double v) {
			return this.acceptValue.test(v);
			}
		}
	
	

	/** Abstract BoxPlot ******************************************************************************************/
	private abstract class AbstractBoxPlot extends AbstractAnalyzer {
		private final Map<String, Counter<Integer>> cat2values= new HashMap<>();
		protected void add(final String cat,int value) {
			Counter<Integer> l = cat2values.get(cat);
			if(l==null) {
				l=new Counter<>();
				cat2values.put(cat, l);
				}	
			l.incr(value);
			}

		@Override
		protected boolean hasAnyDataToSave() {
			return !cat2values.isEmpty();
			}
		
		
		}
	
	private abstract class AbstractFractionOfVariants extends AbstractBoxPlot {
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
	private class VariantTypePerPop extends AbstractAnalyzer {
		private final Map<SamplePopulation.Population,Counter<VariantContext.Type>> pop2data = new HashMap<>();
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData() && h.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY)!=null;
			for(SamplePopulation.Population pop: VcfStats2.this.samplePopulation.getPopulations()) {
				pop2data.put(pop, new Counter<>());
				}
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			for(SamplePopulation.Population pop: VcfStats2.this.samplePopulation.getPopulations()) {
				if(pop.getSampleNames().stream().map(SN->ctx.getGenotype(SN)).filter(G->acceptGenotype(G)).anyMatch(G->G.hasAltAllele())) {
					pop2data.get(pop).incr(ctx.getType());
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
	
	private abstract class AbstractFORMATValue extends AbstractBoxPlot {
		private final Map<String,DataPoint> sample2data = new HashMap<>();
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData() && h.getFormatHeaderLine(getFormatKey())!=null;
			for(String sn:h.getGenotypeSamples()) {
				sample2data.put(sn, new DataPoint());
				}
			}
		protected abstract String getFormatKey();
		
		protected OptionalDouble objectToDouble(Object o) {
			if(o==null) return OptionalDouble.empty();
			if(o.equals(VCFConstants.MISSING_VALUE_v4)) return  OptionalDouble.empty();
			if(o instanceof Number) return OptionalDouble.of(Number.class.cast(o).doubleValue());
			LOG.warn("cannot convert "+o+" to numeric.");
			return OptionalDouble.empty();
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			for(Genotype g:ctx.getGenotypes()) {
				if(!acceptGenotype(g)) continue;
				final OptionalDouble opt= objectToDouble(g.getAnyAttribute(getFormatKey()));
				if(!opt.isPresent()) continue;
				sample2data.get(g.getSampleName()).accept(opt.getAsDouble());
				}
			}
		@Override
		protected void saveAsMultiqcJson(JsonWriter w) throws IOException {
			
			}
		}
	/***************************************************************************/
	private class AFHistogram extends AbstractAnalyzer {
		private final Map<SamplePopulation.Population,Counter<Double>> pop2data = new HashMap<>();
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData() && h.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY)!=null;
			for(SamplePopulation.Population pop: VcfStats2.this.samplePopulation.getPopulations()) {
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
			for(SamplePopulation.Population pop : VcfStats2.this.samplePopulation.getPopulations()) {
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
				final List<Point2D> L = this.pop2data.get(pop).
						entrySet().
						stream().
						map(KV->new Point2D.Double(KV.getKey(),KV.getValue().doubleValue())).
						sorted((A,B)->Double.compare(A.getX(),B.getX())).
						collect(Collectors.toList());
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
	
	
	private static boolean valueInRange(double v,double minV,double maxV) {
		return minV<=v  && v<=maxV;
		}
	
	@Override
	public int doWork(final List<String> args) {
		final List<Analyzer> modules =new ArrayList<>();
		modules.add(new Singleton());
		modules.add(new QUALDensityPlot());
		/*modules.add(new DensityPlot(GATKConstants.FS_KEY,(V)->V,V->true));
		modules.add(new DensityPlot(GATKConstants.MQ_KEY,(V)->V,V->true));
		modules.add(new DensityPlot(GATKConstants.MQRankSum_KEY,(V)->V,V->true).attribute("precision", "10"));
		modules.add(new DensityPlot(GATKConstants.QD_KEY,(V)->V,V->true));
		modules.add(new DensityPlot(GATKConstants.ReadPosRankSum_KEY,(V)->V,V->true));
		modules.add(new DensityPlot(GATKConstants.SOR_KEY,(V)->V,V->true).attribute("precision", "10"));
		*/
		modules.add(new DepthPerSampleBoxPlot());
		modules.add(new GQPerSampleBoxPlot());
		modules.add(new FILTERedPerSampleBoxPlot());
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
		new VcfStats2().instanceMainWithExit(args);
		}
	}
