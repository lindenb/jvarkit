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

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.DoubleConsumer;
import java.util.function.DoublePredicate;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Predicate;


import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gatk.GATKConstants;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.RangeOfDoubles;
import com.github.lindenb.jvarkit.multiqc.MultiqcCustomContent;
import com.github.lindenb.jvarkit.pedigree.SamplePopulation;
import com.github.lindenb.jvarkit.tools.multiqc.MultiqcPostProcessor;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.google.gson.stream.JsonWriter;

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

	@Parameter(names={"-o","--output"},description="output directory")
	private Path outputDirectory = null;
	
	@Parameter(names={"--pipe"},description="write input VCF to stdout")
	private boolean enabled_pipe_to_stdout = false;

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
		protected Predicate<VariantContext> acceptVariant = V->true;
		protected Predicate<Genotype> acceptGT = V->true;

		@Override
		public String getName() {
			return attributes.getOrDefault("name", getClass().getName());
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
	
		AbstractAnalyzer setVariantPredicate(final Predicate<VariantContext> f) {
			this.acceptVariant= f;
			return this;
			}
		AbstractAnalyzer setGenotypePredicate(final Predicate<Genotype> f) {
			this.acceptGT = f;
			return this;
			}
		
		boolean acceptVariant(VariantContext ctx) {
			return this.acceptVariant.test(ctx);
			}
		boolean acceptGenotype(Genotype g) {
			return this.acceptGT.test(g);
			}
		
		protected abstract boolean hasAnyDataToSave();
		
		@Override
		public void saveAsMultiqcJson() throws IOException {
			if(!hasAnyDataToSave()) return;
			final String fname= getName().replaceAll("[^A-Za-z0-9_]", "_").replaceAll("[_]+","_");
			final Path jsonPath = outputDirectory.resolve(fname+"_mqc.json");
			try(JsonWriter w = new JsonWriter(Files.newBufferedWriter(jsonPath))) {
				saveAsMultiqcJson(w);
				w.flush();
				}
			}
		
		protected void boxplot(JsonWriter w,Map<String,List<Number>> content) throws IOException {
			new MultiqcCustomContent().
				attribute(MultiqcCustomContent.KEY_PARENT_ID, StringUtils.md5(StringUtils.ifBlank(VcfStats2.this.main_section_title,MultiqcPostProcessor.class.getName()))).
				attribute(MultiqcCustomContent.KEY_PARENT_NAME, StringUtils.ifBlank(VcfStats2.this.main_section_title,MultiqcPostProcessor.class.getName())).
				attribute(MultiqcCustomContent.KEY_PARENT_DESC,VcfStats2.this.main_section_description).
				attribute(MultiqcCustomContent.KEY_SECTION_ID,String.valueOf(ID_GENERATOR++)).
				attribute(MultiqcCustomContent.KEY_SECTION_NAME,getName()).
				attribute(MultiqcCustomContent.KEY_SECTION_DESC,getDescription()).
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
		private DecimalFormat decimalFormat = null;
		@Override
		public void init(VCFHeader h) {
			String str="#.";
			int p = Integer.parseInt(this.attributes.getOrDefault("precision","1"));
			while(p>1) {
				str+="#";
				p=p/10;
				}
			this.decimalFormat = new DecimalFormat(str);
			this.decimalFormat.setRoundingMode(java.math.RoundingMode.CEILING);
			}
		
		@Override
		protected boolean hasAnyDataToSave() {
			return !counter.isEmpty();
			}
		protected abstract List<Double> toDouble(VariantContext ctx);
		
		protected double toDiscreteValue(double v) {
			return Double.parseDouble(this.decimalFormat.format(v));
			}
		
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
			
			}
		}
	
	private abstract class QUALDensityPlot extends AbstractDensityPlot {
		@Override
		protected List<Double> toDouble(VariantContext ctx) {
			return Arrays.asList(ctx.getPhredScaledQual());
			}
		@Override
		boolean acceptVariant(VariantContext ctx) {
			if(!ctx.hasLog10PError()) return false;
			return super.acceptVariant(ctx);
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
	private class DensityPlot extends AbstractINFODensityPlot {
		final DoubleUnaryOperator converter;
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
		protected final Counter<String> counter = new Counter<>();
		protected long n_variants = 0L;
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData();
			}
		@Override
		protected boolean hasAnyDataToSave() {
			return n_variants>0L;
			}
		protected Map<String,List<Number>> groupByPopulation() {
			Map<String,List<Number>> content=new TreeMap<>();
			for(String sample: this.counter.keySet()) {
				String pop = samplePopulation.getSampleByName(sample).getPopulation().getName();
				List<Number> L = content.get(pop);
				if(!content.containsKey(pop)) {
					L=new ArrayList<>();
					content.put(pop, L);
					}
				L.add(this.counter.count(sample)/(double)n_variants);
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
	private class VariantType extends AbstractFractionOfVariants {
		final Counter<VariantContext.Type> counter=new Counter<>();
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			super.n_variants++;
			counter.incr(ctx.getType());
			}
		@Override
		protected boolean hasAnyDataToSave() {
			return !counter.isEmpty();
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
		modules.add(new DensityPlot(GATKConstants.FS_KEY,(V)->V,V->true));
		modules.add(new DensityPlot(GATKConstants.MQ_KEY,(V)->V,V->true));
		modules.add(new DensityPlot(GATKConstants.MQRankSum_KEY,(V)->V,V->true).attribute("precision", "10"));
		modules.add(new DensityPlot(GATKConstants.QD_KEY,(V)->V,V->true));
		modules.add(new DensityPlot(GATKConstants.ReadPosRankSum_KEY,(V)->V,V->true));
		modules.add(new DensityPlot(GATKConstants.SOR_KEY,(V)->V,V->true).attribute("precision", "10"));
		modules.add(new VariantType());
		
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
