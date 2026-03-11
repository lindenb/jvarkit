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

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.function.DoubleConsumer;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.chart.BarPlot;
import com.github.lindenb.jvarkit.chart.BoxPlotChart;
import com.github.lindenb.jvarkit.chart.Chart;
import com.github.lindenb.jvarkit.chart.DataXY;
import com.github.lindenb.jvarkit.chart.NamedSeries;
import com.github.lindenb.jvarkit.chart.NamedY;
import com.github.lindenb.jvarkit.chart.ScatterXY;
import com.github.lindenb.jvarkit.chart.SeriesXY;
import com.github.lindenb.jvarkit.gatk.GATKConstants;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.DoubleRounder;
import com.github.lindenb.jvarkit.pedigree.SampleToGroup;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

## Example

```
java -jar dist/jvarkit.jar vcfstats src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz |  R --no-save 
```

END_DOC
 */
@Program(name="vcfstats",
	description="Produce VCF statitics",
	keywords={"vcf","stats","R"},
	modificationDate = "20230707",
	creationDate = "20131212",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfStats extends Launcher {
	private static final Logger LOG = Logger.of(VcfStats.class);

	@Parameter(names={"-o","--output"},description="Output directory",required = true)
	private Path outputDirectory = null;
		
	@Parameter(names={"--categories","--phenotypes"},description=SampleToGroup.OPT_DESC,hidden = true)
	private Path sample2catPath = null;

	
	@Parameter(names={"-exclude","--exclude"},description="name of modules to be excluded",hidden=true)
	private String moduleExcludeStr = "";
	@Parameter(names={"--list","-l"},description="list available modules and exit",help = true)
	private boolean list_modules = false;	
	@Parameter(names={"--prefix"},description="file prefix")
	private String prefix = "";
	@Parameter(names={"--title"},description="title as prefix")
	private String extra_title = "";

	@Parameter(names={"--fill-group"},description="assign group 'x' to sample if it's not defined in the file for --categories. Special value = '*' : assign the sample name to it's own private group")
	private String fill_sample_name_unknown_group = "other";

	
	@DynamicParameter(names={"-D"},description="other parameters.",hidden = true)
	private Map<String,String> __dynaParams = new HashMap<>();
	

	private final SampleToGroup sampleToGroup = new SampleToGroup();
	
	private static abstract class DataPoint implements DoubleConsumer,Supplier<OptionalDouble>,Comparable<DataPoint>{
		protected long count=0;
		public long getCount() {
			return this.count;
			}
		public  double getAsDouble() {
			return  get().getAsDouble();
			}
		public abstract boolean isPresent();
		@Override
		public int compareTo(final DataPoint o) {
			if(!isPresent()) {
				if(!o.isPresent()) return 0;
				return -1;
				}
			else if(!o.isPresent()) {
				return 1;
				}
			return Double.compare(this.getAsDouble(), o.getAsDouble());
			}
		}
	private static class DataPointSum extends DataPoint {
		protected double sum=0.0;
		
		@Override
		public boolean isPresent() {
			return true;
			}
		@Override
		public void accept(double v) { sum+=v;super.count++;}
		@Override
		public  OptionalDouble get() {
			return  OptionalDouble.of(sum);
			}
		}
	private static class DataPointAverage extends DataPointSum {
		@Override
		public boolean isPresent() {
			return getCount()>0L;
			}
		@Override
		public  OptionalDouble get() {
			return  isPresent()?OptionalDouble.of(super.sum/(double)getCount()):OptionalDouble.empty();
			}
		}
	
	
	/** Analyzer **/
	private static interface Analyzer {
		void init(VCFHeader h,Map<String,String> properties,final SampleToGroup sn2group);
		void visit(final VariantContext ctx);
		Set<Path> finish(Path outputDir) throws IOException,XMLStreamException;
		public String getName();
		public String getTitle();
		public boolean isEnabled();
		public Analyzer setProperty(String key,String v);
		public String getProperty(String key,String def);
		}
	
	/** AbstractAnalyzer **/
	private static  abstract class AbstractAnalyzer implements Analyzer {
		protected VCFHeader vcfHeader;
		protected boolean enabled=true;
		protected Predicate<VariantContext> _acceptVariant = V->true;
		protected Predicate<Genotype> _acceptGT = V->true;
		protected final Map<String, String> properties = new HashMap<>();
		protected SampleToGroup sampleToGroup = null;
		@Override
		public void init(VCFHeader h, final Map<String, String> properties, final SampleToGroup sampleToGroup) {
			this.vcfHeader = h;
			this.properties.putAll(properties);
			this.sampleToGroup = sampleToGroup;
			}
		@Override
		public String getProperty(final String key,final String def) {
			return this.properties.getOrDefault(key,def);
			}
		
		@Override
		public AbstractAnalyzer setProperty(String key,String v) {
			this.properties.put(key,v);
			return this;
			}
		
		@Override
		public boolean isEnabled() {
			return enabled;
			}
		protected String getLabelForGroup(final String grpName) {
			if(!this.sampleToGroup.hasGroup(grpName)) return grpName;
			final Set<String> sns  = this.sampleToGroup.getSamplesForGroup(grpName);
			if(sns.size()<2) return grpName;
			return grpName+" N="+sns.size();
			}
		
		
		@Override
		public final String getName() {
			return getProperty("name", String.valueOf(getClass().getSimpleName()));
			}
		@Override
		public String getTitle() {
			return getProperty("title",getName());
			}
		AbstractAnalyzer name(final String s) {
			setProperty("name", s);
			return this;
			}
		AbstractAnalyzer description(final String s) {
			return this.setProperty("description", s);
			}
		
		AbstractAnalyzer xlab(final String s) {
			return this.setProperty("xlab", s);
			}
		AbstractAnalyzer ylab(final String s) {
			return this.setProperty("ylab", s);
			}
		AbstractAnalyzer setAcceptVariant(final Predicate<VariantContext> f) {
			this._acceptVariant= f;
			return this;
			}
		
		Predicate<VariantContext> getVariantPredicate() {
			return this._acceptVariant;
			}
		Predicate<Genotype> getGenotypePredicate() {
			return this._acceptGT;
			}
		public boolean acceptVariant(final VariantContext ctx) {
			return this._acceptVariant.test(ctx);
			}
		public boolean acceptGenotype(final Genotype gt) {
			return this._acceptGT.test(gt);
			}
		Set<Path> exportChart(final Path outputDir,Chart chart) throws IOException,XMLStreamException {
			final String basename =  getProperty("filename","file");
			LOG.info("saving "+getName()+" to "+basename);
			final Set<Path> paths = new HashSet<Path>(3);
			Path p = outputDir.resolve(basename+".html");
			chart.savePlotly(p);
			paths.add(p);
			
			p = outputDir.resolve(basename+".xml");
			chart.saveXml(p);
			paths.add(p);
			
			p  = outputDir.resolve(basename+"_mqc.json");
			chart.saveMultiQC(p);
			paths.add(p);
			
			p  = outputDir.resolve(basename+".R");
			chart.saveR(p);
			paths.add(p);
			
			return paths;
			}
		
		}
	
	/** Abstract BoxPlot ******************************************************************************************/
	private abstract static class AbstractBoxPlot extends AbstractAnalyzer {
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
		public Set<Path> finish(final Path dir) {
			if(cat2values.isEmpty()) return Collections.emptySet();
			return Collections.emptySet();
			}
		}
	/***************************************************************************/
	private abstract static class AbstractBarPlot extends AbstractAnalyzer {
		}
	/***************************************************************************/

	private static  class RangeBarPlot extends AbstractAnalyzer {
		
		private final String infoTag;
		private final DoubleRounder doubleRounder;
		private final Map<VariantContext.Type,Counter<Double>> sv2type2ranges = new HashMap<>();
		private final boolean logX;
		
		RangeBarPlot(final String infoTag,int number_of_decimal_after_comma,boolean logX) {
			this.doubleRounder=new DoubleRounder(number_of_decimal_after_comma);
			this.infoTag = infoTag;
			this.logX = logX;
			super.properties.put("filename","info_"+infoTag+"_distribution");
			}
			
		private double round(final double v) {
			return this.doubleRounder.applyAsDouble(v);
			}		
		@Override
		public void init(final VCFHeader h,Map<String,String> prop,SampleToGroup s2g) {
			super.init(h, prop, s2g);
			final VCFInfoHeaderLine info = h.getInfoHeaderLine(this.infoTag);
			this.enabled = info!=null &&
					(info.getType()==VCFHeaderLineType.Float ||info.getType()==VCFHeaderLineType.Integer )
					;
			if(this.enabled) {
				for(VariantContext.Type vt : VariantContext.Type.values()) {
					sv2type2ranges.put(vt, new Counter<>());
					}
				}
			}
		
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			if(!ctx.hasAttribute(this.infoTag)) return;
			try {
				final VariantContext.Type vt = ctx.getType();
				for(String s : ctx.getAttributeAsStringList(this.infoTag,".") ) {
					if(s.equals(".")) continue;
					final double v0;
					try {
						v0 = Double.parseDouble(s);
						}
					catch(NumberFormatException err) {
						continue;
						}
					if(Double.isNaN(v0)) continue;
					if(Double.isInfinite(v0)) continue;
					final double v = round(v0);
					this.sv2type2ranges.get(vt).incr(v);
					}
				}
			catch(final Throwable err) {
				LOG.warn(err);
				}
			}
		

		
		@Override
		public Set<Path> finish(Path outputDir) throws IOException,XMLStreamException{
			
			final List<SeriesXY> L=new ArrayList<SeriesXY>();
			for(VariantContext.Type vt: this.sv2type2ranges.keySet()) {
				final Counter<Double> c = this.sv2type2ranges.get(vt);
				if(c.isEmpty()) continue;
				final SeriesXY series = new SeriesXY(
						vt.name(),
						c.stream().map(KV->new DataXY(KV.getKey(),KV.getValue())).collect(Collectors.toList())
						);
				if(series.isEmpty()) continue;
				series.sort();
				L.add(series);
				}
			if(L.isEmpty()) return Collections.emptySet();
			final ScatterXY chart = new ScatterXY(L);
			chart.setTitle(getTitle());
			chart.setYAxisLabel("Count "+this.infoTag);
			chart.setLogX(this.logX);
			chart.setXAxisLabel(this.logX?"log("+this.infoTag+")":this.infoTag);
			return exportChart(outputDir,chart);
			}
	}

	/***************************************************************************/
	private abstract static class AbstractMultipleBarPlot extends AbstractBarPlot {
		private final Map<String,Counter<String>> horiz2counts = new LinkedHashMap<>();
		private final Set<String> distinct_vertical = new LinkedHashSet<>();
		/** category/count > value */
		//private BiFunction<String,Long,Double> normalizer = (A,L)->L.doubleValue();
		protected void add(final String h,final String v) {
			Counter<String> t = horiz2counts.get(h);
			if(t==null) {
				t=new Counter<>();
				horiz2counts.put(h, t);
				}
			t.incr(v);
			distinct_vertical.add(v);
			}
		
		
		
		
		@Override
		public Set<Path> finish(Path outputDir ) {
			if(horiz2counts.isEmpty()) return Collections.emptySet();
			return Collections.emptySet();
			}
		}
	/***************************************************************************/
	private abstract static class AbstractManhattanPlot extends AbstractAnalyzer {
		
		private Predicate<SAMSequenceRecord> acceptContig ;
		private SAMSequenceDictionary dict;
		private long genomeLength;
		private final Map<String,List<DataPoint>> cat2index = new HashMap<>();
		private final int win_width;
		AbstractManhattanPlot() {
			this.acceptContig = SSR -> SSR.getContig().matches(getProperty("contig.regex","(chr)?[0-9XY][0-9]?"));
			this.win_width = Integer.parseInt(getProperty("manhattan.width","1000"));
			}
		
		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			final SAMSequenceDictionary dict0 = h.getSequenceDictionary();
			if(dict0!=null) {
				this.dict =  new SAMSequenceDictionary(dict0.getSequences().stream().filter(this.acceptContig).collect(Collectors.toList()));
				this.genomeLength = this.dict.getReferenceLength();
				}
			this.enabled = dict!=null && !dict.isEmpty();
			}
		private long pos2genomeindex(String contig,int pos) {
			long n=0;
			for(SAMSequenceRecord ssr:dict.getSequences()) {
				if(ssr.getContig().equals(contig)) {
					return n + pos;
					}
				n+=ssr.getLengthOnReference();
				}
			return -1;
			}
		
		private int pos2pixel(String contig,int pos) {
			return (int)((pos2genomeindex(contig,pos)/(double)this.genomeLength)* this.win_width);
			}
		protected DataPoint createDataPoint(final String category) {
			return new DataPointAverage();
			}
		
		
		
		protected void visit(String category,final Locatable loc,double value) {
			final SAMSequenceRecord ssr= this.dict.getSequence(loc.getContig());
			if(ssr==null || loc.getStart() <1 || loc.getStart() > ssr.getLengthOnReference()) return;
			List<DataPoint> datapoints;
			if(this.cat2index.containsKey(category)) {
				datapoints = this.cat2index.get(category);
				}
			else
				{
				datapoints = new ArrayList<>();
				this.cat2index.put(category,datapoints);
				}
			final int pixl0 = pos2pixel(loc.getContig(),loc.getStart());
			final int pixl1 = pos2pixel(loc.getContig(),loc.getEnd());
			for(int pixl = pixl0;  pixl <= pixl1 ; ++pixl ) {
				while(datapoints.size() <= pixl) {
					datapoints.add(null);
					}
				DataPoint dp  = datapoints.get(pixl);
				if(dp==null) {
					dp = createDataPoint(category);
					datapoints.set(pixl,dp);
					}
				dp.accept(value);
				}
			}
		@Override
		public Set<Path> finish(Path  outdir) {
			if(this.cat2index.isEmpty()) return Collections.emptySet();
			return Collections.emptySet();
			}
		}
	/***************************************************************************/
	private static class SingletonAnalyzer extends AbstractAnalyzer {
		private Counter<String> sample2count = new Counter<>();
		private long n_variants = 0L;
		@Override
		public void init(final VCFHeader h ,Map<String,String> props, SampleToGroup s2g) {
			super.init(h,props,s2g);
			if(this.enabled && h.getNGenotypeSamples()<=1) {
				this.enabled=false;
				}
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			n_variants++;
			Genotype single = null;
			for(Genotype g:ctx.getGenotypes()) {
				if(!g.hasAltAllele()) continue;
				if(!acceptGenotype(g)) continue;
				if(single!=null) return;
				single=g;
				}
			if(single!=null) sample2count.incr(single.getSampleName());
			}
		@Override
		public Set<Path> finish(Path outputDir) throws IOException, XMLStreamException {
			if(sample2count.isEmpty()) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			final List<NamedSeries> series = new ArrayList<NamedSeries>(super.sampleToGroup.getGroupsCount());
			for(String groupName: super.sampleToGroup.getGroups()) {
				final List<NamedY> L2 =new ArrayList<>();
				for(String sn: this.sample2count.keySet()) {
					if(!super.sampleToGroup.hasSampleInGroup(sn,groupName)) continue;
					L2.add(new NamedY(sn,this.sample2count.count(sn)));
					}
				if(L2.isEmpty()) continue;
				
				series.add(new NamedSeries(getLabelForGroup(groupName), L2));
				}
			
			final BoxPlotChart chart = new BoxPlotChart(series);
			chart.sortOnName();
			chart.setTitle(getTitle()+" (n-variants = "+n_variants+")");
			chart.setYAxisLabel("count");
			chart.setXAxisLabel("Group");
			return exportChart(outputDir,chart);
			}
		}
	
	
	
	/***************************************************************************/
		
	private static class CountPerContig extends AbstractAnalyzer {
		private long n_variants = 0L;
		private SAMSequenceDictionary dict = null;
		private final Predicate<SAMSequenceRecord> predicate;
		private final Map<String,Counter<VariantContext.Type>> chromosome2count = new HashMap<>();
		CountPerContig(final Predicate<SAMSequenceRecord> predicate) {
			this.predicate=predicate;
			}
		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			final SAMSequenceDictionary dict0 = h.getSequenceDictionary();
			if(dict0!=null) {
				this.dict = new SAMSequenceDictionary(
					dict0.getSequences().stream().
					filter(this.predicate).
					collect(Collectors.toList())
					);
				
				for(SAMSequenceRecord ssr: this.dict.getSequences()) {
					this.chromosome2count.put(ssr.getSequenceName(), new Counter<>());
					}
				}
			this.enabled = dict!=null && !dict.isEmpty();
			
			}
	
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			final Counter<VariantContext.Type> c=this.chromosome2count.get(ctx.getContig());
			if(c==null) return;
			n_variants++;
			c.incr(ctx.getType());
			}
		
		@Override
		public Set<Path> finish(final Path outputDir) throws IOException, XMLStreamException {
			if(n_variants==0) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			final List<NamedSeries> series=new ArrayList<>();
			for(SAMSequenceRecord ssr: this.dict.getSequences()) {
				final List<NamedY> L2 =new ArrayList<>();
				
				for(VariantContext.Type vtype: VariantContext.Type.values()) {
					L2.add(new NamedY(vtype.name(),this.chromosome2count.get(ssr.getContig()).count(vtype)));
					}
				series.add(new NamedSeries(ssr.getContig(), L2));
				}
			
			final BarPlot chart = new BarPlot(series);
			chart.setTitle(getTitle()+" (n-variants = "+n_variants+")");
			chart.setYAxisLabel("count");
			chart.setXAxisLabel("Chromosome");
			return exportChart(outputDir,chart);
			}
		}
	
	/***************************************************************************/
	private static class GatkDeNovo extends AbstractAnalyzer {
		private final String[] confDeNovos = new String[]{
				GATKConstants.hiConfDeNovo,
				GATKConstants.loConfDeNovo
				};
		private long n_variants=0L;
		private final Map<String,Counter<String>> sample2count=new HashMap<>();
		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			this.enabled = h.getInfoHeaderLine(this.confDeNovos[0])!=null && h.getInfoHeaderLine(this.confDeNovos[1])!=null;
			}
		
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			n_variants++;
			for(String info: this.confDeNovos) {
				if(!ctx.hasAttribute(info)) continue;
				for(final String sn: ctx.getAttributeAsStringList(info, "")) {
					if(StringUtils.isBlank(sn)) continue;
					Counter<String> c= this.sample2count.get(sn);
					if(c==null) {
						c = new Counter<>();
						this.sample2count.put(sn, c);
						}
					c.incr(info);
					}
				}
			}
		@Override
		public Set<Path> finish(final Path outputDir) throws IOException, XMLStreamException {
			if(sample2count.isEmpty()) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			final List<NamedSeries> series=new ArrayList<>();
			for(String sn: this.sample2count.keySet()) {
				final Counter<String> count = this.sample2count.get(sn);
				final List<NamedY> L2 =new ArrayList<>();
				for(String info: confDeNovos) {
					L2.add(new NamedY(info,count.count(info)));
					}
				series.add(new NamedSeries(sn, L2));
				}
			
			final BarPlot chart = new BarPlot(series);
			chart.setTitle(getTitle()+" (n-variants = "+n_variants+")");
			chart.setYAxisLabel("count");
			chart.setXAxisLabel("Sample");
			return exportChart(outputDir,chart);
			}
		}
	/***************************************************************************/
	private static class DragenDeNovo extends AbstractAnalyzer {
		private final String DN = "DN";
		private long n_variants=0L;
		private final Map<String,Counter<String>> sample2count=new HashMap<>();
		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			VCFFormatHeaderLine hdr= h.getFormatHeaderLine(DN);
			this.enabled = hdr!=null && hdr.getType().equals(VCFHeaderLineType.String);
			}
		
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			n_variants++;
			
			for(Genotype g:ctx.getGenotypes()) {
				if(!acceptGenotype(g)) continue;
				if(!g.hasExtendedAttribute(DN)) continue;
				final Object v = g.getExtendedAttribute(DN,"");
				if(v==null) continue;
				final String s= v.toString();
				if(StringUtils.isBlank(s) || !s.equals("DeNovo")) continue;
				Counter<String> c= this.sample2count.get(g.getSampleName());
				if(c==null) {
					c = new Counter<>();
					this.sample2count.put(g.getSampleName(), c);
					}
				c.incr(s + (g.isFiltered() || ctx.isFiltered()?".FILTERED":""));
				}
			}
		@Override
		public Set<Path> finish(final Path outputDir) throws IOException, XMLStreamException {
			if(sample2count.isEmpty()) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			final List<NamedSeries> series=new ArrayList<>();
			for(String sn: this.sample2count.keySet()) {
				final Counter<String> count = this.sample2count.get(sn);
				final List<NamedY> L2 =new ArrayList<>();
				for(String flag: count.keySet()) {
					L2.add(new NamedY(flag,count.count(flag)));
					}
				series.add(new NamedSeries(sn, L2));
				}
			
			final BarPlot chart = new BarPlot(series);
			chart.setTitle(getTitle()+" (n-variants = "+n_variants+")");
			chart.setYAxisLabel("count");
			chart.setXAxisLabel("Sample");
			return exportChart(outputDir,chart);
			}
		}

	/***************************************************************************/
	private static class SVTypeContig extends AbstractMultipleBarPlot {
		SVTypeContig() {
			name("svtype2contig");
			description("SVTYpe per contig");
			}
		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			this.enabled = h.getInfoHeaderLine(VCFConstants.SVTYPE)!=null;
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			final String svType = ctx.getAttributeAsString(VCFConstants.SVTYPE, null);
			if(StringUtils.isBlank(svType)) return;
			super.add(ctx.getContig(), svType);
			}
		}

	/***************************************************************************/
	private static class SampleCount extends AbstractAnalyzer {
		
		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			}
		@Override
		public void visit(final VariantContext ctx) {
			//nothing
			}
		@Override
		public Set<Path> finish(Path outputDir) throws  XMLStreamException,IOException {
			if(super.sampleToGroup.getGroups().stream().mapToInt(N->sampleToGroup.getSamplesForGroup(N).size()).allMatch(C->C<2)) return Collections.emptySet();
			
			final List<NamedSeries> L=new ArrayList<>();
			for(String g : super.sampleToGroup.getGroups()) {
				L.add(new NamedSeries(getLabelForGroup(g), new NamedY("count",super.sampleToGroup.getSamplesForGroup(g).size())));
				}
			if(L.isEmpty()) return Collections.emptySet();
			
			final BarPlot chart = new BarPlot(L);
			chart.setTitle("Samples");
			chart.setXAxisLabel("Group");
			chart.setYAxisLabel("Count");
			return exportChart(outputDir,chart);
			}
		}
	/***************************************************************************/
	private static class SnpEffAnalyzer extends AbstractAnalyzer {
		private long n_variants=0L;
		private Map<String,Counter<String>> group2count=new HashMap<>();
		private SnpEffPredictionParser snpEffParser = null;
		private SnpEffAnalyzer() {
			}
		
		
		
		@Override
		public void init(VCFHeader h, Map<String, String> properties, SampleToGroup sampleToGroup) {
			super.init(h, properties, sampleToGroup);
			if(h.hasGenotypingData()) {
				for(String gn:sampleToGroup.getGroups()) {
					this.group2count.put(gn, new Counter<>());
					}
				}
			else
				{
				this.group2count.put("all", new Counter<>());
				}
			snpEffParser = new SnpEffPredictionParserFactory(h).get();
			this.enabled= h.hasInfoLine(snpEffParser.getTag());
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			this.n_variants++;
			// find first prediction ?
			final String pred = this.snpEffParser.getPredictions(ctx).stream().flatMap(PRED->PRED.getSOTermsStrings().stream()).findFirst().orElse("undefined");
			if(ctx.hasGenotypes()) {
				final Set<String> group_in_variant = new HashSet<>();
				for(final Genotype gt: ctx.getGenotypes()) {
					if(!acceptGenotype(gt)) continue;
					if(!gt.hasAltAllele()) continue;
					for(String groupName: super.sampleToGroup.getGroupsForSample(gt.getSampleName())) {
						group_in_variant.add(groupName);
						}
					}
				for(String groupName: group_in_variant) {
					this.group2count.get(groupName).incr(pred);
					}
				}
			else
				{
				this.group2count.get("all").incr(pred);
				}
			}
			@Override
			public Set<Path> finish(final Path outputDir) throws IOException, XMLStreamException {
				if(n_variants==0L ) {
					LOG.warn("nothing found for "+getTitle());
					return Collections.emptySet();
					}
				final List<NamedSeries> L;
				if(super.vcfHeader.hasGenotypingData()) {
					L =new ArrayList<>(this.group2count.size());
					for(String grpName: this.group2count.keySet()) {
						Counter<String> count = this.group2count.get(grpName);
						final List<NamedY> L2 = count.entrySet().stream()
								.map(KV->new NamedY(KV.getKey(),KV.getValue()) )
								.collect(Collectors.toList());
						L.add(new NamedSeries(getLabelForGroup(grpName),L2));
						}
					}
				else
					{
					L = this.group2count.get("all").entrySet().stream()
							.map(KV->new NamedSeries(KV.getKey(),KV.getValue()) )
							.collect(Collectors.toList());
					}
				if(L.isEmpty()) return Collections.emptySet();
				
				final BarPlot chart = new BarPlot(L);
				chart.setTitle(getProperty("title","")+ " N-variants="+this.n_variants);
				chart.setXAxisLabel(vcfHeader.hasGenotypingData()?"collection":"type");
				chart.setYAxisLabel("Count");
				return exportChart(outputDir,chart);
				}
			}
	
	/***************************************************************************/
	private static class VariantTypeFraction extends AbstractAnalyzer {
		private long n_variants=0L;
		private Map<String,Counter<VariantContext.Type>> group2count=new HashMap<>();
		private final Predicate<VariantContext> variant_filter;
		private boolean do_normalize=false;
		private VariantTypeFraction(final Predicate<VariantContext> variant_filter) {
			this.variant_filter = variant_filter;
			}
		
		VariantTypeFraction setNormalize(boolean b) {
			this.do_normalize = b;
			return this;
			}
		
		@Override
		public void init(VCFHeader h, Map<String, String> properties, SampleToGroup sampleToGroup) {
			super.init(h, properties, sampleToGroup);
			for(String gn:sampleToGroup.getGroups()) {
				this.group2count.put(gn, new Counter<>());
				}
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			this.n_variants++;
			if(!this.variant_filter.test(ctx)) return;
			final Set<String> group_in_variant = new HashSet<String>();
			for(final Genotype gt: ctx.getGenotypes()) {
				if(!acceptGenotype(gt)) continue;
				if(!gt.hasAltAllele()) continue;
				for(String groupName: super.sampleToGroup.getGroupsForSample(gt.getSampleName())) {
					group_in_variant.add(groupName);
					}
				}
			for(String groupName: group_in_variant) {
				this.group2count.get(groupName).incr(ctx.getType());
				}
			}
		@Override
		public Set<Path> finish(final Path outputDir) throws IOException, XMLStreamException {
			if(n_variants==0L ) {
				LOG.warn("nothing found for "+getTitle());
				return Collections.emptySet();
				}
			final List<NamedSeries> L=new ArrayList<>(this.group2count.size());
			
			for(String grpName: this.group2count.keySet()) {
				Counter<VariantContext.Type> count = this.group2count.get(grpName);
				final double divide = do_normalize?this.n_variants:1.0;
				final List<NamedY> L2 = count.stream()
						.map(KV->new NamedY(KV.getKey().name(), KV.getValue()/divide) )
						.collect(Collectors.toList());
				if(L2.isEmpty()) continue;
				L.add(new  NamedSeries(getLabelForGroup(grpName), L2));
				}
			
			
			if(L.isEmpty()) return Collections.emptySet();
			
			final BarPlot chart = new BarPlot(L);
			chart.setTitle(getProperty("title","")+ " N-variants="+this.n_variants);
			chart.setXAxisLabel("collection");
			chart.setYAxisLabel(this.do_normalize?"Percentage":" Count");
			return exportChart(outputDir,chart);
			}
		}
	
	
	
	/***************************************************************************/
	private static class HomozygousPurityFraction extends AbstractAnalyzer {
		private long n_variants=0L;
		private int min_dp = 20;
		private static class CountPure {
			long n_pure_homvar = 0L;
			long n_pure_homref = 0L;
			long n_impure_homvar = 0L;
			long n_impure_homref = 0L;
			long sum() {
				return n_pure_homvar + n_pure_homref + n_impure_homvar + n_impure_homref;
				}
			}
		private Map<String,CountPure> group2count=new HashMap<>();
		private HomozygousPurityFraction() {
			}
		
		@Override
		public void init(VCFHeader h, Map<String, String> properties, SampleToGroup sampleToGroup) {
			super.init(h, properties, sampleToGroup);
			for(String gn:sampleToGroup.getGroups()) {
				this.group2count.put(gn, new CountPure());
				}
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			if(ctx.getNAlleles()!=2) return;
			this.n_variants++;
			for(final Genotype gt: ctx.getGenotypes()) {
				if(!acceptGenotype(gt)) continue;
				if(gt.hasDP() && gt.getDP() < this.min_dp) continue;
				if(!gt.hasAD()) continue;
				final int[] ad = gt.getAD();
				if(ad.length!=2 || (ad[0]==0 && ad[1]==0)) continue;
				final boolean is_pure=(ad[0]==0 || ad[1]==0);
				if(!gt.isHom()) continue;
				for(String groupName : super.sampleToGroup.getGroupsForSample(gt.getSampleName())) {
					final CountPure  c = this.group2count.get(groupName);
					if(gt.isHomRef()) {
						if(is_pure) {
							c.n_pure_homref++;
							}
						else
							{
							c.n_impure_homref++;
							}
						}
					else if(gt.isHomVar()) {
						if(is_pure) {
							c.n_pure_homvar++;
							}
						else
							{
							c.n_impure_homvar++;
							}
						}
					}
				}
		
			}
		@Override
		public Set<Path> finish(final Path outputDir) throws IOException, XMLStreamException {
			if(n_variants==0L && this.group2count.values().stream().allMatch(C->C.sum()==0L)) {
				LOG.warn("nothing found for "+getTitle());
				return Collections.emptySet();
				}
			final List<NamedSeries> L=new ArrayList<>(this.group2count.size());
			
			for(String grpName: this.group2count.keySet()) {
				final CountPure c = this.group2count.get(grpName);
				final List<NamedY> L2 = new ArrayList<>(4);
				L2.add(new NamedY("Impure HOM_REF",c.n_impure_homref));
				L2.add(new NamedY("Impure HOM_VAR",c.n_impure_homvar));
				L2.add(new NamedY("Pure HOM_REF",c.n_pure_homref));
				L2.add(new NamedY("Pure HOM_VAR",c.n_pure_homvar));
				L.add(new  NamedSeries(getLabelForGroup(grpName), L2));
				}
			
			
			if(L.isEmpty()) return Collections.emptySet();
			
			final BarPlot chart = new BarPlot(L);
			chart.setTitle(getProperty("title","Purity of homozygous genotypes with DP>="+this.min_dp+ " for Di-Alleleic variants")+" N-variants="+this.n_variants);
			chart.setXAxisLabel("collection");
			chart.setYAxisLabel("Count");
			return exportChart(outputDir,chart);
			}
		}
	/***************************************************************************/
	private static class FormatFiltersAnalyzer extends AbstractAnalyzer {
		private Map<String,Counter<String>> group2flt= new HashMap<>();
		private long n_variants = 0L;
;		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			if(!h.hasGenotypingData() ||  h.getFormatHeaderLine(VCFConstants.GENOTYPE_FILTER_KEY)==null) {
				this.enabled = false;
				}
			else
				{
				for(String g: super.sampleToGroup.getGroups()) {
					group2flt.put(g, new Counter<>());
					}
				}
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			n_variants++;
			for(final Genotype gt: ctx.getGenotypes()) {				
				if(!acceptGenotype(gt)) continue;
				final String filter = gt.isFiltered()?gt.getFilters():VCFConstants.PASSES_FILTERS_v4;
				for(String grpName: super.sampleToGroup.getGroupsForSample(gt.getSampleName())) {
					group2flt.get(grpName).incr(filter);	
					}
				}
			}
		@Override
		public Set<Path> finish(final Path outputDir) throws IOException, XMLStreamException {
			if(n_variants==0) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			final List<NamedSeries> series=new ArrayList<>();
			for(String groupName: this.group2flt.keySet()) {
				final List<NamedY> L2 =new ArrayList<>();
				final Counter<String> filters = this.group2flt.get(groupName);
				for(String filter : filters.keySet()) {
					L2.add(new NamedY(filter, filters.count( filter)));
					}
				series.add(new NamedSeries(getLabelForGroup(groupName), L2));
				}
			
			final BarPlot chart = new BarPlot(series);
			chart.setTitle(getTitle()+" (n-variants = "+n_variants+")");
			chart.setYAxisLabel("count Filters");
			chart.setXAxisLabel("Group");
			return exportChart(outputDir,chart);
			}
		}
	/***************************************************************************/
	private static class InfoFiltersAnalyzer extends AbstractAnalyzer {
		private final Counter<String> filters = new Counter<String>();
		private long n_pass = 0L;
		private long n_variants = 0L;
;		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			if(h.getFilterLines().isEmpty()) {
				this.enabled = false;
				}
			else
				{
				for(VCFFilterHeaderLine g: h.getFilterLines()) {
					filters.initializeIfNotExists(g.getID());
					}
				}
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			if(!ctx.isFiltered()) {
				n_pass++;
				return;
				}
			n_variants++;
			for(final String flt: ctx.getFilters()) {				
				filters.incr(flt);
				}
			}
		@Override
		public Set<Path> finish(final Path outputDir) throws IOException, XMLStreamException {
			if(n_variants==0) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			final List<NamedSeries> series=new ArrayList<>();
			for(String filter: this.filters.keySet()) {
				series.add(new NamedSeries(filter,this.filters.count(filter)));
				}
			
			final BarPlot chart = new BarPlot(series);
			chart.setTitle(getTitle()+" (n-variants = "+n_variants+") PASS="+n_pass+" ("+((n_pass/(n_pass+n_variants))*100.0)+"%)");
			chart.setYAxisLabel("count Filters");
			chart.setXAxisLabel("Filter");
			return exportChart(outputDir,chart);
			}
		}

	/*********************************************************************/
	private static class SVLen extends AbstractAnalyzer {
		private final String svType;
		private boolean with_genotypes=true;
		private final Map<String,Counter<Integer>> grp2sizes = new HashMap<>();
		SVLen(final String svType) {
			this.svType=svType;
			}
		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			this.enabled = h.getInfoHeaderLine(VCFConstants.SVTYPE)!=null;
			this.with_genotypes = h.hasGenotypingData();
			if(this.with_genotypes && h.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY)==null) {
				with_genotypes=false;
				}
			
			if(with_genotypes) {
				for(String grpName: s2g.getGroups()) {
					grp2sizes.put(grpName, new Counter<>());
					}
				}
			else
				{
				grp2sizes.put("ALL", new Counter<>());
				}
			}
		
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			final String st = ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
			if(StringUtils.isBlank(st) || st.equalsIgnoreCase("BND")) return;
			int svLen;
			if(ctx.hasAttribute("SVLEN")) {
				svLen = ctx.getAttributeAsInt("SVLEN", 0);
				}
			else if(ctx.hasAttribute(VCFConstants.END_KEY)) {
				svLen =  ctx.getLengthOnReference();
				}
			else
				{
				return;
				}
			svLen = Math.abs(svLen);
			if(with_genotypes) {
				final Set<String> groups  = ctx.getGenotypes()
						.stream()
						.filter(G->G.hasAltAllele())
						.flatMap(G->super.sampleToGroup.getGroupsForSample(G.getSampleName()).stream())
						.collect(Collectors.toSet());
				for(String groupName: groups) {
					grp2sizes.get(groupName).incr(svLen);
					}
				}
			else
				{
				grp2sizes.get("ALL").incr(svLen);
				}
			}
	
		@Override
		public Set<Path> finish(final Path outputDir) throws IOException, XMLStreamException {
			final List<SeriesXY> series=new ArrayList<SeriesXY>();
			this.grp2sizes.entrySet().forEach(KV->{
				final SeriesXY L = new SeriesXY(
						(this.with_genotypes? getLabelForGroup( KV.getKey()):"ALL"),
						KV.getValue()
							.entrySet()
							.stream()
							.map(KV2->new DataXY(KV2.getKey(), KV2.getValue()))
							.collect(Collectors.toList())
						);
				L.sort();
				series.add(L);
				});
			
			final ScatterXY chart = new ScatterXY(series);
			chart.setLogY(true);
			chart.setLogX(true);
			chart.setTitle("SVLEN "+this.svType);
			chart.setYAxisLabel("log(Count)");
			chart.setXAxisLabel("log(SVLen)");
			return exportChart(outputDir,chart);
			}
		}
	
	
	/*********************************************************************/
	private static class SampleToSVTypes extends AbstractMultipleBarPlot {
		SampleToSVTypes() {
			xlab("Samples");
			}
		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			this.enabled = h.getInfoHeaderLine(VCFConstants.SVTYPE)!=null && h.getNGenotypeSamples()>0;
			}
	
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			final String st = ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
			if(StringUtils.isBlank(st)) return;
			for(final Genotype gt:ctx.getGenotypes()) {
				if(!getGenotypePredicate().test(gt)) continue;
				if(gt.getAlleles().stream().allMatch(A->A.isReference() || A.isNoCall())) continue;
				super.add(gt.getSampleName(), st);
				}
			}
		}
	/*********************************************************************/
	private static class SampleToSVTypeLen extends AbstractBoxPlot {
		private final String svType;
		SampleToSVTypeLen(final String svType) {
			this.svType = svType;
			name("SVLEN."+svType);
			ylab("SVLEN."+svType);
			}
		@Override
		public void init(final VCFHeader h,Map<String,String> props,SampleToGroup s2g) {
			super.init(h, props, s2g);
			this.enabled = h.getInfoHeaderLine(VCFConstants.SVTYPE)!=null && h.getNGenotypeSamples()>0;
			}
		
		@Override
		public void visit(final VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			final String st = ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
			if(!this.svType.equals(st)) return;
			final int svLen;
			if(!ctx.hasAttribute("SVLEN")) {
				svLen = ctx.getLengthOnReference();
				}
			else {
				svLen = ctx.getAttributeAsInt("SVLEN", 0);
				}
			for(Genotype g:ctx.getGenotypes()) {
				if(!getGenotypePredicate().test(g)) continue;
				if(g.getAlleles().stream().allMatch(G->G.isReference() || G.isNoCall())) continue;
				add(g.getSampleName(), svLen);
				}
			}
		}
	
	/*********************************************************************/
	/**
	 * Base class used to get metrics for FORMAT numeric stuff like DP or GQ
	 */
	private abstract class AbstractFormatNumeric extends AbstractAnalyzer {
		private final String formatTag;
		private long n_variants=0L;
		protected final Map<String,DataPointAverage> sample2count= new HashMap<>();
		protected AbstractFormatNumeric(String formatTag) {
			this.formatTag = formatTag;
			}
		@Override
		public void init(VCFHeader h, Map<String, String> props, SampleToGroup s2g) {
			super.init(h, props, s2g);
			if(!h.hasGenotypingData() || h.getFormatHeaderLine(getFormatKey())==null) {
				this.enabled=false;
				}
			for(String sn: h.getGenotypeSamples()) {
				sample2count.put(sn, new DataPointAverage());
				}
			}
	
		
		protected final  String getFormatKey() {
			return this.formatTag;
			}
		
				
		protected abstract void visit(Genotype gt);
		
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			n_variants++;
			for(Genotype gt:ctx.getGenotypes()) {
				if(!acceptGenotype(gt)) continue;
				visit(gt);
				}
			}
		
		@Override
		public Set<Path> finish(Path outputDir) throws IOException, XMLStreamException {
			if(sample2count.isEmpty()) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			final List<NamedSeries> series = new ArrayList<NamedSeries>(super.sampleToGroup.getGroupsCount());
			for(String groupName: super.sampleToGroup.getGroups()) {
				final List<NamedY> L2 =new ArrayList<>();
				for(String sn: this.sample2count.keySet()) {
					if(!super.sampleToGroup.hasSampleInGroup(sn,groupName)) continue;
					final OptionalDouble od = this.sample2count.get(sn).get();
					if(!od.isPresent()) continue;
					L2.add(new NamedY(sn,od.getAsDouble()));
					}
				if(L2.isEmpty()) continue;
				
				series.add(new NamedSeries(getLabelForGroup(groupName), L2));
				}
			if(series.isEmpty()) return Collections.emptySet();
			final BoxPlotChart chart = new BoxPlotChart(series);
			chart.sortOnName();
			chart.setTitle(getTitle()+" (n-variants = "+n_variants+")");
			chart.setYAxisLabel("count");
			chart.setXAxisLabel("Group");
			return exportChart(outputDir,chart);
			}
		
		}


	/*********************************************************************/
	private class SampleToGQ extends AbstractFormatNumeric {
		SampleToGQ() {
			super(VCFConstants.GENOTYPE_QUALITY_KEY);
			}
		
		@Override
		protected void visit(final Genotype gt) {
			if(!gt.hasGQ()) return;
			super.sample2count.get(gt.getSampleName()).accept(gt.getGQ());
			}
		}
	/*********************************************************************/
	private class SampleToDP extends AbstractFormatNumeric {
		SampleToDP() {
			super(VCFConstants.DEPTH_KEY);
			}
		@Override
		protected void visit(final Genotype gt) {
			if(!gt.hasDP()) return;
			super.sample2count.get(gt.getSampleName()).accept(gt.getDP());
			}
		}
	/*********************************************************************/
	
	/*********************************************************************/
	private static class CountVariantsManhattanPlot extends AbstractManhattanPlot {
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
				visit(ctx.getType().name(),ctx,  1);
			}
		}
	/*********************************************************************/
	private static class CountSamplesVariantsManhattanPlot extends AbstractManhattanPlot {
		@Override
		public void init(VCFHeader h,Map<String,String> props, SampleToGroup s2g) {
			super.init(h,props,s2g);
			if(super.enabled && !h.hasGenotypingData()) {
				super.enabled=false;
				}
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			for(Genotype g:ctx.getGenotypes()) {
				if(!getGenotypePredicate().test(g)) continue;
				if(g.getAlleles().stream().allMatch(G->G.isReference() || G.isNoCall())) continue;
				visit(g.getSampleName(),ctx,  1);
				}
			}
		}
	/***************************************************************************/
	private static class SampleToGenotypeType extends AbstractAnalyzer {
		private final boolean ignore_hom_ref;
		private long n_variants = 0L;
		private boolean do_normalize=false;
		private final Map<String,Counter<GenotypeType>> group2gtype = new HashMap<>();
		SampleToGenotypeType(boolean ignore_hom_ref) {
			this.ignore_hom_ref = ignore_hom_ref;
			}
		
		SampleToGenotypeType setNormalizeFlag(boolean b) {
			this.do_normalize = b;
			return this;
			}
		
		@Override
		public void init(VCFHeader h, Map<String, String> props, SampleToGroup sample2group) {
			super.init(h, props, sample2group);
			if(!h.hasGenotypingData() || h.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY)==null) {
				super.enabled = false;
				}
			for(String groupName: sample2group.getGroups()) {
				this.group2gtype.put(groupName, new Counter<>());
				}
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!ctx.hasGenotypes()) return;
			if(!acceptVariant(ctx)) return;
			n_variants++;
			for(Genotype gt: ctx.getGenotypes()) {
				if(!acceptGenotype(gt)) continue;
				if(ignore_hom_ref && gt.isHomRef()) continue;
				for(String groupname: super.sampleToGroup.getGroupsForSample(gt.getSampleName())) {
					group2gtype.get(groupname).incr(gt.getType());
					}
				}
			}
		@Override
		public Set<Path> finish(Path outputDir) throws IOException, XMLStreamException {
			if(this.group2gtype.values().stream().noneMatch(C->C.getTotal()>0L)) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			boolean use_fraction= group2gtype.keySet().stream().anyMatch(G->super.sampleToGroup.getSamplesForGroup(G).size()>1);
			final List<NamedSeries> series=new ArrayList<>();
			for(String groupName: this.group2gtype.keySet()) {
				final List<NamedY> L2 =new ArrayList<>();
				final double total;
				if(do_normalize) {
						total = Math.max(1.0, (use_fraction?this.group2gtype.get(groupName)
						.entrySet()
						.stream()
						.filter(KV->!ignore_hom_ref || !KV.getKey().equals(GenotypeType.HOM_REF))
						.mapToDouble(KV->KV.getValue())
						.sum():1.0)
						);
						}
					else
						{
						total = 1.0;
						}
				
				for(GenotypeType gtype: GenotypeType.values()) {
					if(this.ignore_hom_ref && gtype.equals(GenotypeType.HOM_REF)) continue;
					L2.add(new NamedY(gtype.name(),this.group2gtype.get(groupName).count(gtype)/total));
					}
				series.add(new NamedSeries(getLabelForGroup(groupName), L2));
				}
			
			final BarPlot chart = new BarPlot(series);
			chart.sortOnName();
			chart.setTitle(getTitle()+" (n-variants = "+n_variants+")");
			chart.setYAxisLabel("count");
			chart.setXAxisLabel("Group");
			return exportChart(outputDir,chart);
			}
		}
	/***************************************************************************/
	private static class ADRatioAnalyzer extends AbstractAnalyzer {
		private final GenotypeType gtype;
		private final AutoMap<String,Counter<Double>,Counter<Double>> group2count= AutoMap.make(SN->new Counter<>());
		private final boolean non_pure_only;
		private final DoubleRounder rounder=new DoubleRounder(2);
		private boolean do_normalize=false;
		ADRatioAnalyzer(final GenotypeType gtype,boolean non_pure_only) {
			this.gtype=gtype;
			this.non_pure_only = non_pure_only;
			}
		ADRatioAnalyzer(final GenotypeType gtype) {
			this(gtype,false);
			}
		ADRatioAnalyzer setNormalizeFlag(boolean b) {
			this.do_normalize= b;
			return this;
			}
		@Override
		public void init(VCFHeader h, Map<String, String> properties, SampleToGroup sampleToGroup) {
			super.init(h, properties, sampleToGroup);
			if(!h.hasGenotypingData()) {
				super.enabled=false;
				}
			if(h.getFormatHeaderLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS)==null) {
				super.enabled=false;
				}
			}
		
		@Override
		public void visit(final VariantContext ctx) {
			if(!ctx.isBiallelic()) {
				return;
				}
			if(!ctx.hasGenotypes()) {
				return ;
				}
			if(!acceptVariant(ctx)) {
				return;
				}
			for(Genotype gt: ctx.getGenotypes()) {
				if(!gt.getType().equals(this.gtype)) {
					continue;
					}
				if(!acceptGenotype(gt)) continue;
				if(!gt.hasAD()) {
					//System.err.println("no AD");
					continue;
					}
				final int[] ad=gt.getAD();
				if(ad==null || ad.length!=2) {
					//System.err.println("no AD.length==2");
					continue;
					}
				if(this.non_pure_only) {
					switch(this.gtype) {
						case HOM_REF: if(ad[1]==0) continue;break;//pure because no ALT
						case HOM_VAR: if(ad[0]==0) continue;break;//putr because no REF
						default: throw new IllegalStateException(this.gtype.name());
						}
					}
				
				final int sum = ad[0]+ad[1];
				if(sum<=0) continue;
				final double f = this.rounder.applyAsDouble(ad[1]/(double)sum);
				for(final String groupName : this.sampleToGroup.getGroupsForSample(gt.getSampleName())) {
					this.group2count.insert(groupName).incr(f);
					}
				}
			}
		@Override
		public Set<Path> finish(Path outputDir) throws IOException,XMLStreamException {
			if(this.group2count.isEmpty()) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			
			final List<SeriesXY> series=new ArrayList<SeriesXY>();
			this.group2count.entrySet().forEach(KV->{
				final SeriesXY L = new SeriesXY(
						getLabelForGroup( KV.getKey()),
						KV.getValue()
							.entrySet()
							.stream()
							.map(KV2->new DataXY(KV2.getKey(), KV2.getValue()))
							.collect(Collectors.toList())
						);
				if(this.do_normalize) L.normalize();	
				L.sort();
				series.add(L);
				});
			
			final ScatterXY chart = new ScatterXY(series);
			chart.setTitle(getTitle());
			chart.setYAxisLabel((do_normalize?"Normalized ":"")+"Count");
			chart.setXAxisLabel("AD Ratio ALT/(REF+ALT)");
			return exportChart(outputDir,chart);
			}
		}
	
	
	/*********************************************************************/
	private static class CrossContaminationAnalyzer extends AbstractAnalyzer {
		private int factor = 100;
		private final AutoMap<String,Counter<Integer>,Counter<Integer>> group2count100 = AutoMap.make(SN->new Counter<>());
		private boolean do_normalize=false;
		@Override
		public void init(VCFHeader h, Map<String, String> properties, final SampleToGroup sampleToGroup) {
			super.init(h, properties, sampleToGroup);
			if(!h.hasGenotypingData()) {
				super.enabled=false;
				}
			if(h.getFormatHeaderLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS)==null) {
				super.enabled=false;
				}
			if(sampleToGroup.getGroups().stream().allMatch(G->sampleToGroup.getSamplesForGroup(G).size()<=1)) {
				super.enabled=false;
				}
			}
		
		CrossContaminationAnalyzer setNormalizeFlag(boolean b) {
			this.do_normalize= b;
			return this;
			}
		
		@Override
		public void visit(VariantContext ctx) {
			if(!acceptVariant(ctx)) return;
			if(!ctx.isBiallelic()) {
				return;
				}
			Genotype singleton = findSingleton(ctx);
			if(singleton==null || !acceptGenotype(singleton) || !singleton.isHet()) return;
			for(String groupName: super.sampleToGroup.getGroupsForSample(singleton.getSampleName())) {
				final Set<String> samples_in_group= super.sampleToGroup.getSamplesForGroup(groupName);
				for(String sn : samples_in_group) {
					final Genotype gt = ctx.getGenotype(sn);
					if(gt==null) continue;
					if(!gt.isHomRef()) continue;
					if(gt.getSampleName().equals(singleton.getSampleName())) continue;// paranoid, useless
					if(!gt.hasAD()) continue;
					final int[] ad= gt.getAD();
					if(ad.length!=2) continue;
					int sum = ad[0]+ad[1];
					if(sum==0) continue;
					final double f = ad[1]/(double)sum;
					final int f_as_int10 = (int)Math.floor(f*(double)this.factor);
					this.group2count100.insert(groupName).incr(f_as_int10);	
					}
				}
			}
		@Override
		public Set<Path> finish(Path outputDir) throws IOException,XMLStreamException {
			if(this.group2count100.isEmpty()) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			
			final List<SeriesXY> series=new ArrayList<SeriesXY>();
			this.group2count100.entrySet().forEach(KV->{
				final SeriesXY L = new SeriesXY(
						getLabelForGroup( KV.getKey()),
						KV.getValue()
							.entrySet()
							.stream()
							.map(KV2->new DataXY(KV2.getKey()/(double)this.factor, KV2.getValue()))
							.collect(Collectors.toList())
						);
				if(this.do_normalize) L.normalize();
				L.sort();
				series.add(L);
				});
			
			final ScatterXY chart = new ScatterXY(series);
			chart.setTitle(getTitle());
			chart.setYAxisLabel((this.do_normalize?"Normalized ":"")+"Count");
			chart.setXAxisLabel("AD Ratio ALT/(REF+ALT)");
			return exportChart(outputDir,chart);
			}
		}
	
	/*********************************************************************/
	private static class CountGenotypesPerVariantAnalyzer extends AbstractAnalyzer {
		private final Map<String,Counter<Integer>> group2count = new HashMap<>();
		private int n_samples=0;
		@Override
		public void init(VCFHeader h, Map<String, String> properties, final SampleToGroup sampleToGroup) {
			super.init(h, properties, sampleToGroup);
			if(!h.hasGenotypingData() || (this.n_samples =h.getNGenotypeSamples())<=1) {
				super.enabled=false;
				}
			if(h.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY)==null) {
				super.enabled=false;
				}
			if(sampleToGroup.getGroups().stream().allMatch(G->sampleToGroup.getSamplesForGroup(G).size()<=1)) {
				super.enabled=false;
				}
			for(String gtpName: sampleToGroup.getGroups()) {
				group2count.put(gtpName, new Counter<>());
				}
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!acceptVariant(ctx)) {
				//System.err.println("skip");
				return;
				}
			for(final String grp:group2count.keySet() ) {
				int n=0;
				for(String sn:super.sampleToGroup.getSamplesForGroup(grp)) {
					Genotype gt = ctx.getGenotype(sn);
					if(gt==null || !gt.hasAltAllele()) continue;
					n++;
					}
				group2count.get(grp).incr(n);
				}
			}
		@Override
		public Set<Path> finish(Path outputDir) throws IOException,XMLStreamException {
			if(this.group2count.values().stream().allMatch(C->C.getMaxCount().orElse(0L)<=0)) {
				LOG.warn("nothing found for "+this.getName()+" "+getTitle());
				return Collections.emptySet();
				}
			final List<SeriesXY> L=new ArrayList<SeriesXY>(this.group2count.size());
			for(String grpName: this.group2count.keySet()) {
				final Counter<Integer> c = this.group2count.get(grpName); 
				final SeriesXY L2 = new SeriesXY(
						getLabelForGroup(grpName),
						c .entrySet()
							.stream()
							.map(KV2->new DataXY(KV2.getKey(),KV2.getValue()))
							.collect(Collectors.toList())
						);
				L2.sort();
				L.add(L2);
				}
	
			
			final ScatterXY chart = new ScatterXY(L);
			chart.setTitle(getTitle());
			chart.setYAxisLabel("log(count Variants)");
			chart.setLogY(true);
			chart.setXAxisLabel("Number of genotypes carrying an ALT allele per variant");
			return exportChart(outputDir,chart);
			}
		}


	private void loadPhenotypes(final VCFHeader header) throws IOException {
		if(!header.hasGenotypingData()) return;
		if(sample2catPath!=null) {
			this.sampleToGroup
				.load(this.sample2catPath)
				.retainSamples(header)
				;
		}
	}
	
	private static Genotype findSingleton(final VariantContext ctx) {
		Genotype single = null;
		for(Genotype g:ctx.getGenotypes()) {
			if(!g.hasAltAllele()) continue;
			if(single!=null) return null;
			single=g;
			}
		return single;
		}
	
	private static boolean isSingletonVariant(final VariantContext ctx) {
		return findSingleton(ctx)!=null;
		}
	
	@Override
	public int doWork(final List<String> args) {
		final List<Analyzer> modules =new ArrayList<>();
		final Map<String,String> properties = new HashMap<>();
		
		for(final String svType: new String[] {"DEL","INV","DUP","INS"}) {
			modules.add(new SampleToSVTypeLen("svType").
					name("SVLEN("+svType+") per sample").
					description("SVLEN("+svType+") per sample").
					xlab("Sample").
					ylab("avg(SVLEN)")
					);
			}
		modules.add(
			new SVTypeContig().
				name("SVTYPE per contig").
				description("SVType per contig").
				xlab("Contig").
				ylab("SVTYPE")
			);
		
		modules.add(new SampleToSVTypes().
				name("Sample To SVTYPE").
				description("Sample To SVTYPE").
				xlab("Samples").
				ylab("SVTYPE")
			);
		
		
		
		modules.add(new RangeBarPlot(GATKConstants.QD_KEY, 1,false).
				setProperty("filename", "info_QD").
				name(GATKConstants.QD_KEY).
				description("Variant Confidence (QUAL) / Quality by Depth.").
				xlab(GATKConstants.QD_KEY).
				ylab("Count Variants")
				);
		

		
		
		/*
		modules.add(new RangeBarPlot(VCFConstants.DEPTH_KEY, 1).
				logX(true).
				name(VCFConstants.DEPTH_KEY).
				description("DEPTH per variant").
				xlab(VCFConstants.DEPTH_KEY).
				ylab("Count Variants")
				);*/

		modules.add(new CountVariantsManhattanPlot().
				name("manhattan count variants").
				description("manhattan count variants").
				xlab("Genome").
				ylab("Count Variants")
				);

		modules.add(new CountSamplesVariantsManhattanPlot().
				name("manhattan count Samples variants").
				description("manhattan count variants per sample").
				xlab("Genome").
				ylab("Count Variants")
				);
		
		

		
		modules.clear();//TODO fix me
		
		for(final String svType: new String[] {"DEL","INV","DUP"}) {
			modules.add( new SVLen(svType) );
		}
		
		modules.add(new GatkDeNovo()
				.setProperty("title","GATK DeNovo")
				.setProperty("filename","gatk_denovo")
				);
		modules.add(new DragenDeNovo()
				.setProperty("title","Dragen DeNovo")
				.setProperty("filename","dragen_denovo")
				);
		
		modules.add(new FormatFiltersAnalyzer()
				.setProperty("title","Sample FILTERs")
				.setProperty("filename","format_flt")
				);
		modules.add(new InfoFiltersAnalyzer()
				.setProperty("title","Variants FILTERs")
				.setProperty("filename","variant_filter")
				);
		modules.add(new SampleToGQ()
					.setProperty("title","Genotype Quality")
					.setProperty("filename","format_genotype_quality")
				);
		modules.add(new SampleToDP()
				.setProperty("title","Sample Depth")
				.setProperty("filename","format_genotype_depth")
				);

		modules.add(new CountGenotypesPerVariantAnalyzer()
				.setProperty("title","Number of Genotypes per Variant")
				.setProperty("filename","count_genotypes_per_variant")
				);
		
		
		modules.add(new RangeBarPlot(GATKConstants.FS_KEY, 1,true).
			setProperty("title","Phred-scaled p-value using Fisher's exact test to detect strand bias").
			setProperty("filename", "info_FS")
			);
		
		modules.add(new RangeBarPlot(GATKConstants.SOR_KEY, 2,false).
			setProperty("title","Symmetric Odds Ratio of 2x2 contingency table to detect strand bias").
			setProperty("filename", "info_SOR")
			);
		
		modules.add(new RangeBarPlot(GATKConstants.MQ_KEY, 1,false).
			setProperty("title","Mean square mapping quality over all the reads at the site").
			setProperty("filename", "info_MQ")
			);

		modules.add(new RangeBarPlot(GATKConstants.MQRankSum_KEY, 2,false).
			setProperty("title","Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities").
			setProperty("filename", "info_MQRankSum")
			);
		
		modules.add(new RangeBarPlot(GATKConstants.ReadPosRankSum_KEY, 2,false)
			.setProperty("title", "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias").
			setProperty("filename", "info_ReadPosRankSum")
			);
		modules.add(new SingletonAnalyzer()
				.setProperty("title","Count singletons")
				.setProperty("filename", "singletons_by_groups")
				);
		
		modules.add(
				new CountPerContig(SSR->SSR.getContig().matches("(chr)?[0-9XY]+"))
					.setProperty("title", "Variant per chromosome")
					.setProperty("filename", "variants_per_chromosome")
				);
		
		modules.add(
				new SampleToGenotypeType(true)
					.setProperty("title", "Genotype Types (without HOM_REF)")
					.setProperty("filename", "genotype_type_ignore_hom_ref")
				);
		
		modules.add(
				new SampleToGenotypeType(false)
					.setProperty("title", "Genotype Types")
					.setProperty("filename", "genotype_type_all")
				);

		
		modules.add(
				new ADRatioAnalyzer(GenotypeType.HET)
					.setProperty("title", "AD Ratio for Singletons HET genotypes, Diallelic Variant")
					.setProperty("filename", "AD_ratio_singleton_HET")
					.setAcceptVariant(VcfStats::isSingletonVariant)
				);
		
		modules.add(
				new ADRatioAnalyzer(GenotypeType.HET)
					.setProperty("title", "AD Ratio for HET genotypes, Diallelic Variant")
					.setProperty("filename", "AD_ratio_HET")
				);
		
		modules.add(
				new ADRatioAnalyzer(GenotypeType.HOM_REF)
					.setProperty("title", "AD Ratio for HOM_REF genotypes, Diallelic Variant")
					.setProperty("filename", "AD_ratio_HOM_REF")
				);
		modules.add(
				new ADRatioAnalyzer(GenotypeType.HOM_VAR)
					.setProperty("title", "AD Ratio for HOM_VAR genotypes, Diallelic Variant")
					.setProperty("filename", "AD_ratio_HOM_VAR")
				);
		
		
		modules.add(
				new ADRatioAnalyzer(GenotypeType.HOM_REF,true)
					.setProperty("title", "AD Ratio for non-pure HOM_REF genotypes, Diallelic Variant")
					.setProperty("filename", "AD_ratio_impure_HOM_REF")
				);
		modules.add(
				new ADRatioAnalyzer(GenotypeType.HOM_VAR,true)
					.setProperty("title", "AD Ratio for non-pure HOM_VAR genotypes, Diallelic Variant")
					.setProperty("filename", "AD_ratio_impure_HOM_VAR")
				);
		
		modules.add(
				new VariantTypeFraction(VC->VC.getNAlleles()>2)
					.setProperty("title", "Multi-allelic Fraction")
					.setProperty("filename", "multi_allelic_fraction")
				);
		
		modules.add(
				new VariantTypeFraction(VC->true)
					.setProperty("title", "Variant Types Fraction")
					.setProperty("filename", "variant_type_fraction")
				);
		
		modules.add(
				new CrossContaminationAnalyzer()
					.setProperty("title", "Cross-contamination: For a singleton HET, what is the AD ratio of the HOM_REF samples in the same group")
					.setProperty("filename", "cross_contamination_singleton_het")
				);
		
		modules.add(
				new HomozygousPurityFraction()
					.setProperty("filename", "homozygous_purity")
				);
		
		modules.add(new SampleCount()
					.setProperty("filename", "sample_count")
				);
		modules.add(new SnpEffAnalyzer()
					.setProperty("filename", "snpeff")
					.setProperty("title", "SNPEFF  predictions")
				);
		
		// update filename with prefix, update title
		for(Analyzer analyzer:modules) {
			analyzer.setProperty("filename", this.prefix+analyzer.getProperty("filename", ""));
			if(!StringUtils.isBlank(this.extra_title)) {
				analyzer.setProperty("title", this.extra_title + " : " + analyzer.getProperty("title", ""));
				}
			}
		
		final String input = oneFileOrNull(args);
		
		// remove modules
		modules.removeIf(M->Arrays.stream(this.moduleExcludeStr.split("[,; \t:]")).anyMatch(S->S.equalsIgnoreCase(M.getName())));
		try {
			if(this.list_modules) {
				for(Analyzer analyzer:modules) {
					System.out.print(analyzer.getName());
					System.out.print("\t");
					System.out.println(analyzer.getTitle());
					}	
					
				return 0;
				}
		
			try(VCFIterator iter= super.openVCFIterator(input)) {	
				final VCFHeader header=iter.getHeader();
				loadPhenotypes(header);
				
				// fill sample without group
				for(String sn:header.getGenotypeSamples()) {
					if(sampleToGroup.hasSample(sn)) continue;
					LOG.info("adding sample "+sn+" to group "+this.fill_sample_name_unknown_group);
					if(this.fill_sample_name_unknown_group.equals("*")) {
						this.sampleToGroup.putSampleGroup(sn, sn);
						}
					else
						{
						this.sampleToGroup.putSampleGroup(sn, this.fill_sample_name_unknown_group);
						}
					}
				LOG.info("n-samples:"+this.sampleToGroup.getSamplesCount()+" -ngroups:"+this.sampleToGroup.getGroupsCount());
				
				
				for(Analyzer analyzer:modules) {
					analyzer.init(header,properties,this.sampleToGroup);
					if(!analyzer.isEnabled()) {
						LOG.warn("module "+analyzer.getName()+" will be disabled. ["+analyzer.getClass().getSimpleName()+"]");
						}
					}
				modules.removeIf(M->!M.isEnabled());
				if(modules.isEmpty()) {
					LOG.warn("no module was enabled");
					}
				final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
				while(iter.hasNext()) {
					final VariantContext ctx = progress.watch(iter.next());
					for(Analyzer analyzer:modules) {
						analyzer.visit(ctx);
						}
					}
				progress.finish();
				
				
				final Set<Path> generated_files = new HashSet<>();
				for(Analyzer analyzer:modules) {
					generated_files.addAll( analyzer.finish(outputDirectory)) ;
					}				
				}
			return 0;
		} catch (final Throwable e) {
			e.printStackTrace();
			LOG.error(e);
			return -1;
			}
		finally
			{
			
			}
		}
	
			
	public static void main(final String[] args)
		{
		new VcfStats().instanceMainWithExit(args);
		}
	}
