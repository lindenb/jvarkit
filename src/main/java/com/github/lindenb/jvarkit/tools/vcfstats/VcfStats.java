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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.DecimalFormat;
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
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.BiFunction;
import java.util.function.DoubleConsumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gatk.GATKConstants;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
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

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
		
	@Parameter(names={"--categories","--phenotypes"},description="A tab delimited file (sample)<tab>(category)",hidden = true)
	private Path sample2catPath = null;

	
	@Parameter(names={"-exclude","--exclude"},description="name of modules to be excluded",hidden=true)
	private String moduleExcludeStr = "";
	@Parameter(names={"--list"},description="list available modules and exit",help = true)
	private boolean list_modules = false;	
	@Parameter(names={"--prefix"},description="file prefix")
	private String prefix = "";
	@DynamicParameter(names={"-D"},description="other parameters.")
	private Map<String,String> __dynaParams = new HashMap<>();
	

	private final Map<String, List<String>> phenotype2samples = new TreeMap<>();
	private final Map<String,String> sample2phenotype = new TreeMap<>();
	private final AttributeMap att = AttributeMap.verbose(AttributeMap.wrap(this.__dynaParams), (S)->{
		LOG.info("undefined parameter \"" + S + "\". Using default.");
		});
	
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
	private interface Analyzer {
		void init(VCFHeader h);
		void visit(final VariantContext ctx);
		void finish(PrintWriter pw);
		public String getName();
		public String getDescription();
		public boolean isEnabled();
		}
	
	/** AbstractAnalyzer **/
	private abstract class AbstractAnalyzer implements Analyzer {
		protected boolean enabled=true;
		private String name="undefined";
		private String description=null;
		private String outputFilename =null;
		private String xlab = "todo_axis_x";
		private String ylab = "todo_axis_y";
		protected Predicate<VariantContext> _acceptVariant = V->true;
		protected Predicate<Genotype> _acceptGT = V->true;

		@Override
		public boolean isEnabled() {
			return enabled;
			}
		public String device() {
			return device(-1);
			}
		public String device(int n) {
			String w = ",width="+String.valueOf((int)Math.max(7.0,n*7/40.0));
			if(n<=0) w="";
			return "pdf("+quote(getOuputFilename()+".pdf")+w+")";
			}
		public String getOuputFilename() {
			return VcfStats.this.prefix+ ( StringUtils.isBlank(outputFilename)?getName().trim().replaceAll("[^A-Za-z0-9]+","_"):outputFilename);
			}
		@Override
		public final String getName() {
			return name;
			}
		@Override
		public String getDescription() {
			return StringUtils.isBlank(description)?getName():this.description;
			}
		AbstractAnalyzer name(final String s) {
			this.name=s;
			return this;
			}
		AbstractAnalyzer description(final String s) {
			this.description= s;
			return this;
			}
		AbstractAnalyzer filename(final String s) {
			this.outputFilename= s;
			return this;
			}
		AbstractAnalyzer xlab(final String s) {
			this.xlab= s;
			return this;
			}
		AbstractAnalyzer ylab(final String s) {
			this.ylab= s;
			return this;
			}
		AbstractAnalyzer acceptVariant(final Predicate<VariantContext> f) {
			this._acceptVariant= f;
			return this;
			}
		AbstractAnalyzer acceptGenotype(final Predicate<Genotype> f) {
			this._acceptGT = f;
			return this;
			}
		Predicate<VariantContext> getVariantPredicate() {
			return this._acceptVariant;
			}
		Predicate<Genotype> getGenotypePredicate() {
			return this._acceptGT;
			}
		
		public String getTitle() { return getName();}
		public String getSubTitle() { return getDescription();}
		public String getYLab() { return ylab;}
		public String getXLab() { return xlab;}
		public String getColor(String category) { return "lightgray";}
		
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
		public void finish(PrintWriter w) {
			if(cat2values.isEmpty()) return;
			w.println(device());
			final List<String> phenotypes = new ArrayList<>(cat2values.keySet());
			Collections.sort(phenotypes, (A,B)->Double.compare(cat2values.get(A).getAverageCount().getAsDouble() ,cat2values.get(B).getAverageCount().getAsDouble()));
			
			
			for(int i=0;i< phenotypes.size();i++) {
				w.println("## "+ phenotypes.get(i));
				w.print("data"+i+" <-c(");
				boolean first = true;
				final Counter<Integer> counter = cat2values.get(phenotypes.get(i));
				for(Integer k:counter.keySet()) {
					long count = counter.count(k);
					for(long u=0;u<count;u++) {
						if(!first) w.print(",");
						first=false;
						w.print(k);
						}
					}
				w.println(")");
				}
			w.print("boxplot(");
			
			for(int i=0;i< phenotypes.size();i++) {
				w.print("data"+i+",");
				}
			
			w.println(
				"main="+quote(getTitle())+
				",sub="+quote(getSubTitle())+
				",ylab="+quote(getYLab())+","+
				"las=2,"+
				"names=c("+ phenotypes.stream().map(S->quote(S)).collect(Collectors.joining(","))+"),"+
				"col= c(" + phenotypes.stream().map(S->quote(getColor(S))).collect(Collectors.joining(",")) + ")" +
				")");
			w.println("dev.off()");
			}
		}
	/***************************************************************************/
	private abstract class AbstractBarPlot extends AbstractAnalyzer {
		public String getColor(String h) { return "lightgray";}
		}
	/***************************************************************************/
	private abstract class AbstractAggregateBarPlot extends AbstractBarPlot {
		protected final Map<String,DataPoint> key2average = new HashMap<>();
		private boolean minY_is_zero = true;
		protected void add(final String key,double v) {
			DataPoint t = key2average.get(key);
			if(t==null) {
				t= createDataPoint();
				key2average.put(key, t);
				}
			t.accept(v);
			}
		
		protected abstract DataPoint createDataPoint();
		
		protected double[] getYLim() {
			final double m;
			if(minY_is_zero) {
				m = 0.0;
				} else {
				m = key2average.values().stream().filter(A->A.isPresent()).
						mapToDouble(A->A.getAsDouble()).
						min().
						orElse(0.0)
						;
				}
			
			final double M= key2average.values().stream().filter(A->A.isPresent()).
					mapToDouble(A->A.getAsDouble()).
					max().
					orElse(1.0)
					;
			return new double[] {m,M};
			}
		
		AbstractAggregateBarPlot setMinYIsZero(boolean b) {
			this.minY_is_zero = b;
			return this;
			}
		
		protected Set<String> getSortedKeys() {
			final Set<String> keys= key2average.keySet().stream().
					sorted((A,B)->key2average.get(A).compareTo(key2average.get(B))).
					collect(Collectors.toCollection(LinkedHashSet::new));
			return keys;
			}
		@Override
		public void finish(PrintWriter w) {
			if(key2average.isEmpty()) return;
			if(key2average.values().stream().allMatch(A->!A.isPresent())) return;
			final Set<String> keys = getSortedKeys();
				
			w.println(device(key2average.size()));
			final double[] ylim = getYLim();
			w.println(
					"plot(c("+
					keys.stream().map(K->String.valueOf(key2average.get(K).get().orElse(0.0))).collect(Collectors.joining(",")) +
					"),main="+quote(getTitle())+","+
					"type="+quote("l")+","+
					"sub="+quote(getSubTitle())+","+
	                "xlab="+quote(getXLab())+","+
	                "ylab="+quote(getYLab())+","+
	                "ylim=c("+ylim[0]+","+ylim[1]+"),"+
	                "xaxt=\"n\")"
	               ); 
			// plot average value
			final double avg = key2average.values().stream().filter(A->A.isPresent()).mapToDouble(A->A.getAsDouble()).average().orElse(0);
			w.println("abline(h="+avg+",col=\"blue\")");
			
			w.println("axis(1,1:"+keys.size()+",c("+keys.stream().map(S->quote(S)).collect(Collectors.joining(","))+"), cex.axis = .7, las=2)");
			
			w.println("dev.off()");
			}
		}
	/***************************************************************************/

	private class RangeBarPlot extends AbstractAnalyzer {
		private  class Range {
			final double lowerBound;
			final double upperBound;
			long count = 0L;
			Range(final double m,final double M) {
				this.lowerBound = m;
				this.upperBound = M;
				}
			boolean contains(final double f) {
				return this.lowerBound <= f && f < this.upperBound;
				}
			@Override
			public String toString() {
				return "("+lowerBound+"/"+upperBound+"(:"+count;
				}
			}
		private final int precision;
		private final String infoTag;
		private final DecimalFormat decimalFormat;
		private final List<Range> ranges = new ArrayList<>();
		private boolean logX=false;
		
		RangeBarPlot(final String infoTag,int precision) {
			this.precision=precision;
			this.infoTag = infoTag;
			String str="#.";
			int p = precision;
			while(p>1) {
				str+="#";
				p=p/10;
				}
			this.decimalFormat = new DecimalFormat(str);
			this.decimalFormat.setRoundingMode(java.math.RoundingMode.CEILING);
			xlab(infoTag);
			xlab("Count Variants");
			}

			
		private double round(final double v) {
			return Double.parseDouble(this.decimalFormat.format(v));
			}		
		@Override
		public void init(final VCFHeader h) {
			final VCFInfoHeaderLine info = h.getInfoHeaderLine(this.infoTag);
			this.enabled = info!=null &&
					(info.getType()==VCFHeaderLineType.Float ||info.getType()==VCFHeaderLineType.Integer )
					;
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			if(!ctx.hasAttribute(this.infoTag)) return;
			try {
				for(String s : ctx.getAttributeAsStringList(this.infoTag,".") ) {
					if(s.equals(".")) continue;
					final double v0 = Double.parseDouble(s);
					if(Double.isNaN(v0)) continue;
					if(Double.isInfinite(v0)) continue;
					final double v = round(v0); 
					int i=0;
					for(i=0;i< this.ranges.size();i++) {
						final Range r = ranges.get(i);
						if(r.contains(v)) {
							r.count++;
							break;
							}
						}
					if(i!=this.ranges.size()) continue;
					final double v1 = v;
					final double v2 = v1 + 1.0/this.precision;
					final Range r = new Range(v1,v2);
					if(!r.contains(v)) throw new IllegalStateException(""+r+" "+v);
					r.count = 1L;	
					this.ranges.add(r);
					Collections.sort(this.ranges,(A,B)->Double.compare(A.lowerBound, B.lowerBound));
					}
				}
			catch(final Throwable err) {
				LOG.warn(err);
				}
			}
		long getCount() {
			return ranges.stream().mapToLong(R->R.count).sum();
			}

		RangeBarPlot logX(boolean b) { this.logX=b; return this;}

		private OptionalDouble getPercentile(double f) {
			final long count = getCount();
			if(count==0L) return OptionalDouble.empty();
			long count2 = (long)(count*f);
			for(Range r: ranges) {
				for(long i=0;i< r.count;i++) {
					if(count2==0L) return OptionalDouble.of(r.lowerBound);
					count2--;
					}
				}
			throw new IllegalStateException();
			}
		
		OptionalDouble getLowPercentile() {
			return getPercentile(0.0025);
			}
		
		OptionalDouble getHighPercentile() {
			return getPercentile(1.0 - 0.0025);
			}

		
		@Override
		public void finish(PrintWriter out) {
			if(this.ranges.isEmpty()) return;
			final long count = getCount();

			out.println(device(ranges.size()));
			
			out.print("x <- c(");
			out.print(ranges.stream().map(R->String.valueOf((R.lowerBound+R.upperBound)/2.0)).collect(Collectors.joining(",")));
			out.println(")");
			
			out.print("y <- c(");
			out.print(ranges.stream().map(R->String.valueOf(R.count/(double)count)).collect(Collectors.joining(",")));
			out.println(")");
			
			double minx = 	this.ranges.stream().mapToDouble(R->R.lowerBound).min().getAsDouble();
			if(minx <=0 && this.logX) minx=1.0;

			out.println("plot(x,y,pch=19,type=\"b\",col=\"blue\"," +
					"log="+quote(this.logX?"x":"")+"," + 
					"xlim=c("+ minx +
						"," +  
						this.ranges.stream().mapToDouble(R->R.upperBound).max().getAsDouble() +
						"),"+
					"ylim=c(0,max(y)),"+
					"main="+ quote(getTitle()) +","+
					"sub="+quote(getSubTitle())+","+
	                "xlab="+quote((this.logX?"log. ":"")+getXLab())+","+
	                "ylab="+quote(getYLab())+","+
	                "las=2)");

			
			OptionalDouble limit = getLowPercentile();
			if(limit.isPresent()) {
				out.println("abline(v="+limit.getAsDouble()+", col=\"green\")");
				}
			limit = getHighPercentile();
			if(limit.isPresent()) {
				out.println("abline(v="+limit.getAsDouble()+", col=\"green\")");
				}

			out.println("dev.off()");
			}
	}

	/***************************************************************************/
	private abstract class AbstractSimpleBarPlot extends AbstractBarPlot {
		protected final Counter<String> counter = new Counter<>();
		protected void add(final String value) {
			counter.incr(value);
			}
		
		@Override
		public void finish(PrintWriter w) {
			if(counter.getCountCategories()==0) return;
			w.println(device(counter.getCountCategories()));
			w.println(
					"barplot(c("+
					counter.keySetDecreasing().stream().map(K->String.valueOf(counter.count(K))).collect(Collectors.joining(",")) +
					"),main="+quote(getTitle())+","+
					"sub="+quote(getSubTitle())+","+
	                "xlab="+quote(getXLab())+","+
	                "ylab="+quote(getYLab())+","+
	                "ylim=c(0,"+counter.getMaxCount().orElse(1L)+"),"+
	                "names=c("+counter.keySetDecreasing().stream().map(S->String.valueOf(S)).map(S->quote(S)).collect(Collectors.joining(","))+"),"+
					"col= c(" + 
					counter.keySetDecreasing().stream().map(K->quote(getColor(K))).collect(Collectors.joining(",")) +
					"),las=2)");
			w.println("dev.off()");
			}
		}

	/***************************************************************************/
	private abstract class AbstractMultipleBarPlot extends AbstractBarPlot {
		private final Map<String,Counter<String>> horiz2counts = new LinkedHashMap<>();
		private final Set<String> distinct_vertical = new LinkedHashSet<>();
		/** category/count > value */
		private BiFunction<String,Long,Double> normalizer = (A,L)->L.doubleValue();
		protected void add(final String h,final String v) {
			Counter<String> t = horiz2counts.get(h);
			if(t==null) {
				t=new Counter<>();
				horiz2counts.put(h, t);
				}
			t.incr(v);
			distinct_vertical.add(v);
			}
		protected boolean isBeside() {
			return false;
			}
		
		protected AbstractMultipleBarPlot normalizer(BiFunction<String,Long,Double> normalizer) {
			this.normalizer = normalizer;
			return this;
			}
		
		@Override
		public void finish(PrintWriter w) {
			if(horiz2counts.isEmpty()) return;
			final List<String> LH = new ArrayList<>(horiz2counts.keySet());
			final List<String> colNames= new ArrayList<>(LH.size());
			for(int i=0;i< LH.size();i++) {
				final String h = LH.get(i);
				final Counter<String> counter = horiz2counts.get(h);
				w.println("#  "+ h);
				final String colName = "var"+i;
				colNames.add(colName);
				w.print(colName+" <- c(");
				w.print(this.distinct_vertical.stream().
					map(V->normalizer.apply(V,counter.count(V)).doubleValue()).
					map(V->String.valueOf(V)).
					collect(Collectors.joining(",")));
				w.println(")");
				}
			w.println("T2 <- as.matrix(data.frame("+String.join(",", colNames) +"))");
			w.println(device(LH.size()));
			w.println("par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)");
			w.println(
				"barplot(T2,main="+quote(getTitle())+","+
				"sub="+quote(getSubTitle())+","+
                "xlab="+quote(getXLab())+","+
                "ylab="+quote(getYLab())+","+
                "col=rainbow("+distinct_vertical.size()+")," +
                "names=c("+
                	horiz2counts.keySet().stream().
                	map(S->quote(String.valueOf(S))).
                	collect(Collectors.joining(","))+
                	"),"+
			"beside="+(isBeside()?"TRUE":"FALSE")+",las=2)");
			if(distinct_vertical.size()>1) {
				w.println("legend("+quote("topright")+", inset=c(-0.2,0), legend = c("+ distinct_vertical.stream().map(S->quote(S.toString())).collect(Collectors.joining(",")) +"),fill=rainbow("+distinct_vertical.size()+"))");
				}
			w.println("dev.off()");
			}
		}
	/***************************************************************************/
	private abstract class AbstractManhattanPlot extends AbstractAnalyzer {
		
		private Predicate<SAMSequenceRecord> acceptContig ;
		private SAMSequenceDictionary dict;
		private long genomeLength;
		private final Map<String,List<DataPoint>> cat2index = new HashMap<>();
		private final int win_width;
		private final int win_height;
		AbstractManhattanPlot() {
			this.acceptContig = SSR -> SSR.getContig().matches(VcfStats.this.att.getAttribute("contig.regex","(chr)?[0-9XY][0-9]?"));
			this.win_width = VcfStats.this.att.getIntAttribute("manhattan.width").orElse(1000);
			this.win_height = VcfStats.this.att.getIntAttribute("manhattan.height").orElse(300);
			}
		
		@Override
		public void init(final VCFHeader h) {
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
		public void finish(PrintWriter pw) {
			if(this.cat2index.isEmpty()) return;
			final double alpha = Math.max(0.1,1.0/this.cat2index.size());
			pw.println("palette <- rainbow("+this.cat2index.size()+",alpha="+alpha+")");
			pw.println("pdf("+quote(getOuputFilename()+".pdf")+",width=40,height=10)");
			
			// write manhattan figure
			pw.println("plot(1,type='n'," +
					"main="+quote(getTitle())+","+
					"sub="+quote(getSubTitle())+","+
					"ylab="+quote(getYLab())+","+
					"xlim = c(0 ,"+this.win_width+"), ylim = c(0,"+this.win_height+"), xaxt='n' )");
			for(int tid=0; tid < this.dict.size(); tid++) {
				final SAMSequenceRecord ssr = this.dict.getSequence(tid);
				pw.println("rect("+ 
						pos2pixel(ssr.getContig(), 1)+",0,"+ 
						pos2pixel(ssr.getContig(), ssr.getLengthOnReference())+","+
						this.win_height+",border=\"NA\",col= rgb("+(tid%2==0?"0,1.0,0.0":"1.0,1,1")+",alpha=0.05))");
				}
			// draw vertical lines
			for(int tid=0; tid +1 < this.dict.size(); tid++) {
				final SAMSequenceRecord ssr = this.dict.getSequence(tid);
				pw.println("abline(v="+ pos2pixel(ssr.getContig(), ssr.getLengthOnReference())+",col=\"gray\")");
				}
			// draw chromosomes names
			for(int tid=0; tid < this.dict.size(); tid++) {
				final SAMSequenceRecord ssr = this.dict.getSequence(tid);
				pw.println("text(x="+pos2pixel(ssr.getContig(), ssr.getLengthOnReference()/2)+",y=-10,labels="+quote(ssr.getSequenceName())+",col=\"black\",cex=0.5)");
				}
			final double maxy = this.cat2index.values().stream().flatMap(L->L.stream()).
					filter(X->X!=null && X.get().isPresent()).
					mapToDouble(X->X.get().getAsDouble()).max().orElse(1.0);
			// plot each cat
			int cat_idx=0;
			for(final String cat: this.cat2index.keySet()) {
				final List<DataPoint> L  = this.cat2index.get(cat);
				final int count= (int)(L.stream().filter(X->X!=null && X.get().isPresent()).count());
				if(count==0) continue;
				final List<Integer> array_x = new ArrayList<>(count);
				final List<Double> array_y = new ArrayList<>(count);
				for(int x=0;x< L.size();++x) {
					final DataPoint pt = L.get(x);
					if(pt==null || !pt.get().isPresent()) continue;
					array_x.add(x);
					array_y.add((pt.get().getAsDouble()/maxy)*win_height);
					}
				pw.println("##" + cat);
				pw.print("lines(x=c(");
				pw.print(array_x.stream().map(X->String.valueOf(X)).collect(Collectors.joining(",")));
				pw.print("),y=c(");
				pw.print(array_y.stream().map(Y->String.valueOf(Y)).collect(Collectors.joining(",")));
				pw.println("),type=\"p\",col=palette["+ (cat_idx+1) +"])");
				cat_idx++;
				}
			//y axis
			//pw.println("axis(2,c(0,"+win_height+"),c(0,"+maxy+"), cex.axis = .7)");

			// Add a legend
			pw.print("legend("+ (win_width+10) +","+maxy+", legend=c("); 
			pw.print( this.cat2index.keySet().stream().map(S->quote(S)).collect(Collectors.joining(",")));
			pw.print("),col=c(");
			cat_idx=0;
			for(final String cat: this.cat2index.keySet()) {
				if(cat_idx>0) pw.print(",");
				pw.print("palette["+(cat_idx+1)+"]");
				cat_idx++;
				}
			
			pw.println("),lty=1:2, cex=0.8)");
			
			pw.println("dev.off()");
			}
		}
	/***************************************************************************/
	private class Singleton extends AbstractSampleToFraction {
		@Override
		public void init(final VCFHeader h) {
			super.init(h);
			if(this.enabled && h.getNGenotypeSamples()<=1) {
				this.enabled=false;
				}
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			super.n_variants++;
			Genotype single = null;
			for(Genotype g:ctx.getGenotypes()) {
				if(g.isNoCall() || g.isHomRef()) continue;
				if(!getGenotypePredicate().test(g)) continue;
				if(single!=null) return;
				single=g;
				}
			if(single!=null) add(single.getSampleName(),1.0);
			}
		}
	
	
	
	/***************************************************************************/
		
	private class NormalizedContigs extends AbstractMultipleBarPlot {
		private SAMSequenceDictionary dict;
		private final Predicate<SAMSequenceRecord> predicate;
		NormalizedContigs(final Predicate<SAMSequenceRecord> predicate) {
			this.predicate=predicate;
			xlab("Samples");
			ylab("Count Variants/length chrom");
			normalizer((CHROM,COUNT)-> COUNT/(double)dict.getSequence(CHROM).getLengthOnReference());
			}
		@Override
		public void init(final VCFHeader h) {
			final SAMSequenceDictionary dict0 = h.getSequenceDictionary();
			if(dict0!=null) {
				this.dict = new SAMSequenceDictionary(
					dict0.getSequences().stream().
					filter(this.predicate).
					collect(Collectors.toList())
					);
				}
			this.enabled = dict!=null && !dict.isEmpty();
			}
		
		@Override
		protected boolean isBeside() {
			return true;
			}
		@Override
		public void visit(VariantContext ctx) {
			final String contig = ctx.getContig();
			if(dict.getSequence(contig)==null) return;
			if(!getVariantPredicate().test(ctx)) return;
			for(Genotype g: ctx.getGenotypes()) {
				if(!getGenotypePredicate().test(g)) continue;
				if(g.getAlleles().stream().allMatch(A->A.isReference() || A.isNoCall())) continue;
				add(g.getSampleName(), contig);
				}
			}
		}
	
	/***************************************************************************/
	private class GatkDeNovo extends AbstractMultipleBarPlot {
		private final String[] confDeNovos = new String[]{
				GATKConstants.hiConfDeNovo,
				GATKConstants.loConfDeNovo
				};

		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.getInfoHeaderLine(this.confDeNovos[0])!=null && h.getInfoHeaderLine(this.confDeNovos[1])!=null;
			}
		
		@Override
		protected boolean isBeside() {
			return true;
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			for(String info: this.confDeNovos) {
				if(!ctx.hasAttribute(info)) continue;
				for(final String s: ctx.getAttributeAsStringList(info, "")) {
					if(StringUtils.isBlank(s)) continue;
					add(s, info+ (ctx.isFiltered()?".FILTER":".PASS"));
					}
				}
			}
		}
	/***************************************************************************/
	private class DragenDeNovo extends AbstractMultipleBarPlot {
		private final String DN = "DN";
		
		@Override
		public void init(final VCFHeader h) {
			VCFFormatHeaderLine hdr= h.getFormatHeaderLine(DN);
			this.enabled = hdr!=null && hdr.getType().equals(VCFHeaderLineType.String);
			}
		
		@Override
		protected boolean isBeside() {
			return true;
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			
			for(Genotype g:ctx.getGenotypes()) {
				if(!getGenotypePredicate().test(g)) continue;
				if(!g.hasExtendedAttribute(DN)) continue;
				final Object v = g.getExtendedAttribute(DN,"");
				if(v==null) continue;
				final String s= v.toString();
				if(StringUtils.isBlank(s) || !s.equals("DeNovo")) continue;
				add(g.getSampleName(), s + (g.isFiltered() || ctx.isFiltered()?".FILTERED":".PASS"));
				}
			}
		}

	/***************************************************************************/
	private class SVTypeContig extends AbstractMultipleBarPlot {
		SVTypeContig() {
			name("svtype2contig");
			description("SVTYpe per contig");
			}
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.getInfoHeaderLine(VCFConstants.SVTYPE)!=null;
			}
		@Override
		public void visit(final VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			final String svType = ctx.getAttributeAsString(VCFConstants.SVTYPE, null);
			if(StringUtils.isBlank(svType)) return;
			super.add(ctx.getContig(), svType);
			}
		}

	/***************************************************************************/
	private class Sample2GTType extends AbstractMultipleBarPlot {
		private final Function<Genotype,String> gt2str;
		Sample2GTType( Function<Genotype,String> gt2str) {
			this.gt2str=gt2str;
			}
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData();
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			for(final Genotype gt: ctx.getGenotypes()) {
				if(!getGenotypePredicate().test(gt)) continue;
				final String value = gt2str.apply(gt);
				if(StringUtils.isBlank(value)) continue;
				super.add(gt.getSampleName(), value);	
				}
			}
		}
	/***************************************************************************/
	private abstract class AbstractSampleToFraction extends AbstractAggregateBarPlot {
		protected long n_variants =0L;
		@Override
		protected DataPoint createDataPoint() {
			return new DataPointSum() {
				@Override
				public OptionalDouble get() {
					if(AbstractSampleToFraction.this.n_variants==0) return OptionalDouble.empty();
					return OptionalDouble.of(getCount()/(double)AbstractSampleToFraction.this.n_variants);
					}
				};
			}
		@Override
		public String getYLab() {
			return super.getYLab()+" (N="+StringUtils.niceInt(this.n_variants)+")";
			}
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData();
			for(final String sn:h.getGenotypeSamples()) {
				key2average.put(sn, createDataPoint());
				}
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			this.n_variants++;
			
			final Predicate<Genotype> filterGT = getGenotypePredicate();
			for(final Genotype gt: ctx.getGenotypes()) {
				if(!filterGT.test(gt)) continue;
				super.add(gt.getSampleName(), 1.0);
				}
			}
		}
	/***************************************************************************/
	private class SampleToGenotypeTypeFraction extends AbstractSampleToFraction {
		private final GenotypeType genotypeType;
		SampleToGenotypeTypeFraction(GenotypeType genotypeType) {
			this.genotypeType = genotypeType;
			name("fraction of "+genotypeType.name());
			ylab("Fraction of "+genotypeType.name()+" genotypes");
			xlab("Sample");
			}
		@Override
		Predicate<Genotype> getGenotypePredicate() {
			return super.getGenotypePredicate().and(GT->GT.getType().equals(this.genotypeType));
			}
		}
	/***************************************************************************/
	private class SampleMultiAllelicFraction extends AbstractSampleToFraction {
		SampleMultiAllelicFraction() {
			name("fraction of multiallelic");
			ylab("fraction of multiallelic");
			xlab("Sample");
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			this.n_variants++;
			if(ctx.getNAlleles()<=2) return;
			
			final Predicate<Genotype> filterGT = getGenotypePredicate();
			for(final Genotype gt: ctx.getGenotypes()) {
				if(!filterGT.test(gt)) continue;
				if(gt.getAlleles().stream().allMatch(A->A.isReference() ||A.isNoCall())) continue;
				super.add(gt.getSampleName(), 1.0);
				}
			}
		}
	/***************************************************************************/
	private class Sample2GTFilters extends AbstractMultipleBarPlot {
		Sample2GTFilters() {
			name("sample2filter");
			description("Genotype FILTERs per Sample");
			xlab("Sample");
			ylab("Filter");
			}
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData();
			if(this.enabled && h.getFormatHeaderLine(VCFConstants.GENOTYPE_FILTER_KEY)==null) {
				this.enabled = false;
				}
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			for(final Genotype gt: ctx.getGenotypes()) {
				if(!gt.isFiltered()) continue;
				if(!getGenotypePredicate().test(gt)) continue;
				super.add(gt.getSampleName(), gt.getFilters());	
				}
			}
		
		}
	
	/*********************************************************************/
	private class SVLen extends AbstractBoxPlot {
		SVLen() {
			name("SVLEN");
			description("SVLEN");
			ylab("median(SVLEN)");
			xlab("SVTYPE = median(SVLEN)");
			}
		@Override
		public void init(VCFHeader h) {
			this.enabled = h.getInfoHeaderLine(VCFConstants.SVTYPE)!=null;
			}
		
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
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
			super.add(st, svLen);
			}
		}
	
	
	/*********************************************************************/
	private class SampleToSVTypes extends AbstractMultipleBarPlot {
		SampleToSVTypes() {
			xlab("Samples");
			}
		@Override
		public void init(VCFHeader h) {
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
	private class SampleToSVTypeLen extends AbstractBoxPlot {
		private final String svType;
		SampleToSVTypeLen(final String svType) {
			this.svType = svType;
			name("SVLEN."+svType);
			ylab("SVLEN."+svType);
			}
		@Override
		public void init(VCFHeader h) {
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
	private abstract class AbstractSampleToAverage extends AbstractAggregateBarPlot {
		protected AbstractSampleToAverage() {
			xlab("sample");
			ylab("avg("+getFormatKey()+")");
			}
		@Override
		protected DataPoint createDataPoint() {
			return new DataPointAverage();
			}
		
		protected abstract String getFormatKey();
		@Override
		public void init(VCFHeader h) {
			this.enabled = h.hasGenotypingData() && h.getFormatHeaderLine(getFormatKey())!=null;
			}
				
		protected abstract void visit(Genotype gt);
		
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			for(Genotype gt:ctx.getGenotypes()) {
				if(!getGenotypePredicate().test(gt)) continue;
				visit(gt);
				}
			}
		}
	/*********************************************************************/
	private abstract class AbstractSampleAvgBoxPlot extends AbstractAnalyzer {
		protected final Map<String,DataPoint> sample2datapoint = new HashMap<>();
		protected AbstractSampleAvgBoxPlot() {
			xlab("phenotype");
			}
	
		protected abstract DataPoint createDataPoint();
		protected abstract OptionalDouble getValueForGenotype(final Genotype gt);
		
		@Override
		public void init(VCFHeader h) {
			this.enabled = VcfStats.this.sample2phenotype!=null &&
					!VcfStats.this.sample2phenotype.isEmpty() &&
					h.hasGenotypingData();
			}
				
		@Override
		public void finish(PrintWriter w) {
			if(sample2datapoint.isEmpty()) return;
			
			final List<String> phenotypes = this.sample2datapoint.keySet().stream().
					map(S->VcfStats.this.sample2phenotype.get(S)).
					collect(Collectors.toCollection(TreeSet::new)).
					stream().
					collect(Collectors.toList());
			if(phenotypes.isEmpty()) return;
			w.println(device());
			
			w.println("palette <- rainbow("+phenotypes.size()+")");
			for(int i=0;i< phenotypes.size();i++) {
				final String phenotype = phenotypes.get(i);
				w.println("## "+ phenotype);
				w.print("data"+i+" <-c(");
				w.print(sample2datapoint.entrySet().stream().
						filter(KV->VcfStats.this.sample2phenotype.get(KV.getKey()).equals(phenotype)).
						map(KV->KV.getValue()).
						filter(KV->KV.isPresent()).
						map(KV->String.valueOf(KV.getAsDouble())).
						collect(Collectors.joining(","))
						);
				w.println(")");
				}
			w.print("boxplot(");
			
			for(int i=0;i< phenotypes.size();i++) {
				w.print("data"+i+",");
				}
			
			w.println(
				"main="+quote(getTitle())+
				",sub="+quote(getSubTitle())+
				",ylab="+quote(getYLab())+","+
				"las=2,"+
				"names=c("+ phenotypes.stream().map(S->quote(S)).collect(Collectors.joining(","))+"),"+
				"col= c(" + IntStream.range(0, phenotypes.size()).mapToObj(S->"palette["+(S+1)+"]").collect(Collectors.joining(",")) + ")" +
				")");
			w.println("dev.off()");
			}

		@Override
		public void visit(final VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
			for(Genotype gt:ctx.getGenotypes()) {
				if(!getGenotypePredicate().test(gt)) continue;
				if(!sample2phenotype.containsKey(gt.getSampleName())) continue;
				final OptionalDouble value =getValueForGenotype(gt);
				if(!value.isPresent()) continue;
				DataPoint dpt = sample2datapoint.get(gt.getSampleName());
				if(dpt==null) {
					dpt = createDataPoint();
					sample2datapoint.put(gt.getSampleName(),dpt);
					}
				dpt.accept(value.getAsDouble());
				}
			}
		}

	/*********************************************************************/
	private class SampleToAvgBoxPlot extends AbstractSampleAvgBoxPlot {
		private String formatKey;
		protected SampleToAvgBoxPlot(final String formatKey) {
			this.formatKey= formatKey; 
			ylab("avg("+getFormatKey()+")");
			}
		@Override
		protected DataPoint createDataPoint() {
			return new DataPointAverage(); 
			}
		
		protected String getFormatKey() {
			return formatKey;
			}
		
		@Override
		protected OptionalDouble getValueForGenotype(final Genotype gt) {
			if(getFormatKey().equals(VCFConstants.DEPTH_KEY)) {
				return gt.hasDP()? OptionalDouble.of(gt.getDP()):OptionalDouble.empty();
				}
			else if(getFormatKey().equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
				return gt.hasGQ()? OptionalDouble.of(gt.getGQ()):OptionalDouble.empty();
				}
			else if(getFormatKey().equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS)) {
				if(!gt.hasAD()) return OptionalDouble.empty();
				final int[] ad = gt.getAD();
				if(ad.length!=2) return  OptionalDouble.empty();
				final int t = ad[0]+ad[1];
				if(t==0) return  OptionalDouble.empty();
				return OptionalDouble.of(ad[1]/(double)t);
				}
			else
				{
				if(!gt.hasAnyAttribute(getFormatKey())) return OptionalDouble.empty();
				final Object o = gt.getAnyAttribute(getFormatKey());
				if(o==null || !(o instanceof Number)) return OptionalDouble.empty();
				return OptionalDouble.of(Number.class.cast(o).doubleValue());
				}
			}
		
		@Override
		public void init(VCFHeader h) {
			super.init(h);
			if(this.enabled) {
				this.enabled =  h.hasGenotypingData() && h.getFormatHeaderLine(getFormatKey())!=null;
				}
			}
		}

	/*********************************************************************/
	private class SampleToGQ extends AbstractSampleToAverage {
		@Override
		public String getFormatKey() {
			return VCFConstants.GENOTYPE_QUALITY_KEY;
			}
		@Override
		protected void visit(final Genotype gt) {
			if(!gt.hasGQ()) return;
			super.add(gt.getSampleName(), gt.getGQ());
			}
		}
	/*********************************************************************/
	private class SampleToDP extends AbstractSampleToAverage {
		@Override
		public String getFormatKey() {
			return VCFConstants.DEPTH_KEY;
			}
		@Override
		protected void visit(final Genotype gt) {
			if(!gt.hasDP()) return;
			super.add(gt.getSampleName(), gt.getDP());
			}
		}
	/*********************************************************************/
	private class SampleToAD extends AbstractSampleToAverage {
		@Override
		public String getFormatKey() {
			return VCFConstants.GENOTYPE_ALLELE_DEPTHS;
			}
		@Override
		protected void visit(final Genotype gt) {
			if(!gt.hasAD()) return;
			final int[] ad = gt.getAD();
			if(ad.length!=2) return;
			final double f =  ad[0]+ad[1]==0 ? 0 : (ad[1]/(double)( ad[0]+ad[1]));
			super.add(gt.getSampleName(), f);
			}
		}
	/*********************************************************************/
	private class CountVariantsManhattanPlot extends AbstractManhattanPlot {
		@Override
		public void visit(VariantContext ctx) {
			if(!getVariantPredicate().test(ctx)) return;
				visit(ctx.getType().name(),ctx,  1);
			}
		}
	/*********************************************************************/
	private class CountSamplesVariantsManhattanPlot extends AbstractManhattanPlot {
		@Override
		public void init(VCFHeader h) {
			super.init(h);
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

	private String quote(final String s) {
		return "\""+StringUtils.escapeC(s)+"\"";
		}
	


	private void loadPhenotypes(VCFHeader header) throws IOException {
		if(!header.hasGenotypingData()) return;
		final Set<String> sampleset = new HashSet<>(header.getGenotypeSamples());
		if(sample2catPath!=null) {
			try(BufferedReader br= IOUtils.openPathForBufferedReading(sample2catPath)) {
				String line;
				while((line=br.readLine())!=null) {
					
					if(line.startsWith("#") || StringUtils.isBlank(line)) {
						continue;
						}
					final String[] tokens = CharSplitter.TAB.split(line);
					if(tokens.length!=2) throw new JvarkitException.TokenErrors(2, tokens);
					if(!sampleset.contains(tokens[0])) {
						continue;
						}
					if(sample2phenotype.containsKey(tokens[0])) {
						throw new IOException("multiple phenotypes for "+tokens[0]+" in "+sample2catPath);
						}
					sample2phenotype.put(tokens[0], tokens[1]);
					List<String> samples = phenotype2samples.get(tokens[1]);
					if(samples==null) {
						samples = new ArrayList<>();
						phenotype2samples.put(tokens[1],samples);
						}
					samples.add(tokens[0]);
				}
			}
		final String other="other";
		for(String sn:sampleset) {
			if(sample2phenotype.containsKey(sn)) continue;
			sample2phenotype.put(sn, other);
			List<String> samples = phenotype2samples.get(other);
			if(samples==null) {
				samples = new ArrayList<>();
				phenotype2samples.put(other,samples);
				}
			samples.add(other);
			}
		}
	}
	
	
	@Override
	public int doWork(final List<String> args) {
		final List<Analyzer> modules =new ArrayList<>();
		for(int side=0;side<2;++side) {
			final String fltstr = side==0?"":" skip FILTEred genotypes";
			final Predicate<Genotype> acceptGT;
			if(side==0) {
				acceptGT = A->true;
				}
			else
				{
				acceptGT = A->!A.isFiltered();
				}
			modules.add(new Singleton().
					acceptGenotype(acceptGT).
					name("Singleton"+fltstr).
					description("Singleton per Sample"+fltstr).
					xlab("Sample").
					ylab("Fraction(Singletons)")
					);
			modules.add(new Sample2GTType(side==0?GT->GT.getType().name():(GT->GT.isFiltered()?VCFConstants.GENOTYPE_FILTER_KEY:GT.getType().name())).
					acceptGenotype(acceptGT).
					name("Genotype-Type / Sample"+fltstr).
					description("Genotype type per Sample"+fltstr).
					xlab("Sample").
					ylab("GT Types")
					);
			modules.add(new SampleToGQ().
					acceptGenotype(acceptGT).
					name("Sample GQ "+fltstr).
					description("Genotype Quality per Sample"+fltstr)
					);
			modules.add(new SampleToDP().
					acceptGenotype(acceptGT).
					name("Sample DP "+fltstr).
					description("Genotype DEPTH per Sample"+fltstr)
					);
			for(final GenotypeType genotype_type: new GenotypeType[] {GenotypeType.HET,GenotypeType.HOM_REF,GenotypeType.HOM_VAR}) {
					modules.add(new SampleToAD().
						setMinYIsZero(false).
						acceptGenotype(acceptGT.and(G->G.getType().equals(genotype_type))).
						name("AD per "+genotype_type.name()+" Sample"+fltstr).
						description("AD per "+genotype_type.name()+" Sample"+fltstr)
						);
					}
			

			}
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
		modules.add(
			new SVLen().
				name("SVLEN per SVTYPE").
				description("SVLEN per SVTYPE").
				xlab("SVTYPE").
				ylab("SVLEN")
				
			);
		modules.add(new SampleToSVTypes().
				name("Sample To SVTYPE").
				description("Sample To SVTYPE").
				xlab("Samples").
				ylab("SVTYPE")
			);
		modules.add(new Sample2GTFilters().
				name("sample Genotype FILTERs").
				description("FILTERs found in Samples").
				xlab("Samples").
				ylab("count(Genotype FILTERs)")
				);
		
		
		modules.add(new RangeBarPlot(GATKConstants.QD_KEY, 1).
				name(GATKConstants.QD_KEY).
				description("Variant Confidence (QUAL) / Quality by Depth.").
				xlab(GATKConstants.QD_KEY).
				ylab("Count Variants")
				);
		
		
		modules.add(new RangeBarPlot(GATKConstants.FS_KEY, 1).
			logX(true).
			name(GATKConstants.FS_KEY).
			description("Phred-scaled p-value using Fisher's exact test to detect strand bias").
			xlab(GATKConstants.FS_KEY).
			ylab("Count Variants")
			);
		
		modules.add(new RangeBarPlot(GATKConstants.SOR_KEY, 10).
			name(GATKConstants.SOR_KEY).
			description("Symmetric Odds Ratio of 2x2 contingency table to detect strand bias").
			xlab(GATKConstants.SOR_KEY).
			ylab("Count Variants")
			);
		
		modules.add(new RangeBarPlot(GATKConstants.MQ_KEY, 1).
			name(GATKConstants.MQ_KEY).
			description("Mean square mapping quality over all the reads at the site").
			xlab(GATKConstants.MQ_KEY).
			ylab("Count Variants")
			);

		modules.add(new RangeBarPlot(GATKConstants.MQRankSum_KEY, 10).
			name(GATKConstants.MQRankSum_KEY).
			description("Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities").
			xlab(GATKConstants.MQRankSum_KEY).
			ylab("Count Variants")
			);
		
		modules.add(new RangeBarPlot(GATKConstants.ReadPosRankSum_KEY, 10).
			name(GATKConstants.ReadPosRankSum_KEY).
			description("Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias").
			xlab(GATKConstants.ReadPosRankSum_KEY).
			ylab("Count Variants")
			);

		
		modules.add(new GatkDeNovo().
				name("GATK DeNovo").
				description("De Novo variants found with GATK").
				xlab("Sample").
				ylab("Count Variants")
				);
		modules.add(new DragenDeNovo().
				name("Dragen DeNovo").
				description("De Novo variants found with Dragen").
				xlab("Sample").
				ylab("Count Variants")
				);
		modules.add(new SampleToAvgBoxPlot(VCFConstants.DEPTH_KEY));
		modules.add(new SampleToAvgBoxPlot(VCFConstants.GENOTYPE_QUALITY_KEY));
		modules.add(new SampleToAvgBoxPlot(VCFConstants.GENOTYPE_ALLELE_DEPTHS));
		
		modules.add(new RangeBarPlot(VCFConstants.DEPTH_KEY, 1).
				logX(true).
				name(VCFConstants.DEPTH_KEY).
				description("DEPTH per variant").
				xlab(VCFConstants.DEPTH_KEY).
				ylab("Count Variants")
				);

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
		
		modules.add(
				new SampleToGenotypeTypeFraction(GenotypeType.NO_CALL).
					setMinYIsZero(false)
				);
		modules.add(
				new NormalizedContigs(SSR->SSR.getContig().matches("(chr)?[XY]")).
					name("xy per sample").
					description("Number of variants per X/Y per Samples")
				);
		modules.add(
				new SampleMultiAllelicFraction().
					setMinYIsZero(false)
				);
		
		
		final String input = oneFileOrNull(args);
		
		// remove modules
		modules.removeIf(M->Arrays.stream(this.moduleExcludeStr.split("[,; \t:]")).anyMatch(S->S.equalsIgnoreCase(M.getName())));
		try {
			if(this.list_modules) {
				try(PrintWriter w = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					for(Analyzer analyzer:modules) {
						w.print(analyzer.getName());
						w.print("\t");
						w.println(analyzer.getDescription());
						}	
					w.flush();
					}
				return 0;
				}
		
		try(VCFIterator iter= super.openVCFIterator(input)) {	
			final VCFHeader header=iter.getHeader();
			loadPhenotypes(header);
			
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
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
			while(iter.hasNext()) {
				final VariantContext ctx = progress.watch(iter.next());
				for(Analyzer analyzer:modules) {
					analyzer.visit(ctx);
					}
				}
			progress.finish();
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(outputFile)) {
				for(Analyzer analyzer:modules) {
					analyzer.finish(pw);
					}
				pw.flush();
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
