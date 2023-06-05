/*
The MIT License (MIT)
Copyright (c) 2023 Pierre Lindenbaum

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
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.Average;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
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
	modificationDate = "20230605",
	creationDate = "20131212",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfStats extends Launcher {
	private static final Logger LOG = Logger.build(VcfStats.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
		
	@Parameter(names={"--categories","--phenotypes"},description="A tab delimited file (sample)<tab>(category)",hidden = true)
	private Path sample2catPath = null;

	
	@Parameter(names={"-exclude","--exclude"},description="name of modules to be excluded")
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
		protected Predicate<VariantContext> acceptVariant = V->true;
		protected Predicate<Genotype> acceptGT = V->true;

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
			this.acceptVariant= f;
			return this;
			}
		AbstractAnalyzer acceptGenotype(final Predicate<Genotype> f) {
			this.acceptGT = f;
			return this;
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
	private abstract class AbstractAverageBarPlot extends AbstractBarPlot {
		protected final Map<String,Average> key2average = new HashMap<>();
		protected void add(final String key,double v) {
			Average t = key2average.get(key);
			if(t==null) {
				t=new Average();
				key2average.put(key, t);
				}
			t.accept(v);
			}
		
		protected double[] getYLim() {
			return new double[] {0.0,key2average.values().stream().mapToDouble(A->A.getAsDouble()).max().orElse(1L)};
		}
		protected Set<String> getSortedKeys() {
			final Set<String> keys= key2average.keySet().stream().
					sorted((A,B)->Double.compare(key2average.get(A).getAsDouble(), key2average.get(B).getAsDouble())).
					collect(Collectors.toCollection(LinkedHashSet::new));
			return keys;
			}
		@Override
		public void finish(PrintWriter w) {
			if(key2average.isEmpty()) return;
			final Set<String> keys = getSortedKeys();
				
			w.println(device(key2average.size()));
			final double[] ylim = getYLim();
			w.println(
					"plot(c("+
					keys.stream().map(K->String.valueOf(key2average.get(K).getAsDouble())).collect(Collectors.joining(",")) +
					"),main="+quote(getTitle())+","+
					"type="+quote("l")+","+
					"sub="+quote(getSubTitle())+","+
	                "xlab="+quote(getXLab())+","+
	                "ylab="+quote(getYLab())+","+
	                "ylim=c("+ylim[0]+","+ylim[1]+"),"+
	                "xaxt=\"n\")"
	               ); 
			// plot average value
			final double avg = key2average.values().stream().filter(A->A.getCount()>0).mapToDouble(A->A.getAsDouble()).average().orElse(0);
			w.println("abline(h="+avg+",col=\"blue\")");
			
			w.println("axis(1,1:"+keys.size()+",c("+keys.stream().map(S->quote(S)).collect(Collectors.joining(","))+"), cex.axis = .7, las=2)");
			
			w.println("dev.off()");
			}
		}
	/***************************************************************************/

	private class RangeBarPlot extends AbstractAverageBarPlot {
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
		public void init(VCFHeader h) {
			VCFInfoHeaderLine info = h.getInfoHeaderLine(this.infoTag);
			this.enabled = info!=null &&
					(info.getType()==VCFHeaderLineType.Float ||info.getType()==VCFHeaderLineType.Integer )
					;
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!super.acceptVariant.test(ctx)) return;
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
					map(V->counter.count(V)).
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
		private abstract class DataPoint {
			public abstract double visit(double v);
			public abstract double get();
			}
		private Predicate<SAMSequenceRecord> acceptContig = SSR -> SSR.getContig().matches("(chr)?[0-9XY][0-9]?");
		private SAMSequenceDictionary dict;
		private long genomeLength;
		private Map<String,List<DataPoint>> cat2index = new HashMap<>();
		private final int win_size=1000;
		@Override
		public void init(VCFHeader h) {
			final SAMSequenceDictionary dict0 = h.getSequenceDictionary();
			if(dict0!=null) {
				this.dict =  new SAMSequenceDictionary(dict0.getSequences().stream().filter(this.acceptContig).collect(Collectors.toList()));
				this.genomeLength = this.dict.getReferenceLength();
				}
			this.enabled = dict!=null && !dict.isEmpty();
			}
		private long pos2index(String contig,int pos) {
			long n=0;
			for(SAMSequenceRecord ssr:dict.getSequences()) {
				if(ssr.getContig().equals(contig)) {
					return n + pos;
					}
				n+=ssr.getLengthOnReference();
				}
			return -1;
			}
		protected abstract DataPoint createDataPoint();
		
		protected void visit(String category,String contig,int pos,double value) {
			final SAMSequenceRecord ssr= this.dict.getSequence(contig);
			if(ssr==null || pos<1 || pos > ssr.getLengthOnReference()) return;
			List<DataPoint> datapoints;
			if(this.cat2index.containsKey(category)) {
				datapoints = this.cat2index.get(category);
				}
			else
				{
				datapoints = new ArrayList<>();
				this.cat2index.put(category,datapoints);
				}
			long genomic_index= pos2index(contig,pos);
			int pixl = (int)((genomic_index/(double)this.genomeLength)*this.win_size);
			while(datapoints.size() <= pixl) {
				datapoints.add(null);
				}
			DataPoint dp  = datapoints.get(pixl);
			if(dp==null) {
				dp = createDataPoint();
				datapoints.set(pixl,dp);
				}
			dp.visit(value);
			}
		}
	/***************************************************************************/
	private class Singleton extends AbstractSimpleBarPlot {
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData() && h.getNGenotypeSamples()>1;
			}
		@Override
		public void visit(VariantContext ctx) {
			if(!super.acceptVariant.test(ctx)) return;
			Genotype single = null;
			for(Genotype g:ctx.getGenotypes()) {
				if(g.isNoCall() || g.isHomRef()) continue;
				if(!acceptGT.test(g)) continue;
				if(single!=null) return;
				single=g;
				}
			if(single!=null) add(single.getSampleName());
			}
		}
	/***************************************************************************/

	private class GatkDeNovo extends AbstractMultipleBarPlot {
		private final String[] confDeNovos = new String[]{"hiConfDeNovo", "loConfDeNovo"};

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
			if(!super.acceptVariant.test(ctx)) return;
			for(String info: this.confDeNovos) {
				if(!ctx.hasAttribute(info)) continue;
				for(final String s: ctx.getAttributeAsStringList(info, "")) {
					if(StringUtils.isBlank(s)) continue;
					add(s, info);
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
			if(!super.acceptVariant.test(ctx)) return;
			for(Genotype g:ctx.getGenotypes()) {
				if(!acceptGT.test(g)) continue;
				if(!g.hasExtendedAttribute(DN)) continue;
				final Object v = g.getExtendedAttribute(DN,"");
				if(v==null) continue;
				final String s= v.toString();
				if(StringUtils.isBlank(s) || s.equals("Inherited")) continue;
				add(g.getSampleName(), s);
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
			if(!super.acceptVariant.test(ctx)) return;
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
			if(!super.acceptVariant.test(ctx)) return;
			for(final Genotype gt: ctx.getGenotypes()) {
				if(!super.acceptGT.test(gt)) continue;
				final String value = gt2str.apply(gt);
				if(StringUtils.isBlank(value)) continue;
				super.add(gt.getSampleName(), value);	
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
			if(!super.acceptVariant.test(ctx)) return;
			for(final Genotype gt: ctx.getGenotypes()) {
				if(!gt.isFiltered()) continue;
				if(!super.acceptGT.test(gt)) continue;
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
			if(!super.acceptVariant.test(ctx)) return;
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
			if(!super.acceptVariant.test(ctx)) return;
			final String st = ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
			if(StringUtils.isBlank(st)) return;
			for(final Genotype gt:ctx.getGenotypes()) {
				if(!super.acceptGT.test(gt)) continue;
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
			if(!super.acceptVariant.test(ctx)) return;
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
				if(!super.acceptGT.test(g)) continue;
				if(g.getAlleles().stream().allMatch(G->G.isReference() || G.isNoCall())) continue;
				add(g.getSampleName(), svLen);
				}
			}
		}
	
	/*********************************************************************/
	private abstract class AbstractSampleToAverage extends AbstractAverageBarPlot {
		protected AbstractSampleToAverage() {
			xlab("sample");
			ylab("avg("+getFormatKey()+")");
			}
		protected abstract String getFormatKey();
		@Override
		public void init(VCFHeader h) {
			this.enabled = h.hasGenotypingData() && h.getFormatHeaderLine(getFormatKey())!=null;
			}
				
		protected abstract void visit(Genotype gt);
		
		@Override
		public void visit(VariantContext ctx) {
			if(!super.acceptVariant.test(ctx)) return;
			for(Genotype gt:ctx.getGenotypes()) {
				if(!super.acceptGT.test(gt)) continue;
				visit(gt);
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
		protected void visit(Genotype gt) {
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
		protected void visit(Genotype gt) {
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
			double f =  ad[0]+ad[1]==0 ? 0 : (ad[1]/(double)( ad[0]+ad[1]));
			super.add(gt.getSampleName(), f);
			}
		}

	
	/*********************************************************************/

	private String quote(final String s) {
		return "\""+StringUtils.escapeC(s)+"\"";
		}
	private String sampleToColor(final String s) {
		return "yellow";
		}


	private int loadPhenotypes() throws IOException {
		if(sample2catPath!=null) {
			try(BufferedReader br= IOUtils.openPathForBufferedReading(sample2catPath)) {
				String line;
				while((line=br.readLine())!=null) {
					
					if(line.startsWith("#") || StringUtils.isBlank(line)) {
						continue;
						}
					final String[] tokens = CharSplitter.TAB.split(line);
					if(tokens.length!=2) throw new JvarkitException.TokenErrors(2, tokens);
					if(sample2phenotype.containsKey(tokens[0])) {
						LOG.error("multiple phenotypes for "+tokens[0]+" in "+sample2catPath);
						return -1;
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
		}
		return 0;
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
					ylab("Count(Singletons)")
					);
			modules.add(new Sample2GTType(side==0?GT->GT.getType().name():(GT->GT.isFiltered()?VCFConstants.GENOTYPE_FILTER_KEY:GT.getType().name())).
					acceptGenotype(acceptGT).
					name("Genotype-Type / Sample"+fltstr).
					description("Singleton per Sample"+fltstr).
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
					name("Sample GQ "+fltstr).
					description("Genotype Quality per Sample"+fltstr)
					);
			for(final GenotypeType genotype_type: new GenotypeType[] {GenotypeType.HET,GenotypeType.HOM_REF,GenotypeType.HOM_VAR}) {
					modules.add(new SampleToAD().
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
		
		
		modules.add(new RangeBarPlot("QD", 1).
				name("QD").
				description("Variant Confidence (QUAL) / Quality by Depth.").
				xlab("QD").
				ylab("Count Variants")
				);
		
		
		modules.add(new RangeBarPlot("FS", 1).
			logX(true).
			name("FQ").
			description("Phred-scaled p-value using Fisher's exact test to detect strand bias").
			xlab("FS").
			ylab("Count Variants")
			);
		
		modules.add(new RangeBarPlot("SOR", 10).
			name("SOR").
			description("Symmetric Odds Ratio of 2x2 contingency table to detect strand bias").
			xlab("SOR").
			ylab("Count Variants")
			);
		
		modules.add(new RangeBarPlot("MQ", 1).
			name("MQ").
			description("Mean square mapping quality over all the reads at the site").
			xlab("MQ").
			ylab("Count Variants")
			);

		modules.add(new RangeBarPlot("MQRankSum", 10).
			name("MQRankSum").
			description("Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities").
			xlab("MQRankSum").
			ylab("Count Variants")
			);
		
		modules.add(new RangeBarPlot("ReadPosRankSum", 10).
			name("ReadPosRankSum").
			description("Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias").
			xlab("ReadPosRankSum").
			ylab("Count Variants")
			);

		modules.add(new RangeBarPlot(VCFConstants.DEPTH_KEY, 1).
				logX(true).
				name(VCFConstants.DEPTH_KEY).
				description("DEPTH per variant").
				xlab(VCFConstants.DEPTH_KEY).
				ylab("Count Variants")
				);
		modules.add(new GatkDeNovo().
				name("GATK DeNovo").
				description("De Novo variants found with GATK").
				xlab("Sample").
				ylab("Count Variants")
				);
		modules.add(new DragenDeNovo().
				name("GATK DeNovo").
				description("De Novo variants found with GATK").
				xlab("Sample").
				ylab("Count Variants")
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
			for(Analyzer analyzer:modules) {
				analyzer.init(header);
				}
			modules.removeIf(M->!M.isEnabled());
			while(iter.hasNext()) {
				final VariantContext ctx = iter.next();
				for(Analyzer analyzer:modules) {
					analyzer.visit(ctx);
					}
				}
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(outputFile)) {
				for(Analyzer analyzer:modules) {
					analyzer.finish(pw);
					}
				pw.flush();
				}
			}
			return 0;
		} catch (final Throwable e) {
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
