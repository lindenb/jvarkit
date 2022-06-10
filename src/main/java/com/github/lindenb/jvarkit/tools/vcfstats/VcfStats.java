/*
The MIT License (MIT)
Copyright (c) 2022 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.Pair;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
/*
BEGIN_DOC

rewritten from scratch 2022

END_DOC
 */
@Program(name="vcfstats",
	description="Produce VCF statitics",
	keywords={"vcf","stats","R"},
	modificationDate = "20220509",
	creationDate = "20131212"
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
	@Parameter(names={"--prefix"},description="file refix")
	private String prefix = "";
	@DynamicParameter(names={"-D"},description="other parameters.")
	private Map<String,String> __dynaParams = new HashMap<>();


	private final Map<String, List<String>> phenotype2samples = new TreeMap<>();
	private final Map<String,String> sample2phenotype = new TreeMap<>();
	private final AttributeMap att = AttributeMap.verbose(AttributeMap.wrap(this.__dynaParams), (S)->{
		LOG.info("undefined parameter \"" + S + "\". Using default.");
		});
	
	private interface Analyzer {
		void init(VCFHeader h);
		void visit(final VariantContext ctx);
		void finish(PrintWriter pw);
		public String getName();
		public String getDescription();
		public boolean isEnabled();
		}
	private abstract class AbstractAnalyzer implements Analyzer {
		protected boolean enabled=true;
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
			return VcfStats.this.prefix+ getName();
			}
		@Override
		public String getDescription() {
			return getName();
			}
		String getInputName() {
			return "<STDIN>";
			}	
		}
	
	private abstract class AbstractBoxPlot<H extends Number> extends AbstractAnalyzer {
		private final Map<String, List<H>> cat2values= new HashMap<>();
		protected void add(final String cat,H value) {
			List<H> l = cat2values.get(cat);
			if(l==null) {
				l=new ArrayList<>();
				cat2values.put(cat, l);
				}	
			l.add(value);
			}
		public String getTitle() { return getName();}
		public String getSubTitle() { return getDescription();}
		public String getYLab() { return "todo axis y";}

		@Override
		public void finish(PrintWriter w) {
			w.println(device());
			final List<String> phenotypes = new ArrayList<>(cat2values.keySet());
			for(int i=0;i< phenotypes.size();i++) {
				final String phenotype = phenotypes.get(i);
				w.println("p"+i+" <- c("+
					cat2values.get(phenotype).
						stream().
						map(X->String.valueOf(X)).
						collect(Collectors.joining(",")) +
					")");
				}
			w.println("boxplot("+
				IntStream.range(0, phenotypes.size()).mapToObj(X->"p"+X).collect(Collectors.joining(","))+","+
				"main="+quote(getTitle())+",sub="+quote(getSubTitle())+",ylab="+quote(getYLab())+","+
				"las=2,"+
				"names=c("+ phenotypes.stream().map(S->quote(S)).collect(Collectors.joining(","))+")"+
				")");
			w.println("dev.off()");
			}
		}
	
	private abstract class AbstractBarPlot<H> extends AbstractAnalyzer {
		public String getXLab() { return "";}
		public String getYLab() { return "";}
		public String getTitle() { return getName();}
		public String getSubTitle() { return getDescription();}
		public String getColor(H h) { return "lightgray";}
		}
	
	private abstract class AbstractSimpleBarPlot<H> extends AbstractBarPlot<H> {
		protected final Counter<H> counter = new Counter<>();
		protected void add(final H value) {
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
	private abstract class AbstractMultipleBarPlot<H,V> extends AbstractBarPlot<H> {
		private final Map<H,Counter<V>> horiz2counts = new LinkedHashMap<>();
		private final Set<V> distinct_vertical = new LinkedHashSet<>();
		protected void add(final H h,final V v) {
			Counter<V> t = horiz2counts.get(h);
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
			final List<H> LH = new ArrayList<>(horiz2counts.keySet());
			final List<String> colNames= new ArrayList<>(LH.size());
			for(int i=0;i< LH.size();i++) {
				final H h = LH.get(i);
				final Counter<V> counter = horiz2counts.get(h);
				w.println("#  "+ h);
				final String colName = "var"+i;
				colNames.add(colName);
				w.print(colName+" <- c(");
				w.print(distinct_vertical.stream().
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
	private class Singleton extends AbstractSimpleBarPlot<String> {
		@Override
		public void init(VCFHeader h) {
			this.enabled = h.hasGenotypingData();
			}
		@Override
		public void visit(VariantContext ctx) {
			Genotype single = null;
			for(Genotype g:ctx.getGenotypes()) {
				if(g.isNoCall() || g.isHomRef()) continue;
				if(single!=null) return;
				single=g;
				}
			if(single!=null) add(single.getSampleName());
			}

		@Override public String getXLab() {return "Samples";}
		@Override public String getYLab() {return "Count Variants";}
		@Override
		public String getName() {
			return "singletons";
			}
		@Override
		public String getDescription() {
			return "singletons per sample";
			}
		}

	/***************************************************************************/
	private class SVTypeContig extends AbstractMultipleBarPlot<String, String> {
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.getInfoHeaderLine(VCFConstants.SVTYPE)!=null;
			}
		@Override
		public void visit(VariantContext ctx) {
			final String svType = ctx.getAttributeAsString(VCFConstants.SVTYPE, null);
			if(StringUtils.isBlank(svType)) return;
			super.add(ctx.getContig(), svType);
			}
		@Override public String getXLab() {return "Contig";}
		@Override public String getYLab() {return "SVTypes";}
		@Override
		public String getName() {
			return "svtype2contig";
			}
		@Override
		public String getDescription() {
			return "SVTYpe per contig";
			}
		}

	/***************************************************************************/
	private class Sample2GTType extends AbstractMultipleBarPlot<String, GenotypeType> {
		@Override
		public void init(final VCFHeader h) {
			this.enabled = h.hasGenotypingData();
			}
		@Override
		public void visit(VariantContext ctx) {
			for(final Genotype gt: ctx.getGenotypes()) {
				super.add(gt.getSampleName(), gt.getType());	
				}
			}
		@Override public String getXLab() {return "Sample";}
		@Override public String getYLab() {return "Genotype Types";}
		@Override
		public String getName() {
			return "gtPerSample";
			}
		@Override
		public String getDescription() {
			return "SVTYpe per contig";
			}
		}
	/*********************************************************************/
	private class SVLen extends AbstractBoxPlot<Integer> {
		@Override
		public void init(VCFHeader h) {
			this.enabled = h.getInfoHeaderLine(VCFConstants.SVTYPE)!=null;
			}
		@Override
		public String getName() {
			return "SVLEN";
			}
		@Override
		public void visit(VariantContext ctx) {
			final String st = ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
			if(StringUtils.isBlank(st) || st.equalsIgnoreCase("BND")) return;
			super.add(st, ctx.getLengthOnReference());
			}
		}
	/*********************************************************************/
	private class SVOverlap1 extends AbstractMultipleBarPlot<String,Integer> {
		private final IntervalTreeMap<Locatable> gnomadTreeMap = new IntervalTreeMap<>();
		private String gnomadSource;
		@Override
		public void init(VCFHeader h) {
			this.gnomadSource = VcfStats.this.att.getAttribute("external.database.vcf").orElse(null);
			if(StringUtils.isBlank(this.gnomadSource)) {
				this.enabled = false;
				return;
				}
			final ContigNameConverter convert =ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(h));
			try(VCFIterator iter = new VCFIteratorBuilder().open(this.gnomadSource)) {
				while(iter.hasNext()) {
					final VariantContext vc = iter.next();
					final String ctg = convert.apply(vc.getContig());
					if(StringUtils.isBlank(ctg)) continue;
					final Interval r = new Interval(ctg,vc.getStart(),vc.getEnd());
					gnomadTreeMap.put(r, r);
					}
				}
			catch(final Throwable err) {
				LOG.error(err);
				this.enabled=false;
				}
			
			
			}
		@Override
		public String getName() {
			return "svoverlap";
			}
		@Override
		public void visit(VariantContext ctx) {
			final String st = ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
			if(StringUtils.isBlank(st) || st.equalsIgnoreCase("BND")) return;
			super.add(st, ctx.getLengthOnReference());
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
		modules.add(new Singleton());
		modules.add(new Sample2GTType());
		modules.add(new SVTypeContig());
		modules.add(new SVLen());
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
		
		
		try(VCFIterator iter= super.openVCFIterator(oneFileOrNull(args))) {	
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
