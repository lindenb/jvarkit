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
package com.github.lindenb.jvarkit.tools.burden;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.date.DurationParser;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.ProcessExecutor;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

@Program(name="optimizervtest",
description="Optimize rvtest",
keywords= {"vcf","burden","fisher"},
modificationDate="20221013",
creationDate="20221013",
generate_doc = false
)
public class OptimizeRVTest extends Launcher {
	private static final Logger LOG = Logger.build( OptimizeRVTest.class).make();
	private static long ID_GENERATOR=0L;
	@Parameter(names={"--pedigree"},description="path to the pedigree",required = true)
	private Path pedPath;
	@Parameter(names={"--rvtest"},description="path rvtest")
	private String rvtestPath="rvtest";
	@Parameter(names={"--n-samples-per-generation","-nsg"},description="Number of samples per generation.")
	private int n_samples_per_generation = 5;
	@Parameter(names={"--tmp-dir","-T"},description="Temporary directory")
	private Path tmpDirectory = IOUtils.getDefaultTempDir();
	@Parameter(names={"--chiasma"},description="Proba chiasma")
	private double proba_chiasma  = 0.01;
	@Parameter(names={"--mutation"},description="Proba mutation")
	private double proba_mutation  = 0.01;
	@Parameter(names={"--duration"},description=DurationParser.OPT_DESC)
	private String durationStr  = "23h";
	@Parameter(names={"--disable"},description="disable factory(ies) by name. Comma separated")
	private String disableFactoriesNames="";

	/*
	private static interface OpCompare {
		public String getName();
		public boolean test(double v1,double v2);
		}

	private static final OpCompare LT_compare = new OpCompare() {
		@Override public String getName() { return "<";};
		@Override public boolean test(double v1, double v2) {return Double.compare(v1, v2)<0;}
		};
	private static final OpCompare LE_compare = new OpCompare() {
		@Override public String getName() { return "<=";};
		@Override public boolean test(double v1, double v2) {return Double.compare(v1, v2)<=0;}
		};
	private static final OpCompare GT_compare = new OpCompare() {
		@Override public String getName() { return ">";};
		@Override public boolean test(double v1, double v2) {return Double.compare(v1, v2)>0;}
		};
	private static final OpCompare GE_compare = new OpCompare() {
		@Override public String getName() { return ">=";};
		@Override public boolean test(double v1, double v2) {return Double.compare(v1, v2)>=0;}
		};
	
	
	private class Range {
		double min;
		double max;
		Range(double x,double y) {
			this.min = Math.min(x, y);
			this.max = Math.max(x, y);
			}
		}*/
	
	private interface Trait {
		public TraitFactory getFactory();
		
		public default String asString() {
			return getFactory().asString(this);
			}
		public default boolean test(final VariantContext ctx) {
			return getFactory().test(ctx, this);
			}
		public default boolean isSame(Trait t) {
			return getFactory() == t.getFactory() && t.getFactory().isSame(this,t);
			}
		public default void export(BufferedWriter w) throws IOException {
			getFactory().export(w,this);
			}
		}
	
	private abstract class TraitFactory {
		boolean enabled = true;
		protected class AbstractTrait implements Trait {
			@Override
			public final TraitFactory getFactory() {
				return TraitFactory.this;
				}
			@Override
			public final String toString() {
				return asString();
				}
			}
		void initialize(VCFHeader v) {
			
			}
		abstract void initialize(VariantContext v);
		boolean isEnabled() {
			return enabled;
			}
		Trait createFirst() {
			return create();
			}

		abstract Trait create();
		Trait mute(Trait t) {
			return create();
			}
		protected Genotype getSingleton(final VariantContext ctx)  {
			Genotype singleton = null;
			for(Genotype g:ctx.getGenotypes()) {
				final String sn = g.getSampleName();
				if(!cases.contains(sn) && !controls.contains(sn)) continue;
				if(g.isHet() || g.isHomVar()) {
					if(singleton==null) {
						singleton  = g;
						}
					else
						{
						return null;
						}
					}
				}
			return singleton;
			}
		abstract String getName();
		abstract String asString(Trait t);
		abstract boolean test(VariantContext ctx,Trait v);
		abstract boolean isSame(Trait a,Trait b);
		abstract void export(BufferedWriter w,Trait v) throws IOException;
		}
	
	
	private abstract class AbstractNumberFactory extends TraitFactory {
		protected class MyTrait extends AbstractTrait {
			double v;
			MyTrait(final double v) {
				this.v = v;
				}
			}
		
		private double minValue=0;
		private double maxValue=0;
		private boolean init_seen=false;
		protected void initialize(final double v) {
			if(!init_seen) {
				minValue=v;	
				maxValue=v;	
				init_seen=true;
				}
			else
				{
				minValue = Math.min(minValue, v);
				maxValue = Math.max(maxValue, v);
				}
			}
		
		protected double getMin() { return this.minValue;}
		protected double getMax() { return this.maxValue;}
		
		private double delta() {
			return getMax() - getMin();
			}
		
		protected double makeNumber() {
			return getMin() + (random.nextDouble() * delta());
			}
		protected double makeFirstNumber() {
			return makeNumber();
			}

		Trait mute(Trait t) {
			final MyTrait t2 = MyTrait.class.cast(t);
			double v = t2.v;
			if(random.nextBoolean()) {
				v+= (getMax() - v)/100.0;
				}
			else
				{
				v-= (v-getMin())/100.0;
				}
			return new MyTrait(v);
			}
		
		
		@Override
		boolean isEnabled() {
			return getMin() < getMax();
			}
		
		@Override
		Trait createFirst() {
			double v = makeFirstNumber();
			if(v<getMin()) {
				throw new IllegalArgumentException("bad init value for "+getName()+" "+v+" < "+getMin());
				//v=getMin();
				}
			if(v>getMax()) {
				throw new IllegalArgumentException("bad init value for "+getName()+" "+v+">"+getMax());
				//v=getMax();
				}
			return new MyTrait(v);
			}

		
		@Override
		Trait create() {
			return new MyTrait(makeNumber());
			}


		@Override
		boolean isSame(Trait a, Trait b) {
			final MyTrait ta = MyTrait.class.cast(a);
			final MyTrait tb = MyTrait.class.cast(b);
			return Double.compare(ta.v,tb.v)==0 ;
			}
		@Override
		void export(BufferedWriter w,Trait t) throws IOException {
			w.write(getName());
			w.write("\t");
			w.write(String.valueOf(MyTrait.class.cast(t).v));
			w.write("\n");
			}

		
		@Override
		public String toString() {
			return getClass().getSimpleName()+" min:"+minValue+" max:"+maxValue+" enabled:"+isEnabled();
			}
		}
	
	
	
	
	private abstract class  SimpleDoubleTagFactory  extends AbstractNumberFactory {
		
		@Override
		void initialize(final VariantContext ctx) {
			if(!ctx.hasAttribute(getTag())) return;
			this.initialize(ctx.getAttributeAsDouble(getTag(), 0));
			}
	
		
		abstract String getTag();
		
		@Override
		String getName() {
			return getTag();
			}
		
		}
	
	
	/** FS **/
	private class  FSFactory  extends SimpleDoubleTagFactory {
		@Override
		String asString(Trait t) {
			return getTag() +" =< " + MyTrait.class.cast(t).v;
			}
		@Override
		protected double makeFirstNumber() {
			return 10.0;
			}
		@Override
		boolean test(final VariantContext ctx,final Trait t) {
			if(!ctx.hasAttribute(getTag())) return true;
			return ctx.getAttributeAsDouble(getTag(), 0.0) <= MyTrait.class.cast(t).v;
			}
		@Override
		String getTag() {
			return "FS";
			}
		}
	/** MQ **/
	private class RMSMappingQualityFactory  extends SimpleDoubleTagFactory {
		@Override
		String asString(Trait t) {
			return getTag() +" >= " + MyTrait.class.cast(t).v;
			}
		@Override
		protected double makeFirstNumber() {
			return 59.0;
			}
		@Override
		boolean test(VariantContext ctx, Trait t) {
			if(!ctx.hasAttribute(getTag())) return true;
			return ctx.getAttributeAsDouble(getTag(), 0.0) >= MyTrait.class.cast(t).v;
			}		@Override
		String getTag() {
			return "MQ";
			}
		}
	/** MQRankSum **/
	private class MQRankSumFactory  extends SimpleDoubleTagFactory {
		@Override
		String asString(Trait t) {
			final double x = Math.abs(MyTrait.class.cast(t).v);
			return ""+(-x)+" <= "+getTag()+" <= "+ x;
			}
		@Override
		protected double makeFirstNumber() {
			return 3.0;
			}
		@Override
		boolean test(VariantContext ctx, Trait t) {
			if(!ctx.hasAttribute(getTag())) return true;
			final double x = Math.abs(MyTrait.class.cast(t).v);
			return (-x) <= ctx.getAttributeAsDouble(getTag(), 0.0) &&
					ctx.getAttributeAsDouble(getTag(), 0.0) < x;
			}	
		@Override
		String getTag() {
			return "MQRankSum";
			}
		}
	
	/** MISSING ./. */
	private class MissingFactory extends AbstractNumberFactory {
		@Override
		void initialize(VariantContext ctx) {
			initialize(getMissing(ctx));
			}
		@Override
		String asString(final Trait t) {
			return getName() +" <= " + MyTrait.class.cast(t).v;
			}
		
		@Override
		protected double makeFirstNumber() {
			return 0.01;
			}
		
		double getMissing(final VariantContext ctx) {
			double n=0.0;
			double miss = 0.0;
			for(Genotype g:ctx.getGenotypes()) {
				final String sn = g.getSampleName();
				if(!cases.contains(sn) && !controls.contains(sn)) continue;
				n++;
				if(g.isNoCall()) miss++;
				}
			if(n==0) throw new IllegalStateException();
			final double f_missing = miss/n;
			return f_missing;
			}
		
		@Override
		boolean test(final VariantContext ctx, final Trait trait) {
			final double f_missing = getMissing(ctx);
			return f_missing <= MyTrait.class.cast(trait).v;
			}
		@Override
		java.lang.String getName() {
			return "F_MISSING";
			}
		}
	
	/** Singleton AD */
	private class SingletonADFactory extends AbstractNumberFactory {
		@Override
		void initialize(VariantContext ctx) {
			initialize(0.1);
			initialize(0.5);
			}
		
		@Override
		String asString(final Trait t) {
			return MyTrait.class.cast(t).v + "<="+ getName() +" <= " + (1.0-MyTrait.class.cast(t).v);
			}
		
		@Override
		protected double makeFirstNumber() {
			return 0.3;
			}
		
		@Override
		boolean test(final VariantContext ctx, final Trait trait) {
			Genotype singleton = getSingleton(ctx);
			if(singleton==null) return true;
			if(!singleton.isHet()) return true;
			if(!singleton.hasAD()) return true;
			int[] ad = singleton.getAD();
			if(ad.length!=2) return true;
			double t = ad[0]+ad[1];
			if(t==0) return true;
			double f=ad[1]/t;
			if(f <  MyTrait.class.cast(trait).v) return false;
			if(f >  (1.0-MyTrait.class.cast(trait).v)) return false;
			return true;
			}
		@Override
		String getName() {
			return "SingletonADRatio";
			}
		}

	/** Singleton GQ */
	private class SingletonGQFactory extends AbstractNumberFactory {
		@Override
		void initialize(VariantContext ctx) {
			final Genotype g = getSingleton(ctx);
			if(g==null || !g.hasGQ()) return;
			initialize((double)g.getGQ());	
			}
		@Override
		String asString(final Trait t) {
			return getName() +" >= " + MyTrait.class.cast(t).v;
			}
		
		@Override
		protected double makeFirstNumber() {
			return 80.0;
			}
		
		@Override
		boolean test(final VariantContext ctx, final Trait trait) {
			final Genotype singleton = getSingleton(ctx);
			if(singleton==null) return true;
			if(!singleton.hasGQ()) return true;
			return singleton.getGQ() >= MyTrait.class.cast(trait).v;
			}
		@Override
		String getName() {
			return "SingletonGQ";
			}
		}
	
	//
	/** Gnomad AF */
	private class GnomadAFFactory extends AbstractNumberFactory {
		final String tag = "gnomad_genome_AF_POPMAX";
		
		private double getAF(VariantContext ctx) {
			return ctx.getAttributeAsDoubleList(tag, 0.0).stream().
					mapToDouble(V->V.doubleValue()).min().orElse(0.0);
			}
		@Override
		void initialize(VariantContext ctx) {
			initialize(getAF(ctx));
			}
		@Override
		String asString(final Trait t) {
			return getName() +" <= " + MyTrait.class.cast(t).v;
			}
		@Override
		protected double makeFirstNumber() {
			return 0.01;
			}
		
		@Override
		boolean test(final VariantContext ctx, final Trait trait) {
			return getAF(ctx) <= MyTrait.class.cast(trait).v ;
			}
		@Override
		String getName() {
			return this.tag;
			}
		}

	private class CaddPhredFactory extends AbstractNumberFactory {
		final String tag = "CADD_PHRED";
		
		private OptionalDouble getScore(VariantContext ctx) {
			return ctx.getAttributeAsStringList(tag,"0").
					stream().
					filter(S->StringUtils.isDouble(S)).
					mapToDouble(S->Double.parseDouble(S)).
					filter(V->V>=0).// unknown value are -999
					min();
			}
		@Override
		void initialize(VariantContext ctx) {
			final OptionalDouble opt = getScore(ctx);
			if(opt.isPresent()) initialize(opt.getAsDouble());
			}
		@Override
		String asString(final Trait t) {
			return getName() +" >= " + MyTrait.class.cast(t).v;
			}
		@Override
		boolean test(final VariantContext ctx, final Trait trait) {
			final OptionalDouble opt = getScore(ctx);
			if(!opt.isPresent()) return true;
			return opt.getAsDouble() >= MyTrait.class.cast(trait).v ;
			}
		@Override
		String getName() {
			return this.tag;
			}
		}
	
	private class PolyXFactory extends AbstractNumberFactory {
		final String tag = "POLYX";
		
		private int getScore(VariantContext ctx) {
			if(!ctx.hasAttribute(tag)) return 0;
			return ctx.getAttributeAsInt(tag,0);
			}
		@Override
		void initialize(VariantContext ctx) {
			initialize(getScore(ctx));
			}
		@Override
		String asString(final Trait t) {
			return getName() +"<= " + MyTrait.class.cast(t).v;
			}
		@Override
		boolean test(final VariantContext ctx, final Trait trait) {
			return getScore(ctx) <= MyTrait.class.cast(trait).v ;
			}
		@Override
		String getName() {
			return this.tag;
			}
		}
	
	private abstract class AbstractPredictionFactory  extends TraitFactory {
		protected final Set<String> all_acns = new HashSet<>();
		protected class MyTrait extends AbstractTrait {
			final Set<String> acns;
			MyTrait(Set<String> acns) {
				this.acns = acns;
				}
			}
		
		@Override
		Trait createFirst() {
			final String[] acns = new String[]{"missense_variant","stop_gained","stop_lost","start_lost","frameshift_variant"};
			final Set<String> set = new HashSet<>(acns.length);
			for(String s: acns) {
				if(!this.all_acns.contains(s)) {
					continue;
					}
				set.add(s);
				}
			if(!set.isEmpty()) {
				return new MyTrait(set);
				}
			else
				{
				return create();
				}
			}
		
		@Override
		Trait create() {
			return new MyTrait(makeSetOfSOTerms());
			}
		
		
		
		@Override
		boolean isSame(Trait a, Trait b) {
			if(a.getFactory()!=b.getFactory()) return false;
			return MyTrait.class.cast(a).acns.equals(MyTrait.class.cast(b).acns);
			}
		abstract String getTag();
		
		Set<String> makeSetOfSOTerms() {
			return makeSetOfSOTerms(this.all_acns);
			}
		
		Set<String> makeSetOfSOTerms(Set<String> collection) {
			final Set<String> set = new HashSet<>();
			final List<String> remains = new ArrayList<>(collection);
			do {
				set.add(remains.get(random.nextInt(remains.size())));
				} while(!remains.isEmpty() && random.nextDouble() <= 0.1);
			return set;
			}

		@Override
		Trait mute(Trait t) {
			MyTrait t2 = MyTrait.class.cast(t);
			double f=random.nextDouble();
			MyTrait t3= new MyTrait(t2.acns);
			if(f<0.5) {
				t3.acns.addAll(makeSetOfSOTerms());
				}
			else 
				{
				Set<String> set2 = makeSetOfSOTerms(t3.acns);
				t3.acns.clear();
				t3.acns.addAll(set2);
				}
			return t3;
			}
		
		@Override
		void export(BufferedWriter w, final Trait t) throws IOException {
			w.write(getName());
			w.write("\t");
			w.write(String.join(",", MyTrait.class.cast(t).acns));
			w.write("\n");
			}
		
		@Override
		String asString(final Trait t) {
			return getTag()+" as SO ("+ String.join(",", MyTrait.class.cast(t).acns)+")";
			}
		@Override
		String getName() {
			return getTag();
			}
		@Override
		public String toString() {
			return getName()+" "+String.join(" ", this.all_acns);
			}
		}

	private class VepPredictionFactory  extends AbstractPredictionFactory {
		private VepPredictionParser vep;
		@Override
		void initialize(final VCFHeader header) {
			super.initialize(header);
			this.vep = new VepPredictionParserFactory(header).get();
			}
		@Override
		void initialize(final VariantContext ctx) {
			vep.getPredictions(ctx).stream().
					flatMap(P->P.getSOTermsStrings().stream()).
					forEach(S->all_acns.add(S));
			}
		
		@Override
		boolean test(final VariantContext ctx,final Trait v) {
			final MyTrait mt = MyTrait.class.cast(v);
			return vep.getPredictions(ctx).stream().
					flatMap(P->P.getSOTermsStrings().stream()).
					anyMatch(S->mt.acns.contains(S));
			}
		@Override
		boolean isEnabled() {
			return this.vep.isValid() && !all_acns.isEmpty();
			}
		@Override
		String getTag() {
			return "CSQ";
			}
		}
	
	private final List<TraitFactory> traitFactories = new ArrayList<>();
	private long n_generations = 0;
	private Random random = new Random(System.currentTimeMillis());
	private final Set<String> cases = new HashSet<>();
	private final Set<String> controls = new HashSet<>();

	
	private class Solution {
		final long id = (++OptimizeRVTest.ID_GENERATOR);
		final long generation = OptimizeRVTest.this.n_generations;
		double pvalue=1.0;
		final List<Trait> traits;
		final List<String> assoc=new ArrayList<>();
		Solution() {
			this(false);
			}
		Solution(boolean is_first) {
			this.traits = new ArrayList<>(traitFactories.size());
			for(TraitFactory t: traitFactories) {
				traits.add(is_first?t.createFirst():t.create());
				}
			}
		
		int size() {
			return traitFactories.size();
			}
		
		final List<Interval> saveVCF(final Path inVCF,final String filename) throws IOException {
			final List<Interval> L = new ArrayList<>();
			try(VCFIterator iter = new VCFIteratorBuilder().open(inVCF)) {
				final VCFHeader hdr = iter.getHeader();
				hdr.addMetaDataLine(new VCFHeaderLine("PVALUE", String.valueOf(this.pvalue)));
				try( VariantContextWriter vw=  new VariantContextWriterBuilder().
						setOutputPath(Paths.get(filename+".vcf.gz")).
						setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).
						setOption(Options.INDEX_ON_THE_FLY).
						build()) {
					vw.writeHeader(hdr);
					while(iter.hasNext()) {
						final VariantContext ctx=iter.next();
						if(!this.traits.stream().allMatch(T->T.test(ctx))) continue;
						vw.add(ctx);
						L.add(new Interval(ctx));
						}
					}
				}
			if(!this.assoc.isEmpty()) {
				Files.write(Paths.get(filename+".assoc"), this.assoc);
				}
			Collections.sort(L);
			int i=0;
			while(i+1 < L.size()) {
				if(L.get(i).overlaps(L.get(i+1))) {
					L.set(i,new Interval(
							L.get(i).getContig(),
							Math.min(L.get(i).getStart(), L.get(i+1).getStart()),
							Math.max(L.get(i).getEnd(), L.get(i+1).getEnd())
							));
					L.remove(i+1);
					}
				else
					{
					i++;
					}
				}
			return L;
			}

		
		void run(final Path inputVCF) throws IOException {
			LOG.info("running "+this.toString());
			final String filename= tmpDirectory.resolve("tmp").toString();
			this.pvalue=1.0;
			final List<Interval> intervals = saveVCF(inputVCF,filename);
			if(intervals.isEmpty()) {
				LOG.info("No interval for this");
				return;
				}
			
			try(BufferedWriter w = Files.newBufferedWriter(Paths.get(filename+".setfile"))) {
				w.write("REGION"+this.id+"\t"+
					intervals.stream().map(R->R.getContig()+":"+R.getStart()+"-"+R.getEnd()).collect(Collectors.joining(","))+"\n"
					);
				w.flush();
				}
			
			
			final StringBuilder commandb = new StringBuilder();
			commandb.append(rvtestPath + " --pheno "+pedPath+"  --setFile "+filename+".setfile --inVcf "+filename +".vcf.gz  --burden exactCMC --out "+filename);
			final String command = commandb.toString();
			final int ret = ProcessExecutor.execute(command);
			if(ret!=0)
				{
				throw new IOException("Error RVTEST returned "+ret+" (!=0).");
				}
			try(BufferedReader br = Files.newBufferedReader(Paths.get(filename+".CMCFisherExact.assoc"))) {
				String line;
				this.assoc.clear();
				while((line=br.readLine())!=null) {
					this.assoc.add(line);
					if(line.startsWith("Range")) continue;
					String[] tokens = line.split("[\\s]+");
					String pvalue = tokens[9];
					try {
						this.pvalue = Double.parseDouble(pvalue);
						LOG.info("ID:"+id+"="+pvalue);
						}
					catch(final Throwable err) {
						err.printStackTrace();
						this.pvalue=1.0;
						}
					}
				}
			
			}
			
		boolean isSame(final Solution other) {
			for(int i=0;i< size();i++) {
				if(!traits.get(i).isSame(other.traits.get(i))) return false;
				}
			return true;
			}
		@Override
		public String toString() {
			return "G"+this.generation+" ID"+id+" P="+pvalue+" traits:["+ traits.stream().map(T->T.toString()).collect(Collectors.joining(" && ")) +"]";
			}
		}
	
	
	private Solution mate(final Solution s1,final Solution s2)  {
		final Solution c = new Solution();
		int side= 0;//always start with same , no when scanning X*Y, we're always starting with the other sample.
		for(int i=0;i< s1.size();i++) {
			c.traits.set(i,(side==0?s1:s2).traits.get(i));
			if(random.nextDouble()< proba_chiasma) {
				side=(side==0?1:0);
				}
			}
		while(this.random.nextDouble()< proba_mutation ) {
			int n = this.random.nextInt(s1.traits.size());
			Trait t = c.traits.get(n);
			c.traits.set(n,traitFactories.get(n).mute(t));
			}
		return c;
		}

	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Path vcfPath = Paths.get(oneAndOnlyOneFile(args));
			IOUtil.assertFileIsReadable(vcfPath);
			IOUtil.assertDirectoryIsWritable(this.tmpDirectory);
			
			
			/* get case controls */
			try(BufferedReader br = Files.newBufferedReader(this.pedPath)) {
				String line;
				while((line=br.readLine())!=null) {
					if(line.startsWith("fid\t")) continue;
					final String[] tokens = line.split("[\\s]+");
					if(tokens[5].equals("2")) {
						this.cases.add(tokens[1]);
						}
					else if(tokens[5].equals("1")) {
						this.controls.add(tokens[1]);
						}
					}
				}

			
			this.traitFactories.addAll(
				Arrays.asList(
					new VepPredictionFactory(),
					new GnomadAFFactory(),
					new FSFactory(),
					new RMSMappingQualityFactory(),
					new MQRankSumFactory(),
					new MissingFactory(),
					new SingletonADFactory(),
					new SingletonGQFactory(),
					new PolyXFactory(),
					new CaddPhredFactory()
					));
			
			
			final Set<String> disable_names = Arrays.stream(disableFactoriesNames.split("[, ;]+")).filter(S->!StringUtil.isBlank(S)).map(S->S.toLowerCase()).collect(Collectors.toSet());
			LOG.info("the following factories will be disabled: "+this.traitFactories.stream().
					filter(F->disable_names.contains(F.getName().toLowerCase())).
					map(F->F.getName()).collect(Collectors.joining(",")));
			this.traitFactories.removeIf(F->disable_names.contains(F.getName().toLowerCase()));
			
			
			LOG.info("initialize factories");
			try(VCFIterator iter = new VCFIteratorBuilder().open(vcfPath)) {
				final VCFHeader header = iter.getHeader();
				//initialize with header
				for(TraitFactory trait:traitFactories) {
					trait.initialize(header);
					}
				while(iter.hasNext()) {
					final VariantContext ctx = iter.next();
					//initialize with current variant
					for(TraitFactory trait:traitFactories) {
						trait.initialize(ctx);
						}
					}
				
				}
			for(TraitFactory trait:traitFactories) {
				LOG.info(trait.toString()+" enabled:"+trait.isEnabled());
			}
			
			traitFactories.removeIf(P->!P.isEnabled());
			
			LOG.info("initialize factories: done");
			if(this.traitFactories.isEmpty()) {
				LOG.info("No trait factory");
				return -1;
			}
			


			LOG.info("cases: "+this.cases.size());
			LOG.info("controls: "+this.controls.size());
			final long stop = System.currentTimeMillis() + new DurationParser().convert(this.durationStr).toMillis();
			
			//write header
			try(BufferedWriter w = Files.newBufferedWriter(Paths.get("best.tsv"),StandardOpenOption.CREATE)) {
				w.write("TIMESTAMP\tDATE\tID\tGENERATION\tPVALUE");
				for(TraitFactory tf:this.traitFactories) w.write("\t"+tf.getName());
				w.write("\n");
				w.flush();
				}

			
			
			Solution best = new Solution(true /* invoke createFirst() */);
			best.run(vcfPath);
			best.saveVCF(vcfPath,"best");

			
			List<Solution> generation = new ArrayList<>();
			generation.add(best);
			while(generation.size()< this.n_samples_per_generation) {
				final Solution c= new Solution();
				c.run(vcfPath);
				generation.add(c);
				}
			while(System.currentTimeMillis()< stop) {
				++n_generations;
				for(Solution sol:generation) {
					if(sol.pvalue < best.pvalue ) {
						LOG.info("new best "+sol);
						best=sol;
						best.saveVCF(vcfPath,"best");

						try(BufferedWriter w = Files.newBufferedWriter(Paths.get("best.tsv"),StandardOpenOption.APPEND,StandardOpenOption.CREATE)) {
							w.write(""+System.currentTimeMillis()+"\t"+StringUtils.now()+"\t"+best.id+"\t"+best.generation+"\t"+best.pvalue);
							for(int i=0;i< best.traits.size();i++) {
								w.write("\t"+ best.traits.get(i));
								}
							w.write("\n");
							w.flush();
							}
						try(BufferedWriter w = Files.newBufferedWriter(Paths.get("best.factory.tsv"),StandardOpenOption.CREATE)) {
							for(int i=0;i< best.traits.size();i++) {
								best.traits.get(i).export(w);
								}
							w.flush();
							}

						}
					}
				LOG.info("Generation "+n_generations+" best:"+best.pvalue+" "+TimeUnit.MILLISECONDS.toMinutes(stop-System.currentTimeMillis())+" minutes");
				final Solution external = new Solution();
				if(generation.stream().noneMatch(S->S.isSame(external)) && !external.isSame(best)) {
					external.run(vcfPath);
					generation.add(external);
					}
				final List<Solution> generation2 = new ArrayList<>(generation.size()*generation.size());
				for(int x=0;x< generation.size();++x) {
					for(int y=0;y< generation.size();++y) {
						final Solution c = mate(generation.get(x),generation.get(y));
						if(generation2.stream().anyMatch(S->S.isSame(c))) continue;
						
						final Solution same = generation.stream().filter(S->S.isSame(c)).findFirst().orElse(null);
						if(same!=null) {
							c.pvalue = same.pvalue;
							}
						else
							{
							c.run(vcfPath);
							}
						if(c.pvalue==1.0) continue;
						generation2.add(c);
						}
					}
				Collections.sort(generation2,(A,B)->Double.compare(A.pvalue,B.pvalue));
				generation = generation2.subList(0, Math.min(generation2.size(),this.n_samples_per_generation));
				while(generation.size()< this.n_samples_per_generation) {
					final Solution x = mate(best, best);
					generation.add(x);
					}
				}
			LOG.info("Done. Best "+best.pvalue);
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new OptimizeRVTest().instanceMainWithExit(args);
	}
}
