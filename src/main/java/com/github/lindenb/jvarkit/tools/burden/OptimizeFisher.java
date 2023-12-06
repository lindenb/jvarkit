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
import java.io.Reader;
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
import java.util.Properties;
import java.util.Random;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.date.DurationParser;
import com.github.lindenb.jvarkit.gatk.GATKConstants;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

## Motivation

loads a (small) VCF and memory and, using a genetic-algorithm, tries, to find the best set of parameters, the best genomic window to find the lowest FisherTest for a burden test.

## Example

```
java -jar TMP/jvarkit.jar optimizefisher --cases TMP/cases.txt --controls TMP/ctrls.txt -o TMP TMP/normalized.vcf.gz
```

END_DOC
*/
@Program(name="optimizefisher",
description="Optimize fisher test on VCF using genetic algo",
keywords= {"vcf","burden","fisher"},
modificationDate="20231205",
creationDate="20221013"
)
public class OptimizeFisher extends Launcher {
	private static final Logger LOG = Logger.build( OptimizeFisher.class).make();
	private static long ID_GENERATOR=0L;
	@Parameter(names={"--cases"},description="path to a file containing the cases samples",required = true)
	private Path casesSamplePath=null;
	@Parameter(names={"--controls"},description="path to a file containing the controls samples",required = true)
	private Path controlsSamplePath=null;
	@Parameter(names={"--properties"},description="optional path to a java property file defining the preferences. Such kind is saved by a run.", hidden = true)
	private Path propertyFile =null;

	
	@Parameter(names={"--n-samples-per-generation","-nsg"},description="Number of samples per generation.")
	private int n_samples_per_generation = 5;
	@Parameter(names={"--output","-o"},description="Output directory",required=true)
	private Path  outputDir = null;
	@Parameter(names={"--chiasma"},description="Proba chiasma")
	private double proba_chiasma  = 0.01;
	@Parameter(names={"--mutation"},description="Proba mutation")
	private double proba_mutation  = 0.01;
	@Parameter(names={"--duration"},description=DurationParser.OPT_DESC)
	private String durationStr  = "23h";
	@Parameter(names={"--disable"},description="disable factory(ies) by name. Comma separated")
	private String disableFactoriesNames="";
	@Parameter(names={"--min-window-size"},description="min sliding window side")
	private int min_sliding_window_size = 100;
	@Parameter(names={"--max-window-size"},description="max sliding window side")
	private int max_sliding_window_size = 100_000;

	
	
	private VCFHeader vcfHeader;
	private final List<VariantContext> variants = new Vector<>(100_00);
	
	
	private class SlidingWindow {
		 int window_size;
		int window_shift;
		
		SlidingWindow() {
			this.window_size =rnd(min_sliding_window_size,max_sliding_window_size);
			this.window_shift = this.window_size/rnd(2,5);
		}
		
		void fill(final Properties props) {
			props.setProperty("window_size", String.valueOf(window_size));
			props.setProperty("window_shift", String.valueOf(window_shift));
		}
		
		void initialize(final Properties props) {
			final String key1 = "window_size";
			final String key2 = "window_shift";
			if(props.containsKey(key1) && props.containsKey(key2)) {
				this.window_size = Integer.parseInt(props.getProperty(key1));
				this.window_shift = Integer.parseInt(props.getProperty(key2));
				}
			}
		SlidingWindow mute() {
			final SlidingWindow sl = new SlidingWindow();
			if(random.nextBoolean()) {
				sl.window_size = Math.max(min_sliding_window_size, (int)(this.window_size* random.nextDouble()));
				}
			else
				{
				sl.window_size = Math.min(max_sliding_window_size, (int)(this.window_size* (1.0 + random.nextDouble())));
				}
			sl.window_shift = sl.window_size/rnd(2,5);
			return sl;
		}	
		
		
		@Override
		public String toString() {
			return "win_size:="+window_size+" win_shift:="+window_shift;
			}
	}
	
	/** an instance of "parameter in the genetic algorithm */
	private interface Trait{
		/** get Owner Factory */
		public TraitFactory getFactory();
		
		public default String asString() {
			return getFactory().asString(this);
			}
		/** compare two traits */
		public default boolean isSame(Trait t) {
			return getFactory() == t.getFactory() && t.getFactory().isSame(this,t);
			}
		}
	
	/** Trait creator */
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
		void initialize(final Properties props) {
			
			}
		void fill(final Trait trait,final Properties props) {
			
			}
		void initialize(final VCFHeader v) {
			
			}
		abstract void initialize(final VariantContext v);
		boolean isEnabled() {
			return enabled;
			}
		/** create the traits the first time */
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
	
	private static short STATE_UNDEFINED=0;
	private static short STATE_SET_BY_USER = 1;
	private static short STATE_FILL = 2;
	
	
	private static interface CompareOp {
		public String getSymbol();
		public boolean test(double v1,double v2);
		}
	private static final CompareOp LowerOrEqual = new  CompareOp() {
		@Override public String getSymbol() { return "<=" ;}
		@Override  public boolean test(double v1,double v2) { return v1 <= v2; }
		};
	private static final CompareOp GreaterOrEqual = new  CompareOp() {
		@Override public String getSymbol() { return ">=" ;}
		@Override  public boolean test(double v1,double v2) { return v1 >= v2; }
		};
	
	private abstract class AbstractNumberFactory extends TraitFactory {
		protected class MyTrait extends AbstractTrait {
			final double v;
			MyTrait(final double v) {
				this.v = v;
				}
			}
		private OptionalDouble startValue = OptionalDouble.empty();
		private double minValue=0;
		private double maxValue=0;
		private short init_state=STATE_UNDEFINED;
		
		@Override
		void initialize(final Properties props) {
			super.initialize(props);
			final String key1 = this.getName()+".min";
			final String key2 = this.getName()+".max";
			final String key3 = this.getName()+".value";
			if(props.containsKey(key1) && props.containsKey(key2)&& props.containsKey(key3)) {
				this.init_state = STATE_SET_BY_USER;
				this.minValue = Double.parseDouble(props.getProperty(key1));
				this.maxValue = Double.parseDouble(props.getProperty(key2));
				this.startValue = OptionalDouble.of(Double.parseDouble(props.getProperty(key3)));
				}
			}
		@Override
		void fill(Trait trait, Properties props) {
			props.setProperty(this.getName()+".min", String.valueOf(getMin()));
			props.setProperty(this.getName()+".max", String.valueOf(getMax()));
			props.setProperty(this.getName()+".value", String.valueOf(MyTrait.class.cast(trait).v));
			}
		@Override
		void initialize(final VariantContext v) {
			final OptionalDouble opt=getValueForVariant(v);
			if(opt.isPresent()) initialize(opt.getAsDouble());
			}
		
		protected void initialize(final double v) {
			if(this.init_state==STATE_SET_BY_USER) {
				//nothing
				}
			else if(this.init_state==STATE_UNDEFINED) {
				minValue=v;	
				maxValue=v;	
				this.init_state = STATE_FILL;
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
			return  startValue.isPresent()?startValue.getAsDouble():makeNumber();
			}

		Trait mute(final Trait t) {
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
		boolean isSame(final Trait a, final Trait b) {
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

		protected abstract CompareOp getCompareOp();
		protected abstract OptionalDouble getValueForVariant(final VariantContext ctx);
		
		
		@Override
		boolean test(VariantContext ctx, Trait t) {
		   final OptionalDouble opt = getValueForVariant(ctx);
		   if(!opt.isPresent()) return true;
		   return getCompareOp().test(opt.getAsDouble(),MyTrait.class.cast(t).v);
		   }
		
		@Override
		final String asString(Trait t) {
			return getName() + " " + getCompareOp().getSymbol()+" "+ MyTrait.class.cast(t).v;
			}
		
		@Override
		public String toString() {
			return getClass().getSimpleName()+" min:"+getMin()+" max:"+getMax()+" enabled:"+isEnabled();
			}
		}
	
	
	private class AFTagFactory extends AbstractNumberFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.LowerOrEqual;
			}
		@Override
		protected double getMax() {
			return 0.01;
			}
	   @Override
	   protected OptionalDouble getValueForVariant(VariantContext ctx) {
			double an=0.0;
			double ac=0.0;
			for(Genotype gt :ctx.getGenotypes()) {
				if(!cases.contains(gt.getSampleName()) && !controls.contains(gt.getSampleName())) continue;
				for(Allele a: gt.getAlleles()) {
					if(a.isNoCall()) continue;
					an++;
					if(!a.isReference()) ac++;
					}
				}
			return OptionalDouble.of( ac/an );
	        }
	   
		   
		@Override
		String getName() {
			return VCFConstants.ALLELE_FREQUENCY_KEY;
			}
		}
	
	
	private abstract class  SimpleDoubleTagFactory  extends AbstractNumberFactory {
		@Override
		protected  OptionalDouble getValueForVariant(VariantContext ctx) {
			if(!ctx.hasAttribute(getTag())) return OptionalDouble.empty();
			return OptionalDouble.of(ctx.getAttributeAsDouble(getTag(), 0.0));
			}
		abstract String getTag();
		
		@Override
		final String getName() {
			return getTag();
			}
		}
	/** QD Quality By Depth **/
	private class  QDFactory  extends SimpleDoubleTagFactory {
		@Override
		protected double getMin() {
			return 0.0;
			}
		@Override
		protected double getMax() {
			return 3.0;
			}
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.GreaterOrEqual;
			}
		
		@Override
		String getTag() {
			return GATKConstants.QD_KEY;
			}
		}
	
	/** FS Fisher Strand **/
	private class  FSFactory  extends SimpleDoubleTagFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.LowerOrEqual;
			}
		@Override
		protected double getMin() {
			return 50.0;
			}
		@Override
		String getTag() {
			return GATKConstants.FS_KEY;
			}
		}
	/** MQ **/
	private class RMSMappingQualityFactory  extends SimpleDoubleTagFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.GreaterOrEqual;
			}
		@Override
		final String getTag() {
			return GATKConstants.MQ_KEY;
			}
		}
	/** MQRankSum **/
	private class MQRankSumFactory  extends SimpleDoubleTagFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.GreaterOrEqual;
			}
		
		@Override
		protected double getMax() {
			return -3.0;
			}
		
		@Override
		final String getTag() {
			return GATKConstants.MQRankSum_KEY;
			}
		}
	
	/** ReadPosRankSumFactory **/
	private class ReadPosRankSumFactory  extends SimpleDoubleTagFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.GreaterOrEqual;
			}
		@Override
		protected double getMax() {
			return -3;
			}
		
		@Override
		final String getTag() {
			return GATKConstants.ReadPosRankSum_KEY;
			}
		}

	/** Strand Odd Ratio **/
	private class SORFactory  extends SimpleDoubleTagFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.LowerOrEqual;
			}
		@Override
		protected double getMin() {
			return 3;
			}
		@Override
		final String getTag() {
			return GATKConstants.SOR_KEY;
			}
		}
	
	/** MISSING ./. */
	private class MissingFactory extends AbstractNumberFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.LowerOrEqual;
			}
		@Override
		protected double getMax() {
			return 0.1;
			}
		@Override
		protected double makeFirstNumber() {
			return 0.01;
			}
		@Override
		protected OptionalDouble getValueForVariant(VariantContext ctx) {
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
			return OptionalDouble.of( f_missing);
			}
		
		@Override
		String getName() {
			return "F_MISSING";
			}
		}
	
	/** horizontal fisher */
	private class FisherH extends AbstractNumberFactory {
		@Override
		void initialize(VariantContext ctx) {
			}
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.GreaterOrEqual;
			}
		@Override
		protected double getMin() {
			return 0;
			}
		@Override
		protected double getMax() {
			return 0.1;
			}
		
		protected OptionalDouble getValueForVariant(VariantContext ctx) {
			int cases_alt = 0;
			int cases_ref = 0;
			int ctrl_alt = 0;
			int ctrl_ref = 0;
			for(Genotype g:ctx.getGenotypes()) {
				if(g.isNoCall()) continue;
				final String sn = g.getSampleName();
				if(cases.contains(sn)) {
					if(g.hasAltAllele()) {
						cases_alt++;
						}
					else
						{
						cases_ref++;
						}
					}
				else if(controls.contains(sn)) {
					if(g.hasAltAllele()) {
						ctrl_alt++;
						}
					else
						{
						ctrl_ref++;
						}
					}
				}
			return OptionalDouble.of( FisherExactTest.compute(
					cases_alt, cases_ref,
					ctrl_alt, ctrl_ref
					).getAsDouble());
			}
		
		
		@Override
		java.lang.String getName() {
			return "FISHERH";
			}
		}
	
	
	/** Singleton AD */
	private class SingletonADFactory extends AbstractNumberFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.GreaterOrEqual;
			}
		
		@Override
		protected double getMin() {
			return 0.25;
			}
		@Override
		protected double getMax() {
			return 0.5;
			}
		@Override
		protected OptionalDouble getValueForVariant(VariantContext ctx) {
			final Genotype singleton = getSingleton(ctx);
			if(singleton==null) return OptionalDouble.empty();
			if(!singleton.isHet()) return OptionalDouble.empty();
			if(!singleton.hasAD()) return OptionalDouble.empty();
			final int[] ad = singleton.getAD();
			if(ad.length!=2) return OptionalDouble.empty();
			final double t = ad[0]+ad[1];
			if(t==0) return OptionalDouble.empty();
			double f=ad[1]/t;
			if(f>0.5) f=1.0-f;
			return OptionalDouble.of(f);
			}
		@Override
		String getName() {
			return "SingletonADRatio";
			}
		}

	/** Singleton GQ */
	private class SingletonGQFactory extends AbstractNumberFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.GreaterOrEqual;
			}
		@Override
		protected double getMin() {
			return 70.0;
			}
		@Override
		protected double makeFirstNumber() {
			return 80.0;
			}
		
		@Override
		protected OptionalDouble getValueForVariant(VariantContext ctx) {
			final Genotype singleton = getSingleton(ctx);
			if(singleton==null) return OptionalDouble.empty();
			if(!singleton.hasGQ()) return OptionalDouble.empty();
			return OptionalDouble.of(singleton.getGQ());
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
		
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.LowerOrEqual;
			}
		
		@Override
		protected OptionalDouble getValueForVariant(VariantContext ctx) {
			return OptionalDouble.of( ctx.getAttributeAsDoubleList(tag, 0.0).stream().
					mapToDouble(V->V.doubleValue()).min().orElse(0.0));
			}
		
		@Override
		protected double makeFirstNumber() {
			return 0.01;
			}
		
		@Override
		String getName() {
			return this.tag;
			}
		}

	private class CaddPhredFactory extends AbstractNumberFactory {
		final String tag = "CADD_PHRED";
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.GreaterOrEqual;
			}
		@Override
		protected OptionalDouble getValueForVariant(VariantContext ctx) {
			return ctx.getAttributeAsStringList(tag,"0").
					stream().
					filter(S->StringUtils.isDouble(S)).
					mapToDouble(S->Double.parseDouble(S)).
					filter(V->V>=0).// unknown value are -999
					min();
			}
		
		@Override
		String getName() {
			return this.tag;
			}
		}
	
	private class PolyXFactory extends AbstractNumberFactory {
		@Override
		protected CompareOp getCompareOp() {
			return OptimizeFisher.LowerOrEqual;
			}
		
		@Override
		protected OptionalDouble getValueForVariant(VariantContext ctx) {
			if(!ctx.hasAttribute(getName())) return OptionalDouble.empty();
			return OptionalDouble.of(ctx.getAttributeAsInt(getName(),0));
			}
		
		@Override
		String getName() {
			return "POLYX";
			}
		}
	
	private abstract class AbstractPredictionFactory  extends TraitFactory {
		protected final Set<String> all_acns = new HashSet<>();
		private String first_create_acns = "missense_variant,stop_gained,stop_lost,start_lost,frameshift_variant";
		protected class MyTrait extends AbstractTrait {
			final Set<String> acns;
			MyTrait(Set<String> acns) {
				this.acns = acns;
				}
			
			}
		
		@Override
		Trait createFirst() {
			final String[] acns = CharSplitter.COMMA.split(first_create_acns);
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
		void initialize(Properties props) {
			final String key=getTag()+".acns";
			if(!props.containsKey(key)) return;
			this.first_create_acns = props.getProperty(key);
			}
		
		@Override
		void fill(Trait trait, Properties props) {
			props.setProperty(getTag()+".acns", String.join(",",MyTrait.class.cast(trait).acns));
			}
		
		@Override
		boolean isSame(final Trait a, final Trait b) {
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
		final long id = (++OptimizeFisher.ID_GENERATOR);
		final long generation = OptimizeFisher.this.n_generations;
		private OptionalDouble optPvalue=OptionalDouble.empty();
		final List<Trait> traits;
		final List<VariantContext> subset = new Vector<>();
		SlidingWindow sliding_window = new SlidingWindow();
		
		Solution() {
			this(false);
			}
		Solution(boolean is_first) {
			this.traits = new ArrayList<>(OptimizeFisher.this.traitFactories.size());
			for(TraitFactory t: OptimizeFisher.this.traitFactories) {
				traits.add(is_first?t.createFirst():t.create());
				}
			}
		
		double getPValue() {
			if(!optPvalue.isPresent()) {
				throw new IllegalStateException("Pvalue was not set");
				}
			return optPvalue.getAsDouble();
			}
		
		int size() {
			return traitFactories.size();
			}
		
		final void save() throws IOException {
				final Path outVcfName = OptimizeFisher.this.outputDir.resolve("best.vcf.gz");
				final VCFHeader hdr = new VCFHeader(OptimizeFisher.this.vcfHeader);
				hdr.addMetaDataLine(new VCFHeaderLine("PVALUE", String.valueOf(getPValue())));
				try( VariantContextWriter vw=  new VariantContextWriterBuilder().
						setOutputPath(outVcfName).
						setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).
						setOption(Options.INDEX_ON_THE_FLY).
						build()) {
					vw.writeHeader(hdr);
					for(VariantContext ctx:this.subset) {
						vw.add(ctx);
						}
					}
				final Properties props = new Properties();
				this.sliding_window.fill(props);

				
				for(Trait t:this.traits) {
					t.getFactory().fill(t, props);
					if(props.containsKey("")) throw new IllegalArgumentException("bug???"+t);
					}
				
				
				final Path propsName = OptimizeFisher.this.outputDir.resolve("best.properties");
				try(BufferedWriter w = Files.newBufferedWriter(propsName,StandardOpenOption.WRITE,StandardOpenOption.TRUNCATE_EXISTING)) {
					props.store(w, "pvalue="+this.getPValue()+" generation="+this.generation);
					w.flush();
					}
			
				final Path tsvName = OptimizeFisher.this.outputDir.resolve("best.tsv");
				try(BufferedWriter w = Files.newBufferedWriter(tsvName,StandardOpenOption.WRITE,StandardOpenOption.APPEND)) {
					w.write(""+System.currentTimeMillis()+"\t"+StringUtils.now()+"\t"+this.id+"\t"+this.generation+"\t"+this.getPValue()+"\t"+getInterval()+"\t"+this.subset.size()+"\t"+
								this.sliding_window.window_size+"\t"+this.sliding_window.window_shift);
					for(int i=0;i< this.traits.size();i++) {
						w.write("\t"+ this.traits.get(i));
						}
					w.write("\n");
					w.flush();
					}
				}
	

		
		void run()  {
			LOG.info("running "+this.toString());
			this.optPvalue = OptionalDouble.empty();
			double pvalue= 1.0;
			this.subset.clear();
			final List<VariantContext> L1 = new Vector<>(OptimizeFisher.this.variants.size());
			final Set<String> contigs = new HashSet<>();
			for(VariantContext ctx: OptimizeFisher.this.variants) {
				boolean keep=true;
				for(Trait t:this.traits) {
					if(!t.getFactory().test(ctx, t)) {
						keep = false;
						break;
						}
					}
				if(!keep) continue;
				L1.add(ctx);
				contigs.add(ctx.getContig());
				}
			if(L1.isEmpty()) {
				LOG.warn("no variant for "+this.toString());
			}
			
			for(final String contig: contigs) {
				final List<VariantContext> L2 = new ArrayList<>(L1);
				L2.removeIf(V->!V.getContig().equals(contig));
				
				int start = L2.get(0).getStart();
				for(;;) {
					final List<VariantContext> L3 = new ArrayList<>();
					for(VariantContext ctx: L2) {
						if(ctx.getStart()< start) continue;
						if(!L3.isEmpty() && CoordMath.getLength(L3.get(0).getStart(), ctx.getStart())> this.sliding_window.window_size) {
							break;
							}
						L3.add(ctx);
						}
					if(L3.isEmpty())
						{
						break;
						}
					final int[]  counts = getFisherCount(L3);
					final double pv = FisherExactTest.compute(counts).getAsDouble();
					if(pv < pvalue) {
						pvalue = pv;
						this.subset.clear();
						this.subset.addAll(L3);
						}
					
					start+= this.sliding_window.window_shift;
					}
				}
					
			LOG.info("value="+pvalue);
			this.optPvalue = OptionalDouble.of(pvalue);
			}
			
		int[] getFisherCount(final List<VariantContext> variants) {
			int case_alt = 0;
			int case_ref = 0;
			int ctrl_alt = 0;
			int ctrl_ref = 0;
			for(String cas : OptimizeFisher.this.cases) {
				boolean has_alt =  variants.stream().map(G->G.getGenotype(cas)).anyMatch(G->G.hasAltAllele());
				if(has_alt) {
					case_alt++;
					}
				else
					{
					case_ref++;
					}
				}
			
			for(String ctrl : OptimizeFisher.this.controls) {
				boolean has_alt = variants.stream().map(G->G.getGenotype(ctrl)).anyMatch(G->G.hasAltAllele());
				if(has_alt) {
					ctrl_alt++;
					}
				else
					{
					ctrl_ref++;
					}
				}
			return new int[] {
					case_alt,case_ref,
					ctrl_alt,ctrl_ref
				};
			}
		
		boolean isSame(final Solution other) {
			if(this.sliding_window.window_size!=other.sliding_window.window_size) return false;
			if(this.sliding_window.window_shift!=other.sliding_window.window_shift) return false;
			for(int i=0;i< size();i++) {
				if(!traits.get(i).isSame(other.traits.get(i))) return false;
				}
			return true;
			}
		
		String getInterval() {
			return this.subset.isEmpty()?
					".":
					this.subset.get(0).getContig()+":"+this.subset.get(0).getStart()+"-"+this.subset.get(this.subset.size()-1).getStart();
			}
		@Override
		public String toString() {
			return "G"+this.generation+" ID"+id+" P="+this.optPvalue+
					" traits:["+ traits.stream().map(T->T.toString()).collect(Collectors.joining(" && ")) +
					"] interval:"+getInterval()+ " N="+ this.subset.size()+" "+
					this.sliding_window;
			}
		}
	
	
	private Solution mate(final Solution s1,final Solution s2)  {
		final Solution c = new Solution();
		c.sliding_window = random.nextBoolean()?s1.sliding_window:s2.sliding_window;
		if(random.nextDouble() < proba_mutation) {
			c.sliding_window = c.sliding_window.mute();
			}
		
		
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

	private int rnd(int m,int M) {
		return m + this.random.nextInt(M-m);
	}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Path vcfPath = Paths.get(oneAndOnlyOneFile(args));
			IOUtil.assertFileIsReadable(vcfPath);
			IOUtil.assertDirectoryIsWritable(this.outputDir);
			
			
			/* get cases */
			try(BufferedReader br = IOUtils.openPathForBufferedReading(this.casesSamplePath)) {
				this.cases.addAll( br.lines().
						filter(S->!StringUtils.isBlank(S) || S.startsWith("#")).
						collect(Collectors.toSet())
						);
				LOG.info("cases N="+this.cases.size());
				
				}
			/* get controls */
			try(BufferedReader br = IOUtils.openPathForBufferedReading(this.controlsSamplePath)) {
				this.controls.addAll( br.lines().
						filter(S->!StringUtils.isBlank(S) || S.startsWith("#")).
						collect(Collectors.toSet())
						);
				LOG.info("controls N="+this.controls.size());
				
				}
			
			try(VCFReader r = VCFReaderFactory.makeDefault().open(vcfPath,false))  {
				this.vcfHeader = r.getHeader();
				this.cases.removeIf(S->!this.vcfHeader.getSampleNameToOffset().containsKey(S));
				this.controls.removeIf(S->!this.vcfHeader.getSampleNameToOffset().containsKey(S));
				
				try(CloseableIterator<VariantContext> iter=r.iterator()) {
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						if(ctx.getAlleles().size()!=2) {
							LOG.error("Use bcftools norm. Expected only two alleles but then I got "+ctx);
							return -1;
							}
						if(ctx.getAlleles().contains(Allele.SPAN_DEL)) {
							continue;
							}
						if(ctx.getGenotypes().stream().
								filter(G->cases.contains(G.getSampleName()) || controls.contains(G.getSampleName())).
								noneMatch(G->G.hasAltAllele())) {
							continue;
							}
						this.variants.add(ctx);
						}
					}
				Collections.sort(this.variants, this.vcfHeader.getVCFRecordComparator());
				}
			LOG.info("N-variants:="+this.variants.size());

			if(this.cases.isEmpty()) {
				LOG.error("No case.");
				return -1;
				}
			if(this.controls.isEmpty()) {
				LOG.error("No control.");
				return -1;
				}
			
			if(cases.stream().anyMatch(S->controls.contains(S))) {
				LOG.error("Common samples between cases and controls.");
				return -1;
				}
			
			if(this.variants.isEmpty()) {
				LOG.error("No variant was found.");
				return -1;
				}
				
			this.traitFactories.addAll(
				Arrays.asList(
					new AFTagFactory(),
					new QDFactory(),
					new VepPredictionFactory(),
					new GnomadAFFactory(),
					new FSFactory(),
					new RMSMappingQualityFactory(),
					new MQRankSumFactory(),
					new ReadPosRankSumFactory(),
					new MissingFactory(),
					new SORFactory(),
					new SingletonADFactory(),
					new SingletonGQFactory(),
					new PolyXFactory(),
					new CaddPhredFactory(),
					new FisherH()
					));
			
			final Properties properties = new Properties();
			if(this.propertyFile!=null) {
				try(Reader r = IOUtils.openPathForBufferedReading(this.propertyFile)) {
					properties.load(r);
					}
				}
			
			final Set<String> disable_names = Arrays.stream(disableFactoriesNames.split("[, ;]+")).filter(S->!StringUtil.isBlank(S)).map(S->S.toLowerCase()).collect(Collectors.toSet());
			LOG.info("the following factories will be disabled: "+this.traitFactories.stream().
					filter(F->disable_names.contains(F.getName().toLowerCase())).
					map(F->F.getName()).collect(Collectors.joining(",")));
			this.traitFactories.removeIf(F->disable_names.contains(F.getName().toLowerCase()));
			
			
			LOG.info("initialize factories");
			for(TraitFactory trait:traitFactories) {
				trait.initialize(vcfHeader);
				trait.initialize(properties);
				}
			for(final VariantContext ctx:this.variants) {
				//initialize with current variant
				for(TraitFactory traitF:traitFactories) {
					traitF.initialize(ctx);
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
			


			final long stop = System.currentTimeMillis() + new DurationParser().convert(this.durationStr).toMillis();
			
			//write header
			try(BufferedWriter w = Files.newBufferedWriter( OptimizeFisher.this.outputDir.resolve("best.tsv"),StandardOpenOption.WRITE,StandardOpenOption.CREATE)) {
				w.append("TIMESTAMP\tDATE\tID\tGENERATION\tPVALUE\tINTERVAL\tN-VARIANTS\tWIN_SIZE\tWIN_SHIFT");
				for(TraitFactory tf:this.traitFactories) w.append("\t"+tf.getName());
				w.append("\n");
				w.flush();
				}

			
			
			Solution best = null;
			
			for(;;) {
				LOG.warning("looking for a first element in the curren config...");
				best = new Solution(true /* invoke createFirst() */);
				if(this.propertyFile!=null) {
					best.sliding_window.initialize(properties);
					}
				best.run();
				if(best.getPValue()>=1.0) {
					continue;
					}
				best.save();
				break;
				}
			
			List<Solution> generation = new ArrayList<>();
			generation.add(best);
			while(generation.size()< this.n_samples_per_generation) {
				final Solution c= new Solution();
				c.run();
				generation.add(c);
				}
			while(System.currentTimeMillis()< stop) {
				++n_generations;
				for(Solution sol:generation) {
					if(sol.getPValue() < best.getPValue() ) {
						LOG.info("new best "+sol);
						best=sol;
						best.save();
						}
					}
				LOG.info("Generation "+n_generations+" best:"+best.getPValue()+" "+TimeUnit.MILLISECONDS.toMinutes(stop-System.currentTimeMillis())+" minutes");
				final Solution external = new Solution();
				if(generation.stream().noneMatch(S->S.isSame(external)) && !external.isSame(best)) {
					external.run();
					generation.add(external);
					}
				final List<Solution> generation2 = new ArrayList<>(generation.size()*generation.size());
				for(int x=0;x< generation.size();++x) {
					for(int y=0;y< generation.size();++y) {
						final Solution c = mate(generation.get(x),generation.get(y));
						if(generation2.stream().anyMatch(S->S.isSame(c))) continue;
						
						final Solution same = generation.stream().filter(S->S.isSame(c)).findFirst().orElse(null);
						if(same!=null) {
							c.optPvalue = same.optPvalue;
							}
						else
							{
							c.run();
							}
						if(c.getPValue()==1.0) continue;
						generation2.add(c);
						}
					}
				if(!generation2.isEmpty()) {
					Collections.sort(generation2,(A,B)->Double.compare(A.getPValue(),B.getPValue()));
					generation = generation2.subList(0, Math.min(generation2.size(),this.n_samples_per_generation));
					}
				while(generation.size()< this.n_samples_per_generation) {
					final Solution x = mate(best, best);
					generation.add(x);
					}
				}
			LOG.info("Done. Best "+best.getPValue());
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new OptimizeFisher().instanceMainWithExit(args);
	}
}
