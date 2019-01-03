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
package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;

public class PedFile implements Iterable<PedFile.Sample>
	{
    private static final Logger LOG= Logger.build(PedFile.class).make();

	static final String EXTENSION=".ped";
	public enum Sex {Male,Female,Unknown};
	public enum Status {Affected,Unaffected,Unknown};
	private final Map<String,Sample> name2sample=new HashMap<>();
	
	public static final Predicate<TrioGenotype> MendelianIncompatibiltyFilter
		= new Predicate<TrioGenotype>()
		{
		@Override
		public boolean test(TrioGenotype trio) {
			return false;
			}
		};
	
    /** class used to define a trio of genotype */
    public static class TrioGenotype
    	{
    	private final Genotype children;
    	private final Genotype father;
    	private final Genotype mother;
    	private TrioGenotype(final Genotype children,final Genotype father,final Genotype mother)
    		{
    		this.children = children;
    		this.father= father;
    		this.mother= mother;
    		}
    	
    	public Genotype getChildren() {
			return children;
		}
    	
    	public Genotype getFather() {
			return father;
		}
    	public Genotype getMother() {
			return mother;
		}
    	
    	private boolean _isMendelianIncompatibility(final Genotype child,final Genotype parent)
    		{
    		if(child==null || parent==null) return false;
    		if(!child.isCalled() || !parent.isCalled()) return false;
    		if(child.getPloidy()!=2 || parent.getPloidy()!=2) return false;
    		for(final Allele childAllele:child.getAlleles())
    			{
    			if(parent.getAlleles().contains(childAllele)) return false;
    			}
    		
    		return true;
    		}
    	public boolean isMendelianIncompatibility()
			{
			if(children==null || !children.isCalled()) return false;
			if(father==null || !father.isCalled()) {
				return this._isMendelianIncompatibility(children,mother);
				}
			if(mother==null || !mother.isCalled()) {
				return this._isMendelianIncompatibility(children,father);
				}
			final Allele alleles[]=new Allele[2];
			for(final Allele af:father.getAlleles())
				{
					alleles[0]=af;
					for(final Allele am:mother.getAlleles())
					{
						alleles[1]=am;
						final Genotype sim = new GenotypeBuilder(children.getSampleName()).alleles(Arrays.asList(alleles)).make();
						if(children.sameGenotype(sim, true)) return false;
					}	
				}
			return true;
			}
    	
    	
    	}

	public Iterable<PedFile.TrioGenotype> getTrios(final VariantContext ctx)
		{
		if(ctx==null || isEmpty()) return Collections.emptyList();
		return new TrioGenotypeIterator(ctx);
		}

	
	private class TrioGenotypeIterator
		implements Iterable<TrioGenotype>
		{
		private final VariantContext ctx;
		TrioGenotypeIterator( final VariantContext ctx)
			{
			this.ctx=ctx;
			}
		private  class MyIterator extends AbstractIterator<TrioGenotype>
			{
			int index=0;
			@Override
			protected TrioGenotype advance() {
				while(this.index<ctx.getNSamples())
					{
					final Genotype gc = ctx.getGenotype(this.index);
					this.index++;
					if(gc==null) continue;
					final PedFile.Sample child = PedFile.this.get(gc.getSampleName());
					if(child==null) continue;
					final PedFile.Sample father= child.getFather();
					final Genotype gf =  (father==null?null:ctx.getGenotype(father.getName()));
					final PedFile.Sample mother= child.getMother();
					final Genotype gm =  (mother==null?null:ctx.getGenotype(mother.getName()));
					if(gf==null && gm==null) continue;
					return new TrioGenotype(gc, gf, gm);
					}
				return null;
				}
			}
		
		@Override
		public Iterator<TrioGenotype> iterator() {
			return new MyIterator();
			}
		}
	
	public interface Sample
		{
		public String getFamily();
		public String getName();
		public default boolean hasFather() { return getFather()!=null;}
		public default boolean hasMother() { return getMother()!=null;}
		public default String getFatherName() { return hasFather()?getFather().getName():null;}
		public default String getMotherName() { return hasMother()?getMother().getName():null;}
		public Sample getFather();
		public Sample getMother();
		public Sex getSex();/*RIP george Michael */
		public Status getStatus();
		public default boolean isAffected() { return getStatus()==Status.Affected;}
		public default boolean isUnaffected() { return getStatus()==Status.Unaffected;}
		public default boolean isMale() { return getSex()==Sex.Male;}
		public default boolean isFemale() { return getSex()==Sex.Female;}
		}
	
	private static final PedFile EMPTY=new PedFile();
	public static PedFile getEmptyInstance()
		{
		return EMPTY;
		}
	
	public boolean isEmpty()
		{
		return this.name2sample.isEmpty();
		}
	
	@Override
	public Iterator<Sample> iterator()
		{
		return name2sample.values().iterator();
		}
	
	public Set<String> getSampleNames()
		{
		return Collections.unmodifiableSet(this.name2sample.keySet());
		}	
	
	public boolean has(final String sample) {
		return this.name2sample.containsKey(sample);
		}
	
	public Sample get(final String sample) {
		return this.name2sample.get(sample);
		}
	
	private class SampleImpl implements PedFile.Sample
		{
		String family;
		String name;
		Sample father=null;
		Sample mother=null;
		Sex sex=Sex.Unknown;
		Status status=Status.Unknown;
		@Override
		public String getFamily() { return family; }
		@Override
		public String getName() { return name; }
		@Override
		public Sample getFather() { return father; }
		@Override
		public Sample getMother() { return mother; }
		@Override
		public Sex getSex() { return sex; }
		@Override
		public Status getStatus() { return status; }
		@Override
		public boolean equals(Object obj)
			{
			if(obj==this) return true;
			if(obj==null || !(obj instanceof SampleImpl)) return false;
			final SampleImpl other=(SampleImpl)obj;
			return this.name.equals(other.name) && 
					this.family.equals(other.family);
			}
		@Override
		public int hashCode()
			{
			return this.name.hashCode();
			}
		@Override
		public String toString()
			{
			return getName();
			}
		}
	
	
	public static PedFile load(final BufferedReader r) throws IOException
		{
		String line;
		final Pattern delim=Pattern.compile("[\t ]+");
		final PedFile ped=new PedFile();
		final Map<SampleImpl,String> sample2fatherid=new HashMap<>();
		final Map<SampleImpl,String> sample2motherid=new HashMap<>();
		while((line=r.readLine())!=null)
			{
			if(line.trim().isEmpty() || line.startsWith("#")) continue;
			final String tokens[]=delim.split(line);
			if(tokens.length<4) throw new IOException("Not enought column in "+line);
			final SampleImpl sample=ped.new SampleImpl();
			sample.family=tokens[0];
			sample.name=tokens[1];
			if(sample.name.isEmpty()) throw new IOException("Sample name cannot be empty");
			if(!(tokens[2].isEmpty() || tokens[2].equals("0"))) {
				sample2fatherid.put(sample,tokens[2]);
				}
			if(!(tokens[3].isEmpty() || tokens[3].equals("0")))
				{
				sample2motherid.put(sample,tokens[3]);
				}
			if(tokens.length>4) {
				if(tokens[4].toLowerCase().equals("M") || tokens[4].equals("1"))
					{
					sample.sex = Sex.Male;
					}
				else if(tokens[4].toLowerCase().equals("F") || tokens[4].equals("2"))
					{
					sample.sex = Sex.Female;
					}
				else
					{
					sample.sex = Sex.Unknown;
					}
				}
			if(tokens.length>5) {
				if(tokens[5].equals("1"))
					{
					sample.status = Status.Affected;
					}
				else if(tokens[5].equals("0"))
					{
					sample.status = Status.Unaffected;
					}
				else
					{
					sample.status = Status.Unknown;
					}
				}
			if(ped.name2sample.containsKey(sample.getName())) {
				throw new IOException("Duplicate sample "+sample);
				}
			ped.name2sample.put(sample.getName(),sample);
			}
		
		for(final SampleImpl c:sample2fatherid.keySet())
			{
			final String parentid = sample2fatherid.get(c);
			final SampleImpl parent = (SampleImpl)ped.name2sample.get(parentid);
			if(parent==null) continue;
			if(parent.getSex()==Sex.Unknown) {
				parent.sex=Sex.Male;
				}
			else if(parent.isFemale()) throw new IOException("father of "+c+" is female");
			c.father=parent;
			}
		
		for(final SampleImpl c:sample2motherid.keySet())
			{
			final String parentid = sample2motherid.get(c);
			final SampleImpl parent = (SampleImpl)ped.name2sample.get(parentid);
			if(parent==null) continue;
			if(parent.getSex()==Sex.Unknown) {
				parent.sex=Sex.Female;
				}
			else if(parent.isMale()) throw new IOException("mother of "+c+" is male");
			c.mother=parent;
			}
		
		return ped;
		}
	
	public static PedFile load(final File file) throws IOException
		{
		FileReader r=null;
		try {
			LOG.info("reading "+file);
			r=new FileReader(file);
			return load(new BufferedReader(r));
		} finally
		{
		CloserUtil.close(r);
		}
		
		}		
	}
