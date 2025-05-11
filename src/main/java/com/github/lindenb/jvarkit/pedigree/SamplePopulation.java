/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.pedigree;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.HasName;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.vcf.VCFHeader;


/**
 * A class parsing sample associated to population
 */
public class SamplePopulation {
	private static final Logger LOG = Logger.of(SamplePopulation.class);
	public static final String OPT_DESC="tab delimited file containing (sample-name)(TAB)(collection-name). Empty lines or starting with '#' are skipped";
	@Parameter(names={"--sample2collection"},description=OPT_DESC)
	private Path filePath = null;

	public interface Population extends HasName,Map<String,Sample> {
		public default Collection<Sample> getSamples() {
			return Collections.unmodifiableCollection(this.values());
			}
		public default boolean contains(final Sample sn) {
			return containsKey(sn.getName());
			}
		
		/** print name followed by number of samples */
		public default String getNameAndCount() {
			return getName()+" (N=" + StringUtils.niceInt(size())+")";
		}
		
		/** get All samples names in that collection, alias of keySet*/
		public default Set<String> getSampleNames() {
			return keySet();
			}
		/** return true if pop contains sample with this name */
		public default boolean hasSample(final String snName) {
			return containsKey(snName);
			}
		}
	public interface Sample extends HasName,Map<String,Population> {
		public default boolean contains(final Population pop) {
			return containsKey(pop.getName());
			}
		public default Collection<Population> getPopulations() {
			return Collections.unmodifiableCollection(this.values());
			}
		public default Population getPopulation() {
			switch(this.size()) {
				case 1: return entrySet().iterator().next().getValue();
				default: throw new IllegalStateException("Expected one and only one population associated to sample "+getName()+" but got "+String.join(",", keySet()));
				}
			}
		}
	
	private static abstract class AbstractNamedMap<T extends HasName> extends AbstractMap<String,T> implements HasName {
		protected final Map<String,T> items = new TreeMap<>();
		final String name;
		AbstractNamedMap(final String name) {
			this.name=name;
			}
		@Override
		public String getName() {
			return name;
			}
		@Override
		public T get(Object key) {
			return items.get(key);
			}
		@Override
		public Set<String> keySet() {
			return Collections.unmodifiableSet(this.items.keySet());
			}
	
		@Override
		public Set<Entry<String, T>> entrySet() {
			return Collections.unmodifiableSet(items.entrySet());
			}
		@Override
		public int hashCode() {
			return getName().hashCode();
			}
		@Override
		public int size() {
			return this.items.size();
			}
		
		@Override
		public String toString() {
			return getName();
			}
		}
	
	private class PopulationImpl extends AbstractNamedMap<Sample> implements Population {
		PopulationImpl(final String name) {
			super(name);
			}
		@Override
		public boolean equals(Object o) {
			if(o==this) return true;
			if(o==null || !(o instanceof PopulationImpl)) return false;
			final PopulationImpl other = PopulationImpl.class.cast(o);
			if(!other.getName().equals(this.getName())) return false;
			return true;
			}
		}
	private class SampleImpl extends AbstractNamedMap<Population> implements Sample {
		SampleImpl(final String name) {
			super(name);
			}
		@Override
		public boolean equals(Object o) {
			if(o==this) return true;
			if(o==null || !(o instanceof SampleImpl)) return false;
			final SampleImpl other = SampleImpl.class.cast(o);
			if(!other.getName().equals(this.getName())) return false;
			return true;
			}
		}
	
	private final AutoMap<String, SampleImpl, SampleImpl> id2sample =  AutoMap.make(S ->new SampleImpl(S));
	private final AutoMap<String, PopulationImpl, PopulationImpl> id2population =  AutoMap.make(S ->new PopulationImpl(S));

	
	public SamplePopulation() {
		}
	/** copy ctor */
	public SamplePopulation(final SamplePopulation cp) {
		for(final Population pop:cp.id2population.values()) {
			for(String sn:pop.keySet()) {
				this.insert(sn, pop.getName());
			}
		}
	}
	
	@Override
	protected SamplePopulation clone() {
		return new SamplePopulation(this);
		}
	
	public void insert(final String sample,final String population) {
		if(StringUtils.isBlank(sample)) throw new IllegalArgumentException("blank sample "+sample);
		if(StringUtils.isBlank(population)) throw new IllegalArgumentException("blank population "+population);
		final PopulationImpl pop = id2population.insert(population);
		final SampleImpl sn = id2sample.insert(sample);
		if(!pop.items.containsKey(sample)) {
			pop.items.put(sample,sn);
			}
		if(!sn.items.containsKey(population)) {
			sn.items.put(population,pop);
			}
		}
	
	public SamplePopulation load(final Path filePath) {
		Objects.requireNonNull(filePath, "path mapping sample to collection is undefined");
		this.filePath=filePath;
		try(BufferedReader br = IOUtils.openPathForBufferedReading(filePath)) {
			String line;
			while((line=br.readLine())!=null) {
				if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
				final String tokens[] = CharSplitter.TAB.split(line);
				if(tokens.length!=2) throw new JvarkitException.TokenErrors(2, tokens);
				insert(tokens[0], tokens[1]);
				}
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		if(id2sample.isEmpty()) LOG.warn("No valid sample was found in "+filePath);
		return this;
		}
	
	/** check that one sample has only one collection */
	public SamplePopulation assertOnePopulationPerSample() {
		for(Sample sn:getSamples()) {
			if(sn.getPopulations().size()!=1) {
				throw new IllegalArgumentException(
				"Expected one and only one population associated to sample "+
					sn.getName()+" but got :"+
					sn.getPopulations().stream().
						map(it->it.getName()).
						collect(Collectors.joining(" and ")
						)
					);
				}
			}
		return this;
		}
	
	public SamplePopulation load() {
		if(this.filePath!=null) {
			return this.load(this.filePath);
			}
		else 
			{
			LOG.warn("No collection was defined.");
			}
		return this;
		}
	
	
	public boolean isEmpty() {
		return this.id2sample.isEmpty();
	}
	
	/** retain cases / controls that are in the VCF header */
	public SamplePopulation retain(final VCFHeader header) {
		if(header==null) throw new NullPointerException("vcf header is null");
		if(!header.hasGenotypingData()) {
			LOG.warn("the vcf header contains no genotype/sample.");
			}
		return retainSamples(header.getGenotypeSamples());
		}
	
	/** retain cases / controls that are in collection */
	public SamplePopulation retainSamples(final Collection<String> keep_samples) {
		if(keep_samples==null) throw new NullPointerException("other header is null");
		final Set<String> sample_to_remove = new HashSet<>(this.id2sample.keySet());
		sample_to_remove.removeAll(keep_samples);
		
		for(final String sn:sample_to_remove) {
			for(final String popid:id2sample.get(sn).keySet()) {
				PopulationImpl pop = this.id2population.get(popid);
				pop.items.remove(sn);
				if(pop.items.isEmpty()) {
					id2population.remove(popid);
					}
				}
			id2sample.remove(sn);
			}
		return this;
		}
	
	/** get All samples */
	public Collection<? extends Sample> getSamples() {
		return id2sample.values();
		}
	/** get All samples */
	public Set<String> getSampleNames() {
		return Collections.unmodifiableSet(this.id2sample.keySet());
		}
	
	public boolean hasSample(final String name) {
		return this.id2sample.containsKey(name);
		}
	
	public boolean hasCollection(final String name) {
		return this.id2population.containsKey(name);
		}
	

	public Sample getSampleByName(final String name) {
		return this.id2sample.get(name);
		}
	/** shortcut to getSampleByName(g.getSampleName()) */
	public Sample getSampleByGenotype(final Genotype g) {
		return getSampleByName(g.getSampleName());
		}

	/** get All population */
	public Collection<? extends Population> getPopulations() {
		return id2population.values();
		}
	/** get All samples */
	public Set<String> getPopulationNames() {
		return Collections.unmodifiableSet(this.id2population.keySet());
		}
	
	public Population getPopulationByName(final String name) {
		return this.id2population.get(name);
		}
	
	/** alias of assertOnePopulationPerSample */
	public SamplePopulation checkOnePopulationPerSample() {
		return assertOnePopulationPerSample();
		}

	}
