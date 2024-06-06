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
package com.github.lindenb.jvarkit.pedigree;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.vcf.VCFHeader;

public class SampleToGroup {
	public static final String OPT_DESC="a tab delimited file with two columns: (sample-name)(TAB)(group-name)";
	
	@SuppressWarnings("serial")
	protected class Sample extends TreeSet<String> {
		final String sample;

		Sample(final String sample) {
			this.sample = sample;
			}
		@Override
		public int hashCode() {
			return sample.hashCode();
			}
		@Override
		public String toString() {
			return sample;
			}
		}
	@SuppressWarnings("serial")
	protected class Group extends TreeSet<String> {
		final String group;
		Group(final String group) {
			this.group=group;
			}
		@Override
		public int hashCode() {
			return group.hashCode();
			}
		@Override
		public String toString() {
			return group;
			}
		}
	
	private final Map<String,Group> groups = new TreeMap<>();
	private final Map<String,Sample> samples = new TreeMap<>();
	
	public boolean hasSample(final String sn) {
		return this.samples.containsKey(sn);
	}
	
	public boolean hasGroup(final String gr) {
		return this.groups.containsKey(gr);
	}

	
	public int getSamplesCount() {
		return this.samples.size();
		}
	public Set<String> getSamples() {
		return Collections.unmodifiableSet(this.samples.keySet());
		}
	public int getGroupsCount() {
		return this.groups.size();
		}
	
	public Set<String> getGroups() {
		return Collections.unmodifiableSet(this.groups.keySet());
		}
	
	public Set<String> getSamplesForGroup(final String groupName) {
		if(!hasGroup(groupName)) throw new IllegalArgumentException("no such sample "+groupName);
		return Collections.unmodifiableSet(this.groups.get(groupName));
		}
	public Set<String> getGroupsForSample(final String sampleName) {
		if(!hasSample(sampleName)) throw new IllegalArgumentException("no such sample "+sampleName);
		return Collections.unmodifiableSet(this.samples.get(sampleName));
		}
	public String getGroupForSample(final String sampleName) {
		final Set<String> groups = getGroupsForSample(sampleName);
		if(groups.size()!=1) throw new IllegalArgumentException("more than one group for sample "+sampleName+" "+groups);
		return groups.iterator().next();
		}
	/** create a default group for those samples if they don't have an affected group */
	public SampleToGroup complete(String groupName,Collection<String> samples) {
		final Set<String> set =new HashSet<>(samples);
		set.removeAll(this.getSamples());
		if(!set.isEmpty()) {
			for(String sn:set) {
				putSampleGroup(sn, groupName);
				}
			}
		return this;
		}
	
	public SampleToGroup assertOneGroupPerSample() {
		Sample opt = this.samples.
				values().
				stream().
				filter(KV->KV.size()!=1).
				findAny().
				orElse(null);
		if(opt!=null) throw new IllegalArgumentException("sample "+opt.sample+" has more than one group "+String.join(" ",opt));
		return this;
		}
	
	
	public SampleToGroup putSampleGroup(final String sample,final String group) {
		Group g=groups.get(group);
		if(g==null) {
			g = new Group(group);
			groups.put(group,g);
			}
		g.add(sample);
		
		Sample s = samples.get(sample);
		if(s==null) {
			s=new Sample(sample);
			samples.put(sample, s);
			}
		s.add(group);
		return this;
		}
	private SampleToGroup cleanup() {
		this.samples.keySet().removeIf(S->samples.get(S).isEmpty());
		this.groups.keySet().removeIf(S->groups.get(S).isEmpty());
		return this;
		}
	public SampleToGroup retainSamples(final VCFHeader h) {
		return this.retainSamples(h.getSampleNamesInOrder());
		}
	public SampleToGroup retainSamples(final Collection<String> samples) {
		final Set<String> toRemove=new HashSet<>(this.samples.keySet());
		toRemove.removeAll(samples);
		return removeSamples(toRemove);
		}
	
	public SampleToGroup removeSamples(final Collection<String> samples_to_be_removed) {
		for(String sn:samples_to_be_removed) {
			this.samples.remove(sn);
			for(final String g: this.groups.keySet()) {
				 this.groups.get(g).remove(sn);
				}
			}
		return cleanup();
		}
	public boolean isEmpty() {
		return this.samples.isEmpty();
		}
	public SampleToGroup load(final Path path) {
		
		try(BufferedReader br=IOUtils.openPathForBufferedReading(path)) {
			String line;
			while((line=br.readLine())!=null) {
				if(StringUtils.isBlank(line)) continue;
				final String[] tokens= CharSplitter.TAB.split(line);
				if(tokens.length!=2) throw new JvarkitException.TokenErrors(2, tokens);
				if(StringUtils.isBlank(tokens[0])) throw new IllegalArgumentException("Empty sample in file "+path);
				if(StringUtils.isBlank(tokens[1])) throw new IllegalArgumentException("Empty group in file "+path);
				this.putSampleGroup(tokens[0], tokens[1]);
				}
			return this;
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	}
