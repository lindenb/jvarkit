/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import htsjdk.variant.vcf.VCFHeader;

/** A collection of samples */
public interface SampleSet {
	/** return the samples in this SampleSet */
	public Set<Sample> getSamples();
	
	/** return the affected samples in this SampleSet */
	public default Set<Sample> getAffectedSamples() {
		return getSamples().stream().filter(S->S.isAffected()).collect(Collectors.toCollection(TreeSet::new));
		}
	/** return the unaffected samples in this SampleSet */
	public default Set<Sample> getUnaffectedSamples() {
		return getSamples().stream().filter(S->S.isUnaffected()).collect(Collectors.toCollection(TreeSet::new));
		}

	
	/** find a sample in this SampleSet . returns the sample with this id or null if not found */
	public default Sample getSampleById(final String id) {
		if(id==null) return null;
		Sample ret = null;
		for(final Sample sample: this.getSamples())
			{
			if(!id.equals(sample.getId()) ) continue;
			if(ret!=null) throw new IllegalStateException("ambigous sample id "+id+" found twice "+ret+" and "+ sample);
			ret = sample;
			}
		return ret;
		}
	/** return all the trios in this SampleSet */
	public default Set<Trio> getTrios() {
		return this. getSamples().
			stream().
			filter(S->S.hasFather() || S.hasMother()).
			map(S->new TrioImpl(S)).
			collect(Collectors.toCollection(TreeSet::new));
		}
	
	/** check all samples-ID in this SampleSet are uniq */
	public default SampleSet checkUniqIds() {
		final Set<Sample> set = this.getSamples();
		final Map<String,Sample> seen = new HashMap<>(set.size());
		for(final Sample sn:set) {
			if(seen.containsKey(sn.getId())) throw new IllegalStateException(
					"Pedigree contains two sample with the same ID: "+sn+" and "+seen.get(sn.getId()));
			seen.put(sn.getId(),sn);
			}
		return this;
		}
	/** retain samples present in Set of ID . check if samples ID are uniq */
	public default Stream<Sample> getSamplesInSet(final Set<String> sampleids) {
		return this.checkUniqIds().getSamples().stream().filter(S->sampleids.contains(S.getId()));
		}
	
	/** retain samples present in VCF header . check if samples ID are uniq */
	public default Stream<Sample> getSamplesInVcfHeader(final VCFHeader header) {
		return this.getSamplesInSet(new HashSet<>(header.getSampleNamesInOrder()));
		}

	}
