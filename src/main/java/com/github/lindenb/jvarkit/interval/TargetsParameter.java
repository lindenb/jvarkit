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
package com.github.lindenb.jvarkit.interval;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Optional;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

public class TargetsParameter {
	@Parameter(names = "--targets",description = "Restrict to comma-separated list of regions. Can be prefixed with \"^\" to request logical complement")
	private String targets_str;
	@Parameter(names = "--targets-file",description = "Restrict to regions listed in BED 'FILE'. Can be prefixed with \"^\" to request logical complement.")
	private String targets_file_str;
	
	private SAMSequenceDictionary dictionary = null;
	public TargetsParameter() {
		}

	public TargetsParameter setDictionary(final SAMSequenceDictionary dictionary) {
		this.dictionary = dictionary;
		return this;
		}
	
	public Predicate<Locatable> makePredicate() {
		final Predicate<Locatable> accept;
		final IntervalTreeMap<Boolean> treeMap = new IntervalTreeMap<>();
		if(targets_file_str!=null && !StringUtils.isBlank(this.targets_str)) {
			throw new IllegalArgumentException("both --targets and --targets-file used");
			}
		else if(targets_file_str!=null) {
			final boolean negate;
			final Path bedFile;
			if(targets_file_str.startsWith("^")) {
				bedFile = Paths.get(targets_file_str.substring(1));
				negate=true;
				}
			else
				{
				bedFile = Paths.get(targets_file_str);
				negate=false;
				}
			try(BedLineReader br=new BedLineReader(bedFile)) {
				if(this.dictionary!=null) br.setContigNameConverter(ContigNameConverter.fromOneDictionary(this.dictionary));
				br.stream().forEach(B->treeMap.put(new Interval(B),Boolean.TRUE));
				}
			if(negate) {
				accept = V->!treeMap.containsOverlapping(V);
				}
			else
				{
				accept = V->treeMap.containsOverlapping(V);
				}
			}
		else if(!StringUtils.isBlank(this.targets_str)) {
			final boolean negate=targets_str.startsWith("^");
			
			final IntervalParser interval_parser=new IntervalParser(this.dictionary);
			for(String rgn:CharSplitter.COMMA.split(
					negate?
					this.targets_str.substring(1):
					this.targets_str
					)) {
				final Optional<SimpleInterval> opt= interval_parser.apply(rgn);
				if(!opt.isPresent()) continue;
				treeMap.put(opt.get().toInterval(),Boolean.TRUE);
				}
			if(negate) {
				accept = V->!treeMap.containsOverlapping(V);
				}
			else
				{
				accept = V->treeMap.containsOverlapping(V);
				}
			}
		else
			{
			accept = V -> true;
			}
		return accept;
		}
}
