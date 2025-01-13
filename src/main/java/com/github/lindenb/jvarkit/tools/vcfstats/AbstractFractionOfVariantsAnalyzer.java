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
package com.github.lindenb.jvarkit.tools.vcfstats;



import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Max;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Min;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.pedigree.SamplePopulation;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.FunctionalMap;

import htsjdk.variant.vcf.VCFHeader;

/**
 * 
 *  Analyzer counting variants per sample
 *
 */
public abstract class AbstractFractionOfVariantsAnalyzer extends AbstractAnalyzer {
	protected long n_variants = 0L;
	protected final Counter<String> counter = new Counter<>();
	
	public void init(VCFHeader vcfHeader, SamplePopulation sample2pop) {
		super.init(vcfHeader,sample2pop);
		for(final String sn: vcfHeader.getSampleNamesInOrder()) {
			counter.initializeIfNotExists(sn);
			}
		}
	
	protected Map<String,List<Number>> groupByPopulation() {
		final AutoMap<String,Number,List<Number>> content= AutoMap.makeList();
		for(String sample: super.vcfHeader.getSampleNamesInOrder()) {
			for(SamplePopulation.Population pop : sample2population.getSampleByName(sample).getPopulations()) {
				content.insert(pop.getName()+ "(N="+pop.size()+")",this.counter.count(sample)/(double)n_variants);
				}
			}
		return content;
		}
	
	
	protected void writeTsv(final Path directory,final FunctionalMap<String, String> options) throws IOException  {
		if(n_variants==0) return;
		String filename="g.tsv";
		try(PrintWriter w = IOUtils.openPathForPrintWriter(directory.resolve(filename))){
			w.println("# "+getName());
			w.println("# "+getDescription());
			w.println("sample\tcount\ttotal\tfraction");
			for(String sample: counter.keySet()) {
				w.print(sample);
				w.print("\t");
				w.print(counter.count(sample));
				w.print("\t");
				w.print(this.n_variants);
				w.print("\t");
				w.print(counter.count(sample)/(double)this.n_variants);
				w.println();
				}
			w.flush();
			}
		final Map<String,List<Number>> bypop=groupByPopulation();
		if(!bypop.isEmpty()) {
			try(PrintWriter w = IOUtils.openPathForPrintWriter(directory.resolve(filename))){
				w.println("# "+getName()+" by group");
				w.println("# "+getDescription()+" by group");
				w.println("group\tmedian\taverage\tmin\tmax");
				for(String groupName: bypop.keySet()) {
					final double[] values=bypop.get(groupName).stream().mapToDouble(V->V.doubleValue()).sorted().toArray();
					w.print(groupName);
					w.print("\t");
					w.print(new Median().evaluate(values));
					w.print("\t");
					w.print(new Mean().evaluate(values));
					w.print("\t");
					w.print(new Min().evaluate(values));
					w.print("\t");
					w.print(new Max().evaluate(values));
					w.println();
					}
				w.flush();
				}
			}
		}
	
	@Override
	public void writeReports(final Path directory,final FunctionalMap<String, String> options) throws IOException {
		writeTsv(directory,options);
		}
	}
