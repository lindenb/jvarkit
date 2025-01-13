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

import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.pedigree.SamplePopulation;
import com.github.lindenb.jvarkit.util.Counter;
import com.google.gson.stream.JsonWriter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public abstract class AbstractDensityPlotAnalyzer extends AbstractAnalyzer {
	private final Map<String,Counter<Double>> pop2counter=new HashMap<>();
	private final String all_names= "ALL";
	@Override
	public void init(VCFHeader vcfHeader, SamplePopulation sample2pop) {
		super.init(vcfHeader,sample2pop);
		if(!vcfHeader.hasGenotypingData()) {
			this.pop2counter.put(all_names, new Counter<>());
			}
		else
			{
			for(SamplePopulation.Population pop: sample2pop.getPopulations()) {
				this.pop2counter.put(pop.getNameAndCount(), new Counter<>());
				}
			}
		}
	
	protected abstract List<Double> toDouble(VariantContext ctx);
	
	protected boolean acceptValue(double v) {
		return true;
		}
	
	@Override
	public void visit(final VariantContext ctx) {
		if(!acceptVariant(ctx)) return;
		final List<Double> values = toDouble(ctx);
		if(ctx.hasGenotypes()) {
			for(SamplePopulation.Population pop: super.sample2population.getPopulations()) {
				if(pop.getSampleNames().stream().
						map(SN->ctx.getGenotype(SN)).
						filter(G->acceptGenotype(G)).
						anyMatch(G->G.hasAltAllele())) {
					for(double v  : values) {
						if(!acceptValue(v)) continue;
						pop2counter.get(pop.getNameAndCount()).incr(v);
						}
					}
				}
			}
		else
			{
			for(double v  : values) {
				if(!acceptValue(v)) continue;
				pop2counter.get(all_names).incr(v);
				}
			}
		}
	
	/*
	protected void saveAsMultiqcJson(JsonWriter w) throws IOException {
		final Map<String,List<Point2D>> hash=new HashMap<>(this.pop2counter.size());
		for(String popName:this.pop2counter.keySet()) {
			final List<Point2D> points = toPoints(this.pop2counter.get(popName));
			if(points.isEmpty()) continue;
				hash.put(popName, points);
			}
		makeMultiqcCustomContent().writeLineGraph(w,hash);
		}*/
	}
