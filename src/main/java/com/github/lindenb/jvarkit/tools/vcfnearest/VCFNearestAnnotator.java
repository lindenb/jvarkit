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
package com.github.lindenb.jvarkit.tools.vcfnearest;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.infotable.VCFInfoTableModelFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VCFNearestAnnotator implements VariantAnnotator {
	private static final Logger LOG = Logger.of(VCFNearestAnnotator.class);
	private String tag = "NEARSET";
	private int max_distance = 0;
	private int limit_count = -1;
	private Path bedPath;
	private final IntervalTreeMap<List<BedLine>> loc2beds = new IntervalTreeMap<>();
	private final List<String> bedColumns = new ArrayList<>();
	
	public VCFNearestAnnotator() {
		}


	
	private void addInterval(BedLine rec) {
		final Interval r=new Interval(rec);
		List<BedLine> genes = loc2beds.get(r);
		if(genes==null) {
			genes=new ArrayList<>();
			loc2beds.put(r, genes);
			}
		genes.add(rec);
		}
	
	public void setBedPath(Path bedPath) {
		this.bedPath = bedPath;
		}
	
	public void setColumns(final List<String> bedcolumns) {
		this.bedColumns.clear();
		this.bedColumns.addAll(bedcolumns);
	}
	
	
	@Override
	public void fillHeader(final VCFHeader header) {
		final List<String> cols2 = new ArrayList<>();
		this.bedColumns.stream().
			filter(S->!StringUtils.isBlank(S)).
			forEach(S->cols2.add(S));
		cols2.add("side");
		cols2.add("distance");
		
		
		
		header.addMetaDataLine(new VCFInfoHeaderLine(
				getTag(), 
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,"Features located less than "+this.max_distance+"bp from variant. "+
						VCFInfoTableModelFactory.encode(cols2)
				));
		
		// load BED
		if(this.bedPath!=null) {
			final SAMSequenceDictionary dict= header.getSequenceDictionary();
			final UnaryOperator<String> converter = dict==null?C->C: ContigNameConverter.fromOneDictionary(dict);
			try(BedLineReader br = new BedLineReader(this.bedPath)) {
				br.setContigNameConverter(converter);
				br.stream().
					filter(B->!StringUtils.isBlank(B.getOrDefault(3,""))).
					forEach(this::addInterval);
				}
			}
	}
	
	private int distanceTo(final VariantContext ctx,final BedLine reg) {
		if(ctx.overlaps(reg)) return 0;
		else if(ctx.getEnd() < reg.getStart()) {
			return CoordMath.getLength(ctx.getEnd(), reg.getStart());
			}
		else
			{
			return CoordMath.getLength(reg.getEnd(), ctx.getStart());
			}
		}
	
	public void setTag(String tag) {
		this.tag = tag;
		}
	
	public String getTag() {
		return tag;
		}
	
	public void setMaxDistance(int max_distance) {
		this.max_distance = max_distance;
		}
	
	public void setLimitCount(int limit_count) {
		this.limit_count = limit_count;
		}
	
	@Override
	public List<VariantContext> annotate(final VariantContext ctx) throws IOException {
		Interval r = new Interval(
				ctx.getContig(), 
				Math.max(1, ctx.getStart()-max_distance),
				ctx.getEnd()+max_distance
				);
		List<BedLine> L1= this.loc2beds.getOverlapping(r).
				stream().
				flatMap(L->L.stream()).
				sorted((A,B)->Integer.compare(distanceTo(ctx,A), distanceTo(ctx,B))).
				collect(Collectors.toList());
		
		if(L1.isEmpty()) return Collections.singletonList(ctx);
		final Set<String> set=new LinkedHashSet<>();
		for(BedLine rec:L1) {
			if(this.limit_count>=0 && set.size() > this.limit_count) break;
			final int side;
			if(ctx.overlaps(rec))  {
				side= 0;
				}
			else if(rec.getEnd() < ctx.getStart() ) {
				side=-1;
				}
			else
				{
				side=1;
				}
			final List<String> row = new ArrayList<>(this.bedColumns.size()+2);
			for(int x=0; x< this.bedColumns.size();++x) {
				if(StringUtils.isBlank(this.bedColumns.get(x))) continue;
				row.add(rec.getOrDefault(x+3,""));
				}
			row.add(String.valueOf(side));
			row.add(String.valueOf(distanceTo(ctx,rec)));
			
			set.add(String.join("|", row));
			}
		return Collections.singletonList(
				new VariantContextBuilder(ctx).
					attribute(getTag(), new ArrayList<String>(set)).
					make()
				);
		}
	
	@Override
	public void close() {
		
		}
}
