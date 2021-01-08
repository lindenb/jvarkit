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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.bedtools;

import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.iterator.LineIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;

/**
BEGIN_DOC

input is bed on standard input or it's a set of interval files (.bed, .interval_list, .gtf, etc... )

output is a BED file:

* 4th column: number of items merged
* 5th column: fraction overlap with previous record
* 6th column an optional label. Concatenation of the distinct label. For a bed record the label is the 4th column. For a VCF record, the label is the name of the samples

## Example

```
$ java -jar dist/bedmergecnv.jar src/test/resources/manta.B00*.vcf.gz | more

chr21	9653162	9653235	1	0.0
chr21	9653169	9653243	1	0.89
chr21	9653357	9653502	1	0.00
chr21	9660834	9660835	3	0.00
chr21	9841718	9842058	2	0.00
chr21	9846433	9846531	1	0.00
chr21	9846456	9846557	1	0.74
chr21	9846481	9846580	1	0.75
chr21	9851247	9854008	1	0.00
chr21	9862002	9862003	1	0.00
chr21	9881949	9882060	1	0.00
chr21	9936419	9936560	1	0.00
chr21	10459126	10459235	3	0.00
chr21	10475193	10475513	1	0.00
chr21	10492877	10493042	3	0.00
chr21	10493815	10493874	1	0.00
chr21	10618665	10618720	1	0.00
chr21	10633657	10633709	1	0.00
chr21	10699831	10700200	2	0.00
(...)
```


END_DOC

 */
@Program(
		name="bedmergecnv",
		description="Merge Bed records if they overlap a fraction of their lengths.",
		keywords={"bed","chromosome","contig"},
		creationDate="20200330",
		modificationDate="20200603"
		)
public class BedMergeCnv
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BedMergeCnv.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile= null;
	@Parameter(names={"-f","--fraction","--overlap"},description="Intervals will be merged if they both overlap this fraction of their lengths. "+FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double fraction = 0.9;

	
	private final IntervalTreeMap<Cluster> clusters = new IntervalTreeMap<>();

	
	private class Cluster implements Locatable{
		private final List<Locatable> intervals = new ArrayList<>();
		Cluster(final Locatable	base) {
			this.intervals.add(base);
			}
		@Override
		public String getContig() {
			return this.intervals.get(0).getContig();
			}
		@Override
		public int getStart() {
			return this.intervals.stream().mapToInt(R->R.getStart()).min().orElse(-1);
			}
		@Override
		public int getEnd() {
			return this.intervals.stream().mapToInt(R->R.getEnd()).max().orElse(-1);
			}
		
		public String getLabel() {
			final String s =  this.intervals.stream().
					map(L->labelFor(L)).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toSet()).
					stream().
					collect(Collectors.joining(","));
			return StringUtils.isBlank(s)?".":s;
			}
		
		private boolean test(final Locatable R1,final Locatable R2) {
			return fractionOverlap(R1,R2) >= BedMergeCnv.this.fraction;
			}

		public boolean add(final Locatable si) {
			if(!this.intervals.stream().allMatch(R->test(R,si))) return false;
			this.intervals.add(si);
			return true;
			}
		}
	
	private String labelFor(final Locatable loc) {
		if(loc instanceof BedLine) {
			final BedLine bed = BedLine.class.cast(loc);
			if(bed.getColumnCount()>3) return bed.get(3);
			}
		if(loc instanceof VariantContext) {
			final VariantContext ctx = VariantContext.class.cast(loc);
			if(ctx.hasGenotypes()) return ctx.getGenotypes().stream().map(G->G.getSampleName()).collect(Collectors.toSet()).stream().collect(Collectors.joining(","));
			}
		return null;
		}
	
	private double fractionOverlap(final Locatable R1,final Locatable R2) {
		if(!R1.overlaps(R2)) return 0.0;
		final int x1 = Math.max(R1.getStart(),R2.getStart());
		final int x2 = Math.min(R1.getEnd(),R2.getEnd());
		final double L = CoordMath.getLength(x1,x2);
		final double L1 = R1.getLengthOnReference();
		final double L2 = R2.getLengthOnReference();
		return Math.min(L/L1 , L/L2);
		}
	private void parseLocatables(final Iterator<? extends Locatable> iter) {
		while(iter.hasNext()) {
			final Locatable rec = iter.next();
			boolean found=false;
			final Interval key = new Interval(rec);
			for(final Cluster bc: this.clusters.getOverlapping(key)) {
				if(bc.add(rec)) { found=true;}
				}
			if(!found) {
				this.clusters.put(key,new Cluster(rec));
				}
			}
		}
	private void parseBed(final Iterator<String> iter) {
		final BedLineCodec codec = new BedLineCodec();
		final Iterator<Locatable> iter2 = new AbstractIterator<Locatable>() {
			@Override
			protected Locatable advance() {
				while(iter.hasNext()) {
					final String line = iter.next();
					if( BedLine.isBedHeader(line) ) continue;
					final BedLine rec = codec.decode(line);
					if(rec==null) continue;
					return rec;
					}
				return null;
				}
			};
		parseLocatables(iter2); 
		}
		
	@Override
	public int doWork(final List<String> args) {
		if(this.fraction<=0) {
			LOG.error("Bad fraction : "+this.fraction);
			return -1;
		}
		PrintWriter w=null;
		try
			{
			if(args.isEmpty()) {
				parseBed(new LineIterator(new InputStreamReader(stdin())));
			} else {
				for(final String fname:args) {
					final IntervalListProvider ivp = IntervalListProvider.from(fname);
					final Iterator<? extends Locatable> iter = ivp.stream().iterator();
					parseLocatables(iter);
				}
			}
			
			w = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			final List<Cluster> ordered = new ArrayList<>(this.clusters.values());
			Collections.sort(ordered,(A,B)->{
				int i= A.getContig().compareTo(B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(),B.getStart());
				if(i!=0) return i;
				return Integer.compare(A.getEnd(),B.getEnd());
			});
			
			Cluster prev = null;
			for(final Cluster bc: ordered) {
				w.print(bc.getContig());
				w.print("\t");
				w.print(bc.getStart()-1);
				w.print("\t");
				w.print(bc.getEnd());
				w.print("\t");
				w.print(bc.intervals.size());
				w.print("\t");
				w.print(prev!=null && prev.contigsMatch(bc)?String.format("%.2f",this.fractionOverlap(prev,bc)):"0.0");
				w.print("\t");
				w.print(bc.getLabel());
				w.println();
				prev = bc;
				}
			w.flush();
			w.close();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			}
		}
	

	public static void main(final String[] args)
		{
		new BedMergeCnv().instanceMainWithExit(args);
		}
	}
