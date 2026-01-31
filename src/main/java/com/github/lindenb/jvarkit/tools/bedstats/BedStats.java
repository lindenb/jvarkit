/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bedstats;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Path;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.PrefixSuffixWriter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.log.Logger;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

/**
BEGIN_DOC

input is bed on standard input or it's a set of interval files (.bed, .interval_list, .gtf, etc... )

output is text file containing multiple chuncks for multiqc


## Example

```
$ java -jar dist/bedstats.jar < in1.bed > out.txt
$ java -jar dist/bedstats.jar in1.bed in2.bed in3.bed > out.txt
```


END_DOC

 */
@Program(
		name="bedstats",
		description="statistics about one or more file",
		keywords={"bed","stats","multiqc"},
		creationDate="20260130",
		modificationDate="20260130",
		jvarkit_amalgamion = true,
		menu="BED Manipulation"
		)
public class BedStats
	extends Launcher
	{
	private static final Logger LOG = Logger.of(BedStats.class);
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile= null;
	@Parameter(names={"--disable-normalize-contig"},description="do not normalize contig name '1'->'chr1'")
	private boolean disable_normalize_contig_name = false;
	@Parameter(names={"-f","--fraction"},description="for overlap between to BED file. Two BED record overlap if they  both share 'x' fraction of overlap compared to their length" +FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double fraction_overlap = 0.0000001;

	

	
	private class BedSource extends AbstractList<BedLine> {
		String name;
		final List<BedLine> _intervals = new ArrayList<BedLine>();
		
		IntervalTreeMap<BedLine> toTreeMap() {
			IntervalTreeMap<BedLine> map = new IntervalTreeMap<BedLine>();
			for(BedLine rec:_intervals) {
				map.put(rec.toInterval(), rec);
				}
			return map;
			}
		void writeBasicStats(JsonWriter w) throws IOException {
			w.name(this.name);
			w.beginObject();
			w.name("number of intervals");
			w.value(this.size());
			if(!this.isEmpty()) {
				w.name("min interval size");
				w.value(this.stream().mapToInt(B->B.getLengthOnReference()).min().getAsInt());
				w.name("max interval size");
				w.value(this.stream().mapToInt(B->B.getLengthOnReference()).max().getAsInt());
				w.name("average interval size");
				w.value(this.stream().mapToInt(B->B.getLengthOnReference()).average().getAsDouble());
				w.name("median interval size");
				w.value(new Median().evaluate(this.stream().mapToDouble(B->B.getLengthOnReference()).toArray()));
				w.name("distinct contigs");
				w.value(this.stream().map(R->R.getContig()).collect(Collectors.toSet()).size());

				if(this.size()>1) {
					final List<Double> dbl = new ArrayList<Double>(this.size());
					for(int i=0;i+1 < this.size();i++) {
						final Locatable rec1 = this.get(i  );
						final Locatable rec2 = this.get(i+1);
						if(!rec1.contigsMatch(rec2)) continue;
						final double distance =  (rec2.getStart()-1) - rec1.getEnd() ;
						dbl.add(distance);
						}
					if(!dbl.isEmpty()) {
						double[] dbla = dbl.stream().mapToDouble(F->F.doubleValue()).sorted().toArray();
						w.name("consecutive mean distance");
						w.value(new Mean().evaluate(dbla));
						w.name("consecutive median distance");
						w.value(new Median().evaluate(dbla));
						}
					}
				}
			w.endObject();
			}
		@Override
		public int size() {
			return this._intervals.size();
			}
		@Override
		public BedLine get(int index) {
			return this._intervals.get(index);
			}
		@Override
		public Stream<BedLine> stream() {
			return  this._intervals.stream();
			}
		}
	
	private String normalizeContig(final String contig) {
		if(disable_normalize_contig_name) return contig;
		if(!contig.startsWith("chr") && (contig.matches("[0-9XY]+") || contig.matches("M(T)?"))) {
			return "chr"+contig;
			}
		else
			{
			return contig;
			}
		}
	
	private BedSource readBedSource(Path pathOrNull) throws IOException {
		try(BedLineReader br = pathOrNull==null?new BedLineReader(stdin(),"stdin"):
					new BedLineReader(pathOrNull)) {
			long N=0;
			final BedSource src = new BedSource();
			src.name = (pathOrNull==null?"stdin":IOUtils.getFilenameWithoutCommonSuffixes(pathOrNull));
			while(br.hasNext()) {
				BedLine rec = br.next();
				rec = rec.renameContig(normalizeContig(rec.getContig()));
				src._intervals.add(rec);
				N+= rec.getLengthOnReference();
				}
			final SmartComparator cmpCtg = new SmartComparator();
			Collections.sort(src._intervals,(A,B)->{
				int i = cmpCtg.compare(A.getContig(),B.getContig());
				if(i!=0) return i;
				i = Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				return Integer.compare(A.getEnd(), B.getEnd());
				});
			LOG.info("COUNT/LENGTH\t" + src.name + "\t"+ src.size()+"\t" + N);
			return src;
			}
		}
	
	
	private boolean overlap(final Locatable R1,final Locatable R2) {
		if(!R1.overlaps(R2)) return false;
		final int x1 = Math.max(R1.getStart(),R2.getStart());
		final int x2 = Math.min(R1.getEnd(),R2.getEnd());
		final double L = CoordMath.getLength(x1,x2);
		final double L1 = R1.getLengthOnReference();
		final double L2 = R2.getLengthOnReference();
		final double f1 = L/L1;
		final double f2 = L/L2;
		return f1>=this.fraction_overlap && f2>=this.fraction_overlap;
		}
		
	@Override
	public int doWork(final List<String> args) {
		final List<BedSource> beds = new ArrayList<>();
		
		try
			{
			List<Path> paths = IOUtils.unrollPaths(args);
			if(paths.isEmpty()) {
				beds.add(readBedSource(null));
			} else {
				for(Path p: paths) {
					beds.add(readBedSource(p));
				}
			}
			
			try( PrintWriter w0 = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				try(PrefixSuffixWriter pw= new  PrefixSuffixWriter(w0) ) {
					
					
					pw.setPrefix("basic_stats\t");
					
					StringWriter sw = new StringWriter();
					try(JsonWriter w = new JsonWriter(sw)) {
						w.setIndent("  ");
						w.beginObject();
						w.name("id");
						w.value("__ID__basic_stats");
						w.name("section_name");
						w.value("Bed stats");
						w.name("section_description");
						w.value("stats about bed files");
						w.name("plot_type");
						w.value("table");
						w.name("pconfig");
						w.beginObject();
							w.name("id");
							w.value("__ID__basic_stats_table");
							w.name("title");
							w.value("Bed Stats");

						w.endObject();
						w.name("data");
						w.beginObject();
						for(BedSource src: beds) {
							src.writeBasicStats(w);
							}
						w.endObject();
					
						w.endObject();
						}
							
					pw.println(sw.toString());
						if(beds.size()>1) {
						pw.setPrefix("heatmap\t");
						sw = new StringWriter();
						try(JsonWriter w = new JsonWriter(sw)) {
							w.setIndent("  ");
							w.beginObject();
							w.name("id");
							w.value("__ID__compare_bed");
							w.name("section_name");
							w.value("Compare bed");
							w.name("section_description");
							w.value("Compare BEDs. Two bed records overlap if they both share "+this.fraction_overlap+" of their length");
							w.name("plot_type");
							w.value("heatmap");
							w.name("pconfig");
							w.beginObject();
								w.name("title");
								w.value("Bed Overlap (f="+this.fraction_overlap+")");
								w.name("min");
								w.value(0.0);
								w.name("max");
								w.value(1.0);
							w.endObject();
							w.name("xcats");
							w.beginArray();
							for(BedSource src:beds) {
								w.value(src.name);
								}
							w.endArray();
							
							w.name("ycats");
							w.beginArray();
							for(BedSource src:beds) {
								w.value(src.name);
								}
							w.endArray();
							
							
							w.name("data");
							w.beginObject();
							for(BedSource srcA:beds) {
								w.name(srcA.name);
								w.beginObject();
								for(BedSource srcB:beds) {
									w.name(srcB.name);
									int count=0;
									final IntervalTreeMap<BedLine> map = srcB.toTreeMap();
									for(BedLine recA: srcA) {
										if(map.getOverlapping(recA).stream().anyMatch(RECB->overlap(recA,RECB))) {
											count++;
											}
										}
									w.value(srcA.isEmpty()?0:count/(double)srcA.size());
									}
								w.endObject();
								}
							w.endObject();
						
							w.endObject();
							}
						pw.println(sw.toString());
						}
					
					pw.flush();
					}
				}
			
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		
		}
	

	public static void main(final String[] args)
		{
		new BedStats().instanceMainWithExit(args);
		}
	}
