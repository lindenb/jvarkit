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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

/**
BEGIN_DOC

# Motivation

For WGS sequencing, a tool I used to split the genome and parallelize things...

# Example

```
$ wget -q -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz" | gunzip -c | cut -f 2,3,4 > jeter.gaps.txt

$ wget -q -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz"  | gunzip -c | cut -f 3,5,6 > jeter.genes.txt

$ java -jar dist/faidxsplitter.jar -gap jeter.gaps.txt -gene jeter.genes.txt -R src/test/resources/human_b37.dict

1	10000	177417
1	227417	267719
1	317719	471368
1	521368	1521369
1	1520369	2634220
1	2684220	3689209
1	3688209	3845268
1	3995268	4995269
1	4994269	6241184
1	6240184	7830767
(...)

```

 

END_DOC
 */

@Program(
	name="faidxsplitter",
	description="Split bed of reference genome into overlapping parts",
	keywords={"vcf","reference","bed"},
	modificationDate="20190307"
	)
public class FaidxSplitter  extends Launcher
	{
	private final Logger LOG=Logger.build(FaidxSplitter.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected File outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	protected Path faidx = null;
	@Parameter(names={"-gap","--gaps","--gap"},description="gap file. A bed file containing the known gaps in the genome. "
			+ "E.g: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz ")
	protected String gapSource = null;
	@Parameter(names={"-gene","--gene","--genes"},description="gene file. You shouldn't break a record in this file.")
	protected String geneSource = null;
	@Parameter(names={"-w","--size"},description="BED Size. The genome should be split into BED of 'w' size.")
	protected int window_size = 1_000_000;
	@Parameter(names={"-x","--overlap"},description="Overlap  Size. The resulting BED region should overlap with 'x' bases.")
	protected int overlap_size = 1_000;
	@Parameter(names={"-s","--small"},description="If it remains 's' bases in the BED split to the end of the chromosome, extends the current BED.")
	protected int small_size = 1_000;
	@Parameter(names={"--exclude"},description="A regular expression to exclude some chromosomes from the dictionary.")
	protected String excludeChromStr = "(chr)?(M|MT)$";

	
	private long count_echoed = 0L;
	
	private void echo(final PrintWriter pw,final Interval interval) {
		if(interval.length()<1) {
			LOG.warn("Ignoring empty interval "+interval);
			return;
			}
		pw.println(interval.getContig()+ "\t" + (interval.getStart()-1) + "\t" + interval.getEnd());
		pw.flush();
		++count_echoed;
		}
	
	private List<Interval> split(final Interval contig,final Interval gap)
		{
		if(!contig.overlaps(gap)) throw new IllegalStateException();
		if(contig.length()==0) return Collections.emptyList();
		if(gap.contains(contig)) return Collections.emptyList();
		final List<Interval> intervals1 = new ArrayList<>(2);
		
		if(contig.getStart() < gap.getStart())
			{
			final Interval left = new Interval(contig.getContig(),contig.getStart(),gap.getStart()-1);
			if(left.length()>0) intervals1.add(left);
			}
		if(contig.getEnd() > gap.getEnd())
			{
			final Interval right = new Interval(contig.getContig(),gap.getEnd()+1,contig.getEnd());
			if(right.length()>0) intervals1.add(right);
			}
		
		return intervals1;
		}
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter pw=null;
		try {
			if(this.faidx==null) {
				LOG.error("reference undefined");
				return -1;
				}
			if(window_size<=0 || overlap_size >= window_size || this.small_size<0) {
				LOG.info("bad parameters");
				return -1;
				}
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.faidx);
			
			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(dict);
			final BedLineCodec codec = new BedLineCodec();
			final Pattern excludePat = Pattern.compile(this.excludeChromStr);
			final IntervalTreeMap<Interval> intervals1 = new IntervalTreeMap<>();
			dict.getSequences().
				stream().
				filter(R->!excludePat.matcher(R.getSequenceName()).matches()).
				map(SSR->new Interval(SSR.getSequenceName(),1,SSR.getSequenceLength())).
				forEach(I->intervals1.put(I, I));		
			
			if(intervals1.isEmpty()) {
				LOG.error("No sequence in dict after filtering.");
				return -1;
				}
			
			final IntervalTreeMap<Interval> geneMap = new IntervalTreeMap<>();
			final Function<Interval, Interval> findGene = I->{
				final Collection<Interval> coll = geneMap.getOverlapping(I);
				if(coll.isEmpty()) return null;
				final int min = coll.stream().mapToInt(G->G.getStart()).min().getAsInt();
				final int max = coll.stream().mapToInt(G->G.getEnd()).max().getAsInt();
				return new Interval(I.getContig(),min,max+this.overlap_size+1);
				};
			if(this.geneSource!=null) {
				LOG.info("loading genes "+this.geneSource);
				String line;
				final Set<String> notFound = new HashSet<>();
				final BufferedReader br=IOUtils.openURIForBufferedReading(this.geneSource);
				while((line=br.readLine())!=null)
					{
					final BedLine gene = codec.decode(line);
					if(gene==null) continue;
					final String contig2 = contigNameConverter.apply(gene.getContig());
					if(StringUtil.isBlank(contig2)) {
						if(notFound.add(gene.getContig())) LOG.warn("gene contig not found :"+gene.getContig());
						continue;
					}
					if(excludePat.matcher(contig2).matches()) continue;
					
					final Interval g = new Interval(contig2, gene.getStart(), gene.getEnd());
					
					if(geneMap.getOverlapping(g).stream().anyMatch(X->X.contains(g))) continue;
					
					geneMap.put(g, g);
					}
				br.close();
				}
			
			if(this.gapSource!=null) {
				LOG.info("loading gaps "+this.gapSource);
				final Set<String> notFound = new HashSet<>();
				String line;
				final BufferedReader br=IOUtils.openURIForBufferedReading(this.gapSource);
				while((line=br.readLine())!=null)
					{
					final BedLine gap0 = codec.decode(line);
					if(gap0==null) continue;
					final String contig2 = contigNameConverter.apply(gap0.getContig());
					if(StringUtil.isBlank(contig2)) {
						if(notFound.add(gap0.getContig())) LOG.warn("gap contig not found :"+gap0.getContig());
						continue;
					}
					if(excludePat.matcher(contig2).matches()) continue;
					final Interval gap = new Interval(contig2,gap0.getStart(),gap0.getEnd());
					
					
					final Collection<Interval> genes = geneMap.getOverlapping(gap);
					for(final Interval gene_interval :genes) {
						if(!gap.overlaps(gene_interval)) throw new IllegalStateException();
						LOG.warn("gene "+gene_interval +" overlap gap "+gap+ " "+this.split(gene_interval, gap));
						geneMap.remove(gene_interval);
						this.split(gene_interval,gap).forEach(N->geneMap.put(N, N));
						}
					if(geneMap.containsOverlapping(gap)) throw new IllegalStateException();
						
					
					final Collection<Interval> list = intervals1.getOverlapping(gap);
					for(final Interval i:list) {
						if(!gap.overlaps(i)) throw new IllegalStateException();
						intervals1.remove(i);
						this.split(i,gap).forEach(N->intervals1.put(N, N));
						}
						
					}
				br.close();
			}
			
			final Comparator<String> contigCmp = new ContigDictComparator(dict);
			final List<Interval> intervals2 = new ArrayList<>(intervals1.values());
			intervals2.sort((A,B)->{
				final int i= contigCmp.compare(A.getContig(),B.getContig());
				if(i!=0) return i;
				return A.getStart()-B.getStart();
				});
			/* start writing output */
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			for(final Interval interval : intervals2) {
				if(interval.length()<= this.window_size)
					{
					echo(pw,interval);
					continue;
					}
				int start =interval.getStart();
				while(start < interval.getEnd())
					{
					int chromEnd  = Math.min(interval.getEnd(),start+this.window_size);
					
					// extend interval while there is a gene that shouldn't be broken */
					for(;;)
						{
						final Interval record = new Interval(
								interval.getContig(),
								start,
								chromEnd
								);


						final Interval gene = findGene.apply(record);
						// no gene 
						if( gene==null ||record.contains(gene) ) break;
						
						if(gene.getStart() < record.getStart()) throw new IllegalStateException(
								"gene start "+gene.getStart()+" < "+start+" "+
								"\ngene:\t" +gene+"\ninterval\t"+interval+"\nrecord\t"+record
								);
						
						chromEnd = Math.min(interval.getEnd(),gene.getEnd() /* not here, already added above  + this.overlap_size */);
						if(chromEnd>=interval.getEnd()) break;
						}
					
					if(interval.getEnd() - chromEnd<= this.small_size)
						{
						chromEnd = interval.getEnd();
						}
					
					echo(pw,new Interval( interval.getContig(), start, chromEnd ));
					if(chromEnd >=interval.getEnd()) break;
					start = chromEnd - this.overlap_size +1;
					}
				}
			pw.flush();
			pw.close();
			pw = null;
			LOG.info("Done N="+this.count_echoed);
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(pw);
			}
		}
	
	public static void main(final String[] args) {
		new FaidxSplitter().instanceMainWithExit(args);
		}
	}
