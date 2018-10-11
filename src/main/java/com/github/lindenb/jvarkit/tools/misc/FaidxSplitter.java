/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

/**
BEGIN_DOC

# Motivation

For WGS sequencing, find 

END_DOC
 */

@Program(
	name="faidxsplitter",
	description="Split bed of reference genome into overlapping parts",
	keywords={"vcf"}
	)
public class FaidxSplitter  extends Launcher
	{
	private final Logger LOG=Logger.build(FaidxSplitter.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected File outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	protected File faidx = null;
	@Parameter(names={"-g","--gaps"},description="gap file. A bed file containing the known gaps in the genome. "
			+ "E.g: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz ")
	protected String gapSource = null;
	@Parameter(names={"-G","--gene"},description="gene file. You shouldn't break a record in this file.")
	protected String geneSource = null;
	@Parameter(names={"-w","--size"},description="BED Size. The genome should be split into BED of 'w' size.")
	protected int window_size = 1_000_000;
	@Parameter(names={"-x","--overlap"},description="Overlap  Size. The resulting BED region should overlap with 'x' bases.")
	protected int overlap_size = 1_000;
	
	private void echo(final PrintWriter pw,final Interval interval) {
		pw.println(interval.getContig()+ "\t" + (interval.getStart()-1) + "\t" + interval.getEnd());
		}
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter pw=null;
		try {
			if(this.faidx==null) {
				LOG.error("reference undefined");
				return -1;
				}
			if(window_size<=0 || overlap_size >= window_size) {
				LOG.info("bad parameters");
				return -1;
				}
			final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(this.faidx);
			if(dict==null)
				{
				LOG.error("no dictionary found");
				return -1;
				}
			final BedLineCodec codec = new BedLineCodec();
			final IntervalTreeMap<Interval> intervals1 = new IntervalTreeMap<>();
			dict.getSequences().
				stream().
				map(SSR->new Interval(SSR.getSequenceName(),1,SSR.getSequenceLength())).
				forEach(I->intervals1.put(I, I));		
			
			if(this.gapSource!=null) {
				String line;
				final BufferedReader br=IOUtils.openURIForBufferedReading(this.gapSource);
				while((line=br.readLine())!=null)
					{
					final BedLine gap = codec.decode(line);
					if(gap==null) continue;
					final Collection<Interval> list = intervals1.getOverlapping(gap);
					for(Interval i:list) {
						if(!gap.overlaps(i)) throw new IllegalStateException();
						if(gap.contains(i))
							{
							intervals1.remove(i);
							}
						else
							{
							intervals1.remove(i);
							if(i.getStart() < gap.getStart())
								{
								Interval left = new Interval(i.getContig(),i.getStart(),gap.getStart()-1);
								if(left.length()>0) intervals1.put(left, left);
								}
							if(i.getEnd() > gap.getEnd())
								{
								Interval right = new Interval(i.getContig(),gap.getEnd()+1,i.getEnd());
								if(right.length()>0) intervals1.put(right, right);
								}
							}
						}
					}
				br.close();
			}
			
			final IntervalTreeMap<Interval> geneMap = new IntervalTreeMap<>();
			final Function<Interval, Interval> findGene = I->{
				final Collection<Interval> coll = geneMap.getOverlapping(I);
				if(coll.isEmpty()) return null;
				int min = coll.stream().mapToInt(G->G.getStart()).min().getAsInt();
				int max = coll.stream().mapToInt(G->G.getEnd()).max().getAsInt();
				return new Interval(I.getContig(),min,max);
				};
			if(this.geneSource!=null) {
				String line;
				final BufferedReader br=IOUtils.openURIForBufferedReading(this.geneSource);
				while((line=br.readLine())!=null)
					{
					final BedLine gene = codec.decode(line);
					if(gene==null) continue;
					final Interval g = new Interval(gene.getContig(), gene.getStart(), gene.getEnd());
					geneMap.put(g, g);
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
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			for(Interval interval : intervals2) {
				if(interval.length()<= this.window_size)
					{
					echo(pw,interval);
					continue;
					}
				int start =interval.getStart();
				while(start+window_size < interval.getEnd())
					{
					Interval record = new Interval(interval.getContig(),start, start+window_size);
					Interval gene = findGene.apply(record);
					if(gene==null || record.contains(gene)) {
						echo(pw,record);
						}
					else
						{
						if(gene.getStart() < start) throw new IllegalStateException();
						record = new Interval(interval.getContig(),start,gene.getEnd() + this.overlap_size + 10);
						echo(pw,record);
						}
					start = record.getEnd()-overlap_size;
					}
				}
			pw.flush();
			pw.close();
			pw = null;
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
