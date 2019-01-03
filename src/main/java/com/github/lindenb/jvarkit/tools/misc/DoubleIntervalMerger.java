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
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.GenomicIndexUtil;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;

public class DoubleIntervalMerger extends Launcher {
	private static final Logger LOG = Logger.build(DoubleIntervalMerger.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-d","--distance"},description="Max distance between two read to test if they both end at the same ~ position.")
	private int fuzzy_distance = 100;

	private static long ID_GENERATOR = 0;
	
	private static class DoubleInterval
		{
		final long id = ++ID_GENERATOR;
		final String contig1;
		int start1;
		int end1;
		final String contig2;
		int start2;
		int end2;
		int count = 0;
		DoubleInterval(final Interval interval1,final Interval interval2) {
			this.contig1 = interval1.getContig();
			this.start1 = interval1.getStart();
			this.end1 = interval1.getEnd();
			this.contig2 = interval2.getContig();
			this.start2 = interval2.getStart();
			this.end2 = interval2.getEnd();
			}
		@Override
		public int hashCode() {
			return Long.hashCode(this.id);
			}
		@Override
		public boolean equals(Object obj) {
			return this.id == DoubleInterval.class.cast(obj).id;
			}
		}
	
	
	
	private final Map<Integer, Set<DoubleInterval>> bin2set = new HashMap<>(GenomicIndexUtil.MAX_BINS);
	
	private int bin(final DoubleInterval bp)
		{
		return GenomicIndexUtil.regionToBin(
				start1(bp),
				end1(bp)
				);
		}
	
	private int start1(final DoubleInterval bp)
		{
		return Math.max(1,(bp.start1-1)-fuzzy_distance);
		}
	private int end1(final DoubleInterval bp) {
		return (bp.end1)+fuzzy_distance;
		}
	
	private List<DoubleInterval> toList() {
		return this.bin2set.
				values().
				stream().
				flatMap(S->S.stream()).
				collect(Collectors.toList());
		}
	
	private void insert(final DoubleInterval bp)
		{
		boolean done=false;
		while(!done)
			{
			done = true;
			final int b= bin(bp);
			if(b>=GenomicIndexUtil.MAX_BINS) throw new IllegalStateException();
			final BitSet bitset = GenomicIndexUtil.regionToBins(start1(bp), end1(bp));
			for(int x=0;x<GenomicIndexUtil.MAX_BINS;++x)
				{
				if(!bitset.get(x)) continue;
				final Set<DoubleInterval> set  = this.bin2set.get(b);
				if(set==null) continue;
				DoubleInterval overlapper = null;
				for(DoubleInterval other : set) {
					if(bp.tid1!=other.tid1) throw new IllegalStateException();
					if(bp.tid2!=other.tid2) continue;
					if( overlap( bp.start1  ,bp.end1 , other.start1, other.end1 ) &&
						overlap( bp.start2  ,bp.end2 , other.start2, other.end2 ))
						{
						overlapper = other;
						break;
						}
					}
				if(overlapper!=null)
					{
					//System.err.println("merging "+bp+" "+overlapper);
					set.remove(overlapper);
					//
					bp.start1 = Math.min(bp.start1, overlapper.start1);
					bp.end1 = Math.max(bp.end1, overlapper.end1);
					//
					bp.start2 = Math.min(bp.start2, overlapper.start2);
					bp.end2 = Math.max(bp.end2, overlapper.end2);
					//
					bp.count += overlapper.count;
					done=false;
					break;
					}
				}
			}
		final int b= bin(bp);
		Set<DoubleInterval> set  = this.bin2set.get(b);
		if(set==null) {
			set = new HashSet<>();
			this.bin2set.put(b,set);
			}
		set.add(bp);
		}
	
	private Interval parseInterval(final String tokens[],int off)
		{
		String contig = tokens[off+0];
		int start = Integer.parseInt(tokens[off+1]);
		int end = Integer.parseInt(tokens[off+2]);
		return new Interval(contig,start,end);
		}
	
	private DoubleInterval makeDoubleInterval(final Interval i1,final Interval i2) {
		int d = i1.getContig().compareTo(i2.getContig());
		if(d>0){
			return makeDoubleInterval(i2,i1);
			}
		if(d==0 )
			{
			d= Integer.compare(i1.getStart(),i2.getStart());
			if(d>0) return makeDoubleInterval(i2,i1);
			
			}
		return new DoubleInterval(i1, i2);
		}
	
	@Override
	public int doWork(final List<String> args) {
		BufferedReader br = null;
		PrintWriter pw = null;
		try
			{
			String line;
			br = super.openBufferedReader(oneFileOrNull(args));
			while((line=br.readLine())!=null) {
				if(StringUtil.isBlank(line)) continue;
				if(line.startsWith("#")) continue;
				final String tokens[] = CharSplitter.TAB.split(line);
				final Interval i1 = parseInterval(tokens,0);
				final Interval i2 = parseInterval(tokens,4);
				DoubleInterval di = makeDoubleInterval(i1, i2);
				insert(di);
				
			}
			
			br.close();
			br=null;
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			
			pw.flush();
			pw.close();
			pw=null;
			return 0;
			}
		catch(final Exception err)
			{
			return -1;
			}
		finally
			{
			CloserUtil.close(pw);
			CloserUtil.close(br);
			}
		
		}
	
	
	public static void main(final String[] args) {
		new DoubleIntervalMerger().instanceMainWithExit(args);

	}

}
