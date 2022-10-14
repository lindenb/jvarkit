/*
The MIT License (MIT)

Copyright (c) 2022 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.ibd;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SortingCollection;

@Program(
		name="ibdreformater01",
		keywords={"ibd"},
		description="reformater ibd data for Fabien Laporte",
		generate_doc = false,
		creationDate = "20221014",
		modificationDate = "20221014"
		)

public class IBDReformater01 extends Launcher{
	private static final Logger LOG = Logger.build(IBDReformater01.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outFile=null;
	@Parameter(names={"--bim"},description="Bim File",required = true)
	private Path bimFile =null;
	@Parameter(names={"--ibd"},description="IBD File",required = true)
	private Path ibdFile =null;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();

	
public static class BimRecord implements Locatable, Comparable<BimRecord> {
	final String contig;
	final int pos;
	boolean keep =false;
	BimRecord(final String contig,int pos) {
		this.contig = contig;
		this.pos = pos;
		}
	@Override
	public String getContig() {
		return contig;
		}
	@Override
	public int getStart() {
		return pos;
		}
	@Override
	public int getEnd() {
		return pos;
		}
	@Override
	public int compareTo(BimRecord o) {
		int i = contig.compareTo(o.contig);
		if(i!=0) return i;
		return Integer.compare(pos, o.pos);
		}
	@Override
	public String toString() {
		return contig+":"+pos;
		}
	}

private static class IBDRecord implements Locatable{
	String sn1;
	int haplo1;
	String sn2;
	int haplo2;
	
	String contig;
	int start;
	int end;
	@Override
	public String getContig() {
		return contig;
		}
	@Override
	public int getStart() {
		return start;
		}
	@Override
	public int getEnd() {
		return end;
		}
	public int compareSample(IBDRecord o) {
		int i = sn1.compareTo(o.sn1);
		if(i!=0) return i;
		return sn2.compareTo(o.sn2);
		}
	public int compareSamplePos(IBDRecord o) {
		int i = compareSample(o);
		if(i!=0) return i;
		i = contig.compareTo(o.contig);
		if(i!=0) return i;
		i = Integer.compare(start, o.start);
		if(i!=0) return i;
		return Integer.compare(end, o.end);
		}
	}

private static class IBDRecordCodec extends AbstractDataCodec<IBDRecord> {
	@Override
	public IBDRecord decode(DataInputStream dis) throws IOException {
		return null;
		}
	@Override
	public void encode(DataOutputStream dos, IBDRecord object) throws IOException {
		
		}
	@Override
	public IBDRecordCodec clone() {
		return new IBDRecordCodec();
		}
}



@Override
public int doWork(List<String> args) {
	try {
		final IntervalTreeMap<BimRecord> bimRecordsIntervalTreeMap = new IntervalTreeMap<>();

		final List<BimRecord> bims = new ArrayList<>(100_000);
		final Pattern ws = Pattern.compile("[\\s]+");
		try(BufferedReader br = IOUtils.openPathForBufferedReading(this.bimFile)) {
			String line;
			while((line=br.readLine())!=null) {
				final String[] tokens = ws.split(line);
				final BimRecord rec = new BimRecord(tokens[0], Integer.parseInt(tokens[3]));
				bimRecordsIntervalTreeMap.put(new Interval(rec), rec);
				bims.add(rec);
				}
			}
		Collections.sort(bims);
		LOG.info("number of record in "+this.bimFile+" "+bims.size());
		
		SortingCollection<IBDRecord> sorter = 
				SortingCollection.newInstance(
						IBDRecord.class,
						new IBDRecordCodec(),
						(A,B)->A.compareSamplePos(B),
						writingSortingCollection.getMaxRecordsInRam(),
						writingSortingCollection.getTmpPaths()
						);
		sorter.setDestructiveIteration(true);
		try(BufferedReader br = IOUtils.openPathForBufferedReading(this.ibdFile)) {
			String line;
			while((line=br.readLine())!=null) {
				final String[] tokens = ws.split(line);
				final IBDRecord rec = new IBDRecord();
				if(tokens[0].compareTo(tokens[2]) <= 0) {
					rec.sn1 = tokens[0];
					rec.haplo1 = Integer.parseInt(tokens[1]);
					rec.sn2 = tokens[2];
					rec.haplo2 = Integer.parseInt(tokens[3]);
					}
				else
					{
					rec.sn2 = tokens[0];
					rec.haplo2 = Integer.parseInt(tokens[1]);
					rec.sn1 = tokens[2];
					rec.haplo1 = Integer.parseInt(tokens[3]);
					}
				rec.contig = tokens[4];
				
				rec.start = Integer.parseInt(tokens[5]);
				rec.end = Integer.parseInt(tokens[6]);
				boolean found=false;
				for(BimRecord bim:bimRecordsIntervalTreeMap.getOverlapping(rec)) {
					bim.keep=true;
					found=true;
					}
				if(found) {
					sorter.add(rec);
					}
				}
			
			}
		sorter.doneAdding();
		bims.removeIf(B->!B.keep);
		if(bims.isEmpty()) {
			LOG.error("no overlapping record from BIM was found");
			return -1;
			}
		try(PrintWriter pw =super.openPathOrStdoutAsPrintWriter(this.outFile)) {
			pw.append("Couple");
			for(BimRecord rec:bims) {
				pw.append("\t");
				pw.append(rec.toString());
				}
			pw.append("\n");
			try(CloseableIterator<IBDRecord> iter0 = sorter.iterator()) {
				final EqualIterator<IBDRecord> iter = new EqualIterator<>(iter0, (A,B)->A.compareSample(B));
				while(iter.hasNext()) {
					final List<IBDRecord> row = iter.next();
					final IBDRecord first = row.get(0);
					pw.print(first.sn1);
					pw.print("/");
					pw.print(first.sn2);
					for(BimRecord rec:bims) {
						pw.append("\t");
						final Set<Integer> ibdset = new TreeSet<>();
						row.stream().
								filter(R->R.overlaps(rec)).
								forEach(R->{ibdset.add(R.haplo1);ibdset.add(R.haplo2);});
						if(ibdset.isEmpty()) ibdset.add(0);
						pw.append(ibdset.stream().map(I->I.toString()).collect(Collectors.joining(",")));
						}
					pw.println();
					}
				iter.close();
				}
			pw.flush();
			}
		sorter.cleanup();
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	}

public static void main(String[] args) {
	new IBDReformater01().instanceMainWithExit(args);
	}
}
