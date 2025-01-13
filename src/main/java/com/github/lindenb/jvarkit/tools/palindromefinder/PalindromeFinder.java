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
package com.github.lindenb.jvarkit.tools.palindromefinder;

import java.io.Closeable;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.LongList;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequence;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequenceReader;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;

@Program(
		name="palindromefinder",
		description="palindrome finder using suffix array",
		keywords={"palindrome","fasta"},
		creationDate = "20240424",
		modificationDate = "20240424",
		generate_doc = false
		)
public class PalindromeFinder  extends Launcher {
private static final Logger LOG = Logger.build(PalindromeFinder.class).make();

@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private Path outfile=null;
@Parameter(names={"-m"},description="min size")
private int min_palindrome_size = 10;
@Parameter(names={"-d"},description="min distance bewteen segment 1 and segment 2" + DistanceParser.OPT_DESCRIPTION ,splitter = NoSplitter.class,converter = DistanceParser.StringConverter.class)
private int min_distance = 0;
@Parameter(names={"-D"},description="max distance bewteen segment 1 and segment 2" + DistanceParser.OPT_DESCRIPTION ,splitter = NoSplitter.class,converter = DistanceParser.StringConverter.class)
private int max_distance = 300_000_000;
@Parameter(names={"-n"},description="skip lowercase letters in the fasta")
private boolean skip_lowercase = false;
@Parameter(names={"-p"},description="palindrome only")
private boolean palindrome_only = false;
@Parameter(names={"-w"},description="run genome wide")
private boolean genome_wide = false;


@Parameter(names={"-g"},description="min GC%." + FractionConverter.OPT_DESC ,splitter = NoSplitter.class,converter = FractionConverter.class)
private double min_gc=0.0;
@Parameter(names={"-G"},description="max GC%." + FractionConverter.OPT_DESC ,splitter = NoSplitter.class,converter = FractionConverter.class)
private double max_gc=1.0;
@ParametersDelegate
private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();

private static final char EOF='\0';
private List<FastaSequence> current_genome = Collections.emptyList();


private static class Position {
	int tid;
	int pos;
	boolean negativeStrand;
	Position(byte tid,int pos, boolean neg) {
		this.tid=tid;
		this.pos=pos;
		this.negativeStrand= neg;
		}
	}

private class PositionCodec extends AbstractDataCodec<Position> {
	@Override
	public Position decode(DataInputStream dis) throws IOException {
		try {
			byte tid = (current_genome.size()==1?0:dis.readByte());
			int p=dis.readInt();
			boolean negativeStrand=dis.readBoolean();
			return new Position(tid,p,negativeStrand);
			}
		catch(EOFException err) {
			return null;
			}
		}
	@Override
	public void encode(DataOutputStream dos, Position p) throws IOException {
		if(current_genome.size()>1) dos.writeByte(p.tid);
		dos.writeInt(p.pos);
		dos.writeBoolean(p.negativeStrand);
		}
	@Override
	public PositionCodec clone() {
		return new PositionCodec();
		}
	}

private class PositionComparator implements Comparator<Position> {
	
	@Override
	public int compare(Position pos1, Position pos2) {
		int extend=0;
		for(;;) {
			char c1 = charAt(pos1,extend);
			char c2 = charAt(pos2,extend);
			if(c1==EOF) {
				if(c2==EOF) return 0;
				return -1;
				}
			else if(c1!=c2) {
				if(c2==EOF) {
					return 1;
					}
				int d = Character.compare(c1, c2);
				if(d!=0) return d;
				}
			extend++;
			}
		}
	}

private class StoredDatabase implements Closeable,LongList<Position>{
	final File db;
	final long nItems;
	final RandomAccessFile fp;
	long cache_index=-1L;
	final List<Position> cache=new ArrayList<>(10_000);
	StoredDatabase(final File db, long nItems)  throws IOException {
		this.db=db;
		this.nItems=nItems;
		this.fp= new RandomAccessFile(db, "r");
		}
	@Override
	public long size() {
		return nItems;
		}
	@Override
	public Position get(long idx) {
		if(cache_index==-1L || idx < cache_index || idx >= cache_index+cache.size()) {
			final int sizeOf = Integer.BYTES+1+(current_genome.size()==1?0:Byte.BYTES);
			cache_index=idx;
			cache.clear();
			try {
				this.fp.seek(sizeOf*idx);
				while(cache_index+cache.size() < size() && cache.size()< 10_000) {
					byte tid = (current_genome.size()==1?0:fp.readByte());
					int p=this.fp.readInt();
					boolean n= this.fp.readBoolean();
					cache.add(new Position(tid,p, n));
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		return cache.get((int)(idx-this.cache_index));
		}
	
	@Override
	public void close() throws IOException {
		db.delete();
		}
	}

private char charAt(Position p,int extend) {
	char c;
	final FastaSequence current_seq = current_genome.get(p.tid);
	if(p.negativeStrand) {
		int i=p.pos-extend;
		if(i<0) return EOF;
		c = AcidNucleics.complement(current_seq.charAt(i));
		}
	else
		{
		int i=p.pos+extend;
		if(i>=current_seq.length()) return EOF;
		c = current_seq.charAt(i);
		}
	if(skip_lowercase && Character.isLowerCase(c)) return EOF;
	if(!AcidNucleics.isATGC(c)) return EOF;
	c= Character.toUpperCase(c);
	return c;
	}

private StoredDatabase buildStoredDatabase() throws IOException {
	SortingCollection<Position> sorter=null;
	
	try {
		
		sorter = SortingCollection.newInstance(
				Position.class,
				new PositionCodec(),
				new PositionComparator(),
				this.writingSortingCollection.getMaxRecordsInRam(),
				this.writingSortingCollection.getTmpPaths()
				);
		sorter.setDestructiveIteration(true);
		for(int k=0;k< current_genome.size();++k) {
			final FastaSequence current_seq = current_genome.get(k);
			for(int i=0;i< current_seq.length();i++) {
				char c = current_seq.charAt(i);
				if(skip_lowercase && Character.isLowerCase(c)) continue;
				if(!AcidNucleics.isATGC(c)) continue;
				sorter.add(new Position((byte)k,i,true));
				sorter.add(new Position((byte)k,i,false));
				}
			}
			sorter.doneAdding();
		File db= File.createTempFile("tmp.", ".bin");
		long nItems=0L;
		try(CloseableIterator<Position> iter=sorter.iterator()) {
			try(DataOutputStream dos=new DataOutputStream(new FileOutputStream(db))) {
				while(iter.hasNext()) {
					final Position p = iter.next();
					if(current_genome.size()>1) dos.writeByte(p.tid);
					dos.writeInt(p.pos);
					dos.writeBoolean(p.negativeStrand);
					nItems++;
					}
				dos.flush();
				}
			}
		return new StoredDatabase(db,nItems);
		}
	catch(IOException err) {
		throw err;
		}
	finally
		{
		if(sorter!=null) sorter.cleanup();
		}
	}

static boolean is_reverse_complement(char c1,char c2) {
	return Character.toUpperCase(AcidNucleics.complement(c1)) == Character.toUpperCase(c2);
	}

private float toGCPercent(Position p,int len) {
	float gc=0;
	for(int x=0;x< len;++x) {
		switch(charAt(p,x)) {
			case 'c': case 'C':
			case 'G': case 'g':
			case 's': case 'S': gc++;break;
			default:break;
			}
		}
	return gc/(len);
	}

private void scan(final PrintWriter out,final List<FastaSequence> sequences) throws IOException {
	if(sequences.isEmpty()) return;
	if(sequences.size()>=Byte.MAX_VALUE) {
		throw new IOException("Cannot process more than "+Byte.MAX_VALUE+" sequences but got "+sequences.size());
		}
	this.current_genome= Collections.unmodifiableList(sequences);
	SimpleInterval prevLoc1=null;
	SimpleInterval prevLoc2=null;
	try(StoredDatabase db=buildStoredDatabase()) {
		for(long i=0;i+1< db.size();i++) {
			final Position pos1 = db.get(i);
			for(long j=i+1;j< db.size();j++) {
				final Position pos2 = db.get(j);
				if(this.palindrome_only && (pos1.tid!=pos2.tid || pos1.negativeStrand==pos2.negativeStrand)) continue;
				int extend=0;
				for(;;) {
					char c1 = charAt(pos1,extend);
					char c2 = charAt(pos2,extend);
					if(c1!=c2 || c1==EOF || c2==EOF) break;
					extend++;
					}
				if(extend<this.min_palindrome_size)  break;
				
				
				float gc = toGCPercent(pos1,extend);
				if(gc<this.min_gc || gc> this.max_gc) continue;
				int start1 = pos1.negativeStrand?pos1.pos-(extend-1):pos1.pos;
				int end1  = pos1.negativeStrand?pos1.pos+1:pos1.pos+extend;
				final SimpleInterval  loc1=new SimpleInterval(current_genome.get(pos1.tid).getName(),start1+1,end1);
				
				int start2 = pos2.negativeStrand?pos2.pos-(extend-1):pos2.pos;
				int end2  = pos2.negativeStrand?pos2.pos+1:pos2.pos+extend;
				final SimpleInterval  loc2=new SimpleInterval(current_genome.get(pos2.tid).getName(),start2+1,end2);
				
				int distance;
				if(loc1.contigsMatch(loc2)) {
					distance = loc1.getDistanceTo(loc2);
					if(distance < this.min_distance) continue;
					if(distance > this.max_distance) continue;
					}
				else
					{
					distance=-1;
					}
				
				if(prevLoc1!=null && prevLoc2!=null &&
						(
						(prevLoc1.overlaps(loc1) && prevLoc2.overlaps(loc2)) ||
						(prevLoc1.overlaps(loc2) && prevLoc2.overlaps(loc1))
						)	
					) continue;
				prevLoc1=loc1;
				prevLoc2=loc2;
				
				
				out.print(loc1.toBed3());
				out.print('\t');
				out.print(pos1.negativeStrand?'-':'+');
				out.print('\t');
				out.print(loc2.toBed3());
				out.print('\t');
				out.print(pos2.negativeStrand?'-':'+');
				out.print('\t');

				if(loc1.withinDistanceOf(loc2, 1)) {
					out.print('.');
					out.print('\t');
					out.print('.');
					out.print('\t');
					out.print(current_genome.get(pos1.tid).subSequence(Math.min(start1, start2), Math.max(end1, end2)));
					}
				else
					{
					out.print(current_genome.get(pos1.tid).subSequence(start1,end1));
					out.print('\t');
					out.print(current_genome.get(pos2.tid).subSequence(start2,end2));
					out.print('\t');
					out.print('.');
					}
				out.print('\t');
				out.print(extend);
				out.print('\t');
				out.print(gc);
				out.print('\t');
				out.print(distance);
				out.println();
				}
			}
		}
	
	}


private void scan(PrintWriter out,InputStream in) throws IOException {
	try(CloseableIterator<FastaSequence> r=new FastaSequenceReader().iterator(in)) {
		if(genome_wide) {
			scan(out,r.toList());
			}
		else
			{
			while(r.hasNext()) {
				scan(out,Collections.singletonList(r.next()));
				}
			}
		}
	}

@Override
public int doWork(List<String> args) {
	if(min_palindrome_size<4) {
		LOG.error("min size is too low");
		return -1;
		}
	try {
		final List<Path> paths = IOUtils.unrollPaths(args);
		try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outfile)) {
			if(paths.isEmpty()) {
				scan(out,System.in);
				}
			else
				{
				for(Path p:paths) {
					try(InputStream in = IOUtils.openPathForReading(p)) {
						scan(out,in);
						}
					}
				}
			out.flush();
			}
		return 0;
		}
	catch(Throwable err) {
		LOG.error(err);
		return -1;
		}
	}

public static void main(String[] args) {
	new PalindromeFinder().instanceMainWithExit(args);
	}
}
