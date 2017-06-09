package com.github.lindenb.jvarkit.tools.genscan;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
@Program(
		name="genscan2",
		description="Paint a Genome Scan picture from a Tab delimited file (CHROM/POS/VALUE1/VALUE2/....).",
		keywords={"chromosome","reference","chart","visualization"})
public class GenScan2 extends Launcher {
	
	private static final Logger LOG = Logger.build(GenScan2.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidxFile = null;
	@Parameter(names={"-r","--region"},description="One or more Specific region to observe. empty string is whole genome. Or the name of a chromosome. Or an interval.")
	private Set<String> regions = new HashSet<>();
	@Parameter(names={"-miny","--miny","--ymin","-ymin"},description="min Y value")
	private double minY=0;
	@Parameter(names={"-maxy","--maxy","--ymax","-ymax"},description="max Y value")
	private double maxY=100.0;
	@Parameter(names={"-W","--width"},description="width")
	private int WIDTH=10000;
	@Parameter(names={"-H","--height"},description="height")
	private int HEIGHT=700;
	@Parameter(names={"-dbc","--distancebetweencontigs"},description="number of pixels between contigs")
	private double distanceBetweenContigs = 1;
	
	private SAMSequenceDictionary dict=null;
	private List<DisplayRange> displayRanges = new ArrayList<>();
	private long genomeViewLength=0L;
	
	private interface Data extends Locatable
		{
		public double getValue();
		}

	
	private abstract class DisplayRange implements Locatable
		{
		double startX;
		double width=0;
		
		
		public abstract String getLabel();
		public abstract int length();
		public boolean contains(final Data data) {
			if(!data.getContig().equals(this.getContig())) return false;
			if(data.getEnd() < this.getStart()) return false;
			if(data.getStart() > this.getEnd()) return false;
			
			return true;
			}
		}
	private class DisplayInterval extends DisplayRange
		{
		private final Interval interval;
		DisplayInterval(final String s)
			{
			final IntervalParser intervalParser = new IntervalParser(GenScan2.this.dict);
			this.interval = intervalParser.parse(s);
			if(this.interval==null) {
				throw new IllegalArgumentException("cannot parse interval "+s);
				}
			}
		@Override
		public String getLabel() {
			return this.getContig()+":"+this.getStart()+"-"+this.getEnd();
			}
		
		@Override
		public String getContig() {
			return this.interval.getContig();
			}
		@Override
		public int getStart() {
			return this.interval.getStart();
			}
		
		@Override
		public int getEnd() {
			return this.interval.getEnd();
			}
		
		@Override
		public int length() {
			return this.interval.length();
			}
		}
	
	private class DisplayContig extends DisplayRange
		{
		private final SAMSequenceRecord ssr;
		DisplayContig(final String contig) {
			this.ssr = GenScan2.this.dict.getSequence(contig);
			if(this.ssr==null) {
				throw new IllegalArgumentException("cannot find contig in dictionary "+contig);
				}
			}
		DisplayContig(final SAMSequenceRecord ssr) {
			this.ssr = ssr;
			if(this.ssr==null) {
				throw new IllegalArgumentException("SAMSequenceRecord==null");
				}
			}
		
		@Override
		public String getLabel() {
			return getContig();
			}
		
		@Override
		public String getContig() {
			return this.ssr.getSequenceName();
			}
		
		@Override
		public int getStart() {
			return 1;
			}
		
		@Override
		public int getEnd() {
			return this.ssr.getSequenceLength();
			}
		
		@Override
		public int length() {
			return this.ssr.getSequenceLength();
			}
		}
	
	
	private Data parseData(final String line) {
		return null;
	}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.faidxFile==null) {
			LOG.error("Reference missing");
			return -1;
			}
		
		BufferedReader r=null;
		try {
			this.dict = SAMSequenceDictionaryExtractor.extractDictionary(this.faidxFile);
			for(final String rgnstr: this.regions) {
				if(rgnstr.trim().isEmpty()) continue;
				final DisplayRange contig;
				if(rgnstr.contains(":")) {
					contig = new DisplayInterval(rgnstr);
				} else
					{
					contig = new DisplayContig(rgnstr);
					}
				if(contig.length()==0) continue;
				this.displayRanges.add(contig);
				}
			// no region, add whole genome
			for(final SAMSequenceRecord ssr:this.dict.getSequences()) {
				this.displayRanges.add(new DisplayContig(ssr));
			}
			// sort on chrom/start
			Collections.sort(this.displayRanges, (DR1,DR2)->{
				final int tid1 = dict.getSequenceIndex(DR1.getContig());
				final int tid2 = dict.getSequenceIndex(DR2.getContig());
				int i=tid1-tid2;
				if(i!=0) return i;
				i = DR1.getStart() - DR2.getStart();
				if(i!=0) return i;
				i = DR1.getEnd() - DR2.getEnd();
				return i;
			});
			this.genomeViewLength = this.displayRanges.stream().mapToLong(DR->DR.length()).sum();
			double imgWidth = this.WIDTH - distanceBetweenContigs*(this.displayRanges.size()-1);
			
			final AtomicInteger countLowY = new AtomicInteger(0);
			final AtomicInteger countHighY = new AtomicInteger(0);
			
			r = super.openBufferedReader(oneFileOrNull(args));
			r.lines().
				filter(L->!(L.trim().isEmpty() || L.startsWith("#"))).
				map(L->parseData(L)).
				filter(D->D!=null).
				forEach(data->{
					if(data.getValue()< GenScan2.this.minY) {countLowY.incrementAndGet() ; return;}
					if(data.getValue()> GenScan2.this.maxY) {countHighY.incrementAndGet() ; return;}
					for(final DisplayRange dr: GenScan2.this.displayRanges) {
						if(!dr.contains(data)) continue;
						
						}
					
				});
			
			r.close();r=null;
			return 0;
			} 
		catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		finally {
			CloserUtil.close(r);
			this.dict =null;
			}
		}
	
	
	public static void main(String[] args) {
		new GenScan2().instanceMainWithExit(args);
	}

}
