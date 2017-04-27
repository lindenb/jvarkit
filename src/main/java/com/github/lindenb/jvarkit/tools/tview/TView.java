package com.github.lindenb.jvarkit.tools.tview;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
//import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


public class TView implements Closeable
	{
	private static final Logger LOG = Logger.build(TView.class).make();
    private final static int INSERTION_OPCODE = -1; 
	private boolean showClip=false;
	private Interval interval=null;
	@Parameter(names={"-R","--reference"},description="Indexed Reference file.")
	private File referenceFile=null;

	
	private int distance_between_reads=1;
	
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;

	
	public TView() {
		
	}
	
	public int initialize() throws IOException
		{
		if(this.referenceFile!=null) {
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.referenceFile);
			}
		
		return 0;
		}
	
	@Override
	public void close() {
		CloserUtil.close(this.indexedFastaSequenceFile);
		this.indexedFastaSequenceFile=null;
		}
	
	public Predicate<SAMRecord> getRecordFilter() {
		return R->true;
	}
	
	public Function<SAMRecord,String> getRecordGroup() {
		return R->"ALL";
	}
	
	
	public Function<SAMRecord,Integer> left() {
		return R->showClip?
			R.getUnclippedStart():
			R.getAlignmentStart();
	}
	public Function<SAMRecord,Integer> right() {
		return R->showClip?
			R.getUnclippedEnd():
			R.getAlignmentEnd();
	}

	
	void paint(PrintStream out,List<String> args) {
		GenomicSequence genomicSequence=new GenomicSequence(indexedFastaSequenceFile, interval.getContig());

		
		final Map<String, List<SAMRecord>> group2record=new TreeMap<>();
		
		final SamReaderFactory srf=SamReaderFactory.makeDefault().
				referenceSequence(this.referenceFile).
				validationStringency(ValidationStringency.LENIENT)
				; 
		for(String s:IOUtils.unrollFiles(args))
			{
			SamReader samReader=  null;
			samReader= srf.open(SamInputResource.of(s));
			SAMRecordIterator iter = samReader.query(
					this.interval.getContig(),
					this.interval.getStart(),
					this.interval.getEnd(),
					false);
			while(iter.hasNext())
				{
				final SAMRecord rec = iter.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(!getRecordFilter().test(rec)) continue;
				if( !rec.getContig().equals(interval.getContig())) continue;
				if(right().apply(rec) < this.interval.getStart()) continue;
				if(this.interval.getEnd() < left().apply(rec) ) continue;
				final String group=getRecordGroup().apply(rec);
				if(group==null || group.isEmpty()) continue;
				List<SAMRecord> records = group2record.get(group);
				if( records == null) {
					records = new ArrayList<>();
					group2record.put(group,records);
					}
				records.add(rec);
				}
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			}
		final Map<Integer,Integer> posToInsertLength =new HashMap<>(this.interval.length());
		for(int i=interval.getStart();i<=interval.getEnd();++i)
			{
			posToInsertLength.put(i,0);
			}
		
		final Predicate<Integer> inInterval = new Predicate<Integer>() {
			@Override
			public boolean test(Integer pos) {
				return interval.getStart() <= pos && pos <= interval.getEnd();
			}
		};
		for(final String groupName : group2record.keySet())
			{
			for(final SAMRecord rec: group2record.get(groupName)) {
				int ref = rec.getAlignmentStart();
				for(CigarElement ce:rec.getCigar()) {
					final CigarOperator op = ce.getOperator();
					if(op.equals(CigarOperator.I) && inInterval.test(ref)) 
						{
						Integer maxL=posToInsertLength.get(ref);
						if(maxL==null) maxL=0;
						posToInsertLength.put(ref,Math.max(maxL,ce.getLength()));
						}					
					if(op.consumesReferenceBases())
						{
						ref += ce.getLength();
						}
					}
				}
			}
		int ref= this.interval.getStart();
		final int pixelToRef[]=new int[this.interval.length()];
		for(int i=0;i< pixelToRef.length ;++i )
			{
			pixelToRef[i] = ref ;
			}
		
		final Function<Integer,Integer> position2column = new Function<Integer, Integer>() {
			@Override
			public Integer apply(Integer pos) {
				if(!inInterval.test(pos)) return -1;
				for(int x=0;x < pixelToRef.length;++x) {
					int p= pixelToRef[x];
					if(p==INSERTION_OPCODE) continue;
					if(p==pos) return x;
					if(p>pos) return -1;
				}
				return -1;
			}
		};

		/* paint reference */
		 ref = this.interval.getStart();
		int x=0;
		while(x< pixelToRef.length)
		{
			if(ref%15==0) {
			String s= String.format("%d", ref);
			for(int i=0;i< s.length() && x< pixelToRef.length;++x) {
				out.print(s.charAt(i));
				x++;
				ref++;
				}
			}
			else
			{
			out.print(" ");
			++x;
			++ref;
			}
		}
		out.println();
		
		
		/* paint base */
		ref = this.interval.getStart();
		x=0;
		while(x< pixelToRef.length)
			{
			out.print(genomicSequence.charAt(ref-1));
			++x;
			}
		out.println();
		
		
		for(final String groupName : group2record.keySet())
			{
			final List<List<SAMRecord>> rows = new ArrayList<>();
			for(final SAMRecord rec: group2record.get(groupName)) {
				int y=0;
				for(y=0;y< rows.size();++y)
					{
					final List<SAMRecord> row = rows.get(y);
					final SAMRecord last = row.get(row.size()-1);
					if(right().apply(last)+1 < left().apply(rec)) {
						row.add(rec);
						break;
						}
					}
				if(y==rows.size()) {
					List<SAMRecord> row = new ArrayList<>();
					row.add(rec);
					rows.add(row);
					}
				}
			for(List<SAMRecord> row : rows)
				{
				x=0;
				ref=interval.getStart();
				while(x < pixelToRef.length && ref < interval.getEnd())
					{
					char c=' ';
					++x;
					}
				out.println();
				}
			
			}
		
		}
	
    
    private int left(final SAMRecord rec) {
    	return(showClip?rec.getUnclippedStart():rec.getAlignmentStart());
    }
    
    private int right(final SAMRecord rec) {
    	return(showClip?rec.getUnclippedEnd():rec.getAlignmentEnd());
    }

    
	}
