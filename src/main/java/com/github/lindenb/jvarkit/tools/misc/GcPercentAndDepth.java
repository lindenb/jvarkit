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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceContig;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceGenome;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceGenomeFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

/**
BEGIN_DOC

## History

* 2017-12-09 : currently rewriting everything... do not use...

## Example

```bash
$ java -jar dist/gcanddepth.jar -R ref.fasta -b capture.bed 1.bam 2.bam ... > result.tsv
```
END_DOC
 */
@Program(name="gcpercentanddepth",
description="Extracts GC% and depth for multiple bam using a sliding window",
keywords={"gc%","depth","coverage"})
public class GcPercentAndDepth extends Launcher
	{
	private static Logger LOG=Logger.build(GcPercentAndDepth.class).make();

	
	@Parameter(names="-o",description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outPutFile=null;
	@Parameter(names="-w",description=" (window size)")
	private int windowSize=100;
	@Parameter(names="-s",description=" (window shift) ")
	private int windowStep=50;
	@Parameter(names="-m",description=" min depth ")
	private int min_depth=0;
	@Parameter(names="-R",description=ReferenceGenomeFactory.OPT_DESCRIPTION)
	private File refFile=null;
	@Parameter(names="-B",description=" (file.bed) (optional). If not defined: use whole genome. Warning memory consumming: must alloc sizeof(int)*win.size()*num(samples).")
	private File bedFile=null;
	@Parameter(names="-n",description=" skip window if Reference contains one 'N'.")
	private boolean skip_if_contains_N=false;
	@Parameter(names="-x",description=" don't print genomic index.")
	private boolean hide_genomic_index=false;
	@Parameter(names={"-filter","--filter"},description="[20171219]"+SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-partition","--partition"},description="[20171219]"+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	//@Parameter(names={"-percentile","--percentile"},description="[20171219] data percentile method")
	//private Percentile percentile = Percentile.average();
	
	
	private SAMSequenceDictionary samSequenceDictionary=null;

	
	/** A bed segment from the catpure */
	private class RegionCaptured
		implements Locatable,
			Comparable<RegionCaptured>,Iterable<RegionCaptured.SlidingWindow>
		{
		final SAMSequenceRecord ssr;
		final int start0;
		final int end0;
		
		RegionCaptured(final SAMSequenceRecord ssr,int start0,int end0)
			{
			this.ssr = ssr;
			this.start0 = start0;
			this.end0 = end0;
			}
		
		public SAMSequenceRecord getSAMSequenceRecord()
			{
			return this.ssr;
			}
		
		@Override
		public String getContig()
			{
			return getSAMSequenceRecord().getSequenceName();
			}
		
		@Override
		public int getStart()
			{
			return this.start0+1;
			}
		@Override
		public int getEnd()
			{
			return this.end0+1;
			}

		public int length()
			{
			return this.getEnd()-this.getStart();
			}
		
		@Override
		public int compareTo(final RegionCaptured o)
			{
			int i=this.ssr.getSequenceIndex()-o.ssr.getSequenceIndex();
			if(i!=0) return i;
			i=start0-o.start0;
			if(i!=0) return i;
			i=end0-o.end0;
			if(i!=0) return i;
			return 0;
			}
		
		@Override
		public String toString() {
			return "roi:"+getContig()+":"+getStart()+"-"+getEnd();
			}
		
		Iterator<SlidingWindow> windows()
			{
			return new MyIter();
			}
		
		@Override
		public Iterator<SlidingWindow> iterator()
			{
			return windows();
			}
		
		private class MyIter
			extends AbstractIterator<SlidingWindow>
			{
			private int index_in_roi=0;
			
			@Override
			protected SlidingWindow advance() {
				final SlidingWindow sw =  new SlidingWindow(index_in_roi);
				if(!sw.isValid()) return null;
				++index_in_roi;
				return sw;
				}
			}
		
		/** a Sliding window from the RegionCaptured */
		public class SlidingWindow implements Locatable
			{
			int index_in_roi;
			private SlidingWindow(int index_in_roi)
				{
				this.index_in_roi=index_in_roi;
				}
			
			public long getGenomicIndex()
				{
				long n=0;
				for(final SAMSequenceRecord x: GcPercentAndDepth.this.samSequenceDictionary.getSequences())
					{
					if(x.getSequenceName().equals(this.getContig()))
						{
						break;
						}
					n+= x.getSequenceLength();
					}
				n+= this.getStart();
				return n;
				}
			@Override
			public String getContig()
				{
				return RegionCaptured.this.getContig();
				}
			@Override
			public int getStart()
				{
				return 1 + index_in_roi*windowStep + RegionCaptured.this.start0;
				}
			@Override
			public int getEnd()
				{
				return Math.min(
					getStart()+windowSize ,
					RegionCaptured.this.getEnd()
					);
				}
			public int length()
				{
				return getEnd()-getStart();
				}
			boolean isValid()
				{
				if(getStart() > RegionCaptured.this.getEnd())
					{
					return false;
					}
				if(getStart()+windowSize <=   RegionCaptured.this.getEnd())
					{
					return true;
					}
				/** avoid side effect */
				if(getStart()+windowSize >= RegionCaptured.this.getSAMSequenceRecord().getSequenceLength())
					{
					return false;
					}
				int overhang=getStart()+windowSize-RegionCaptured.this.getEnd();
				return overhang>windowSize/2;
				}
			@Override
			public String toString() {
				return "win:"+getContig()+":"+getStart()+"-"+getEnd()+" "+this.index_in_roi;
				}
			}
		}
	
	
	/** constructor */
	public GcPercentAndDepth()
		{
		}
	
	
	@Override
	public int doWork(final List<String> args)
		{
		
		if(this.windowSize<=0)
			{
			LOG.error("Bad window size.");
			return -1;
			}
		if(this.windowStep<=0)
			{
			LOG.error("Bad window step.");
			return -1;
			}
		
		if(this.refFile==null)
			{
			LOG.error("Undefined REF File");
			return -1;
			}
		
		if(args.isEmpty())
			{
			LOG.error("Illegal Number of arguments.");
			return -1;
			}
		ReferenceGenome indexedFastaSequenceFile= null;
				
		List<SamReader> readers=new ArrayList<SamReader>();

		PrintWriter out=null;
		try
			{
			LOG.info("Loading "+this.refFile);
			indexedFastaSequenceFile=	new ReferenceGenomeFactory().
					openFastaFile(this.refFile);
			this.samSequenceDictionary = indexedFastaSequenceFile.getDictionary();
			if(this.samSequenceDictionary==null)
				{
				LOG.error("Cannot get sequence dictionary for "+this.refFile);
				return -1;
				}
			
			out= super.openFileOrStdoutAsPrintWriter(outPutFile);
			
			Set<String> all_samples=new TreeSet<String>();
			/* create input, collect sample names */
			for(int optind=0;
					optind< args.size();
					++optind)
				{
				LOG.info("Opening "+args.get(optind));
				final SamReader samFileReaderScan= super.openSamReader(args.get(optind));
				readers.add(samFileReaderScan);
				
				final SAMFileHeader header= samFileReaderScan.getFileHeader();
				if(!SequenceUtil.areSequenceDictionariesEqual(this.samSequenceDictionary, header.getSequenceDictionary()))
					{
					LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(this.samSequenceDictionary, header.getSequenceDictionary()));
					return -1;
					}
				
				for(final SAMReadGroupRecord g: header.getReadGroups())
					{
					final String sample = this.partition.apply(g);
					if(StringUtil.isBlank(sample))
						{
						LOG.warning("Read group "+g.getId()+" has no sample in merged dictionary");
						continue;
						}
					all_samples.add(sample);
					}
				}
			
			
			LOG.info("N "+this.partition.name()+"="+all_samples.size());
			
			/* print header */
			out.print("#");
			if( !this.hide_genomic_index)
				{
				out.print("id");
				out.print("\t");
				}
			out.print("chrom");
			out.print("\t");
			out.print("start");
			out.print("\t");
			out.print("end");
			out.print("\t");
			out.print("GCPercent");
			for(final String sample:all_samples)
				{
				out.print("\t");
				out.print(sample);
				}
			out.println();
			
			
			final List<RegionCaptured> regionsCaptured=new ArrayList<RegionCaptured>();
			if(bedFile!=null)
				{
				LOG.info("Reading BED:" +bedFile);
				final BedLineCodec bedLineCodec =new BedLineCodec();
				
				BufferedReader r=IOUtils.openFileForBufferedReading(bedFile);
				r.lines().
					filter(L->!L.startsWith("#")).
					filter(L->!StringUtil.isBlank(L)).
					map(L->bedLineCodec.decode(L)).
					filter(B->B!=null).
					forEach(B->{
						final SAMSequenceRecord ssr= this.samSequenceDictionary.getSequence(B.getContig());
						if(ssr==null)
							{
							LOG.warning("Cannot resolve "+B.getContig());
							return;
							}
							
						final RegionCaptured roi=new RegionCaptured(
								ssr,
								B.getStart()-1,
								B.getEnd()
								);
						regionsCaptured.add(roi);
						});
				CloserUtil.close(r);
				LOG.info("end Reading BED:" +bedFile);
				Collections.sort(regionsCaptured);
				}
			else
				{
				LOG.info("No capture, peeking everything");
				for(final SAMSequenceRecord ssr:this.samSequenceDictionary.getSequences())
					{
					final RegionCaptured roi=new RegionCaptured(ssr,0,ssr.getSequenceLength());
					regionsCaptured.add(roi);
					}
				}
			
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(this.samSequenceDictionary).logger(LOG);
			ReferenceContig genomicSequence=null;
			for(final RegionCaptured roi:regionsCaptured)
				{
				if(genomicSequence==null || !genomicSequence.hasName(roi.getContig()))
					{
					genomicSequence= indexedFastaSequenceFile.getContig(roi.getContig());
					if(genomicSequence==null)
						{
						LOG.error(JvarkitException.ContigNotFoundInDictionary.getMessage(roi.getContig(), this.samSequenceDictionary));
						return -1;
						}
					}
				Map<String,int[]> sample2depth=new HashMap<String,int[]>();
				Map<String,Double> sample2meanDepth=new HashMap<String,Double>();
				for(final String sample:all_samples)
					{
					int depth[]=new int[roi.length()];
					Arrays.fill(depth, 0);
					sample2depth.put(sample, depth);
					}
				List<CloseableIterator<SAMRecord>> iterators = new ArrayList<CloseableIterator<SAMRecord>>();
				for(final SamReader r:readers)
					{
					iterators.add(r.query(roi.getContig(), roi.getStart(), roi.getEnd(), false));
					}
				
				final MergingIterator<SAMRecord> merginIter=
						new MergingIterator<>(new SAMRecordCoordinateComparator(),iterators);
					
				while(merginIter.hasNext())
					{
					final SAMRecord rec=merginIter.next();
					if(rec.getReadUnmappedFlag()) continue;
					if(this.filter.filterOut(rec)) continue;
										
					final String sample= this.partition.getPartion(rec,null);
					if(sample==null ) continue;
					
					final int depth[]=sample2depth.get(sample);
					if(depth==null) continue;
					final Cigar cigar=rec.getCigar();
					if(cigar==null) continue;
					
					int refpos1=rec.getAlignmentStart();
					
					for(final CigarElement ce: cigar.getCigarElements())
						{
						final CigarOperator op = ce.getOperator();
						if(!op.consumesReferenceBases() ) continue;
						if(op.consumesReadBases())
							{
							for(int i=0;i< ce.getLength();++i)
								{
								if(refpos1+i < roi.getStart()) continue;
								if(refpos1+i > roi.getEnd()) break;
								depth[refpos1+i-roi.getStart()]++;
								}
							}
						refpos1 += ce.getLength();
						}
					}
				merginIter.close();
				
				for(final RegionCaptured.SlidingWindow win: roi)
					{
					double total=0f;
					int countN=0;
					for(int pos1=win.getStart();pos1<=win.getEnd();++pos1)
						{
						switch(genomicSequence.charAt(pos1-1))
							{
							case 'c':case 'C':
							case 'g':case 'G':		
							case 's':case 'S':
								{
								total++;
								break;
								}
							case 'n':case 'N':countN++;break;
							default:break;
							}
						}
					if(skip_if_contains_N && countN>0) continue;
 					double GCPercent=total/(double)win.length();
					
					int max_depth_for_win=0;
					sample2meanDepth.clear();
					for(final String sample:all_samples)
						{
						int depth[]=sample2depth.get(sample);
						double sum=0;
						for(int pos=win.getStart();
								pos<win.getEnd() && (pos-roi.getStart())< depth.length;
								++pos)
							{
							sum+=depth[pos-roi.getStart()];
							}		
						double mean= (sum/(double)depth.length);
						max_depth_for_win=Math.max(max_depth_for_win, (int)mean);
						sample2meanDepth.put(sample,mean);
						}
					if(max_depth_for_win< this.min_depth) continue;
					if(!this.hide_genomic_index)
						{
						out.print(win.getGenomicIndex());
						out.print("\t");
						}
					out.print(win.getContig());
					out.print("\t");
					out.print(win.getStart()-1);
					out.print("\t");
					out.print(win.getEnd());
					out.print("\t");
					out.printf("%.2f",GCPercent);
					
					for(String sample:all_samples)
						{
						out.print("\t");
						out.printf("%.2f",(double)sample2meanDepth.get(sample));
						}
					out.println();
					}
				}
			progress.finish();
			out.flush();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			for(SamReader r:readers) CloserUtil.close(r);
			CloserUtil.close(indexedFastaSequenceFile);
			CloserUtil.close(out);
			}	
		}
	
	public static void main(final String[] args) {
		new GcPercentAndDepth().instanceMainWithExit(args);
		}
	}
