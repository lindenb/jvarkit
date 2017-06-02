/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import java.util.regex.Pattern;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.MergingSamRecordIterator;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

/**
BEGIN_DOC
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
	private SAMSequenceDictionary firstSamDict=null;
	@Parameter(names="-R",description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File refFile=null;
	@Parameter(names="-B",description=" (file.bed) (optional). If not defined: use whole genome. Warning memory consumming: must alloc sizeof(int)*win.size()*num(samples).")
	private File bedFile=null;
	@Parameter(names="-n",description=" skip window if Reference contains one 'N'.")
	private boolean skip_if_contains_N=false;
	@Parameter(names="-N",description=" (file) . chrom.name.helper .")
	private String chromNameFile=null;
	@Parameter(names="-x",description=" don't print genomic index.")
	private boolean hide_genomic_index=false;
	@Parameter(names={"-filter","--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter filter  = SamFilterParser.buildDefault();

	
	/** A bed segment from the catpure */
	private class RegionCaptured
		implements Comparable<RegionCaptured>,Iterable<RegionCaptured.SlidingWindow>
		{
		int tid;
		int start0;
		int end0;
		
		public SAMSequenceRecord getSAMSequenceRecord()
			{
			return firstSamDict.getSequence(this.tid);
			}
		
		public String getChromosome()
			{
			return getSAMSequenceRecord().getSequenceName();
			}
		
		public int getStart()
			{
			return this.start0;
			}
		public int getEnd()
			{
			return this.end0;
			}

		public int length()
			{
			return this.getEnd()-this.getStart();
			}
		
		@Override
		public int compareTo(RegionCaptured o)
			{
			int i=tid-o.tid;
			if(i!=0) return i;
			i=start0-o.start0;
			if(i!=0) return i;
			i=end0-o.end0;
			if(i!=0) return i;
			return 0;
			}
		
		@Override
		public String toString() {
			return "roi:"+getChromosome()+":"+getStart()+"-"+getEnd();
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
			implements Iterator<SlidingWindow>
			{
			int index_in_roi=0;
			private SlidingWindow make()
				{
				return new SlidingWindow(index_in_roi);
				}
			@Override
			public boolean hasNext()
				{
				return make().isValid();
				}
			@Override
			public SlidingWindow next()
				{
				SlidingWindow w=make();
				if(!w.isValid()) throw new IllegalStateException();
				this.index_in_roi++;
				return w;
				}
			
			@Override
			public void remove()
				{
				throw new UnsupportedOperationException();
				}
			}
		
		/** a Sliding window from the RegionCaptured */
		public class SlidingWindow
			{
			int index_in_roi;
			private SlidingWindow(int index_in_roi)
				{
				this.index_in_roi=index_in_roi;
				}
			
			public long getGenomicIndex()
				{
				long n=0;
				for(int t=0;t< RegionCaptured.this.tid;++t)
					{
					n+= firstSamDict.getSequence(t).getSequenceLength();
					}
				n+= this.getStart();
				return n;
				}
			
			public String getChromosome()
				{
				return RegionCaptured.this.getChromosome();
				}
			public int getStart()
				{
				return index_in_roi*windowStep + RegionCaptured.this.start0;
				}
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
				if(getStart() >= RegionCaptured.this.getEnd())
					{
					return false;
					}
				if(getStart()+windowSize<   RegionCaptured.this.getEnd())
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
				return "win:"+getChromosome()+":"+getStart()+"-"+getEnd()+" "+this.index_in_roi;
				}
			}
		}
	
	
	/** constructor */
	private GcPercentAndDepth()
		{
		}
	
	
	
	@Override
	public int doWork(List<String> args)
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
		
		if(refFile==null)
			{
			LOG.error("Undefined REF file");
			return -1;
			}
		
		if(args.isEmpty())
			{
			LOG.error("Illegal Number of arguments.");
			return -1;
			}
		IndexedFastaSequenceFile indexedFastaSequenceFile=null;
		List<SamReader> readers=new ArrayList<SamReader>();

		PrintWriter out=null;
		try
			{
			Map<String,String> resolveChromName=new HashMap<String, String>();

			LOG.info("Loading "+refFile);
			indexedFastaSequenceFile=new IndexedFastaSequenceFile(refFile);
			SAMSequenceDictionary dict= indexedFastaSequenceFile.getSequenceDictionary();
			if(dict==null)
				{
				LOG.error("Cannot get sequence dictionary for "+refFile+". "+
						"picard.dictionary.needed");
				return -1;
				}
			/* load chrom aliases */
			if(chromNameFile!=null)
				{
				LOG.info("Reading "+chromNameFile);
				LineIterator r=IOUtils.openURIForLineIterator(chromNameFile);
				while(r.hasNext())
					{
					String line=r.next();
					if(line.isEmpty()) continue;
					String tokens[]=line.split("[\t]");
					if(tokens.length<2)
						{
						LOG.warning("Bad conversion line (column count<2) in "+line);
						continue;
						}
					resolveChromName.put(tokens[0], tokens[1]);
					resolveChromName.put(tokens[1], tokens[0]);
					}
				CloserUtil.close(r);
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
				
				SAMFileHeader header= samFileReaderScan.getFileHeader();
				if(firstSamDict==null)
					{
					firstSamDict=header.getSequenceDictionary();
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(firstSamDict, header.getSequenceDictionary()))
					{
					LOG.error("not.the.same.sequence.dictionaries");
					return -1;
					}
				
				for(SAMReadGroupRecord g: header.getReadGroups())
					{
					if(g.getSample()==null || g.getSample().isEmpty())
						{
						LOG.warning("Read group "+g.getId()+" has no sample in merged dictionary");
						continue;
						}
					all_samples.add(g.getSample());
					}
				}
			
			LOG.info("NSample:"+all_samples.size());
			
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
			for(String sample:all_samples)
				{
				out.print("\t");
				out.print(sample);
				}
			out.println();
			
			
			List<RegionCaptured> regionsCaptured=new ArrayList<RegionCaptured>();
			if(bedFile!=null)
				{
				int N=0;
				LOG.info("Reading BED:" +bedFile);
				Pattern tab=Pattern.compile("[\t]");
				LineIterator r=IOUtils.openFileForLineIterator(bedFile);
				while(r.hasNext())
					{
					String line=r.next();
					if(line.startsWith("#") || line.startsWith("browser") || line.startsWith("track") || line.isEmpty()) continue;
					String tokens[]=tab.split(line,4);
					if(tokens.length<3)
						{
						LOG.warning("No enough column in "+line+" "+bedFile);
						continue;
						}
					RegionCaptured roi=new RegionCaptured();
					roi.tid=dict.getSequenceIndex(tokens[0]);
					if(roi.tid==-1)
						{
						String altName= resolveChromName.get(tokens[0]);
						if(altName!=null)
								{
							LOG.info(tokens[0]+" was resolved to "+altName);
								roi.tid=firstSamDict.getSequenceIndex(altName);
								}
						else
								{
							LOG.warning("Cannot resolve "+tokens[0]+ " from "+new ArrayList<String>(resolveChromName.keySet()));
								}
						}
					
					
					if(roi.tid==-1)
						{
						LOG.warning("not in reference: chromosome "+tokens[0]+" in \""+line+"\" "+bedFile);
						continue;
						}
					roi.start0=Integer.parseInt(tokens[1]);
					roi.end0=Integer.parseInt(tokens[2]);
					regionsCaptured.add(roi);
					++N;
					}
				CloserUtil.close(r);
				LOG.info("end Reading BED:" +bedFile+"  N="+N);
				Collections.sort(regionsCaptured);
				}
			else
				{
				LOG.info("No capture, peeking everything");
				for(SAMSequenceRecord ssr:this.firstSamDict.getSequences())
					{
					RegionCaptured roi=new RegionCaptured();
					roi.tid=ssr.getSequenceIndex();
					roi.start0=0;
					roi.end0=ssr.getSequenceLength();
					regionsCaptured.add(roi);
					}
				}
			
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(firstSamDict);
			GenomicSequence genomicSequence=null;
			for(RegionCaptured roi:regionsCaptured)
				{
				if(out.checkError()) break;
				String chromName=roi.getChromosome();
				
				if(genomicSequence==null || !genomicSequence.getChrom().equals(chromName))
					{
					String altChromName=resolveChromName.get(chromName);
					if(altChromName!=null && genomicSequence!=null && genomicSequence.getChrom().equals(altChromName))
						{
						//ok, good chrom
						}
					else if(indexedFastaSequenceFile.getSequenceDictionary().getSequence(chromName)!=null)
						{
						genomicSequence=new GenomicSequence(indexedFastaSequenceFile,chromName);
						}
					else if(altChromName!=null  && indexedFastaSequenceFile.getSequenceDictionary().getSequence(altChromName)!=null)
						{
						genomicSequence=new GenomicSequence(indexedFastaSequenceFile,altChromName);
						}
					else
						{
						
						LOG.error(
								"when looking for genomic sequence "+ "chrom.missing.in.sequence.dictionary" +" "+
										chromName+
								" conversions available:"+resolveChromName);
						return -1;
						}
					}
				Map<String,int[]> sample2depth=new HashMap<String,int[]>();
				Map<String,Double> sample2meanDepth=new HashMap<String,Double>();
				for(String sample:all_samples)
					{
					int depth[]=new int[roi.length()];
					Arrays.fill(depth, 0);
					sample2depth.put(sample, depth);
					}
				List<CloseableIterator<SAMRecord>> iterators = new ArrayList<CloseableIterator<SAMRecord>>();
				for(SamReader r:readers)
					{
					iterators.add(r.query(roi.getChromosome(), roi.getStart()+1, roi.getEnd(), false));
					}
				
				Iterator<SAMRecord> merginIter=null;
				
				merginIter=null;
				if(iterators.isEmpty())
					{
					merginIter=Collections.emptyIterator();
					}
				else if(iterators.size()==1)
					{
					merginIter=iterators.get(0);
					}
				else
					{
					merginIter=new MergingSamRecordIterator(iterators);
					}
				while(merginIter.hasNext())
					{
					final SAMRecord rec=merginIter.next();
					if(rec.getReadUnmappedFlag()) continue;
					if(this.filter.filterOut(rec)) continue;
										
					SAMReadGroupRecord g=rec.getReadGroup();
					if(g==null) continue;
					progress.watch(rec);
					String sample=g.getSample();
					if(sample==null ) continue;
					int depth[]=sample2depth.get(sample);
					if(depth==null) continue;
					Cigar cigar=rec.getCigar();
					if(cigar==null) continue;
					
					int refpos1=rec.getAlignmentStart();
					
					for(CigarElement ce: cigar.getCigarElements())
						{
						CigarOperator op = ce.getOperator();
						if(!op.consumesReferenceBases() ) continue;
						if(op.consumesReadBases())
							{
							for(int i=0;i< ce.getLength();++i)
								{
								int pos0 = (refpos1+i-1);
								if(pos0< roi.getStart()) continue;
								if(pos0>= roi.getEnd()) continue;
								depth[pos0-roi.getStart()]++;
								}
							}
						refpos1 += ce.getLength();
						}
					}
				CloserUtil.close(merginIter);
				
				for(RegionCaptured.SlidingWindow win: roi)
					{
					double total=0f;
					int countN=0;
					for(int pos=win.getStart();pos<win.getEnd();++pos)
						{
						switch(genomicSequence.charAt(pos))
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
					for(String sample:all_samples)
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
					out.print(win.getChromosome());
					out.print("\t");
					out.print(win.getStart());
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
	
	public static void main(String[] args) {
		new GcPercentAndDepth().instanceMainWithExit(args);
		}
	}
