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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.bamstats04;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;

/**
BEGIN_DOC

## Example

```

$ head genes.bed
1	179655424	179655582	ZORG
1	179656788	179656934	ZORG

$ java -jar  dist/bamstats05.jar -B genes.bed --mincoverage 10  > out.txt

$ head out.txt
#chrom	start	end	gene	sample	length	mincov	maxcov	avg	nocoverage.bp	percentcovered
1	179655424	179656934	ZORG	SAMPLE1	304	27	405	216.80921052631578	0	100

```

END_DOC


*/
@Program(name="bamstats05",
description="Coverage statistics for a BED file, group by gene",
keywords={"bam","coverage","statistics","bed"}
)
public class BamStats05 extends Launcher
	{
	private static final Logger LOG = Logger.build(BamStats05.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-m","--mincoverage"},description="min coverage to say the position is not covered")
	private int MIN_COVERAGE = 0 ;

	@Parameter(names={"-B","--bed"},description="bed file (columns: chrom start end GENE)",required=true)
	private File BEDILE = null;
	
	@Parameter(names={"--groupby"},description="Group Reads by")
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;

	@Parameter(names={"-f","--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter filter  = SamFilterParser.buildDefault();
	
	private Map<String, List<Interval>> readBedFile(final File bedFile) throws IOException
    	{
    	final Map<String, List<Interval>> gene2interval=new TreeMap<String, List<Interval>>();
    	
    	BufferedReader bedIn=null;
    	try
    		{
    		bedIn=IOUtils.openFileForBufferedReading(bedFile);
    		final BedLineCodec codec = new BedLineCodec();
    		String line=null;
			while((line=bedIn.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				final BedLine bedLine = codec.decode(line);
				if(bedLine==null) continue;
				if(bedLine.getColumnCount()<4)
					{
					throw new IOException("bad bed line in "+line+" "+bedFile);
					}
				final String chrom = bedLine.getContig();
				final int chromStart1= bedLine.getStart();
				final int chromEnd1= bedLine.getEnd();
				final String gene = bedLine.get(3);
				if(gene.isEmpty())  throw new IOException("bad bed gene in "+line+" "+bedFile);
				 List<Interval> intervals = gene2interval.get(gene);
				 if(intervals==null)
				 	{
					 intervals=new ArrayList<>();
					 gene2interval.put(gene,intervals);
				 	} 
				 else if(!intervals.get(0).getContig().equals(chrom))
				 	{
					throw new IOException("more than one chromosome for gene:"+gene);
				 	}
				else
				 	{
					for(final Interval interval:intervals)
						{
						if(interval.getEnd()<chromStart1) continue;
						if(interval.getStart()>chromEnd1) continue;
						throw new IOException("overlapping region: "+line+" and "+interval);
						}
				 	}
				intervals.add(new Interval(chrom, chromStart1, chromEnd1));
				}
    		bedIn.close();
    		return gene2interval;
    		}
    	finally
    		{
    		CloserUtil.close(bedIn);
    		}
    	}
	
	protected  int doWork(
			final PrintWriter pw,
			final Map<String, List<Interval>> gene2interval,
			final String filename,
			final SamReader IN) throws Exception
		{
		
		try
			{
			LOG.info("Scanning "+filename);
			final SAMFileHeader header = IN.getFileHeader();
			final List<SAMReadGroupRecord> rgs = header.getReadGroups();
			if(rgs==null || rgs.isEmpty())
				throw new IOException("No read groups in "+filename);
			final Set<String> groupNames = this.groupBy.getPartitions(rgs);
			
			for(final String partition : groupNames)
				{
				if(partition.isEmpty()) throw new IOException("Empty "+groupBy.name());
				for(final String gene: gene2interval.keySet())
					{
					int geneStart = Integer.MAX_VALUE;
					int geneEnd = 0;
					
					final List<Integer> counts = new ArrayList<>();
					
					
					for(final Interval interval:gene2interval.get(gene))
						{
						geneStart = Math.min(geneStart, interval.getStart()-1);
						geneEnd = Math.max(geneEnd, interval.getEnd());
		
						/* picard javadoc:  - Sequence name - Start position (1-based) - End position (1-based, end inclusive)  */
						int interval_counts[]=new int[interval.getEnd()-interval.getStart()+1];
						if(interval_counts.length==0) continue;
						Arrays.fill(interval_counts, 0);
						
						if(IN.getFileHeader().getSequenceIndex(interval.getContig())==-1)
							{
							throw new IllegalArgumentException("NO DICT FOR \""+interval.getContig()+"\"");
							}
						
						/**
						 *     start - 1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
			    		*	   end - 1-based, inclusive end of interval of interest. Zero implies end of the reference sequence. 
						 */
					
						SAMRecordIterator r=IN.query(new QueryInterval[]{
								new QueryInterval(
								header.getSequenceIndex(interval.getContig()),
								interval.getStart(),
								interval.getEnd()
								)},false)
								;
						while(r.hasNext())
							{
							final SAMRecord rec=r.next();
							if(rec.getReadUnmappedFlag()) continue;
							if(filter.filterOut(rec)) continue;
							
							if(!rec.getReferenceName().equals(interval.getContig())) continue;
							
							final SAMReadGroupRecord rg = rec.getReadGroup();
							if(rg==null || !partition.equals(this.groupBy.apply(rg))) continue;
							final Cigar cigar=rec.getCigar();
							if(cigar==null) continue;
				    		int refpos1=rec.getAlignmentStart();
				    		for(final CigarElement ce:cigar.getCigarElements())
				    			{
				    			final CigarOperator op=ce.getOperator();
				    			if(!op.consumesReferenceBases()) continue;
				    			if(op.consumesReadBases())
				    				{
				    				for(int i=0;i< ce.getLength();++i)
			    		    			{
										if(refpos1+i>= interval.getStart() && refpos1+i<=interval.getEnd())
											{
											interval_counts[refpos1+i-interval.getStart()]++;
											}
					    				}
				    				}
				    			refpos1+=ce.getLength();
				    			}
							}/* end while r */
						r.close();
						for(int d: interval_counts)
							{
							counts.add(d);
							}
						}/* end interval */
						
						
						Collections.sort(counts);
						
						int count_no_coverage=0;
						double mean=0;
						for(int cov:counts)
							{
							if(cov<=MIN_COVERAGE) ++count_no_coverage;
							mean+=cov;
							}
						mean/=counts.size();
						
					pw.println(
							gene2interval.get(gene).get(0).getContig()+"\t"+
							geneStart+"\t"+geneEnd+"\t"+gene+"\t"+partition+"\t"+
							counts.size()+"\t"+
							counts.get(0)+"\t"+
							counts.get(counts.size()-1)+"\t"+
							mean+"\t"+
							count_no_coverage+"\t"+
							(int)(((counts.size()-count_no_coverage)/(double)counts.size())*100.0)
							);
					}//end gene
				}//end sample
		return RETURN_OK;
		}
	catch(final Exception err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(IN);
		}
	}
	
	@Override
	public int doWork(List<String> args) {
		if(BEDILE==null)
			{
			LOG.error( "missing bed file");
			return -1;
			}
		SamReader in=null;
		BufferedReader r=null;
		PrintWriter pw=null;
		try
			{
			Map<String, List<Interval>> gene2interval = readBedFile(BEDILE);
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			pw.println("#chrom\tstart\tend\tgene\tsample\tlength\tmincov\tmaxcov\tavg\tnocoverage.bp\tpercentcovered");
			final SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			final Set<String> files = new  HashSet<>();
			if(args.isEmpty())
				{
				LOG.info("reading BAM path from stdin");
				r = new BufferedReader(new InputStreamReader(stdin()));
				String line;
				while((line=r.readLine())!=null)
					{
					if(line.startsWith("#") || line.trim().isEmpty()) continue;
					if(!line.endsWith(".bam"))
						{
						LOG.error("line should end with .bam :"+line);
						return -1;
						}
					}
				CloserUtil.close(r);
				}
			
			else
				{
				for(final String fname:args)
					{
					files.addAll(IOUtils.unrollFiles(
							Collections.singletonList(fname)
							));
					}
				}
			
			for(final String f:files)
				{
				in = srf.open(new File(f));
				int tl =doWork(pw,gene2interval,f,in);
				CloserUtil.close(in);
				in=null;
				if(tl!=0) return tl;
				}
			pw.flush();
			pw.close();
			pw=null;
			return RETURN_OK;
			}
		catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(r);
			CloserUtil.close(pw);
			}
		}

	
	public static void main(String[] args) throws IOException
		{
		new BamStats05().instanceMainWithExit(args);
		}

	}

