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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BamStats05 extends AbstractCommandLineProgram
	{
	private File BEDILE=null;
	private boolean PROPERLY_PAIRED=true;
	private int MMQ=0;
	private int MIN_COVERAGE=0;
	
    @SuppressWarnings("resource")
	private Map<String, List<Interval>> readBedFile(File bedFile) throws IOException
    	{
    	final Map<String, List<Interval>> gene2interval=new TreeMap<String, List<Interval>>();
    	
    	BufferedReader bedIn=null;
    	Pattern tab=Pattern.compile("[\t]");
    	String tokens[];
    	try
    		{
    		bedIn=IOUtils.openFileForBufferedReading(bedFile);
    		String line=null;
			while((line=bedIn.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				tokens=tab.split(line,6);
				if(tokens.length<4)
					{
					throw new IOException("bad bed line in "+line+" "+bedFile);
					}
				String chrom=tokens[0];
				int chromStart1= 1+Integer.parseInt(tokens[1]);//add one
				int chromEnd1= Integer.parseInt(tokens[2]);
				String gene = tokens[3];
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
					for(Interval interval:intervals)
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
	
	protected void doWork(
			final Map<String, List<Interval>> gene2interval,
			final String filename,
			final SamReader IN) throws Exception
		{
		try
			{
			this.info("Scanning "+filename);
			final SAMFileHeader header = IN.getFileHeader();
			List<SAMReadGroupRecord> rgs = header.getReadGroups();
			if(rgs==null || rgs.isEmpty())
				throw new IOException("No read groups in "+filename);
			Set<String> samples = new TreeSet<>();
			for(SAMReadGroupRecord rg:rgs)
				{
				String sample = rg.getSample();
				if(sample==null || sample.trim().isEmpty())
					{
					throw new IOException("Empty sample in "+rg);
					}
				samples.add(sample);
				}
			for(String sample : samples)
				{
				for(String gene: gene2interval.keySet())
					{
					int geneStart = Integer.MAX_VALUE;
					int geneEnd = 0;
					
					final List<Integer> counts = new ArrayList<>();
					
					for(Interval interval:gene2interval.get(gene))
						{
						geneStart = Math.min(geneStart, interval.getStart()-1);
						geneEnd = Math.max(geneEnd, interval.getEnd());
		
						/* picard javadoc:  - Sequence name - Start position (1-based) - End position (1-based, end inclusive)  */
						int interval_counts[]=new int[interval.getEnd()-interval.getStart()+1];
						if(interval_counts.length==0) continue;
						Arrays.fill(interval_counts, 0);
						
						
						/**
						 *     start - 1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
			    		*	   end - 1-based, inclusive end of interval of interest. Zero implies end of the reference sequence. 
						 */
						SAMRecordIterator r=IN.queryOverlapping(
								interval.getContig(),
								interval.getStart(),
								interval.getEnd()
								);
						while(r.hasNext())
							{
							SAMRecord rec=r.next();
							if(rec.getReadUnmappedFlag()) continue;
							if(rec.getReadFailsVendorQualityCheckFlag()) continue;
							if(rec.getDuplicateReadFlag() ) continue;
							if(rec.getReadPairedFlag())
								{
								if(PROPERLY_PAIRED && !rec.getProperPairFlag()) continue;
								}
							if(rec.isSecondaryOrSupplementary()) continue;
							if(!rec.getReferenceName().equals(interval.getContig())) continue;
							if(rec.getMappingQuality()==255 ||
								rec.getMappingQuality()==0 ||
								rec.getMappingQuality()< this.MMQ)
								{
								continue;
								}
							SAMReadGroupRecord rg = rec.getReadGroup();
							if(rg==null || !sample.equals(rg.getSample())) continue;
							Cigar cigar=rec.getCigar();
							if(cigar==null) continue;
				    		int refpos1=rec.getAlignmentStart();
				    		for(CigarElement ce:cigar.getCigarElements())
				    			{
				    			CigarOperator op=ce.getOperator();
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
						for(int d: interval_counts)
							{
							counts.add(d);
							}
						}/* end interval */
						
						IN.close();
						
						Collections.sort(counts);
						
						int count_no_coverage=0;
						double mean=0;
						for(int cov:counts)
							{
							if(cov<=MIN_COVERAGE) ++count_no_coverage;
							mean+=cov;
							}
						mean/=counts.size();
						
					System.out.println(
							gene2interval.get(gene).get(0).getContig()+"\t"+
							geneStart+"\t"+geneEnd+"\t"+gene+"\t"+sample+"\t"+
							counts.size()+"\t"+
							counts.get(0)+"\t"+
							counts.get(counts.size()-1)+"\t"+
							mean+"\t"+
							count_no_coverage+"\t"+
							(int)(((counts.size()-count_no_coverage)/(double)counts.size())*100.0)
							);
					}//end gene
				}//end sample
		}
	finally
		{
		CloserUtil.close(IN);
		}
	}
	
	@Override
	public String getProgramDescription() {
		return "Coverage statistics for a BED file, group by gene. It uses the Cigar string instead of the start/end to get the voverage";
		}
	
	@Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"BamStats05";
    }
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -m min coverage to say the position is not covered");
		out.println(" -B (file) bed file (chrom start end GENE)");
		out.println(" -p use orphan reads (not only properly paired)");
		out.println(" -q (int) min mapping quality:"+MMQ);
		super.printOptions(out);
		}
	@Override
	public int doWork(String[] args) {
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"m:B:pq:"))!=-1)
			{
			switch(c)
				{
				case 'm': this.MIN_COVERAGE=Integer.parseInt(opt.getOptArg());break;
				case 'B': this.BEDILE=new File(opt.getOptArg());break;
				case 'p': this.PROPERLY_PAIRED = false; break;
				case 'q': this.MMQ=Integer.parseInt(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(opt.getOptInd()==args.length)
			{
			error("Illegal number of args");
			return -1;
			}
		if(BEDILE==null)
			{
			error("missing bed file");
			return -1;
			}
		SamReader in=null;
		try
			{
			Map<String, List<Interval>> gene2interval = readBedFile(BEDILE);
			System.out.println("#chrom\tstart\tend\tgene\tsample\tlength\tmincov\tmaxcov\tmean\tnocoverage.bp\tpercentcovered");
			SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			for(int i=opt.getOptInd();i< args.length;++i)
				{
				in = srf.open(new File(args[i]));
				doWork(gene2interval,args[i],in);
				CloserUtil.close(in);
				in=null;
				}
			return 0;
			}
		catch (Exception e) {
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}

	
	public static void main(String[] args) throws IOException
		{
		new BamStats05().instanceMainWithExit(args);
		}

	}

