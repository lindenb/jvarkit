/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.pcr;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class PcrSliceReads extends AbstractCommandLineProgram
	{
	private boolean clip_reads=false;
	private IntervalTreeMap<Interval> bedIntervals=new IntervalTreeMap<Interval>();
	private File fileout = null;
	private boolean binary=false;
	private enum AmbiguityStrategy {zero,random,closer};
	private AmbiguityStrategy ambiguityStrategy= AmbiguityStrategy.closer; 
	private Random random=new Random(0L);//0L for reproductive calculations
	@Override
	public String getProgramDescription() {
		return "Mark PCR reads to their PCR amplicon https://www.biostars.org/p/149687/";
		}
	
	@Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"PcrSliceReads";
    	}
	
	private List<Interval> findIntervals(String chrom,int start,int end)
		{
		Interval i= new Interval(chrom,start,end);
		return  new ArrayList<Interval>(this.bedIntervals.getOverlapping(i));
		}
	
	
	@SuppressWarnings("resource")
	private int run(SamReader reader)
		{
		ReadClipper readClipper = new ReadClipper();
		SAMFileHeader header1= reader.getFileHeader();
		SAMFileHeader header2 = header1.clone();
		header2.addComment(getProgramName()+" "+getVersion()+": Processed with "+getProgramCommandLine());
		header2.setSortOrder(SortOrder.unsorted);
		
		for(SAMReadGroupRecord srg:header1.getReadGroups())
			{
			
			for(Interval i:this.bedIntervals.keySet())
				{
				//create new read group for this interval
				SAMReadGroupRecord rg=new SAMReadGroupRecord(srg.getId()+"_"+this.bedIntervals.get(i).getName(), srg);
				header2.addReadGroup(rg);
				}
			}
		
		SAMFileWriter sw=null;
		SAMRecordIterator iter = null;
		try
			{
			SAMFileWriterFactory sfw=new SAMFileWriterFactory();
			
			if( this.fileout == null )
				{
				if( this.binary)
					{
					sw = sfw.makeBAMWriter(header2, false, System.out);
					}
				else
					{
					sw = sfw.makeSAMWriter(header2, false, System.out);
					}
				}
			else
				{
				sw = sfw.makeSAMOrBAMWriter(header2, false, this.fileout);
				}
			SAMSequenceDictionaryProgress progress =new SAMSequenceDictionaryProgress(header1);
			iter =  reader.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec= progress.watch(iter.next());
				if(rec.getReadUnmappedFlag())
					{
					sw.addAlignment(rec);
					continue;
					}
				
				if(!rec.getReadPairedFlag())
					{
					//@doc if the read is not a paired-end read ,  then the quality of the read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				
				if(rec.getMateUnmappedFlag())
					{
					//@doc if the MATE is not mapped ,  then the quality of the current read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				if(!rec.getProperPairFlag())
					{
					//@doc if the properly-paired flag is not set,  then the quality of the current read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				
				if(rec.getMateReferenceIndex()!=rec.getReferenceIndex())
					{
					//@doc if the read and the mate are not mapped on the same chromosome,  then the quality of the current read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				
				if(rec.getReadNegativeStrandFlag()==rec.getMateNegativeStrandFlag())
					{
					//@doc if the read and the mate are mapped on the same strand,  then the quality of the current read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				int chromStart;
				int chromEnd;
				if(rec.getAlignmentStart() < rec.getMateAlignmentStart())
					{
					if(rec.getReadNegativeStrandFlag())
						{
						//@doc if the read(POS) < mate(POS) and read is mapped on the negative strand, then the quality of the current read is set to zero
						rec.setMappingQuality(0);
						sw.addAlignment(rec);
						continue;
						}
					else
						{
						chromStart = rec.getAlignmentStart();
						chromEnd = rec.getMateAlignmentStart();
						}
					}
				else 
					{
					if(!rec.getReadNegativeStrandFlag())
						{
						//@doc if the read(POS) > mate(POS) and read is mapped on the positive strand, then the quality of the current read is set to zero
						rec.setMappingQuality(0);
						sw.addAlignment(rec);
						continue;
						}
					else
						{
						chromStart = rec.getMateAlignmentStart();
						chromEnd = rec.getAlignmentStart();//don't use getAlignmentEnd, to be consistent with mate data
						}
					}
				
				
				
				
				List<Interval> fragments = findIntervals(rec.getContig(),chromStart,chromEnd);
				if(fragments.isEmpty())
					{
					//@doc if no BED fragment is found overlapping the read, then the quality of the read is set to zero
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				Interval best=null;
				if(fragments.size()==1)
					{
					best = fragments.get(0);
					}
				else switch(this.ambiguityStrategy)
					{
					case random:
						{
						best = fragments.get(this.random.nextInt(fragments.size()));
						break;
						}
					case zero:
						{
						rec.setMappingQuality(0);
						sw.addAlignment(rec);
						continue;
						}
					case closer:
						{
						int best_distance=Integer.MAX_VALUE;
						for(Interval frag : fragments)
							{
							int distance= Math.abs(chromStart-frag.getStart()) + Math.abs(frag.getEnd() - chromEnd);
							if(best==null || distance < best_distance)
								{
								best = frag;
								best_distance = distance;
								}
							}
						break;
						}
					default: throw new IllegalStateException();
					}
				
				
				
				if(clip_reads)
					{
					rec = readClipper.clip(rec, best);
					if(rec.getMappingQuality()==0)
						{
						sw.addAlignment(rec);
						continue;
						}
					}
				SAMReadGroupRecord rg = rec.getReadGroup();
				if(rg == null )
					{
					throw new IOException("Read "+rec+" is missing a Read-Group ID . See http://broadinstitute.github.io/picard/ http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups");
					}
				rec.setAttribute("RG",rg.getId()+"_"+best.getName());
				sw.addAlignment(rec);
				}
			progress.finish();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sw);
			}
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -a (strategy) if a read is mapped on multiple PCR fragments, how to resolve ambiguity ? default:"+this.ambiguityStrategy.name()+" . Where strategy is"); 
		out.println("     "+AmbiguityStrategy.random.name() + " : choose a random fragment"); 
		out.println("     "+AmbiguityStrategy.zero.name() + " : set MAPQ to zero and ignore."); 
		out.println("     "+AmbiguityStrategy.closer.name() + " : choose PCR fragment closest to NGS-fragment boundaries."); 
		out.println(" -o (file) output file (default stdout)"); 
		out.println(" -b force binary for stdout (optional)"); 
		out.println(" -B (file) bed file containing non-overlapping PCR fragments. Column name is required."); 
		out.println(" -c clip read to their PCR fragments. see https://github.com/lindenb/jvarkit/wiki/PcrClipReads"); 
		super.printOptions(out);
		}

	
	@SuppressWarnings("resource")
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		File bedFile=null;
		while((c=opt.getopt(args,getGetOptDefault()+"o:bB:ca:"))!=-1)
			{
			switch(c)
				{
				case 'a': this.ambiguityStrategy = AmbiguityStrategy.valueOf(opt.getOptArg());break;
				case 'c': this.clip_reads=true; break;
				case 'B': bedFile =new File(opt.getOptArg());break;
				case 'b': binary=true;break;
				case 'o': fileout = new File(opt.getOptArg());break;				
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
		if(bedFile==null)
			{
			error("undefined bed file");
			return -1;
			}
		BufferedReader r=null;
		SamReader samReader=null;
		try {
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			if(opt.getOptInd()==args.length)
				{
				samReader = srf.open(SamInputResource.of(System.in));
				}
			else if(opt.getOptInd()+1==args.length)
				{
				samReader = srf.open(SamInputResource.of(args[opt.getOptInd()]));
				}
			else
				{
				error("illegal number of args");
				return -1;
				}
			
			Pattern tab= Pattern.compile("[\t]");
			 r= IOUtils.openFileForBufferedReading(bedFile);
			String line;
			while((line=r.readLine())!=null)
				{
				String tokens[]=tab.split(line);
				if(tokens.length<4)
					{
					error("Bad bed line "+line);
					return -1;
					}
				String chrom = tokens[0];
				int chromStart1 = Integer.parseInt(tokens[1])+1;
				int chromEnd1 = Integer.parseInt(tokens[2])+0;
				if(chromStart1<1 || chromStart1>chromEnd1)
					{
					error("Bad bed line "+line);
					return -1;
					}
				String name = tokens[3].trim();
				if(name.isEmpty())
					{
					error("Bad bed line (name missing) in  "+line);
					return -1;
					}
				Interval i =new Interval(chrom, chromStart1, chromEnd1,false,name);
				this.bedIntervals.put(i, i);
				}
			return run(samReader);
			}
		catch (Exception e) {
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(samReader);
			}
		}

	
	public static void main(String[] args) {
		new PcrSliceReads().instanceMain(args);
		}

}
