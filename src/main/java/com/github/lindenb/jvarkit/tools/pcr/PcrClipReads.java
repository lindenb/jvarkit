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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
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

public class PcrClipReads extends AbstractCommandLineProgram
	{
	private IntervalTreeMap<Interval> bedIntervals=new IntervalTreeMap<Interval>();
	private File fileout = null;
	private boolean binary=false;
	@Override
	public String getProgramDescription() {
		return "Soft clip bam files based on PCR target regions https://www.biostars.org/p/147136/";
		}
	
	@Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"PcrClipReads";
    	}
	
	private Interval findInterval(final SAMRecord rec)
		{
		if(rec.getReadUnmappedFlag()) return null;
		return findInterval(rec.getContig(), rec.getAlignmentStart(), rec.getAlignmentEnd());
		}
	private Interval findInterval(String chrom,int start,int end)
		{
		Interval i= new Interval(chrom,start,end);
		Collection<Interval> L=this.bedIntervals.getOverlapping(i);
		Iterator<Interval> iter = L.iterator();
		if(iter.hasNext())
			{
			Interval j = iter.next();
			if(iter.hasNext()  ) throw new IllegalStateException("Overlapping PCR intervals : "+j+" "+iter.next());
			return j;
			}
		return null;
		}
	
	
	private int run(SamReader reader)
		{
		SAMFileHeader header1= reader.getFileHeader();
		SAMFileHeader header2 = header1.clone();
		header2.addComment(getProgramName()+" "+getVersion()+": Processed with "+getProgramCommandLine());
		header2.setSortOrder(SortOrder.unsorted);
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
				Interval fragment = findInterval(rec);
				if(fragment==null)
					{
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				// strand is '-' and overap in 5' of PCR fragment
				if( rec.getReadNegativeStrandFlag() &&
					fragment.getStart()< rec.getAlignmentStart() &&
					rec.getAlignmentStart()< fragment.getEnd())
					{
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					continue;
					}
				// strand is '+' and overap in 3' of PCR fragment
				if( !rec.getReadNegativeStrandFlag() &&
					fragment.getStart()< rec.getAlignmentEnd() &&
					rec.getAlignmentEnd()< fragment.getEnd())
					{
					rec.setMappingQuality(0);
					sw.addAlignment(rec);
					
					continue;
					}
				
				// contained int PCR fragment
				if(rec.getAlignmentStart()>= fragment.getStart() && rec.getAlignmentEnd()<=fragment.getEnd())
					{
					sw.addAlignment(rec);
					
					continue;
					}
				
				int newStart = rec.getAlignmentStart();
				
				Cigar cigar = rec.getCigar();
				if(cigar==null)
					{
					warning("cigar missing in "+rec);
					sw.addAlignment(rec);
					continue;
					}
				
				List<CigarOperator> operators = new ArrayList<>();
				//expand cigar	
				for(CigarElement ce:cigar.getCigarElements())
					{
					CigarOperator op = ce.getOperator();
					for(int x=0;x < ce.getLength();++x)
						{
						operators.add(op);
						}
					}
				
				/* 5' side */				
				if(rec.getAlignmentStart() < fragment.getStart())
					{
					int operator_index =0;
					newStart = fragment.getStart();
					int refPos1 = rec.getUnclippedStart();
					
					while(	operator_index < operators.size() &&
							refPos1< fragment.getStart())
						{
						CigarOperator op = operators.get(operator_index);
						if(op.consumesReferenceBases() )
							{
							refPos1++;
							}
						if(op.consumesReadBases())
							{
							operators.set(operator_index, CigarOperator.S);
							}
						operator_index++;
						}
					//shouln't start with a problem
					while(operator_index < operators.size())
						{
						CigarOperator op = operators.get(operator_index);
						if(op.equals(CigarOperator.M) || op.equals(CigarOperator.EQ))
							{
							break;
							}
						if(op.consumesReadBases())//insertion...
							{
							operators.set(operator_index, CigarOperator.S);
							newStart++;
							}
						operator_index++;
						}
					int x=0;
					int y=0;
					for(x=0;x< operator_index;++x)
						{
						CigarOperator op = operators.get(y);
						if(!(op.equals(CigarOperator.S) || op.equals(CigarOperator.H)))
							{
							operators.remove(y);
							}
						else
							{
							++y;
							}
						}
					
					}

				
				/* 3' side */				
				if(rec.getAlignmentEnd() > fragment.getEnd())
					{
					int operator_index = operators.size()-1;
					int refPos1 = rec.getUnclippedEnd();

					while(operator_index >=0 &&
						refPos1 > fragment.getEnd())
						{
						CigarOperator op = operators.get(operator_index);
						if(op.consumesReferenceBases() )
							{
							refPos1--;
							}
						if(op.consumesReadBases())
							{
							operators.set(operator_index, CigarOperator.S);
							}
						operator_index--;
						}
					
					//shouln't end with a problem
					while(operator_index >=0 )
						{
						CigarOperator op = operators.get(operator_index);
						if(op.equals(CigarOperator.M) || op.equals(CigarOperator.EQ))
							{
							break;
							}
						if(op.consumesReadBases())//insertion...
							{
							operators.set(operator_index, CigarOperator.S);
							}
						operator_index--;
						}
					int y=operators.size()-1;
					int len = (operators.size()-1)-operator_index;
					for(int x=0;x< len;++x)
						{
						CigarOperator op = operators.get(y);
						if(!(op.equals(CigarOperator.S) || op.equals(CigarOperator.H)))
							{
							operators.remove(y);
							}
						else
							{
							--y;
							}
						}
					
					}
				
				
				
				//build new cigar
				boolean found_M=false;
				List<CigarElement> newcigarlist=new ArrayList<>();
				int c1 = 0;
				while(c1 < operators.size())
					{
					CigarOperator op =  operators.get(c1);
					if(op.equals(CigarOperator.M) || op.equals(CigarOperator.EQ))
						{
						found_M=true;
						}
					int c2=c1;
					while(c2< operators.size() && op.equals(operators.get(c2)))
						{
						++c2;
						}
					newcigarlist.add(new CigarElement(c2-c1,op));
					c1=c2;
					}
				
				if(!found_M)
					{
					
					rec.setReadUnmappedFlag(true);
					rec.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
					rec.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
					rec.setMappingQuality(0);
					//rec.setCigar(null);
					//rec.setCigarString(null);
					sw.addAlignment(rec);
					continue;
					}
				cigar = new Cigar(newcigarlist);				
				rec.setCigar(cigar);
				rec.setAlignmentStart(newStart);
			
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
		out.println(" -o (file) output file (default stdout)"); 
		out.println(" -b force binary for stdout (optional)"); 
		out.println(" -B (file) bed file containing non-overlapping PCR fragments."); 
		super.printOptions(out);
		}

	
	@SuppressWarnings("resource")
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		File bedFile=null;
		while((c=opt.getopt(args,getGetOptDefault()+"o:bB:"))!=-1)
			{
			switch(c)
				{
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
				if(tokens.length<3)
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
				Interval i =new Interval(chrom, chromStart1, chromEnd1);
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
		new PcrClipReads().instanceMain(args);
		}

}
