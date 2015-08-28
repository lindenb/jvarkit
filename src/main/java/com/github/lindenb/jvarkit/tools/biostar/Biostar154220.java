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
package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.AbstractBamWriterProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class Biostar154220 extends AbstractBamWriterProgram
	{
	
	private int capDepth=20;
	public void setCapDepth(int capDepth) {
		this.capDepth = capDepth;
	}
	
	public int getCapDepth() {
		return capDepth;
	}
	
	
	@Override
	public String getProgramDescription()
		{
		return "Cap BAM to a given coverage. see https://www.biostars.org/p/154220";
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return DEFAULT_WIKI_PREFIX+"Biostar154220";
		}

	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -n (int) number of reads. default:"+getCapDepth()); 
		super.printOptions(out);
		}

	@SuppressWarnings("resource")
	private int doWork( SamReader in) throws IOException
		{
		SAMFileHeader header= in.getFileHeader();
		if(header.getSortOrder()!=SAMFileHeader.SortOrder.unsorted)
			{
			error("input should be unsorted, reads sorted on REF/query-name e.g: see https://github.com/lindenb/jvarkit/wiki/SortSamRefName");
			return -1;
			}
		SAMSequenceDictionary dict=header.getSequenceDictionary();
		if(dict==null)
			{
			error("no dict !");
			return -1;
			}
		SAMFileWriter out=null;
		SAMRecordIterator iter=null;
		int prev_tid=-1;
		int depth_array[]=null;
		try
			{
			SAMFileHeader header2=header.clone();
			header2.addComment(getProgramName()+" "+getVersion()+" "+getProgramCommandLine());
			out = openSAMFileWriter(header2, true);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			iter = in.iterator();
			List<SAMRecord> buffer=new ArrayList<>();
			for(;;)
				{
				SAMRecord rec =null;
				
				if(iter.hasNext())
					{
					rec = progress.watch(iter.next());
					}
				
				if(rec!=null && rec.getReadUnmappedFlag())
					{
					out.addAlignment(rec);
					continue;
					}
				//no more record or 
				if(!buffer.isEmpty() &&
					rec!=null &&
					buffer.get(0).getReadName().equals(rec.getReadName()) &&
					buffer.get(0).getReferenceIndex().equals(rec.getReferenceIndex())
					)
					{
					buffer.add(progress.watch(rec));
					}
				else if(buffer.isEmpty() && rec!=null)
					{
					buffer.add(progress.watch(rec));
					}
				else //dump buffer
					{
					if(!buffer.isEmpty())
						{
						final int tid = buffer.get(0).getReferenceIndex();
						if(prev_tid==-1 || prev_tid!=tid)
							{
							SAMSequenceRecord ssr=dict.getSequence(tid);
							prev_tid=tid;
							depth_array=null;
							System.gc();
							info("Alloc memory for contig "+ssr.getSequenceName()+" N="+ssr.getSequenceLength()+"*sizeof(int)");
							depth_array=new int[ssr.getSequenceLength()+1];//use a +1 pos
							Arrays.fill(depth_array, 0);
							}
						//position->coverage for this set of reads
						Counter<Integer> readposition2coverage=new Counter<Integer>();
						
						boolean dump_this_buffer=true;
						for(SAMRecord sr:buffer)
							{
							if(!dump_this_buffer) break;
							if(sr.isSecondaryOrSupplementary()) continue;
							if(sr.getDuplicateReadFlag()) continue;
							if(sr.getMappingQuality()==0) continue;
							
							
							Cigar cigar=sr.getCigar();
							if(cigar==null)
								{
								throw new IOException("Cigar missing in "+rec.getSAMString());
								}
							int refPos1=sr.getAlignmentStart();
							for(CigarElement ce:cigar.getCigarElements())
								{
								final CigarOperator op =ce.getOperator();
								if(!op.consumesReferenceBases()) continue;
								if(op.consumesReadBases())
									{
									for(int x=0;x<ce.getLength() && refPos1+x< depth_array.length;++x)
										{
										int cov = (int)readposition2coverage.incr(refPos1+x);
										if( depth_array[refPos1+x]+cov > this.getCapDepth())
											{
											dump_this_buffer=false;
											break;
											}
										}
									}
								if(!dump_this_buffer) break;
								refPos1+=ce.getLength();
								}
							}
						if(dump_this_buffer)
							{
							//consumme this coverage
							for(Integer pos:readposition2coverage.keySet())
								{
								depth_array[pos]+= (int)readposition2coverage.count(pos);
								}
							for(SAMRecord sr:buffer)
								{
								out.addAlignment(sr);
								}
							}
						
						buffer.clear();
						}
					if(rec==null) break;
					buffer.add(rec);
					}
				}
			depth_array=null;
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
			CloserUtil.close(out);
			}
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "n:"))!=-1)
			{
			switch(c)
				{
				case 'n': this.setCapDepth(Integer.parseInt(opt.getOptArg()));break;				
				case 'o': this.setOutputFile(new File(opt.getOptArg()));break;				
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		if(getCapDepth()<0) // -1 == infinite
			{
			error("Bad depth:"+getCapDepth());
			return -1;
			}
		SamReader in=null;
		try
			{
			SamReaderFactory srf= SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			if(opt.getOptInd()==args.length)
				{
				in = srf.open(SamInputResource.of(System.in));
				}
			else if(opt.getOptInd()+1==args.length)
				{
				in = srf.open(SamInputResource.of(args[opt.getOptInd()]));
				}
			else
				{
				error("Illegal number of arguments");
				return -1;
				}
			return doWork(in); 
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}

	
	public static void main(String[] args) throws IOException
		{
		new Biostar154220().instanceMainWithExit(args);
		}
		

	}
