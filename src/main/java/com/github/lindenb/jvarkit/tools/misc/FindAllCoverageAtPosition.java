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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Map;
import java.util.TreeMap;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;

public class FindAllCoverageAtPosition extends AbstractCommandLineProgram
	{
	private static class Mutation
		{
		String chrom;
		int pos;
		Mutation(String chrom,int pos)
			{
			this.chrom=chrom;
			this.pos=pos;
			}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + chrom.hashCode();
			result = prime * result + pos;
			return result;
			}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)return true;
			Mutation other = (Mutation) obj;
			if (pos != other.pos) return false;
			 if (!chrom.equals(other.chrom))
				return false;
			
			return true;
		}
		
		@Override
		public String toString() {
			return chrom+":"+pos;
			}
		
		}
	private Mutation mutation=null;
	private PrintWriter out=null;
	private SamReaderFactory samReaderFactory;
    private FindAllCoverageAtPosition()
    	{
    	samReaderFactory=  SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
    	}		
    @Override
    protected String getOnlineDocUrl() {
    	return "https://github.com/lindenb/jvarkit/wiki/FindAllCoverageAtPosition";
    	}
    
    @Override
    public String getProgramDescription() {
    	return "Find depth at specific position in a list of BAM files.";
    	}
    
  
    
    private Mutation convertFromSamHeader(File f,SAMFileHeader h)
    	{
    	SAMSequenceDictionary dict=h.getSequenceDictionary();
    	if(dict==null)
    		{
    		warning("No dictionary in "+h);
    		return null;
    		}
    	SAMSequenceRecord rec=dict.getSequence(this.mutation.chrom);
    	if(rec!=null) return this.mutation;
    	String chromName=this.mutation.chrom;
		if(chromName.startsWith("chr"))
			{
			rec=dict.getSequence(chromName.substring(3));
			if(rec!=null) return new Mutation(rec.getSequenceName(),this.mutation.pos);
			}
		else
			{
			rec=dict.getSequence("chr"+chromName);
			if(rec!=null) return new Mutation(rec.getSequenceName(),this.mutation.pos);
			}
	
		if(chromName.equals("MT") &&  dict.getSequence("chrM")!=null)
			{
			return new Mutation("chrM",this.mutation.pos);
			}
		if(chromName.equals("chrM") &&  dict.getSequence("MT")!=null)
			{
			return new Mutation("MT",this.mutation.pos);
			}
    	return null;
    	}

    private void scan(BufferedReader in) throws IOException
    	{
    	final String DEFAULT_SAMPLE_NAME="(undefined)";
    	String line;
    	while((line=in.readLine())!=null)
			{
    		if(out.checkError()) break;
			if(line.isEmpty() || line.startsWith("#")) continue;
			File f=new File(line);
			if(!f.isFile()) continue;
			if(!f.canRead()) continue;
			String filename=f.getName();
			if(!filename.endsWith(".bam")) continue;
			
    			
			SamReader samReader=null;
			SAMRecordIterator iter=null;
			try
				{
				samReader = this.samReaderFactory.open(f);
				if(!samReader.hasIndex())
					{
					warning("no index for "+f);
					continue;
					}
				SAMFileHeader header=samReader.getFileHeader();
				Map<String, Counter<CigarOperator>> sample2count=new TreeMap<String, Counter<CigarOperator>>();
				for(SAMReadGroupRecord rg:header.getReadGroups())
					{
					if(rg!=null)
						{
						String sn= rg.getSample();
						if(sn!=null && !sn.trim().isEmpty())
							{
							sample2count.put(sn, new Counter<CigarOperator>());
							}
						}
					}
				
				if(sample2count.isEmpty())
					{
					sample2count.put(DEFAULT_SAMPLE_NAME, new Counter<CigarOperator>());
					}
				
				Mutation m = convertFromSamHeader(f,header);
				if(m==null) continue;
				iter=samReader.query(m.chrom, m.pos-1, m.pos+1,	false);
				while(iter.hasNext())
					{
					SAMRecord rec=iter.next();
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getReadFailsVendorQualityCheckFlag()) continue;
					if(rec.getNotPrimaryAlignmentFlag()) continue;
					if(rec.isSecondaryOrSupplementary()) continue;
					if(rec.getDuplicateReadFlag()) continue;
					Cigar cigar=rec.getCigar();
					if(cigar==null) continue;
					String sampleName=DEFAULT_SAMPLE_NAME;
					SAMReadGroupRecord rg=rec.getReadGroup();
					if(rg!=null)
						{
						String sn= rg.getSample();
						if(sn!=null && !sn.trim().isEmpty())
							{
							sampleName=sn;
							}
						}
					Counter<CigarOperator> counter= sample2count.get(sampleName);
					if(counter==null)
						{
						counter=new Counter<CigarOperator>();
						sample2count.put(sampleName, counter);
						}	
					int ref= rec.getUnclippedStart();
					for(int k=0;k<cigar.numCigarElements() && ref< m.pos+1;++k)
						{
						CigarElement ce=cigar.getCigarElement(k);
						CigarOperator op=ce.getOperator();
						switch(op)
							{
							case P: break;
							case I: 
								{
								if(ref==m.pos) counter.incr(op);
								break;
								}
							case D:case N:
							case M: case X: case EQ: 
							case H:
							case S:
								{
								for(int i=0;i< ce.getLength();++i )
									{
									if(ref==m.pos)
										{
										counter.incr(op);
										break;
										}	
									ref++;
									}
								break;
								}
							default: throw new RuntimeException("unknown operator:"+op);
							}
						}
					}
				
				
				
				for(String sample:sample2count.keySet())
					{
					Counter<CigarOperator> counter= sample2count.get(sample);
					
					out.print(f);
					out.print('\t');
					out.print(m.chrom);
					out.print('\t');
					out.print(m.pos);
					out.print('\t');
					out.print(sample);
					out.print('\t');
					out.print(
							counter.count(CigarOperator.M)+
							counter.count(CigarOperator.EQ)+
							counter.count(CigarOperator.X)
							);
					for(CigarOperator op:CigarOperator.values())
						{
						out.print('\t');
						out.print(counter.count(op));
						}
					out.println();
					}
				}
			catch(Exception err)
				{
				error(err);
				}
			finally
				{
				CloserUtil.close(iter);
				CloserUtil.close(samReader);
				}    				
			}		
	    
    	}
    
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -p chrom:pos . Add this chrom/position. Required");
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"p:"))!=-1)
			{
			switch(c)
				{
				case 'p':
					{
					String s=opt.getOptArg();
					int colon=s.indexOf(':');
					if(colon==-1 || colon+1==s.length())
						{
						error("Bad chrom:pos "+s);
						return -1;
						}
					
					String chrom=s.substring(0,colon).trim();
					if(chrom.isEmpty())
						{
						error("Bad chrom:pos "+s);
						return -1;
						}
					Mutation m=new Mutation(chrom, Integer.parseInt(s.substring(colon+1)));
					info("adding "+m);
					this.mutation = m;
					break;
					}
				default:
					{
					switch(super.handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS:return 0;
						case OK:break;
						}
					}
				}
			}
		if(this.mutation==null)
			{
			error("position not defined");
			return -1;
			}
		try
			{
			this.out=new PrintWriter(System.out);
			
			out.print("#File");
			out.print('\t');
			out.print("CHROM");
			out.print('\t');
			out.print("POS");
			out.print('\t');
			out.print("SAMPLE");
			out.print('\t');
			out.print("DEPTH");
			for(CigarOperator op:CigarOperator.values())
				{
				out.print('\t');
				out.print(op.name());
				}
			out.println();

			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				scan(new BufferedReader(new InputStreamReader(System.in)));
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("Reading from "+filename);
					BufferedReader r=IOUtils.openURIForBufferedReading(filename);
					scan(r);
					r.close();
					}
				}
			this.out.flush();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{

			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FindAllCoverageAtPosition().instanceMainWithExit(args);

	}

}
