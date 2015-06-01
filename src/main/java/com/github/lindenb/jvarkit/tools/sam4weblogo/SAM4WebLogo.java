package com.github.lindenb.jvarkit.tools.sam4weblogo;

import java.io.PrintStream;
import java.io.PrintWriter;


import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SAM4WebLogo extends AbstractCommandLineProgram
	{
    private SAM4WebLogo()
    	{
    	
    	}
    @Override
    protected String getOnlineDocUrl()
    	{
    	return "https://github.com/lindenb/jvarkit/wiki/SAM4WebLogo";
    	}
    @Override
    public String getProgramDescription()
    	{
    	return "Sequence logo for different alleles or generated from SAM/BAM http://www.biostars.org/p/73021";
    	}
    
    
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -R Region to observe: chrom:start-end . REQUIRED.");
		out.println(" -c use clipped bases.");
		super.printOptions(out);
		}
    @Override
    public int doWork(String[] args)
    	{
    	boolean useClip=false;
    	Interval interval=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt getopt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=getopt.getopt(args,getGetOptDefault()+ "r:c"))!=-1)
			{
			switch(c)
				{
				case 'c': useClip=true;break;
				case 'r':
					{
					interval=parseInterval(getopt.getOptArg());
					if(interval==null)
						{
						error("Bad interval "+getopt.getOptArg());
						return -1;
						}
					break;
					}
				default: 
					{
					switch(handleOtherOptions(c, getopt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}

		if(interval==null)
			{
			error("Undefined interval.");
			return -1;
			}
		
		PrintWriter out=new PrintWriter(System.out);
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		boolean warningInsertion=false;
		try {
			if(getopt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				samReader=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(getopt.getOptInd()+1==args.length)
				{
				samReader=SamFileReaderFactory.mewInstance().open(args[getopt.getOptInd()]);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
		
			if(samReader.hasIndex())
					{
					iter=samReader.queryOverlapping(
							interval.getContig(),
							interval.getStart(),
							interval.getEnd()
							);
					}
			else
					{
					iter=samReader.iterator();
					}
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samReader.getFileHeader().getSequenceDictionary());
	       while(iter.hasNext())
                {
                SAMRecord rec=iter.next();
                progress.watch(rec);
                if(rec.getReadUnmappedFlag()) continue;
                if(!rec.getReferenceName().equals(interval.getContig())) continue;
                if(rec.getAlignmentEnd() < interval.getStart() ) continue;
                if(rec.getAlignmentStart() > interval.getEnd() ) continue;
                Cigar cigar=rec.getCigar();
                if(cigar==null) continue;
                byte bases[]=rec.getReadBases();
                
                StringBuilder seq=new StringBuilder(interval.length());
                int readPos=0;
                int refPos=rec.getUnclippedStart();
                for(int i=0;i< cigar.numCigarElements();++i)
                	{
                	CigarElement ce=cigar.getCigarElement(i);
            		CigarOperator op=ce.getOperator();
                	switch(op)
                		{
                		case P:break;
                		case I:
                			{
                			warningInsertion=true;
                			break;
                			}
                		case D: case N:
                			{
		        			for(int j=0;j< ce.getLength() && refPos <=interval.getEnd()  ;++j)
		        				{
		        				if(refPos>= interval.getStart())
		        					{
		        					seq.append('-');
		        					}
		        				refPos++;
		        				}
		        			break;
                			}
                		case H:
		        			{
		        			for(int j=0;j< ce.getLength() && refPos <=interval.getEnd()  ;++j)
		        				{
		        				if(refPos>= interval.getStart() && useClip)
		        					{
		        					seq.append('N');
		        					}
		        				refPos++;
		        				}
		        			break;
		        			}
                		case S:
                			{
                			for(int j=0;j< ce.getLength() && refPos <=interval.getEnd()  ;++j)
                				{
                				if(refPos>= interval.getStart() && useClip)
                					{
                					seq.append((char)bases[readPos]);
                					}
                				readPos++;
                				refPos++;
                				}
                			break;
                			}
                		case M:case X: case EQ:
                			{
                			for(int j=0;j< ce.getLength() && refPos <=interval.getEnd() ;++j)
                				{
                				if(refPos>= interval.getStart())
                					{
                					seq.append((char)bases[readPos]);
                					}
                				readPos++;
                				refPos++;
                				}
                			break;
                			}
                		default:throw new IllegalStateException("Not handled. op:"+op);
                		}	
                	}
                if(seq.length()==0) continue;
                for(int i= interval.getStart();
                	i< (useClip?rec.getUnclippedStart():rec.getAlignmentStart());
                	++i)
                	{
                	seq.insert(0, '-');
                	}
                while(seq.length()< interval.length())
	            	{
	            	seq.append('-');
	            	}
            	out.print(">"+rec.getReadName());
            	if(rec.getReadPairedFlag())
	            	{
	            	if(rec.getFirstOfPairFlag()) out.print("/1");
	            	if(rec.getSecondOfPairFlag()) out.print("/2");
	            	}
            	out.println();
            	out.println(seq);
                }
	       progress.finish();
	        if(warningInsertion)
	        	{
	        	warning("Some reads contained insertions.");
	        	}
	        
			} 
		catch (Exception e) {
			e.printStackTrace();
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			out.flush();
			}
		return 0;
		}

private static Interval parseInterval(String reg)
	{
	try
			{	 
				int colon = reg.indexOf(':');
				if(colon<1) throw new IllegalArgumentException("bad region "+reg);
				int hyphen = reg.indexOf('-');
		
				String s=reg.substring(0,colon);
				int start= Integer.parseInt(reg.substring(colon+1,hyphen));
				int end=Integer.parseInt(reg.substring(hyphen+1));
				return new Interval(s, start, end);
			}
			catch(Exception err)
			{
				System.err.println("bad interval "+reg);
				return null;
			}
		}
public static void main(String[] args)
	{
	new SAM4WebLogo().instanceMainWithExit(args);
	}
}
