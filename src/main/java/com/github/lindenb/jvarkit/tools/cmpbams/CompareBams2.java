package com.github.lindenb.jvarkit.tools.cmpbams;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;


public class CompareBams2  extends AbstractCommandLineProgram
	{
	private static final Logger LOG=Logger.getLogger(CompareBams2.class.getName());
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Compare two or more BAM files.";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM files to process.",
    		minElements=2,optional=false)
	public List<File> IN=new ArrayList<File>();
    @Option(shortName= "L", doc="restrict to that region (chr:start-end)",optional=true)
	public String REGION=null;
    
    @Option(shortName= "FLG", doc="use SAM Flag when comparing.",optional=true)
    public boolean useSamFlag=false;

    @Option(shortName= "CIGAR", doc="use CIGAR when comparing.",optional=true)
    public boolean useCigar=false;

	
	private class MatchComparator
		implements Comparator<Match>
		{
		@Override
		public int compare(Match m0, Match m1)
			{
			int i=m0.readName.compareTo(m1.readName);
			if(i!=0) return i;
			i=m0.num_in_pair-m1.num_in_pair;
			if(i!=0) return i;
			//i= m0.bamIndex - m1.bamIndex;//NO ! (when comparing two Set<Match>)
			//if(i!=0) return i;
			i= m0.tid - m1.tid;
			if(i!=0) return i;
			i= m0.pos - m1.pos;
			if(i!=0) return i;
			i= m0.flag - m1.flag;
			if(i!=0) return i;
			i= m0.cigar.compareTo(m1.cigar);
			return 0;
			}
		}
	
	private class MatchCodec
		extends AbstractDataCodec<Match>
		{
		@Override
		public MatchCodec clone()
			{
			return new MatchCodec();
			}
		@Override
		public Match decode(DataInputStream dis) throws IOException
			{
			Match m=new Match();
			try {
				m.readName=dis.readUTF();
				}
			catch(IOException err)
				{
				return null;
				}
			m.bamIndex=dis.readInt();
			m.tid=dis.readInt();
			m.pos=dis.readInt();
			m.num_in_pair=dis.readInt();
			if(useSamFlag) m.flag=dis.readInt();
			if(useCigar) m.cigar=dis.readUTF();
			return m;
			}
		@Override
		public void encode(DataOutputStream dos, Match match)
				throws IOException
			{
			dos.writeUTF(match.readName);
			dos.writeInt(match.bamIndex);
			dos.writeInt(match.tid);
			dos.writeInt(match.pos);
			dos.writeInt(match.num_in_pair);
			if(useSamFlag) dos.writeInt(match.flag);
			if(useCigar) dos.writeUTF(match.cigar);
			}
		
		}
	
	private static class Match
		{
		String readName;
		int num_in_pair=0;
		int tid=-1;
		int bamIndex=-1;
		int pos=-1;
		int flag=0;
		String cigar="";

		@Override
		public int hashCode()
			{
			int result = 1;
			result = 31 * result + num_in_pair;
			result = 31 * result + pos;
			result = 31 * result + tid;
			result = 31 * result + readName.hashCode();
			result = 31 * result + bamIndex;
			result = 31 * result + cigar.hashCode();
			return result;
			}
		@Override
		public boolean equals(Object obj)
			{
			if (this == obj) { return true; }
			if (obj == null) { return false; }
			Match other = (Match) obj;
			if (num_in_pair != other.num_in_pair) { return false; }
			if (tid != other.tid) { return false; }
			if(tid==-1) return true;
			if (pos != other.pos) { return false; }
			if (bamIndex != other.bamIndex) { return false; }
			if(!readName.equals(other.readName)) return false;
			return true;
			}
		}
	
	
	
	private void print(final Set<Match> set,final SAMSequenceDictionary dict)
		{
		boolean first=true;
		for(Match m:set)
			{
			if(!first)System.out.print(',');
			first=false;
			if(m.tid<0){ System.out.print("unmapped"); continue;}
			SAMSequenceRecord ssr=(dict==null?null:dict.getSequence(m.tid));
			String seqName=(ssr==null?null:ssr.getSequenceName());
			if(seqName==null) seqName="tid"+m.tid;
			System.out.print(String.valueOf(seqName+":"+(m.pos)));
			if(this.useSamFlag) System.out.print("="+m.flag);
			if(this.useCigar) System.out.print("/"+m.cigar);
			}
		if(first) System.out.print("(empty)");
		}
	
	@Override
	protected int doWork()
		{
		SAMFileReader samFileReader=null;
		try
			{
			if(this.IN.size() <2)
				{
				System.err.println("Need more bams please");
				return -1;
				}
			
			
			final MatchComparator matchComparator=new MatchComparator();
			SortingCollection<Match> database=SortingCollection.newInstance(
					Match.class,
					new MatchCodec(),
					matchComparator,
					super.MAX_RECORDS_IN_RAM
					);
			database.setDestructiveIteration(true);
	
			List<SAMSequenceDictionary> sequenceDictionaries=new ArrayList<SAMSequenceDictionary>(this.IN.size());
			
			for(int currentSamFileIndex=0;
					currentSamFileIndex<this.IN.size();
					currentSamFileIndex++ )
				{
				long nReads=0L;
				File samFile=this.IN.get(currentSamFileIndex);
				LOG.info("Opening "+samFile);
				samFileReader=new SAMFileReader(samFile);
				samFileReader.setValidationStringency(super.VALIDATION_STRINGENCY);
				SAMSequenceDictionary dict=samFileReader.getFileHeader().getSequenceDictionary();
				sequenceDictionaries.add(dict);
				
				Interval interval=null;
				if(REGION!=null)
					{
					interval=IntervalUtils.parseOne(dict, REGION);
					if(interval==null)
						{
						System.err.println("Cannot parse "+REGION+" (bad syntax or not in dictionary");
						return -1;
						}
					}
				
				Iterator<SAMRecord> iter=null;
				if(interval==null)
					{
					iter=samFileReader.iterator();
					}
				else
					{
					iter=samFileReader.queryOverlapping(interval.getSequence(), interval.getStart(), interval.getEnd());
					}
				
				while(iter.hasNext() )
					{
					if(nReads++%10000000==0) LOG.info("in "+samFile+" count:"+nReads);
					SAMRecord rec=iter.next();
					Match m=new Match();
					if(rec.getReadPairedFlag())
						{
						m.num_in_pair=(rec.getFirstOfPairFlag()?1:2);
						}
					else
						{
						m.num_in_pair=0;
						}
					m.readName=rec.getReadName();
					m.bamIndex=currentSamFileIndex;
					m.flag=rec.getFlags();
					m.cigar=rec.getCigarString();
					if(m.cigar==null) m.cigar="";
					if(rec.getReadUnmappedFlag())
						{
						m.tid=-1;
						m.pos=-1;
						}
					else
						{
						m.tid=rec.getReferenceIndex();
						m.pos=rec.getAlignmentStart();
						}
					database.add(m);
					}
				samFileReader.close();
				samFileReader=null;
				LOG.info("Close "+samFile);
				}
			database.doneAdding();
			LOG.info("Writing results....");
			
			//compute the differences for each read
			System.out.print("#READ-Name\t");
			for(int x=0;x<this.IN.size();++x)
				{
				for(int y=x+1;y<this.IN.size();++y)
					{
					if(!(x==0 && y==1)) System.out.print("|");
					System.out.print(IN.get(x));
					System.out.print(" ");
					System.out.print(IN.get(y));
					}
				}
			for(int x=0;x<this.IN.size();++x)
				{
				System.out.print("\t"+IN.get(x));
				}
			System.out.println();
			
			/* create an array of set<Match> */
			List<Set<Match>> matches=new ArrayList<Set<CompareBams2.Match>>(this.IN.size());
			while(matches.size() < this.IN.size())
				{
				matches.add(new TreeSet<CompareBams2.Match>(matchComparator));
				}
			
			CloseableIterator<Match> iter=database.iterator();
			String currReadName=null;
			int curr_num_in_pair=-1;
			for(;;)
				{
				Match nextMatch = null;
				if(iter.hasNext())
					{
					nextMatch = iter.next();
					}
				if(nextMatch==null ||
					(currReadName!=null && !currReadName.equals(nextMatch.readName)) ||
					(curr_num_in_pair!=-1 && curr_num_in_pair!=nextMatch.num_in_pair))
					{
					if(currReadName!=null)
						{
						System.out.print(currReadName);
						if(curr_num_in_pair>0)
							{
							System.out.print("/");
							System.out.print(curr_num_in_pair);
							}
						System.out.print("\t");
						
						
						for(int x=0;x<this.IN.size();++x)
							{
							Set<Match> first=matches.get(x);
							for(int y=x+1;y<this.IN.size();++y)
								{
								if(!(x==0 && y==1)) System.out.print("|");
								Set<Match> second=matches.get(y);
								if(first.size()==second.size() && first.containsAll(second))
									{
									System.out.print("EQ");
									}
								else
									{
									System.out.print("NE");
									}
								}
							}
	
						for(int x=0;x<this.IN.size();++x)
							{
							System.out.print("\t");
							print(matches.get(x),sequenceDictionaries.get(x));
							}
						
						System.out.println();
						}
					if(nextMatch==null) break;
					for(Set<Match> set:matches) set.clear();
					}
				currReadName=nextMatch.readName;
				curr_num_in_pair=nextMatch.num_in_pair;
				matches.get(nextMatch.bamIndex).add(nextMatch);
				}
			
			iter.close();
			}
		catch(Exception err)
			{
			err.printStackTrace();
			return -1;
			}
		finally
			{
			if(samFileReader!=null) samFileReader.close();
			}
		return 0;
		}
		
	public static void main(String[] args) throws Exception
		{
		new CompareBams2().instanceMainWithExit(args);
		}
}
