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
package com.github.lindenb.jvarkit.tools.cmpbams;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;
import com.github.lindenb.jvarkit.util.picard.SortingCollectionFactory;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;


public class CompareBams2  extends AbstractCompareBams2
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(CompareBams2.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractCompareBams2.AbstractCompareBams2Command
		{    

	
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
			i= compareTid(m0.bamIndex,m0.tid ,m1.bamIndex, m1.tid);
			if(i!=0) return i;
			i= m0.pos - m1.pos;
			if(i!=0) return i;
			i= m0.flag - m1.flag;
			if(i!=0) return i;
			i= m0.cigar.compareTo(m1.cigar);
			return 0;
			}
		}
	
	private class MatchOrderer
	implements Comparator<Match>
		{
		@Override
		public int compare(Match m0, Match m1)
			{
			int i=m0.readName.compareTo(m1.readName);
			if(i!=0) return i;
			i=m0.num_in_pair-m1.num_in_pair;
			return i;
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
	
	private class Match
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
			if (compareTid(this.bamIndex,tid,other.bamIndex,other.tid)!=0) { return false; }
			if(tid==-1) return true;
			if (pos != other.pos) { return false; }
			if (bamIndex != other.bamIndex) { return false; }
			if(!readName.equals(other.readName)) return false;
			return true;
			}
		}
	
	private String norm(String s1)
		{
		if(s1.startsWith("chr")) s1=s1.substring(3);
		if(s1.startsWith("0")) s1=s1.substring(1);
		if(s1.equals("MT")) s1="M";
		return s1;
		}
	
	private int compare(String s1,String s2)
		{
		return norm(s1).compareToIgnoreCase(norm(s2));
		}
	
	private int compareTid(int file_id1,int tid1,int file_id2,int tid2)
		{
		if(samSequenceDictAreTheSame) return tid1-tid2;
		if(tid1==-1)
			{
			return tid2==-1?0:-1;
			}
		if(tid2==-1)
			{
			return 1;
			}
		String chrom1=this.sequenceDictionaries.get(file_id1).getSequence(tid1).getSequenceName();
		String chrom2=this.sequenceDictionaries.get(file_id2).getSequence(tid2).getSequenceName();
		if(chrom1==null)
			{
			return chrom2==null?0:-1;
			}
		if(chrom2==null) return 1;
		return compare(chrom1,chrom2);
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
	
	
	private boolean samSequenceDictAreTheSame=true;
	private List<SAMSequenceDictionary> sequenceDictionaries=new ArrayList<SAMSequenceDictionary>();
	private List<File> IN=new ArrayList<File>();
    private SortingCollectionFactory<Match> sortingFactory=new SortingCollectionFactory<Match>();
	
    private boolean same(Set<Match> set1,Set<Match> set2)
    	{
    	for(Match m0:set1)
    		{
    		for(Match m1:set2)
	    		{
    			int i=m0.readName.compareTo(m1.readName);
    			if(i!=0) continue;
    			i=m0.num_in_pair-m1.num_in_pair;
    			if(i!=0) continue;
    			//i= m0.bamIndex - m1.bamIndex;//NO ! (when comparing two Set<Match>)
    			//if(i!=0) return i;
    			i= m0.tid - m1.tid;
    			if(i!=0) continue;
    			i= Math.abs(m0.pos - m1.pos);
    			if(i>this.distance_tolerance) continue;
    			i= m0.flag - m1.flag;
    			if(i!=0) continue;
    			i= m0.cigar.compareTo(m1.cigar);
    			if(i!=0) continue;
    			return true;
	    		}
    		}
    	return false;
    	}
    
    @Override
    public Collection<Throwable> call() throws Exception
    	{
    	List<String> args = super.getInputFiles();	
		this.sortingFactory.setTmpDirs(this.getTmpDirectories());
		for(String arg: args)
			{
			this.IN.add(new File(arg));
			}
		
		try
			{
			return doWork();
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			
			}
		}

	private Collection<Throwable> doWork()
			{
			SamReader samFileReader=null;
			try
				{
				if(this.IN.size() <2)
					{
					wrapException("Need more bams please");
					}
				
				
				
				this.sortingFactory.setComparator(new MatchOrderer());
				this.sortingFactory.setCodec(new MatchCodec());
				this.sortingFactory.setComponentType(Match.class);
				this.samSequenceDictAreTheSame=true;
				SortingCollection<Match> database=this.sortingFactory.make();
				database.setDestructiveIteration(true);
		
				
				for(int currentSamFileIndex=0;
						currentSamFileIndex<this.IN.size();
						currentSamFileIndex++ )
					{
					File samFile=this.IN.get(currentSamFileIndex);
					LOG.info("Opening "+samFile);
					samFileReader=SamFileReaderFactory.mewInstance().open(samFile);
					SAMSequenceDictionary dict=samFileReader.getFileHeader().getSequenceDictionary();
					if(dict.isEmpty())
						{
						return wrapException("Empty Dict  in "+samFile);
						}
					
					if(!this.sequenceDictionaries.isEmpty() && !SequenceUtil.areSequenceDictionariesEqual(this.sequenceDictionaries.get(0), dict))
						{
						this.samSequenceDictAreTheSame=false;
						LOG.warn("FOOL !! THE SEQUENCE DICTIONARIES ARE **NOT** THE SAME. I will try to compare anyway but it will be slower.");
						}
					sequenceDictionaries.add(dict);
					
					
					Interval interval=null;
					if(REGION!=null)
						{
						interval=IntervalUtils.parseOne(dict, REGION);
						if(interval==null && REGION.startsWith("chr"))
							interval=IntervalUtils.parseOne(dict, REGION.substring(3));
						if(interval==null)
							interval=IntervalUtils.parseOne(dict, "chr"+REGION);
						if(interval==null)
							{
							return wrapException("Cannot parse "+REGION+" (bad syntax or not in dictionary)");
							}
						}
					
					SAMRecordIterator iter=null;
					if(interval==null)
						{
						iter=samFileReader.iterator();
						}
					else
						{
						iter=samFileReader.queryOverlapping(interval.getContig(), interval.getStart(), interval.getEnd());
						}
					SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
					while(iter.hasNext() )
						{
						SAMRecord rec=iter.next();
						progress.watch(rec);
						if(rec.getReadUnmappedFlag())
							{
							if(rec.getMappingQuality() < this.min_mapq) continue;
							if(rec.isSecondaryOrSupplementary()) continue;
							}
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
					iter.close();
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
				final MatchComparator match_comparator=new MatchComparator();
				List<Set<Match>> matches=new ArrayList<Set<MyCommand.Match>>(this.IN.size());
				while(matches.size() < this.IN.size())
					{
					matches.add(new TreeSet<MyCommand.Match>(match_comparator));
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
									if(same(first,second))
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
					if(System.out.checkError()) break;
					}
				
				iter.close();
				return RETURN_OK;
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(samFileReader);
				}
			}
		}
	
	public static void main(String[] args) throws Exception
		{
		new CompareBams2().instanceMainWithExit(args);
		}
}
