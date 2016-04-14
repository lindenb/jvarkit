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
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMSequenceDictionaryHelper;

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
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(CompareBams2.class);

	private boolean samSequenceDictAreTheSame=true;
	private List<SAMSequenceDictionary> sequenceDictionaries=new ArrayList<SAMSequenceDictionary>();
	private PrintWriter out;
	
	private class MatchComparator
		implements Comparator<Match>
		{
		@Override
		public int compare(final Match m0, final Match m1)
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
		public int compare(final Match m0, final Match m1)
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
		public Match decode(final DataInputStream dis) throws IOException
			{
			final Match m=new Match();
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
		public void encode(final DataOutputStream dos, final Match match)
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
		public boolean equals(final Object obj)
			{
			if (this == obj) { return true; }
			if (obj == null) { return false; }
			final Match other = (Match) obj;
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
	
	private int compare(final String s1,final String s2)
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
		final String chrom1=this.sequenceDictionaries.get(file_id1).getSequence(tid1).getSequenceName();
		final String chrom2=this.sequenceDictionaries.get(file_id2).getSequence(tid2).getSequenceName();
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
			if(!first)this.out.print(',');
			first=false;
			if(m.tid<0){ this.out.print("unmapped"); continue;}
			final SAMSequenceRecord ssr=(dict==null?null:dict.getSequence(m.tid));
			String seqName=(ssr==null?null:ssr.getSequenceName());
			if(seqName==null) seqName="tid"+m.tid;
			this.out.print(String.valueOf(seqName+":"+(m.pos)));
			if(super.useSamFlag) this.out.print("="+m.flag);
			if(super.useCigar) this.out.print("/"+m.cigar);
			}
		if(first) this.out.print("(empty)");
		}
	
	
	
	
	
	private final List<File> IN=new ArrayList<File>();
	
    private boolean same(final Set<Match> set1,final Set<Match> set2)
    	{
    	for(final Match m0:set1)
    		{
    		for(final Match m1:set2)
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
    public Collection<Throwable> call() throws Exception {
    	this.IN.clear();
		for(final String s:getInputFiles())
			{
			this.IN.add(new File(s));
			}
		
		try
			{
			return doWork();
			}
		catch(final Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			this.IN.clear();
			}
		}

	private Collection<Throwable> doWork()
		{
		SortingCollection<Match> database = null;
		SamReader samFileReader=null;
		CloseableIterator<Match> iter=null;
		try
			{
			if(this.IN.size() <2)
				{
				return wrapException("Need more bams please");
				}
			
			database = SortingCollection.newInstance(
					Match.class,
					new MatchCodec(),
					new MatchOrderer(),
					super.getMaxRecordsInRam(),
					super.getTmpDirectories()
					);
			this.samSequenceDictAreTheSame=true;
			database.setDestructiveIteration(true);
	
			
			for(int currentSamFileIndex=0;
					currentSamFileIndex<this.IN.size();
					currentSamFileIndex++ )
				{
				File samFile=this.IN.get(currentSamFileIndex);
				LOG.info("Opening "+samFile);
				samFileReader= super.createSamReaderFactory().open(samFile);
				final SAMSequenceDictionary dict=samFileReader.getFileHeader().getSequenceDictionary();
				if(dict==null || dict.isEmpty())
					{
					return wrapException("Empty Dict  in "+samFile);
					}
				
				if(!this.sequenceDictionaries.isEmpty() && !SequenceUtil.areSequenceDictionariesEqual(this.sequenceDictionaries.get(0), dict))
					{
					this.samSequenceDictAreTheSame=false;
					LOG.warn("FOOL !! THE SEQUENCE DICTIONARIES ARE **NOT** THE SAME. I will try to compare anyway but it will be slower.");
					}
				this.sequenceDictionaries.add(dict);
				
				
				final Optional<Interval> interval;
				if(REGION!=null && !REGION.trim().isEmpty())
					{
					final SAMSequenceDictionaryHelper dix = new SAMSequenceDictionaryHelper();
					interval = dix.parseInterval(REGION);
					
					if(!interval.isPresent())
						{
						return wrapException("Cannot parse "+REGION+" (bad syntax or not in dictionary)");
						}
					}
				else
					{
					interval = Optional.empty();
					}
				
				
				SAMRecordIterator it=null;
				if(!interval.isPresent())
					{
					it=samFileReader.iterator();
					}
				else
					{
					it=samFileReader.queryOverlapping(
							interval.get().getContig(),
							interval.get().getStart(),
							interval.get().getEnd()
							);
					}
				final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
				while(it.hasNext() )
					{
					final SAMRecord rec=progress.watch(it.next());
					if(!rec.getReadUnmappedFlag())
						{
						if(rec.getMappingQuality() < this.min_mapq) continue;
						if(rec.isSecondaryOrSupplementary()) continue;
						}
					final Match m=new Match();
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
					if(m.cigar==null ) m.cigar="";
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
				it.close();
				samFileReader.close();
				samFileReader=null;
				LOG.info("Close "+samFile);
				}
			database.doneAdding();
			LOG.info("Writing results....");
			
			this.out = super.openFileOrStdoutAsPrintWriter();
			
			//compute the differences for each read
			this.out.print("#READ-Name\t");
			for(int x=0;x<this.IN.size();++x)
				{
				for(int y=x+1;y<this.IN.size();++y)
					{
					if(!(x==0 && y==1)) this.out.print("|");
					this.out.print(IN.get(x));
					this.out.print(" ");
					this.out.print(IN.get(y));
					}
				}
			for(int x=0;x<this.IN.size();++x)
				{
				this.out.print("\t"+IN.get(x));
				}
			this.out.println();
			
			/* create an array of set<Match> */
			final MatchComparator match_comparator=new MatchComparator();
			final List<Set<Match>> matches=new ArrayList<Set<CompareBams2.Match>>(this.IN.size());
			while(matches.size() < this.IN.size())
				{
				matches.add(new TreeSet<CompareBams2.Match>(match_comparator));
				}
			
			iter = database.iterator();
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
						this.out.print(currReadName);
						if(curr_num_in_pair>0)
							{
							this.out.print("/");
							this.out.print(curr_num_in_pair);
							}
						this.out.print("\t");
						
						
						for(int x=0;x<this.IN.size();++x)
							{
							final Set<Match> first=matches.get(x);
							for(int y=x+1;y<this.IN.size();++y)
								{
								if(!(x==0 && y==1)) this.out.print("|");
								Set<Match> second=matches.get(y);
								if(same(first,second))
									{
									this.out.print("EQ");
									}
								else
									{
									this.out.print("NE");
									}
								}
							}
	
						for(int x=0;x<this.IN.size();++x)
							{
							this.out.print("\t");
							print(matches.get(x),sequenceDictionaries.get(x));
							}
						
						this.out.println();
						}
					if(nextMatch==null) break;
					for(Set<Match> set:matches) set.clear();
					}
				currReadName=nextMatch.readName;
				curr_num_in_pair=nextMatch.num_in_pair;
				matches.get(nextMatch.bamIndex).add(nextMatch);
				if(this.out.checkError()) break;
				}
			
			iter.close();
			this.out.flush();
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			if(database!=null) database.cleanup();
			CloserUtil.close(samFileReader);
			CloserUtil.close(this.out);this.out=null;
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		new CompareBams2().instanceMainWithExit(args);
		}
}
