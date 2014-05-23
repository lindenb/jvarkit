package com.github.lindenb.jvarkit.tools.cmpbams;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;
import com.github.lindenb.jvarkit.util.picard.SortingCollectionFactory;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;


public class CompareBamAndBuild  extends AbstractCommandLineProgram
	{
	private File bamFiles[]=new File[2];
	private SAMSequenceDictionary sequenceDictionaries[]=new SAMSequenceDictionary[2];
	private LiftOver liftOver=null;
    private String REGION=null;
    private SortingCollectionFactory<Match> sortingFactory=new SortingCollectionFactory<Match>();
	private int distance_tolerance=10;

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
			if(m0.firstBamFile!=m1.firstBamFile) throw new IllegalStateException();
			i=m0.tid-m1.tid;
			if(i!=0) return i;
			i=m0.pos-m1.pos;
			if(i!=0) return i;
			return 0;
			}
		}
	
	private class MatchOrdererInSortingCollection
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
			m.firstBamFile=dis.readBoolean();
			m.tid=dis.readInt();
			m.pos=dis.readInt();
			m.num_in_pair=dis.readInt();
			return m;
			}
		@Override
		public void encode(DataOutputStream dos, Match match)
				throws IOException
			{
			dos.writeUTF(match.readName);
			dos.writeBoolean(match.firstBamFile);
			dos.writeInt(match.tid);
			dos.writeInt(match.pos);
			dos.writeInt(match.num_in_pair);
			}
		
		}
	
	private class Match
		{
		String readName;
		int num_in_pair=0;
		int tid=-1;
		boolean firstBamFile=true;
		int pos=-1;

		@Override
		public int hashCode()
			{
			int result = 1;
			result = 31 * result + num_in_pair;
			result = 31 * result + pos;
			result = 31 * result + tid;
			result = 31 * result + readName.hashCode();
			result = 31 * result + (firstBamFile?0:1);
			return result;
			}
		
		String getChrom()
			{	
			if(tid==-1) return null;
			return sequenceDictionaries[firstBamFile?0:1].getSequence(tid).getSequenceName();
			}
		
		
		Interval getLiftOver()
			{
			if(tid==-1) return null;
			Interval src=new Interval(getChrom(),pos,pos);
			if(!firstBamFile) return src;
			return liftOver.liftOver(src);
			}
		
		@Override
		public boolean equals(Object obj)
			{
			if (this == obj) { return true; }
			if (obj == null) { return false; }
			Match other = (Match) obj;
			if (num_in_pair != other.num_in_pair) { return false; }
			if (compareStrings(getChrom(), other.getChrom())!=0) { return false; }
			if(tid==-1) return true;
			if (pos != other.pos) { return false; }
			if (firstBamFile != other.firstBamFile) { return false; }
			if(!readName.equals(other.readName)) return false;
			return true;
			}
		}
	
	private int compareStrings(String s1,String s2)
		{
		if(s1==null)
			{
			return s2==null?0:-1;
			}
		if(s2==null) return 1;
		return 0;
		}
	
	/*
	private int compareTid(boolean is_first1,int tid1,boolean is_first2,int tid2)
		{
		if(tid1==-1)
			{
			return tid2==-1?0:-1;
			}
		if(tid2==-1)
			{
			return 1;
			}
		String chrom1=this.sequenceDictionaries[is_first1?0:1].getSequence(tid1).getSequenceName();
		String chrom2=this.sequenceDictionaries[is_first2?0:1].getSequence(tid2).getSequenceName();
		if(chrom1==null)
			{
			return chrom2==null?0:-1;
			}
		if(chrom2==null) return 1;
		return compare(chrom1,chrom2);
		}*/
	
	private void print(final Set<Match> set)
		{
		boolean first=true;
		for(Match m:set)
			{
			if(!first)System.out.print(',');
			first=false;
			if(m.tid<0){ System.out.print("unmapped"); continue;}
			if(m.firstBamFile)
				{
				System.out.print(String.valueOf(m.getChrom()+":"+m.pos+"->"));
				}
			
			Interval interval=m.getLiftOver();
			if(interval==null)
				{
				System.out.print("liftoverfail");
				}
			else
				{
				System.out.print(String.valueOf(interval.getSequence()+":"+interval.getStart()));
				}
			}
		if(first) System.out.print("(empty)");
		}
	
	@Override
	public String getProgramDescription() {
		return "Compare two  BAM files mapped on two different builds. Requires a liftover chain file.";
		}
	
	
	
	
    private boolean same(Set<Match> set1,Set<Match> set2)
    	{
    	for(Match m0:set1)
    		{
    		Interval i0=m0.getLiftOver();
    		if(i0==null || i0.getSequence()==null) continue;
    		for(Match m1:set2)
	    		{
    			int i=m0.readName.compareTo(m1.readName);
    			if(i!=0) continue;
    			i=m0.num_in_pair-m1.num_in_pair;
    			if(i!=0) continue;
    			//i= m0.bamIndex - m1.bamIndex;//NO ! (when comparing two Set<Match>)
    			//if(i!=0) return i;
    			Interval i1=m1.getLiftOver();
        		if(i1==null || i1.getSequence()==null) continue;
    			
    			i= compareStrings(i0.getSequence(),i1.getSequence());
    			if(i!=0) continue;
    			i= Math.abs(i0.getStart() -i1.getStart());
    			if(i>this.distance_tolerance) continue;
    			return true;
	    		}
    		}
    	return false;
    	}
    
    @Override
    protected String getOnlineDocUrl() {
    	return "https://github.com/lindenb/jvarkit/wiki/CmpBamAndBuild";
    	}
    
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -c (chain file) Lift Over file from bam1 to bam2. REQUIRED.");
		out.println(" -d (dist) distance tolerance between two alignments.");
		out.println(" -r (region) restrict to that region chr:start-end");
		out.println(" -n (int) "+getMessageBundle("max.records.in.ram")+" Optional.");
		out.println(" -T (dir) "+getMessageBundle("add.tmp.dir")+" Optional.");
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"r:c:n:d:T:"))!=-1)
			{
			switch(c)
				{
				case 'd':distance_tolerance=Integer.parseInt(opt.getOptArg());break;
				case 'r':REGION=opt.getOptArg();break;
				case 'c':
					{
					info("Loading lift over file");
					this.liftOver=new LiftOver(new File(opt.getOptArg()));
					break;
					}
				case 'n': sortingFactory.setMaxRecordsInRAM(Integer.parseInt(opt.getOptArg()));break;
				case 'T': this.addTmpDirectory(new File(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(this.liftOver==null)
			{
			error("Undefined lift Over file");
			return -1;
			}
		if(opt.getOptInd()+2!=args.length)
			{
			error("Illegal number of arguments. Expected two indexed BAMS.");
			return -1;
			}
		this.sortingFactory.setTmpDirs(this.getTmpDirectories());
		
		try
			{

			
			
			
			this.sortingFactory.setComparator(new MatchOrdererInSortingCollection());
			this.sortingFactory.setCodec(new MatchCodec());
			this.sortingFactory.setComponentType(Match.class);
			SortingCollection<Match> database=this.sortingFactory.make();
			database.setDestructiveIteration(true);
	
			
			for(int currentSamFileIndex=0;
					currentSamFileIndex<2;
					currentSamFileIndex++ )
				{
				File samFile=new File(args[opt.getOptInd()+currentSamFileIndex]);
				this.bamFiles[currentSamFileIndex]=samFile;
				SamReader samFileReader=SamFileReaderFactory.mewInstance().open(samFile);
				SAMSequenceDictionary dict=samFileReader.getFileHeader().getSequenceDictionary();
				this.sequenceDictionaries[currentSamFileIndex]=dict;
				if(dict.isEmpty())
					{
					error("Empty Dict  in "+samFile);
					samFileReader.close();
					return -1;
					}
				
			
				Interval interval=null;
				if(REGION!=null)
					{
					interval=IntervalUtils.parseOne(dict, REGION);
					if(interval==null)
						{
						error("Cannot parse "+REGION+" (bad syntax or not in dictionary");
						samFileReader.close();
						return -1;
						}
					}
				
				SAMRecordIterator iter=null;
				if(interval==null)
					{
					iter=samFileReader.iterator();
					}
				else
					{
					iter=samFileReader.queryOverlapping(interval.getSequence(), interval.getStart(), interval.getEnd());
					}
				SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
				while(iter.hasNext() )
					{
					SAMRecord rec=iter.next();
					progress.watch(rec);
					if(rec.isSecondaryOrSupplementary()) continue;
					
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
					m.firstBamFile=currentSamFileIndex==0;
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
				info("Close "+samFile);
				}
			database.doneAdding();
			info("Writing results....");
			
			//compute the differences for each read
			System.out.print("#READ-Name\tCOMPARE");
			for(File f:this.bamFiles)
				{
				System.out.print("\t"+f);
				}
			System.out.println();
			
			/* create an array of set<Match> */
			final MatchComparator match_comparator=new MatchComparator();
			List<Set<Match>> matches=new ArrayList<Set<CompareBamAndBuild.Match>>(2);
			while(matches.size() < 2)
				{
				matches.add(new TreeSet<CompareBamAndBuild.Match>(match_comparator));
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
						
						if(same(matches.get(0),matches.get(1)))
							{
							System.out.print("EQ");
							}
						else
							{
							System.out.print("NE");
							}
						
	
						for(int x=0;x<2;++x)
							{
							System.out.print("\t");
							print(matches.get(x));
							}
						
						System.out.println();
						}
					if(nextMatch==null) break;
					for(Set<Match> set:matches) set.clear();
					}
				currReadName=nextMatch.readName;
				curr_num_in_pair=nextMatch.num_in_pair;
				matches.get(nextMatch.firstBamFile?0:1).add(nextMatch);
				}
			
			iter.close();
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

	public static void main(String[] args) throws Exception
		{
		new CompareBamAndBuild().instanceMainWithExit(args);
		}
}
