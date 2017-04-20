package com.github.lindenb.jvarkit.tools.cmpbams;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.PeekableIterator;

@Program(name="commbams",description="Equivalent of unix 'comm' for bams")
public class CommBams extends Launcher {
	private static final Logger LOG=Logger.build(CommBams.class).make();
	@Parameter(names={"-o","--out"},description="output file . Default:stdout")
	private File outputFile = null;
	@Parameter(names={"-1","--hide1"},description="suppress read unique to file 1")
	private boolean hide1=false;
	@Parameter(names={"-2","--hide2"},description="suppress read unique to file 2")
	private boolean hide2=false;
	@Parameter(names={"-3","--hide2"},description="suppress reads present in both files")
	private boolean hide3=false;
	@Parameter(names={"-st","--samtools"},description="Data was sorted using samtools sort -n algorithm (!= picard) see https://github.com/samtools/hts-specs/issues/5")
	private boolean samtoolsquerysort = false;

	private static int side(final SAMRecord rec)
		{
		return CompareBams4.side(rec);
		}
	
	private static String readName(final SAMRecord rec) {
		if(rec.getReadPairedFlag()) {
			if(rec.getFirstOfPairFlag()) return rec.getReadName()+"/1";
			if(rec.getSecondOfPairFlag()) return rec.getReadName()+"/2";
			throw new RuntimeException("Side for flag "+rec.getReadName()+":"+rec.getFlags()+"?");
		}
		else {
			return rec.getReadName();
		}
	}
	
	private void dump(PrintWriter out,final String name,int col)
		{
		boolean first=false;
		for(int i=1;i<=3;++i)
			{
			if(i==1 && hide1) continue;
			if(i==2 && hide2) continue;
			if(i==3 && hide3) continue;
			if(!first) out.print("\t");
			if(i==col)
				{
				out.print(name);
				}
			else
				{
				out.print(".");
				}
			first=false;
			}
		out.println();
		}
	
	@Override
	public int doWork(final List<String> args) {
		final Comparator<SAMRecord> comparator;
		if(args.size() !=2)
			{
			LOG.info("Expected two and only two bams please, but got "+args.size());
			return -1;
			}

		if( this.samtoolsquerysort) {
			LOG.info("using the samtools sort -n comparator");
			comparator = new CompareBams4.SamToolsReadNameComparator();		
			}
		else
			{
			comparator = new CompareBams4.SimpleReadNameComparator();
			}
		final SamReader samFileReaders[]={null,null};
		@SuppressWarnings("unchecked")
		final PeekableIterator<SAMRecord> iters[]=new PeekableIterator[]{null,null};
		PrintWriter out=null;
		final List<List<SAMRecord>> recordLists=new ArrayList<>();
		
		if(hide1 && hide2 && hide3) {
			LOG.error("all flags hide** are on");
			return -1;
		}
		try
			{
			for(int i=0;i< args.size() && i< samFileReaders.length;++i)
				{
				final String samFile=args.get(i);
				LOG.info("opening "+samFile);
				samFileReaders[i]=super.openSamReader(samFile);
				final SAMFileHeader header = samFileReaders[i].getFileHeader();
				if(header.getSortOrder()!=SAMFileHeader.SortOrder.queryname) {
					LOG.error("Expected "+samFile+" to be sorted on "+SAMFileHeader.SortOrder.queryname+" but got "+header.getSortOrder());
					return -1;
					}
				
				iters[i] = new PeekableIterator<>(samFileReaders[i].iterator());
				recordLists.add(new ArrayList<>());
				}
			out= super.openFileOrStdoutAsPrintWriter(outputFile);
			
			for(;;) {
				for(int i=0;i< 2;++i) {
					if(recordLists.get(i).isEmpty())
						{
						while(iters[i].hasNext()) {
							final SAMRecord rec = iters[i].peek();
							
							if(!recordLists.get(i).isEmpty() && comparator.compare(recordLists.get(i).get(0),rec)>0)
								{
								LOG.error("Something is wrong in sort order of "+args.get(i)+" : got\n\t"
										+rec+" "+side(rec)+"\nafter\n\t"+ recordLists.get(i).get(0)+" "+side(recordLists.get(i).get(0))+"\nSee also option (samtools querysort)"
										);
								return -1;
								}
							else if( recordLists.get(i).isEmpty() ||
								comparator.compare(recordLists.get(i).get(0),rec)==0)
								{
								recordLists.get(i).add(iters[i].next());
								}
							else
								{	
								break;
								}
							}
						}
				}
				
				final SAMRecord rec0=(recordLists.get(0).isEmpty()?null:recordLists.get(0).get(0));
				final SAMRecord rec1=(recordLists.get(1).isEmpty()?null:recordLists.get(1).get(0));
				
				if(rec0==null && rec1==null) break;
	
						
				if((rec0==null && rec1!=null) ||
				   (rec0!=null && rec1!=null && comparator.compare(rec0,rec1)>0)
					)
					{
					dump(out,readName(rec1),2);
					recordLists.get(1).clear();
					}
				else if((rec0!=null && rec1==null) ||
						(rec0!=null && rec1!=null && comparator.compare(rec0, rec1)<0))
					{
					dump(out,readName(rec0),1);
					recordLists.get(0).clear();
					}
				else
					{
					dump(out,readName(rec0),3);
					recordLists.get(1).clear();
					recordLists.get(0).clear();
					}
				}
			out.flush();
			out.close();out=null;
			return 0;
			}
		catch(Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			for(int i=0;i< samFileReaders.length;++i){
				CloserUtil.close(iters[i]);
				CloserUtil.close(samFileReaders[i]);
				}
			CloserUtil.close(out);
		}
		}
	
	public static void main(final String[] args) {
		new CommBams().instanceMainWithExit(args);

	}

}
