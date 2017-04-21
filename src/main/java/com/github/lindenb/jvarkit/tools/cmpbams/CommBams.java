package com.github.lindenb.jvarkit.tools.cmpbams;

import java.io.File;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.PeekableIterator;
/**
 BEGIN_DOC
 
 END_DOC
 */
@Program(name="commbams",description="Equivalent of unix 'comm' for bams sorted on queryname")
public class CommBams extends Launcher {
	private static final Logger LOG=Logger.build(CommBams.class).make();
	@Parameter(names={"-o","--out"},description="output file . Default:stdout")
	private File outputFile = null;
	@Parameter(names={"-1","--hide1"},description="suppress read unique to file 1")
	private boolean hide1=false;
	@Parameter(names={"-2","--hide2"},description="suppress read unique to file 2")
	private boolean hide2=false;
	@Parameter(names={"-3","--hide3"},description="suppress reads present in both files")
	private boolean hide3=false;
	@Parameter(names={"-st","--samtools"},description="Data was sorted using samtools sort -n algorithm (!= picard) see https://github.com/samtools/hts-specs/issues/5")
	private boolean samtoolsquerysort = false;
	@Parameter(names={"-delim","--delimiter"},description="Output delimiter")
	private String delim = "\t";
	@Parameter(names={"-empty","--empty"},description="Empty content symbol")
	private String emptySymbol = ".";

	
	
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
		boolean first=true;
		for(int i=1;i<=3;++i)
			{
			if(i==1 && hide1) continue;
			if(i==2 && hide2) continue;
			if(i==3 && hide3) continue;
			if(!first) out.print(this.delim);
			if(i==col)
				{
				out.print(name);
				}
			else
				{
				out.print(this.emptySymbol);
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
		Optional<SAMRecord> rec0  = Optional.empty();
		Optional<SAMRecord> rec1  = Optional.empty();
		
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
				}
			out= super.openFileOrStdoutAsPrintWriter(outputFile);
			
			for(;;) {
				for(int i=0;i< 2;++i) {
					Optional<SAMRecord> optRec = (i==0?rec0:rec1);
					if(!optRec.isPresent())
						{
						while(iters[i].hasNext()) {
							final SAMRecord rec = iters[i].peek();
							
							if(optRec.isPresent() && comparator.compare(optRec.get(),rec)>0)
								{
								LOG.error("Something is wrong in sort order of "+args.get(i)+" : got\n\t"
										+rec+" "+side(rec)+"\nafter\n\t"+ optRec.get()+" "+side(optRec.get())+"\nSee also option (samtools querysort)"
										);
								return -1;
								}
							// equals
							else if( optRec.isPresent()&& comparator.compare(optRec.get(),rec)==0)
								{
								iters[i].next();//consumme
								}
							//it's a new
							else if(!optRec.isPresent()){
								optRec = Optional.of(iters[i].next());//consumme
								if(i==0)
									{
									rec0=optRec;
									}
								else
									{
									rec1=optRec;
									}
								}
							// compare <0
							else
								{	
								break;
								}
							}
						}
				}
				
				
				if(!rec0.isPresent() && !rec1.isPresent()) break;
	
						
				if((!rec0.isPresent() && rec1.isPresent()) ||
				   (rec0.isPresent() && rec1.isPresent() && comparator.compare(rec0.get(),rec1.get())>0)
					)
					{
					if(!hide2) dump(out,readName(rec1.get()),2);
					rec1 =Optional.empty();
					}
				else if((rec0.isPresent() && !rec1.isPresent()) ||
						(rec0.isPresent() && rec1.isPresent() && comparator.compare(rec0.get(), rec1.get())<0))
					{
					if(!hide1) dump(out,readName(rec0.get()),1);
					rec0 =Optional.empty();
					}
				else
					{
					if(!hide3) dump(out,readName(rec0.get()),3);
					rec0 =Optional.empty();
					rec1 =Optional.empty();
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
