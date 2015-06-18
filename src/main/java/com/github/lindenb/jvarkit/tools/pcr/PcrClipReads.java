package com.github.lindenb.jvarkit.tools.pcr;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
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
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.readers.LineIterator;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.picard.CigarIterator;
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
					sw.addAlignment(rec);
					continue;
					}
				Cigar cigar = rec.getCigar();
				if(cigar==null)
					{
					warning("cigar missing in "+rec);
					sw.addAlignment(rec);
					continue;
					}

				List<CigarElement> cigarlist = new ArrayList<>();
				//expand cigar	
				for(CigarElement ce:cigar.getCigarElements())
					{
					CigarOperator op = ce.getOperator();
					for(int x=0;x < ce.getLength();++x)
						{
						cigarlist.add(new CigarElement(1, op));
						}
					}
				if(rec.getAlignmentStart()< fragment.getStart())
					{
					int i=0;
					int refPos1 = rec.getUnclippedStart();
					while(i< cigarlist.size() && refPos1< fragment.getStart())
						{
						CigarElement ce = cigarlist.get(i);
						if(!ce.getOperator().consumesReadBases())
							{
							cigarlist.remove(i);
							continue;
							}

						
						++i;
						}
					}
				int read5 =0;
				
				
						
					
				
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
		super.printOptions(out);
		}

	
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
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
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
