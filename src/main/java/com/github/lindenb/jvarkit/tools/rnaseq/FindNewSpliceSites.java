package com.github.lindenb.jvarkit.tools.rnaseq;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.tribble.readers.LineIterator;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene.Exon;

public class FindNewSpliceSites extends AbstractCommandLineProgram
	{
	private IntervalTreeMap<KnownGene> knownGenesMap=new IntervalTreeMap<KnownGene>();
	private int max_distance=10;
	private SAMFileWriter sfw=null;
	private SAMFileWriter weird=null;

	private FindNewSpliceSites()
		{
		}
	
	private boolean is_close_to(int d1,int d2)
		{
		return Math.abs(d1-d2)<=this.max_distance;
		}
	


	private boolean findJunction(
			Collection<KnownGene> genes,
			int start1,
			int end1
			)
		{
		for(KnownGene g:genes)
			{
			for(int k=0;k+1< g.getExonCount();++k)
				{
				Exon ex0=g.getExon(k);
				Exon ex1=g.getExon(k+1);
				if( is_close_to(ex0.getEnd()+1,start1) &&
					is_close_to(ex1.getStart()+1,end1))
					{
					return true;
					}
				}
			}
		return false;
		}
	private static boolean isMatch(CigarElement e)
		{
		switch(e.getOperator())
			{
			case X:case EQ: case M: return true;
			default: return false;
			}
		}
	private void scanRead(
			SAMRecord rec,
			SAMSequenceDictionary dict
			)
			{
			Cigar cigar=rec.getCigar();
			//if(cigar==null || !rec.getCigarString().contains("N")) return; //aleady checked

		
			
			Interval interval=new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
			Collection<KnownGene> genes=this.knownGenesMap.getOverlapping(interval);
			if(genes.isEmpty())
				{
				return;
				}

			
			int refPos1=rec.getAlignmentStart();
			
			
			for(int cIdx=0;cIdx< cigar.numCigarElements();++cIdx)
				{
				CigarElement ce=cigar.getCigarElement(cIdx);
				switch(ce.getOperator())
					{
					case S: break;
					case I: break;
					case N:
						{
						if(cIdx+1<cigar.numCigarElements() &&
							isMatch(cigar.getCigarElement(cIdx+1)) &&	
							!findJunction(genes,refPos1-1,refPos1+ce.getLength()))
							{
							this.sfw.addAlignment(rec);//unknown junction
							return;
							}
						refPos1+=ce.getLength();	
						break;
						}
					case D:
					case X:
					case EQ:
					case M:
						{
						refPos1+=ce.getLength();
						break;
						}
					case H:case P: break;//ignore
					default:throw new RuntimeException("operator not handled. ops.");
					}
				}
			
			
			}

	
	
	private static boolean isWeird(SAMRecord rec,SAMSequenceDictionary dict)
		{
		if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag() && 
				rec.getReferenceIndex().equals(rec.getMateReferenceIndex()) &&
				(
				rec.getReadNegativeStrandFlag()==rec.getMateNegativeStrandFlag() ||
				(rec.getReadNegativeStrandFlag()&& !rec.getMateNegativeStrandFlag() && rec.getAlignmentStart() < rec.getMateAlignmentStart()) ||
				(!rec.getReadNegativeStrandFlag() && rec.getMateNegativeStrandFlag() && rec.getAlignmentStart() > rec.getMateAlignmentStart())
				))
			{
			if(rec.getAlignmentStart() < rec.getMateAlignmentStart())
				{
				if(rec.getMateAlignmentStart() < rec.getAlignmentEnd()) return true;
				}
			if(rec.getAlignmentStart() > rec.getMateAlignmentStart())
				{
				if(rec.getMateAlignmentStart() + Math.abs(rec.getInferredInsertSize() ) > rec.getAlignmentStart()) return true;
				}
			return true;
			}
		return false;
		}
	
	private void scan(SamReader in) 
		{
		SAMSequenceDictionary dict=in.getFileHeader().getSequenceDictionary();
		if(dict==null) throw new PicardException("Sequence dictionary missing");
		SAMRecordIterator iter=in.iterator();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
		
		while(iter.hasNext())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.isSecondaryOrSupplementary()) continue;
			progress.watch(rec);
			
			if(isWeird(rec,dict))
				{
				this.weird.addAlignment(rec);
				continue;
				}
			
			for(CigarElement ce:rec.getCigar().getCigarElements())
				{
				if(ce.getOperator().equals(CigarOperator.N))
					{
					scanRead(rec,dict);
					break;
					}
				}	
			}
		iter.close();
		progress.finish();
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FindNewSpliceSites";
		}

	@Override
	public String getProgramDescription()
		{
		return "use the 'N' operator in the cigar string to find unknown splice sites";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-k (uri) UCSC Known Gene URI. Required. May be specified multiple times");
		out.println("-d (int) max distance between known splice site and cigar end default:"+this.max_distance);
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		SAMFileWriterFactory swf=new SAMFileWriterFactory();
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		Set<String> kgUris=new HashSet<String>();
		while((c=opt.getopt(args,getGetOptDefault()+"k:d:"))!=-1)
			{
			switch(c)
				{
				case 'k': kgUris.add( opt.getOptArg() );break;
				case 'd': max_distance=Integer.parseInt(opt.getOptArg());break;
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
		
		if(kgUris.isEmpty())
			{
			error("known Gene file undefined");
			return -1;
			}
		
		SamReader sfr=null;
		try
			{

			Pattern tab=Pattern.compile("[\t]");
			for(String kgUri: kgUris)
				{
				info("Opening "+kgUri);
				LineIterator r=IOUtils.openURIForLineIterator(kgUri);
				while(r.hasNext())
					{
					KnownGene g=new KnownGene(tab.split(r.next()));
					if(g.getExonCount()==1) continue;//need spliced one
					this.knownGenesMap.put(new Interval(g.getContig(), g.getTxStart()+1, g.getTxEnd()), g);
					}
				info("Done reading: "+kgUri);
				}
			if(opt.getOptInd()==args.length)
				{
				sfr=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				sfr=SamFileReaderFactory.mewInstance().open(filename);
				}
			else
				{
				error("Illegal number of args");
				return -1;
				}
			SAMFileHeader header=sfr.getFileHeader().clone();
			SAMProgramRecord p=header.createProgramRecord();
			p.setCommandLine(getProgramCommandLine());
			p.setProgramVersion(getVersion());
			p.setProgramName(getProgramName());
			this.sfw=swf.makeSAMWriter(header, true,System.out);
			
			header=sfr.getFileHeader().clone();
			p=header.createProgramRecord();
			p.setCommandLine(getProgramCommandLine());
			p.setProgramVersion(getVersion());
			p.setProgramName(getProgramName());
			this.weird=swf.makeSAMWriter(header,true, new NullOuputStream());
			
			scan(sfr);
			sfr.close();
			info("Done");
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
			CloserUtil.close(this.sfw);
			CloserUtil.close(this.weird);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FindNewSpliceSites().instanceMainWithExit(args);

	}

}
