package com.github.lindenb.jvarkit.tools.rnaseq;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.regex.Pattern;

import org.broad.tribble.readers.LineIterator;

import net.sf.picard.PicardException;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlign;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlignFactory;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamWriterFactory;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene.Exon;

public class FindNewSpliceSites extends AbstractCommandLineProgram
	{
	private IntervalTreeMap<KnownGene> knownGenesMap=new IntervalTreeMap<KnownGene>();
	private int max_distance=10;
	private int max_extend_gene=2000;
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
	
	private void topHat(
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
					case D:
					case N:
						{
						refPos1+=ce.getLength();	
						break;
						}
					case X:
					case EQ:
					case M:
						{
						if(cIdx>0 &&
							cigar.getCigarElement(cIdx-1).getOperator().equals(CigarOperator.N) &&	
							findJunction(genes,refPos1,refPos1+ce.getLength()))
							{
							return;//known transcript
							}
						refPos1+=ce.getLength();
						break;
						}
					case H:case P: break;//ignore
					default:throw new RuntimeException("operator not handled. ops.");
					}
				}
			
			this.sfw.addAlignment(rec);
			}

	
	private void bwaMem(
			SAMRecord rec,
			SAMSequenceDictionary dict,
			ArrayList<OtherCanonicalAlign> xpAligns)
		{
		if(xpAligns.isEmpty()) return;
		
		
		Interval interval=new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
		
		Collection<KnownGene> genes=this.knownGenesMap.getOverlapping(interval);
		
		if(genes.isEmpty())
			{
			return;
			}
		
		
		
		//remove XP aligns if no overlap/ transcripts
		int i=0;
		List<OtherCanonicalAlign>  weirdPslAlignments=new ArrayList<OtherCanonicalAlign>();
		while(i< xpAligns.size())
			{
			OtherCanonicalAlign xp=xpAligns.get(i);
			
			boolean found=false;
			for(KnownGene g:genes)
				{
				if(!rec.getReferenceIndex().equals(xp.getChromIndex())) continue;
				if(g.getTxEnd()+this.max_extend_gene < xp.getPos()) continue;
				if(g.getTxStart()>xp.getAlignmentEnd()+this.max_extend_gene) continue;
				found=true;
				
				if((xp.getStrand()=='-')!=rec.getReadNegativeStrandFlag())
					{
					weirdPslAlignments.add(xp);
					}
				break;
				}
			
			//xp overlap read
			if(found && !(xp.getAlignmentEnd()<rec.getAlignmentStart()||xp.getPos()>rec.getAlignmentEnd()) )
				{
				found=false;
				}
			
			if(found)
				{
				++i;
				}
			else
				{
				xpAligns.remove(i);
				}
			}
		if(xpAligns.isEmpty())
			{
			return;
			}
		
		if(weirdPslAlignments.size()==xpAligns.size())
			{
			this.weird.addAlignment(rec);
			return;
			}
		
		boolean found_known_junction=false;
		for(OtherCanonicalAlign xp:xpAligns)
			{
			if((xp.getStrand()=='-')!=rec.getReadNegativeStrandFlag()) continue;
			if( findJunction(genes,rec.getAlignmentEnd(),xp.getPos()) ||
				findJunction(genes,xp.getAlignmentEnd(),rec.getAlignmentStart())
				)
				{
				found_known_junction=true;
				break;
				}
			
			}
		if(found_known_junction) return;
		this.sfw.addAlignment(rec);
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
	
	private void scan(SAMFileReader in) 
		{
		in.setValidationStringency(ValidationStringency.LENIENT);
		SAMSequenceDictionary dict=in.getFileHeader().getSequenceDictionary();
		if(dict==null) throw new PicardException("Sequence dictionary missing");
		SAMRecordIterator iter=in.iterator();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
		OtherCanonicalAlignFactory xpFactory=new OtherCanonicalAlignFactory(in.getFileHeader());;
		
		
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
			
			boolean has_N=false;
			for(CigarElement ce:rec.getCigar().getCigarElements())
				{
				if(ce.getOperator().equals(CigarOperator.N))
					{
					has_N=true;
					break;
					}
				}	
			if(has_N)
				{
				topHat(rec, dict);
				}
			else
				{
				ArrayList<OtherCanonicalAlign> xpAligns=new ArrayList<OtherCanonicalAlign>(xpFactory.getXPAligns(rec));
				if(!xpAligns.isEmpty())
					{
					bwaMem(rec, dict, xpAligns);
					}
				}
			}
		iter.close();
		progress.finish();
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Find new splice sites";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-k (uri) UCSC Known Gene URI. Required.");
		out.println("-d (int) max distance between known splice site and cigar end default:"+this.max_distance);
		out.println("-g (int) max pb to extend gene. default:"+this.max_extend_gene);
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		SamWriterFactory swf=SamWriterFactory.newInstance();
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		String kgUri=null;
		
		while((c=opt.getopt(args,getGetOptDefault()+"k:d:g:"))!=-1)
			{
			switch(c)
				{
				case 'k': kgUri=opt.getOptArg();break;
				case 'd': max_distance=Integer.parseInt(opt.getOptArg());break;
				case 'g': max_extend_gene=Integer.parseInt(opt.getOptArg());break;
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
		
		if(kgUri==null)
			{
			error("known Gene file undefined");
			return -1;
			}
		
		SAMFileReader sfr=null;
		try
			{

			Pattern tab=Pattern.compile("[\t]");
			info("Opening "+kgUri);
			LineIterator r=IOUtils.openURIForLineIterator(kgUri);
			while(r.hasNext())
				{
				KnownGene g=new KnownGene(tab.split(r.next()));
				if(g.getExonCount()==1) continue;//need spliced one
				this.knownGenesMap.put(new Interval(g.getChr(), g.getTxStart()+1, g.getTxEnd()), g);
				}
			info("Done reading: "+kgUri);
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				sfr=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				sfr=new SAMFileReader(new File(filename));
				}
			else
				{
				error("Illegal number of args");
				return -1;
				}
			sfr.setValidationStringency(ValidationStringency.SILENT);
			SAMFileHeader header=sfr.getFileHeader().clone();
			SAMProgramRecord p=header.createProgramRecord();
			p.setCommandLine(getProgramCommandLine());
			p.setProgramVersion(getVersion());
			p.setProgramName(getProgramName());
			this.sfw=swf.make(header, System.out);
			
			header=sfr.getFileHeader().clone();
			p=header.createProgramRecord();
			p.setCommandLine(getProgramCommandLine());
			p.setProgramVersion(getVersion());
			p.setProgramName(getProgramName());
			this.weird=swf.make(header, new NullOuputStream());
			
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
