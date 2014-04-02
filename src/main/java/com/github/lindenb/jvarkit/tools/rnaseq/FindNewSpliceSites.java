package com.github.lindenb.jvarkit.tools.rnaseq;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.regex.Pattern;

import org.broad.tribble.readers.LineIterator;

import net.sf.picard.PicardException;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
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
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene.Exon;
import com.github.lindenb.jvarkit.util.ucsc.PslAlign;

public class FindNewSpliceSites extends AbstractCommandLineProgram
	{
	private IntervalTreeMap<KnownGene> knownGenesMap=new IntervalTreeMap<KnownGene>();
	private int max_distance=10;
	private int max_extend_gene=2000;
	private PrintWriter pslWriter=new PrintWriter(System.out);
	private PrintWriter weirdPslWriter=new PrintWriter(new NullOuputStream());
	private FindNewSpliceSites()
		{
		}
	
	private boolean is_close_to(int d1,int d2)
		{
		return Math.abs(d1-d2)<=this.max_distance;
		}
	
	
	private PslAlign createPslAlign(SAMRecord rec,SAMFileReader in)
		{
		PslAlign a=new PslAlign();
		
		a.setMatches(rec.getCigar().getReferenceLength()*2);
		
		a.setStrand('+');
		a.setQName(rec.getReadName());
		a.setQSize(rec.getReadLength()*2);
		a.setQStart(0);
		a.setQEnd(rec.getReadLength()*2);
		
		a.setTName(rec.getReferenceName());
		a.setTStart(Math.min( rec.getAlignmentStart(),rec.getMateAlignmentStart())-1);
		a.setTEnd(Math.max( rec.getAlignmentEnd(),rec.getMateAlignmentStart()+rec.getReadLength()));
		a.setTSize(in.getFileHeader().getSequenceDictionary().getSequence(rec.getReferenceIndex()).getSequenceLength());
		
		if(rec.getAlignmentStart() < rec.getMateAlignmentStart())
			{
			a.addBlock(0,rec.getAlignmentStart()-1,rec.getReadLength());
			a.addBlock(rec.getReadLength(),rec.getMateAlignmentStart()-1,rec.getReadLength());			
			}
		else
			{
			a.addBlock(0,rec.getMateAlignmentStart()-1,rec.getReadLength());
			a.addBlock(rec.getReadLength(),rec.getAlignmentStart()-1,rec.getReadLength());			
			}
		
		return a;
		}

	
	private PslAlign createPslAlign(SAMRecord rec,OtherCanonicalAlign xp,SAMFileReader in)
		{
		PslAlign a=new PslAlign();
		
		a.setMatches(rec.getCigar().getReferenceLength()+xp.getCigar().getReferenceLength());
		
		a.setStrand('+');
		a.setQName(rec.getReadName());
		a.setQSize(rec.getCigar().getReferenceLength()+xp.getCigar().getReferenceLength());
		a.setQStart(0);
		a.setQEnd(rec.getCigar().getReferenceLength()+xp.getCigar().getReferenceLength());
		
		a.setTName(rec.getReferenceName());
		a.setTStart(Math.min(rec.getAlignmentStart(),xp.getPos())-1);
		a.setTEnd(Math.max(rec.getAlignmentEnd(),xp.getAlignmentEnd())-1);
		a.setTSize(in.getFileHeader().getSequenceDictionary().getSequence(rec.getReferenceIndex()).getSequenceLength());
		
		if(rec.getAlignmentStart() < xp.getPos())
			{
			a.addBlock(0,rec.getAlignmentStart()-1,rec.getCigar().getReferenceLength());
			a.addBlock(rec.getCigar().getReferenceLength(),xp.getPos()-1,xp.getCigar().getReferenceLength());			
			}
		else
			{
			a.addBlock(0,xp.getPos()-1,xp.getCigar().getReferenceLength());
			a.addBlock(xp.getCigar().getReferenceLength(),rec.getAlignmentStart()-1,rec.getCigar().getReferenceLength());
			}
		
		return a;
		}

	
	
	private void scan(SAMFileReader in) 
		{
		long count_intergenic=0;
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
			List<OtherCanonicalAlign> xpAligns=new ArrayList<OtherCanonicalAlign>(xpFactory.getXPAligns(rec));
			if(xpAligns.isEmpty()) continue;
			Interval interval=new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
			List<PslAlign> weirdPslAlignments=new ArrayList<PslAlign>();
			
			Collection<KnownGene> genes=this.knownGenesMap.getOverlapping(interval);
			
			if(genes.isEmpty())
				{
				count_intergenic++;
				continue;
				}
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
					if(rec.getMateAlignmentStart() < rec.getAlignmentEnd()) continue;
					}
				if(rec.getAlignmentStart() > rec.getMateAlignmentStart())
					{
					if(rec.getMateAlignmentStart() + Math.abs(rec.getInferredInsertSize() ) > rec.getAlignmentStart()) continue;
					}
				weirdPslWriter.println(	createPslAlign(rec,in));
				continue;
				}
			
			
			//remove XP aligns if no overlap/ transcripts
			int i=0;
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
						weirdPslAlignments.add(createPslAlign(rec, xp,in));
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
				continue;
				}
			
			if(weirdPslAlignments.size()==xpAligns.size())
				{
				for(PslAlign a:weirdPslAlignments) 
					this.weirdPslWriter.println(a);
				continue;
				}
			List<PslAlign> novoPslAlignments=new ArrayList<PslAlign>();
			
			boolean found_known_junction=false;
			for(OtherCanonicalAlign xp:xpAligns)
				{
				if((xp.getStrand()=='-')!=rec.getReadNegativeStrandFlag()) continue;
				for(KnownGene g:genes)
					{
					for(int k=0;k+1< g.getExonCount();++k)
						{
						Exon ex0=g.getExon(k);
						Exon ex1=g.getExon(k+1);
						if( is_close_to(ex0.getEnd()+1,rec.getAlignmentEnd()) &&
							is_close_to(ex1.getStart()+1,xp.getPos()))
							{
							found_known_junction=true;
							}
						else if( is_close_to(ex0.getEnd()+1,xp.getAlignmentEnd()) &&
								 is_close_to(ex1.getStart()+1,rec.getAlignmentStart()))
							{
							found_known_junction=true;
							}
						if(found_known_junction) break;
						}
					if(found_known_junction) break;
					}
				if(found_known_junction) break;
				novoPslAlignments.add(createPslAlign(rec, xp,in));
				}
			if(found_known_junction) continue;
			for(PslAlign a:novoPslAlignments)
				{
				this.pslWriter.println(a);
				}
			}
		iter.close();
		progress.finish();
		info("Intergenic: "+count_intergenic);
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
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		String kgUri=null;
		
		while((c=opt.getopt(args,getGetOptDefault()+"k:d:w:"))!=-1)
			{
			switch(c)
				{
				case 'k': kgUri=opt.getOptArg();break;
				case 'd': max_distance=Integer.parseInt(opt.getOptArg());break;
				case 'w':
					{
					try
						{
						weirdPslWriter=new PrintWriter(new File(opt.getOptArg()));break;
						}
					catch(IOException err)
						{
						error(err);
						return -1;
						}
					}
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
			
			scan(sfr);
			sfr.close();
			info("Done");
			pslWriter.flush();
			weirdPslWriter.flush();
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
			CloserUtil.close(pslWriter);
			CloserUtil.close(weirdPslWriter);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FindNewSpliceSites().instanceMainWithExit(args);

	}

}
