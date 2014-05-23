package com.github.lindenb.jvarkit.tools.misc;

import java.awt.Color;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.ucsc.PslAlign;

public class SamToPsl extends AbstractCommandLineProgram
	{
	private PrintWriter out=new PrintWriter(System.out);
	private boolean handle_paired_reads=true;
	private boolean output_bed12=false;
	private SamToPsl()
		{
		}
	
	
	
	private static class Align
		{
		int read0;
		int ref0;
		int len;
		Align(int read0,int ref0,int len)
			{
			this.read0=read0;
			this.ref0=ref0;
			this.len=len;
			}
		}
	private static List<Align> scancigar(int read0,int ref0,SAMRecord rec)
		{
		List<Align> L=new ArrayList<Align>();
		for(CigarElement ce:rec.getCigar().getCigarElements())
			{
			switch(ce.getOperator())
				{
				case I:
					{
					read0+=ce.getLength();
					break;
					}
				case N:
				case D:
					{
					ref0+=ce.getLength();
					break;
					}
				case M:
				case X:
				case EQ:
					{
					L.add(new Align(read0, ref0, ce.getLength()));
					read0+=ce.getLength();
					ref0+=ce.getLength();
					break;
					}
				case P:
				case H:
				case S:
					{
					break;
					}
				default:
					{
					throw new IllegalArgumentException("No handled:"+rec.getCigarString());
					}
				}
			}
		return L;
		}
	
	
	private List<PslAlign> makePslAlign(SAMRecord rec,SAMSequenceDictionary dict)
		{
		List<PslAlign> aligns=new ArrayList<PslAlign>();

		PslAlign a=new PslAlign();
		
		
		
		a.setStrand(rec.getReadNegativeStrandFlag()?'-':'+');
		a.setQName(rec.getReadName()+(rec.getReadPairedFlag()?(rec.getFirstOfPairFlag()?"/1":"/2"):"")+"_"+rec.getFlags());
		
		a.setTName(rec.getReferenceName());
		SAMSequenceRecord ssr=dict.getSequence(rec.getReferenceIndex());
		if(ssr==null)
			{
			throw new PicardException("Cannot get SAMSequenceRecord dict for "+rec.getReferenceName());
			}
		a.setTSize( ssr.getSequenceLength());
		
		int qBaseInsert=0;
		int tBaseInsert=0;
		int readLength=0;
		for(CigarElement ce:rec.getCigar().getCigarElements())
			{
			switch(ce.getOperator())
				{
				case I: qBaseInsert+=ce.getLength();readLength+=ce.getLength();break;
				case N:
				case D: tBaseInsert+=ce.getLength();break;
				case M://threw
				case X:
				case EQ:
				case H:
				case S:
					{
					readLength+=ce.getLength();
					break;
					}
				default:break;
				}
			}
		a.setQBaseInsert(qBaseInsert);
		a.setTBaseInsert(tBaseInsert);
		
		boolean treat_as_paired=this.handle_paired_reads;
		if(rec.getSupplementaryAlignmentFlag()) treat_as_paired=false;

		if(rec.getReadPairedFlag())
			{
			if(rec.getMateUnmappedFlag())
				{
				treat_as_paired=false;
				}
			else
				{
				if(!rec.getMateReferenceIndex().equals(rec.getReferenceIndex()))
					{
					treat_as_paired=false;
					}
				else
					{
					if(rec.getReadNegativeStrandFlag()==rec.getMateNegativeStrandFlag())
						treat_as_paired=false;
					if(rec.getAlignmentStart() <= rec.getMateAlignmentStart())
						{
						if(rec.getMateAlignmentStart() <= rec.getAlignmentEnd())
							treat_as_paired=false;
						}
					}
				}
			}
		else
			{
			treat_as_paired=false;
			}
		
		
		if(!treat_as_paired)
			{
			a.setMatches(rec.getCigar().getReferenceLength());
			a.setQSize(readLength);
			a.setQEnd(readLength-(rec.getUnclippedEnd()-rec.getAlignmentEnd()));
			a.setTEnd(rec.getAlignmentEnd());
			
			int readPos=rec.getAlignmentStart()-rec.getUnclippedStart();
			int refPos=rec.getAlignmentStart()-1;
			a.setQStart(readPos);
			a.setTStart(refPos);
			
			for(Align ctg:scancigar(readPos,refPos,rec))
				{
				a.addBlock(ctg.read0, ctg.ref0, ctg.len);
				}
			
			}
		else 
			{
			a.setMatches(rec.getCigar().getReferenceLength()+1);
			a.setQSize(readLength+1);
			
			
			
			if(rec.getAlignmentStart()< rec.getMateAlignmentStart())
				{
				
				a.setQEnd(readLength+1);
				a.setTEnd(rec.getMateAlignmentStart());
				
				int readPos=rec.getAlignmentStart()-rec.getUnclippedStart();
				int refPos=rec.getAlignmentStart()-1;
				a.setQStart(readPos);
				a.setTStart(refPos);
								
				for(Align ctg:scancigar(readPos,refPos,rec))
					{
					a.addBlock(ctg.read0, ctg.ref0, ctg.len);
					}
				
				a.addBlock(
						readLength,
						rec.getMateAlignmentStart()-1,
						1
						);
				}
			else
				{
				a.setQStart(0);
				a.setQEnd((readLength+1)-(rec.getUnclippedEnd()-rec.getAlignmentEnd()));
				a.setTStart(rec.getMateAlignmentStart()-1);
				a.setTEnd(rec.getAlignmentEnd());
				
				
				a.addBlock(
						0,
						rec.getMateAlignmentStart()-1,
						1
						);
				
				for(Align ctg:scancigar(1,rec.getAlignmentStart(),rec))
					{
					a.addBlock(ctg.read0, ctg.ref0, ctg.len);
					}
				}
			}
		
		aligns.add(a);
		
		return aligns;
		}
	
	private void scan(SamReader in) 
		{
		SAMSequenceDictionary dict=in.getFileHeader().getSequenceDictionary();
		if(dict==null) throw new PicardException("Sequence dictionary missing...");
		SAMRecordIterator iter=in.iterator();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);		
		while(iter.hasNext() && !this.out.checkError())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			progress.watch(rec);
			for(PslAlign a:makePslAlign(rec,dict))
				{
				out.println(toString(a,rec));
				}
			
			
			}
		progress.finish();
		iter.close();
		}
	
	private String toString(PslAlign a, SAMRecord rec)
		{
		if(output_bed12)
			{
			StringBuilder b=new StringBuilder();
		    //chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
			b.append(a.getTName());
			//chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
			b.append('\t');
			b.append(a.getTStart());
			//chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 
			b.append('\t');
			b.append(a.getTEnd());
			//name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
			b.append('\t');
			b.append(a.getTName());
			//score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray: shade 	  	  	  	  	  	  	  	  	 
			b.append('\t');
			b.append((int)(((rec.getMappingQuality()>=255?0:rec.getMappingQuality())/255.0))*1000);
			//strand - Defines the strand - either '+' or '-'.
			b.append('\t');
			b.append((char)a.getStrand());
			//thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
			b.append('\t');
			b.append(a.getTStart());
			//thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
			b.append('\t');
			b.append(a.getTEnd());
			//itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
			b.append('\t');
			Color c=Color.BLACK;
			if(rec.getReadPairedFlag())
				{
				c=Color.ORANGE;
				if(rec.getProperPairFlag()) c=Color.GREEN;
				}
			b.append(c.getRed()).append(",").append(c.getGreen()).append(",").append(c.getBlue());
			//blockCount - The number of blocks (exons) in the BED line.
			b.append('\t');
			b.append(a.getBlockCount());
			//blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
			b.append('\t');
			for(PslAlign.Block block:a.getBlocks())
				{
				b.append(block.size());
				b.append(',');
				}
			//blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 
			b.append('\t');
			for(PslAlign.Block block:a.getBlocks())
				{
				b.append(block.getTStart()-a.getTStart());
				b.append(',');
				}
			return b.toString();
			}
		else
			{
			return a.toString();
			}
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Convert SAM/BAM to PSL http://genome.ucsc.edu/FAQ/FAQformat.html#format2 or BED12";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/SamToPsl ";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.print("-s treat all reads as single end");
		out.print("-B export as BED 12");
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"sB"))!=-1)
			{
			switch(c)
				{
				case 's': this.handle_paired_reads=false;break;
				case 'B': this.output_bed12=true;break;
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
		
		SamReader sfr=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				sfr=SamFileReaderFactory.mewInstance().openStdin();
				scan(sfr);
				sfr.close();
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					sfr=SamFileReaderFactory.mewInstance().open(filename);
					scan(sfr);
					sfr.close();
					}
				}
			out.flush();
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
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamToPsl().instanceMainWithExit(args);
		}
	}
