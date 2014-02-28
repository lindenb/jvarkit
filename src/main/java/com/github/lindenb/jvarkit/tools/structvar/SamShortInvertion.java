package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.PrintStream;
import java.util.LinkedList;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
//import com.github.lindenb.jvarkit.util.picard.SamWriterFactory;

public class SamShortInvertion extends AbstractCommandLineProgram
	{
	private int min_clip_length=10;
	private boolean ignore_poly_x=false;
	private int max_size_inversion=2000;
	private float pct_identity=1.0f;
	private GenomicSequence genomicSequence=null;
	
	private class ShortRead
		{
		SAMRecord rec;
	
		
		ShortRead(SAMRecord rec)
			{
			this.rec=rec;
			}
		
		public boolean isClipped()
			{
			return isClipped(0) || isClipped(1);
			}
		private CigarElement getCigarElement(int side)
			{
			Cigar cigar=rec.getCigar();
			CigarElement ce=cigar.getCigarElement(side==0?0:cigar.numCigarElements()-1);
			return ce;
			}
		public boolean isClipped(int side)
			{
			CigarElement ce=getCigarElement(side);
			return ce.getOperator().equals(CigarOperator.S) && ce.getLength() >= min_clip_length;
			}
		
		public String getClippedSequence(int side)
			{
			CigarElement ce=getCigarElement(side);
			if(side==0)
				{
				return rec.getReadString().substring(0, ce.getLength());
				}
			else
				{
				return rec.getReadString().substring(rec.getReadLength()-ce.getLength());
				}
			}
		public String getGenomicSequence(int side)
			{
			CigarElement ce=getCigarElement(side);
			int start,end;
			if(side==0)
				{
				start=rec.getUnclippedStart()-1;
				end=start+ce.getLength();
				}
			else
				{
				
				end=rec.getUnclippedEnd();
				start=end-ce.getLength();
				}
			return genomicSequence.subSequence(start, end).toString();
			}
		
		}
	
	@Override
	public String getProgramDescription() {
		return "";
		}
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -R (ref) "+getMessageBundle("reference.faidx"));
		out.println(" -m (int) max size of inversion . Default: "+max_size_inversion);
		out.println(" -c (int) min size of clipped read. default:"+min_clip_length);
		out.println(" -x ignore poly-X");
		out.println(" -f (0<f<1.0) fraction identity clipped seq/reference.");
		super.printOptions(out);
		}
	
	private static boolean same(char c1,char c2)
		{
		if(c1=='N' || c2=='N') return false;
		return Character.toUpperCase(c1)==Character.toUpperCase(c2);
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		//String samTag="XI";
		//SamWriterFactory swf=SamWriterFactory.newInstance();
		File faidx=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:m:c:xf:"))!=-1)
			{
			switch(c)
				{
				case 'f':pct_identity=Float.parseFloat(opt.getOptArg());break;
				case 'x': ignore_poly_x=true;break;
				case 'R':faidx=new File(opt.getOptArg());break;
				case 'm':max_size_inversion= Integer.parseInt(opt.getOptArg());break;
				case 'c':min_clip_length=Integer.parseInt(opt.getOptArg());break;
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
		if(faidx==null)
			{
			error(getMessageBundle("reference.undefined"));
			return -1;
			}
		IndexedFastaSequenceFile indexedFastaSequenceFile=null;
		SAMFileReader r=null;
		PrintStream out=System.out;
		//SAMFileWriter w=null;
		try
			{
			indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				r=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				r=new SAMFileReader(new File(filename));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			r.setValidationStringency(ValidationStringency.SILENT);
			SAMFileHeader header=r.getFileHeader();
			SAMProgramRecord spr=header.createProgramRecord();
			spr.setProgramName(getProgramName());
			spr.setProgramVersion(getVersion());
			spr.setCommandLine(getProgramCommandLine());
			
			//w=swf.make(header, System.out);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			SAMRecordIterator it= r.iterator();
			LinkedList<ShortRead> buffer=new LinkedList<ShortRead>();
			while(it.hasNext())
				{
				SAMRecord rec=it.next();
				progress.watch(rec);
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				
				Cigar cigar=rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				ShortRead shortRead=new ShortRead(rec);
				
				//this read must be soft clipped in 5' or 3'
				if(!shortRead.isClipped()) continue;
				
				if(!buffer.isEmpty() && !(buffer.getFirst().rec.getReferenceIndex().equals(rec.getReferenceIndex())))
					{
					buffer.clear();
					}
				while(!buffer.isEmpty() && buffer.getFirst().rec.getAlignmentEnd()+this.max_size_inversion < rec.getAlignmentStart())
					{
					buffer.removeFirst();
					}
				
				boolean changed=false;
				for(int side=0;side<2;++side)
					{
					if(!shortRead.isClipped(side))continue;
					String clippedSeq=shortRead.getClippedSequence(side);
					
					
					if(ignore_poly_x)
						{
						int y=1;
						while(y< clippedSeq.length() &&
								same(clippedSeq.charAt(0),clippedSeq.charAt(y)))
							{
							++y;
							}
						if(y==clippedSeq.length() ) continue;
						
						Counter<String> dinucl=new Counter<String>();
						for(y=0;y+1< clippedSeq.length();++y)
							{
							dinucl.incr(clippedSeq.substring(y, y+2));
							}
						if(dinucl.getCountCategories()<=2) continue;
						/*
						if(dinucl.count(dinucl.getMostFrequent())*2.0 >= dinucl.getTotal())
							{
							continue;
							}*/
						}					
					clippedSeq=AcidNucleics.reverseComplement(clippedSeq);
					

					if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(rec.getReferenceName()))
						{
						this.genomicSequence=new GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
						}
					int genomeStart;
					int genomeEnd;
					/*
					if(side==0)
						{
						genomeStart=rec.getAlignmentEnd();
						genomeEnd=Math.min(
								genomicSequence.length()-clippedSeq.length(),
								rec.getAlignmentEnd()+this.max_size_inversion
								);
						}
					else
						{
						genomeStart=Math.max(1,rec.getAlignmentStart()-this.max_size_inversion);
						genomeEnd=rec.getAlignmentStart()-clippedSeq.length();

						}*/
					
					genomeStart=Math.max(1,rec.getAlignmentStart()-this.max_size_inversion);
					genomeEnd=Math.min(
							genomicSequence.length()-clippedSeq.length(),
							rec.getAlignmentEnd()+this.max_size_inversion
							);
					
					
					for(int x=genomeStart;
							x< genomeEnd && !changed;
							++x)
						{
						if(!(x+clippedSeq.length()<rec.getAlignmentStart() || x>rec.getAlignmentEnd()))
							{
							continue;
							}
						
						int y=0;
						if(pct_identity>=1)
							{
							while(y< clippedSeq.length() &&
								same(clippedSeq.charAt(y),genomicSequence.charAt(x+y))
								)	
								{
								++y;
								}
							if(y!=clippedSeq.length()) continue;
							}
						else
							{
							float nmatch=0;
							for(y=0; y< clippedSeq.length();++y)
								{
								if(same(clippedSeq.charAt(y),genomicSequence.charAt(x+y)))
									{
									nmatch++;
									}
								}
							if(nmatch/(float)clippedSeq.length() < this.pct_identity)
								{
								continue;
								}
							}
						changed=true;
						
						boolean found_clip_in_buffer=false;
						String genomicClip=shortRead.getGenomicSequence(side);
						for(ShortRead prev:buffer)
							{
							if(!prev.isClipped(side==0?1:0)) continue;
							String complRead=prev.getClippedSequence(side==0?1:0);
							if(complRead.contains(genomicClip) || complRead.contains(AcidNucleics.reverseComplement(genomicClip)))
								{
								found_clip_in_buffer=true;
								}
							}
						
						//String msg= ""+(side==0?rec.getAlignmentStart():rec.getAlignmentEnd())+
						//		","+x+","+clippedSeq+","+(side==0?genomeEnd-x:x-genomeStart)+","+shortRead.getGenomicSequence(side);
						//rec.setAttribute(samTag,msg);
						//info("Found "+rec+" ["+rec.getReferenceName()+":"+(side==0?rec.getAlignmentStart():rec.getAlignmentEnd())+"]" +"  "+msg+" "+buffer.size());
						
						
						
						out.print(rec.getReferenceName());
						out.print("\t");
						out.print(rec.getAlignmentStart());
						out.print("\t");
						out.print(rec.getAlignmentEnd());
						out.print("\t");
						out.print(rec.getReadName());
						out.print("\t");
						out.print(rec.getFlags());
						out.print("\t");
						out.print(rec.getCigarString());
						out.print("\t");
						out.print(rec.getReadString());
						out.print("\t");
						out.print(side==0?"5_prime":"3_prime");
						out.print("\t");
						out.print(clippedSeq);
						out.print("\t");
						out.print(clippedSeq.length());
						out.print("\t");
						out.print(10*((int)(x/10.0)));
						out.print("\t");
						out.print(x<rec.getAlignmentStart()?
								rec.getAlignmentStart()-x:
								x-rec.getAlignmentEnd());
						out.print("\t");
						out.print(found_clip_in_buffer);
						out.println();
						}
					}
				buffer.add(shortRead);
				
				if(!changed) continue;
				
				
				
				
				
				
				//w.addAlignment(rec);
				}
			
			it.close();
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
			CloserUtil.close(indexedFastaSequenceFile);
			CloserUtil.close(r);
			//CloserUtil.close(w);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamShortInvertion().instanceMainWithExit(args);
		}
	}
