package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
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
import com.github.lindenb.jvarkit.util.picard.SamWriterFactory;

public class SamShortInvertion extends AbstractCommandLineProgram
	{
	private int min_clip_length=10;
	private boolean ignore_poly_x=false;
	private int max_size_inversion=2000;
	private float pct_identity=1.0f;
	
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
		String samTag="XI";
		SamWriterFactory swf=SamWriterFactory.newInstance();
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
				case 'm':faidx=new File(opt.getOptArg());break;
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
		SAMFileWriter w=null;
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
			GenomicSequence genomicSequence=null;
			w=swf.make(header, System.out);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			SAMRecordIterator it= r.iterator();
			while(it.hasNext())
				{
				SAMRecord rec=it.next();
				progress.watch(rec);
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				Cigar cigar=rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				boolean changed=false;
				for(int side=0;side<2;++side)
					{
					CigarElement ce=cigar.getCigarElement(side==0?0:cigar.numCigarElements()-1);
					if(ce.getOperator()!=CigarOperator.S) continue;
					if((float)ce.getLength() < min_clip_length) continue;
					String clippedSeq;
					if(side==0)
						{
						clippedSeq=rec.getReadString().substring(0, ce.getLength());
						}
					else
						{
						clippedSeq=rec.getReadString().substring(rec.getReadLength()-ce.getLength());
						}
					
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
					
					if(clippedSeq.contains("N")) continue;
					
					
					
					if(genomicSequence==null || !genomicSequence.getChrom().equals(rec.getReferenceName()))
						{
						genomicSequence=new GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
						}
					int genomeStart;
					int genomeEnd;
					
					if(side==0)
						{
						genomeStart=Math.max(1,rec.getAlignmentStart()-this.max_size_inversion);
						genomeEnd=rec.getAlignmentStart()-clippedSeq.length();
						}
					else
						{
						genomeStart=rec.getAlignmentEnd();
						genomeEnd=Math.min(
								genomicSequence.length()-clippedSeq.length(),
								rec.getAlignmentEnd()+this.max_size_inversion
								);
						}
					
					for(int x=genomeStart;
							x< genomeEnd && !changed;
							++x)
						{
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
						String msg= ""+(side==0?rec.getAlignmentStart():rec.getAlignmentEnd())+
								","+x+","+clippedSeq+","+(side==0?genomeEnd-x:x-genomeStart);
						rec.setAttribute(samTag,msg);
						info("Found "+rec+" ["+rec.getReferenceName()+":"+(side==0?rec.getAlignmentStart():rec.getAlignmentEnd())+"]" +"  "+msg);
						}
					
					}
				if(!changed) continue;
				w.addAlignment(rec);
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
			CloserUtil.close(w);
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
