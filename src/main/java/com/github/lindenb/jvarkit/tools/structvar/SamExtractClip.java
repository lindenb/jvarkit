package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class SamExtractClip extends AbstractCommandLineProgram
	{
	private int min_clip_length=5;
	private boolean print_clipped_read=false;
	private boolean print_original_read=false;
	
	
	@Override
	public String getProgramDescription() {
		return "Extract Clipped Sequences from a SAM. Ouput is a FASTQ";
		}
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -m (int) Min size of clipped read. default:"+min_clip_length);
		out.println(" -c print Clipped read");
		out.println(" -r print original Read");
		out.println(" -o (file.fastq) file out. Default:stdout");
		super.printOptions(out);
		}
	
	
	
	
	@Override
	public int doWork(String[] args)
		{
		File fatqsout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"m:cro:"))!=-1)
			{
			switch(c)
				{
				case 'r': print_original_read=true;break;
				case 'c': print_clipped_read=true;break;
				case 'm':min_clip_length=Integer.parseInt(opt.getOptArg());break;
				case 'o':fatqsout=new File(opt.getOptArg());break;
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
	
		SAMFileReader r=null;
		BasicFastqWriter out=null;
		try
				{
				if(fatqsout!=null)
					{
					info("writing to "+fatqsout);
					out=new BasicFastqWriter(fatqsout);
					}
				else
					{
					info("writing to stdout");
					out=new BasicFastqWriter(System.out);
					}
				if(opt.getOptInd()==args.length)
					{
					info("Reading from stdin");
					r=new SAMFileReader(System.in);
					run(r,out);
					r.close();
					}
				else 
					{
					for(int optind=opt.getOptInd();optind<args.length;++optind)
						{
						String filename=args[optind];
						info("Reading from "+filename);
						r=new SAMFileReader(new File(filename));
						run(r,out);
						r.close();
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
				CloserUtil.close(r);
				CloserUtil.close(out);
				}
		}
		
		private void run(SAMFileReader r,FastqWriter out)
			{
			int startend[]=new int[2];
			r.setValidationStringency(ValidationStringency.SILENT);
			SAMFileHeader header=r.getFileHeader();
			SAMProgramRecord spr=header.createProgramRecord();
			spr.setProgramName(getProgramName());
			spr.setProgramVersion(getVersion());
			spr.setCommandLine(getProgramCommandLine());
			
			//w=swf.make(header, System.out);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			SAMRecordIterator it= r.iterator();
			while(it.hasNext())
				{
				SAMRecord rec=it.next();
				progress.watch(rec);
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				
				Cigar cigar=rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				
				String suffix="";
				if(rec.getReadPairedFlag())
					{
					suffix=(rec.getFirstOfPairFlag()?"/1":"/2");
					}
				
			
				
				
				startend[0]=0;
				startend[1]=rec.getReadLength();
				boolean found=false;
				for(int side=0;side<2;++side)
					{
					CigarElement ce=cigar.getCigarElement(side==0?0:cigar.numCigarElements()-1);
					if(!ce.getOperator().equals(CigarOperator.S)) continue;
					if(ce.getLength() < min_clip_length) continue;
					
					found=true;
					String clippedSeq;
					String clippedQual;
					if(side==0)
						{
						startend[0]=ce.getLength();
						clippedSeq= rec.getReadString().substring(0, startend[0]);
						clippedQual=rec.getBaseQualityString().substring(0, startend[0]);
						}
					else
						{
						startend[1]=rec.getReadLength()-ce.getLength();
						clippedSeq= rec.getReadString().substring(startend[1]);
						clippedQual=rec.getBaseQualityString().substring(startend[1]);
						}
					
					out.write(new FastqRecord(
							rec.getReadName()+suffix+":"+side+":"+rec.getReferenceName()+":"+rec.getAlignmentStart()+":"+rec.getFlags(),
							clippedSeq,
							"",
							clippedQual
							));
					}
				if(!found) continue;
				
				String bases=rec.getReadString();
				String qual=rec.getBaseQualityString();
				if( rec.getReadNegativeStrandFlag())
					{
					bases=AcidNucleics.reverseComplement(bases);
					qual=new StringBuilder(qual).reverse().toString();
					}
				
				if(print_original_read)
					{
					out.write(new FastqRecord(
							rec.getReadName()+suffix,
							bases,
							"",
							qual
							));
					}
				
				if(print_clipped_read)
					{
					out.write(new FastqRecord(
							rec.getReadName()+suffix+":clipped",
							bases.substring(startend[0], startend[1]),
							"",
							qual.substring(startend[0], startend[1])
							));
					}
				}
			
			it.close();
			progress.finish();
			}

		
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamExtractClip().instanceMainWithExit(args);
		}
	}
