package com.github.lindenb.jvarkit.tools.mem;

import java.io.File;
import java.util.ArrayList;
import java.util.List;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlign;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlignFactory;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

public class BWAMemNOp extends AbstractCommandLineProgram
	{
	private int min_soft_clip_length=10;
	private boolean print_only_spit_read=false;
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BWAMemNOp";
		}
	
	@Override
	public String getProgramDescription() {
		return "Merge the other BWA-MEM alignements with its SA:Z:* attributes to an alignment containing a cigar string with 'N' (  Skipped region from the reference.)";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-S only print converted reads.");
		out.println("-o (file.bam) output file. Default: stdout.");
		out.println("-m (int) min soft clip length Default: "+min_soft_clip_length);
		super.printOptions(out);
		}
	static class CigarEvt
		{
		CigarOperator op;
		int read0;
		int ref1;
		}
	private static List<CigarEvt> cigarEvents(int read0,int ref1,Cigar cigar)
		{
		List<CigarEvt> L=new ArrayList<CigarEvt>();
		for(CigarElement ce:cigar.getCigarElements())
			{
			for(int i=0;i< ce.getLength();++i)
				{
				CigarEvt evt=new CigarEvt();
				evt.op=ce.getOperator();
				evt.read0=read0;
				evt.ref1=ref1;
				L.add(evt);
				switch(ce.getOperator())
					{
					case P:break;
					case H:break;
					case I:
					case S:
						{
						read0++;
						break;
						}
					case D:
					case N:
						{
						ref1++;
						break;
						}
					case X:case EQ:case M:
						{
						read0++;
						ref1++;
						break;
						}
					default:throw new IllegalStateException();
					}
				}
			}
		return L;
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		File bamOut=null;
		SAMFileWriterFactory swf=new SAMFileWriterFactory();
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:Sm:"))!=-1)
			{
			switch(c)
				{
				case 'm':min_soft_clip_length=Math.min(1, Integer.parseInt(opt.getOptArg()));break;
				case 'o': bamOut=new File(opt.getOptArg());break;
				case 'S': print_only_spit_read=true;break;
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
		
		SamReader r=null;
		SAMFileWriter w=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				r=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				r=SamFileReaderFactory.mewInstance().open(new File(filename));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			SAMFileHeader header=r.getFileHeader();
			OtherCanonicalAlignFactory ocaf=new OtherCanonicalAlignFactory(header);
			SAMProgramRecord prg=header.createProgramRecord();
			prg.setCommandLine(getProgramCommandLine());
			prg.setProgramName(getProgramName());
			prg.setProgramVersion(getVersion());
			if(bamOut==null)
				{
				w=swf.makeSAMWriter(header, true, System.out);
				}
			else
				{
				w=swf.makeSAMOrBAMWriter(header, true, bamOut);
				}
			SAMRecordFactory samRecordFactory=new DefaultSAMRecordFactory();
			SAMRecordIterator iter=r.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				if(rec.getSupplementaryAlignmentFlag())
					{
					continue;
					}
				if(rec.getReadUnmappedFlag())
					{
					if(!print_only_spit_read) w.addAlignment(rec);
					continue;
					}
				Cigar cigar1=rec.getCigar();
				if(cigar1==null ||
					cigar1.isEmpty() ||
					!(cigar1.getCigarElement(cigar1.numCigarElements()-1).getOperator().equals(CigarOperator.S) ||
					  cigar1.getCigarElement(0).getOperator().equals(CigarOperator.S)
					 )
					) //last or first is soft clipping
					{
					if(!print_only_spit_read) w.addAlignment(rec);
					continue;
					}
				rec.getAlignmentStart();
				List<OtherCanonicalAlign> xps=ocaf.getXPAligns(rec);
				if(xps.isEmpty())
					{
					if(!print_only_spit_read) w.addAlignment(rec);
					continue;
					}
				
				boolean found_one=false;
				for(OtherCanonicalAlign  xp:xps)
					{
					if(!rec.getReferenceName().equals(xp.getReferenceName())) continue;
					if(xp.getReadNegativeStrandFlag()!=rec.getReadNegativeStrandFlag() ) continue;
					
					Cigar cigar2=xp.getCigar();
					if(cigar2==null || cigar2.isEmpty() )
						{
						continue;
						}
					
					SAMRecord newrec=null;
					List<CigarEvt> L1=null;
					List<CigarEvt> L2=null;
					if( cigar1.getCigarElement(cigar1.numCigarElements()-1).getOperator().equals(CigarOperator.S) &&
						cigar1.getCigarElement(cigar1.numCigarElements()-1).getLength() >= this.min_soft_clip_length &&
						cigar2.getCigarElement(0).getOperator().equals(CigarOperator.S)  &&
						cigar2.getCigarElement(0).getLength() >= this.min_soft_clip_length &&
					   rec.getAlignmentEnd()< xp.getAlignmentStart()
					   )
						{
						newrec=samRecordFactory.createSAMRecord(header);
						int ref1=rec.getAlignmentStart();
						newrec.setAlignmentStart(ref1);
						L1=cigarEvents(0, ref1, cigar1);
						L2=cigarEvents(0, xp.getAlignmentStart(), cigar2);
						}
					else  if(
							cigar2.getCigarElement(cigar2.numCigarElements()-1).getOperator().equals(CigarOperator.S) &&
							cigar2.getCigarElement(cigar2.numCigarElements()-1).getLength() >= this.min_soft_clip_length &&
							cigar1.getCigarElement(0).getOperator().equals(CigarOperator.S) &&
							cigar1.getCigarElement(0).getLength() >= this.min_soft_clip_length &&
							xp.getAlignmentEnd()< rec.getAlignmentStart()
							)
						{
						newrec=samRecordFactory.createSAMRecord(header);
						int ref1=xp.getAlignmentStart();
						newrec.setAlignmentStart(ref1);
						L1=cigarEvents(0, ref1, cigar2);
						L2=cigarEvents(0, rec.getAlignmentStart(), cigar1);
						}
					
					if(newrec==null) continue;
					
					newrec.setFlags(rec.getFlags());
					newrec.setReadName(rec.getReadName());
					newrec.setReadBases(rec.getReadBases());
					newrec.setMappingQuality(rec.getMappingQuality());
					newrec.setReferenceIndex(rec.getReferenceIndex());
					newrec.setBaseQualities(rec.getBaseQualities());
					if(found_one)
						{
						newrec.setNotPrimaryAlignmentFlag(true);
						}
					
					found_one=true;
					
					for(SAMTagAndValue tav: rec.getAttributes())
						{
						if(tav.tag.equals(ocaf.getAttributeKey())) continue;
						if(tav.tag.equals("NM")) continue;
						newrec.setAttribute(tav.tag, tav.value);
						}
					if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag())
						{
						newrec.setMateAlignmentStart(rec.getMateAlignmentStart());
						newrec.setMateReferenceIndex(rec.getMateReferenceIndex());
						newrec.setInferredInsertSize(rec.getInferredInsertSize());
						}
					while(!L1.isEmpty() &&
						(L1.get(L1.size()-1).op.equals(CigarOperator.S) ||
						 L1.get(L1.size()-1).op.equals(CigarOperator.D) ||
						 L1.get(L1.size()-1).op.equals(CigarOperator.H)))
						{
						L1.remove(L1.size()-1);
						}
					while(!L2.isEmpty() && L2.get(0).read0<=L1.get(L1.size()-1).read0)
						{
						L2.remove(0);
						}
					
					List<CigarElement> cigarElements=new ArrayList<CigarElement>();
					int i=0;
					while(i< L1.size())
						{
						int j=i+1;
						while(j< L1.size() && L1.get(i).op.equals(L1.get(j).op))
							{
							j++;
							}
						cigarElements.add(new CigarElement(j-i, L1.get(i).op));
						i=j;
						}
					
					//add 'N'
					cigarElements.add(
							new CigarElement((L2.get(0).ref1-L1.get(L1.size()-1).ref1)-1,
							CigarOperator.N)
							);
					
					
					i=0;
					while(i< L2.size())
						{
						int j=i+1;
						while(j< L2.size() && L2.get(i).op.equals(L2.get(j).op))
							{
							j++;
							}
						cigarElements.add(new CigarElement(j-i, L2.get(i).op));
						i=j;
						}
					
					//cleanup : case where  'S' is close to 'N'
					i=0;
					while(i+1< cigarElements.size())
						{
						CigarElement ce1=cigarElements.get(i);
						CigarElement ce2=cigarElements.get(i+1);
						
						if( i>0 &&
							ce1.getOperator().equals(CigarOperator.S) && 
							ce2.getOperator().equals(CigarOperator.N) )
							{
							cigarElements.set(i, new CigarElement(
								ce1.getLength(),
								CigarOperator.X));
							}
						else if(i+2 < cigarElements.size() &&
								ce1.getOperator().equals(CigarOperator.N) && 
								ce2.getOperator().equals(CigarOperator.S) )
							{
							cigarElements.set(i+1, new CigarElement(
								ce2.getLength(),
								CigarOperator.X));
							}
						i++;
							
						}

					
					
					newrec.setCigar(new Cigar(cigarElements));
					List<SAMValidationError> validations=newrec.isValid();
					if(validations!=null)
						{
						for(SAMValidationError err:validations)
							{
							warning(err.getType()+":"+err.getMessage());
							}
						}
					w.addAlignment(newrec);	
					}
				if(!found_one)
					{
					if(!print_only_spit_read)  w.addAlignment(rec);	
					}
				}
			iter.close();
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
			CloserUtil.close(w);
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new BWAMemNOp().instanceMainWithExit(args);
	}

}
