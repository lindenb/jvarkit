package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

public class SamChangeReference extends AbstractCommandLineProgram
	{
	private SAMSequenceDictionary dict=null;
	
	private SamChangeReference()
		{
		
		}
	
	private static class RefName
		{
		String ref;
		int start1;
		}
	
	private RefName convertName(String name)
		{
		RefName r=new RefName();
		r.ref=name;
		r.start1=1;
		int colon=name.indexOf(name);
		if(colon!=-1)
			{
			r.ref=name.substring(0,colon);
			int dash=name.indexOf('-',colon+1);
			if(dash==-1)
				{
				r.start1=Integer.parseInt(name.substring(colon+1,dash));
				}
			else
				{
				r.start1=Integer.parseInt(name.substring(colon+1));
				}
			}
		if(dict.getSequence(r.ref)==null)
			{
			throw new PicardException("The reference sequence '"+r.ref+"' is not declared in the dictionary");
			}
		return r;
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/SamChangeReference";
		}
	
	@Override
	public String getProgramDescription() {
		return " BAM file.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (filenameout.bam) default: SAM as stdout");
		out.println(" -R (file) "+ getMessageBundle("reference.faidx")+".Required");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File faidx=null;
		File bamout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:R:"))!=-1)
			{
			switch(c)
				{
				case 'R': faidx=new File(opt.getOptArg());break;
				case 'o': bamout=new File(opt.getOptArg());break;
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
		
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			info("loading dict");
			this.dict=new SAMSequenceDictionaryFactory().load(faidx);
			
			
			if(opt.getOptInd()==args.length)
				{
				sfr=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File fin=new File(args[opt.getOptInd()]);
				sfr=SamFileReaderFactory.mewInstance().open(fin);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			
			SAMFileHeader header1=sfr.getFileHeader();
			SAMFileHeader header2=header1.clone();
			header2.setSequenceDictionary(this.dict);			
			SAMProgramRecord prog=header2.createProgramRecord();
			prog.setCommandLine(this.getProgramCommandLine());
			prog.setProgramName(getProgramName());
			prog.setProgramVersion(getVersion());
			
			
			
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1.getSequenceDictionary());

			SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
			boolean presorted=false;
			if(bamout!=null)
				{
				info("saving to "+bamout);
				sfw=sfwf.makeSAMOrBAMWriter(header2, presorted, bamout);
				}
			else
				{
				sfw=sfwf.makeSAMWriter(header2, presorted, System.out);
				}
			
			
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec1=iter.next();
				progress.watch(rec1);
				RefName newName1=null;
				RefName newName2=null;
				
				
				if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
					{
					newName1=convertName(rec1.getReferenceName());
					}
				if(rec1.getReadPairedFlag() &&
					!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getMateReferenceName()))
					{
					newName2=convertName(rec1.getMateReferenceName());
					}
				rec1.setHeader(header2);
				
				if(newName1!=null)
					{
					if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
						{
						rec1.setReferenceName(newName1.ref);
						}
					int pos1=rec1.getAlignmentStart();
					if(pos1!=0)
						{
						rec1.setAlignmentStart(pos1 + (newName1.start1-1) );
						}
					}
				if(rec1.getReadPairedFlag())
					{
					if(newName2!=null)
						{
						if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getMateReferenceName()))
							{
							rec1.setMateReferenceName(newName2.ref);
							}
						int pos1=rec1.getMateAlignmentStart();
						if(pos1!=0)
							{
							rec1.setMateAlignmentStart(pos1 + (newName2.start1-1) );
							}
						}
					}
				sfw.addAlignment(rec1);
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
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamChangeReference().instanceMainWithExit(args);
		}
	}
