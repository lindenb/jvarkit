package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.File;
import java.io.PrintStream;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import net.sf.picard.PicardException;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloserUtil;

/**
 * https://github.com/lindenb/jvarkit/wiki/SAM2Tsv
 */
public class Sam2Tsv
	extends AbstractCommandLineProgram
	{
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	private boolean printAlignment=false;
	/** lines for alignments */
	private StringBuilder L1=null;
	private StringBuilder L2=null;
	private StringBuilder L3=null;

	
	
	private void printAln(
			final SAMRecord rec,
			Integer readPos,
			Character readChar,
			Integer refPos,
			Character refChar,
			CigarOperator op
			)
			{
			System.out.print(rec.getReadName());
			System.out.print("\t");
			System.out.print(rec.getReadUnmappedFlag()?".":rec.getReferenceName());
			System.out.print("\t");
			if(readPos!=null && readChar!=null)
				{
				System.out.print(readPos);
				System.out.print("\t");
				System.out.print(readChar);
				System.out.print("\t");
				}
			else
				{
				System.out.print(".\t.\t");
				}
			if(refPos!=null && refChar!=null)
				{
				System.out.print(refPos);
				System.out.print("\t");
				System.out.print(refChar);
				System.out.print("\t");
				}
			else
				{
				System.out.print(".\t.\t");
				}
			System.out.print(op==null?".":op.name());
			System.out.println();
			if(this.printAlignment)
				{
				L1.append(readChar==null?'-':readChar);
				L3.append(refChar==null?'-':refChar);
				if(refChar!=null && readChar!=null && Character.toUpperCase(refChar)==Character.toUpperCase(readChar))
					{
					L2.append('|');
					}
				else
					{
					L2.append(' ');
					}
				}
			}
	
	private void printAln(final SAMRecord rec)
		{
		if(rec==null) return;
		byte readbases[]=rec.getReadBases();
		if(readbases==null || rec.getReadUnmappedFlag())
			{
			printAln(rec,null,null,null,null,null);
			return;
			}
		
		if(genomicSequence==null || !genomicSequence.getChrom().equals(rec.getReferenceName()))
			{
			genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, rec.getReferenceName());
			}
		

		
		
		
		 int readIndex = 0;
		 int refIndex = rec.getAlignmentStart();
		 				 
		 for (final CigarElement e : rec.getCigar().getCigarElements())
			 {
			 switch (e.getOperator())
				 {
				 case H : break; // ignore hard clips
				 case P : break; // ignore pads
				 case I : //cont.
				 case S :
				 		{
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			if(readIndex>=0 && readIndex< readbases.length)
				 				{
				 				printAln(rec,
				 						readIndex,
				 						(char)(readbases[readIndex]),
				 						null,
				 						null,
				 						e.getOperator()
				 						);
				 				}
				 			readIndex++;
				 			}
				 		break;
				 		}
				 case N :  //cont. -- reference skip
				 case D :
				 		{
				 		for(int i=0;i<e.getLength();++i)
				 			{
				 			if(refIndex>=1 && refIndex<=this.genomicSequence.length())
				 				{
				 				printAln(rec,
				 						null,
				 						null,
				 						refIndex,
				 						genomicSequence.charAt(refIndex-1),
				 						e.getOperator()
				 						);
				 				}
				 			refIndex++;
				 			}
				 		break;
				 		}
				 case M :
				 case EQ :
				 case X :
			 			{
				 		final int length = e.getLength();
				 		for(int i=0;i<length;++i)
				 			{
				 			char baseRead='*';
				 			char baseRef='*';
				 			
				 			if(readIndex>=0 && readIndex< readbases.length)
				 				{
					 			baseRead=(char)(rec.getReadBases()[readIndex]);
				 				}
				 			if(refIndex>=1 && refIndex<= genomicSequence.length())
					 			{
				 				baseRef=genomicSequence.charAt(refIndex-1);
					 			}
				 			printAln(rec,
				 					readIndex,
				 					baseRead,
			 						refIndex,
			 						baseRef,
			 						e.getOperator()
			 						);
				 			
				 			refIndex++;
				 			readIndex++;
				 			}
				 		break;
			 			}
					
				 default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
				 }

			 }
	
		
		
		 if(printAlignment)
				{
				
				int len=Math.max(rec.getReadNameLength(), rec.getReferenceName().length())+2;
				
				System.out.printf(":%"+len+"s %8d %s %-8d\n",
						rec.getReferenceName(),
						rec.getUnclippedStart(),
						L1.toString(),
						rec.getUnclippedEnd()
						);
				System.out.printf(":%"+len+"s %8s %s\n",
						"",
						"",
						L2.toString()
						);

				System.out.printf(":%"+len+"s %8d %s %-8d\n",
						rec.getReadName(),
						1,
						L3.toString(),
						rec.getReadLength()
						);

				L1.setLength(0);
				L2.setLength(0);
				L3.setLength(0);
				}

		}
	
	
	
	private void scan(SAMFileReader r) 
		{
		r.setValidationStringency(ValidationStringency.LENIENT);
		SAMRecordIterator iter=null;
		try{
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(r.getFileHeader().getSequenceDictionary());
			iter=r.iterator();	
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				progress.watch(rec);
				printAln(rec);
				if(System.out.checkError()) break;
				}
			}
		catch(Exception err)
			{
			error(err);
			throw new PicardException(String.valueOf(err.getMessage()),err);
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/SAM2Tsv";
		}
	
	@Override
	public String getProgramDescription() {
		return "Prints the SAM alignments as a TAB delimited file.";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -r (reference) Fasta Reference indexed file . REQUIRED.");
		out.println(" -A display alignment.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File refFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt getopt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=getopt.getopt(args,getGetOptDefault()+ "r:R:A"))!=-1)
			{
			switch(c)
				{
				case 'A': this.printAlignment=true;break;
				case 'R': case 'r': refFile=new File(getopt.getOptArg());break;
				default: 
					{
					switch(handleOtherOptions(c, getopt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		
		if(refFile==null)
			{
			error("Undefined REF file");
			return -1;
			}
		
		if(printAlignment)
			{
			L1=new StringBuilder();
			L2=new StringBuilder();
			L3=new StringBuilder();
			}
		
		SAMFileReader samFileReader=null;
		try
			{
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(refFile);
			samFileReader=null;
			if(getopt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				samFileReader=new SAMFileReader(System.in);
				scan(samFileReader);
				samFileReader.close();
				}
			else 
				{
				for(int optind=getopt.getOptInd();optind< args.length;++optind)
					{
					File bamFile=new File(args[optind]);
					info("Reading "+bamFile);
					samFileReader=new SAMFileReader(bamFile);
					scan(samFileReader);
					samFileReader.close();
					}
				}
			return 0;
			}
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(indexedFastaSequenceFile);
			CloserUtil.close(samFileReader);
			}
		}
	
	public static void main(String[] args)
		{
		new Sam2Tsv().instanceMainWithExit(args);
		}
}
