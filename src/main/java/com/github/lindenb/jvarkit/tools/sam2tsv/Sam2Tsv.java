package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.CigarIterator;
import com.github.lindenb.jvarkit.util.picard.IntervalUtils;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class Sam2Tsv extends CommandLineProgram
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Prints the SAM alignments as a TAB delimited file. ";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM files to process.",minElements=0)
	public List<File> IN=new ArrayList<File>();

    @Option(shortName= "L", doc="restrict to that region (chr:start-end)",optional=true)
	public String REGION=null;
    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Indexex reference",optional=false)
	public File  REF=null;
    @Option(shortName= "A", doc="Use Alignment output.",optional=false)
	public boolean ALN=false;

    private Interval interval=null;
	private PrintWriter pw=new PrintWriter(System.out);
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	
	private void print(final SAMRecord rec)
		{
		if(ALN)
			{
			printAln(rec);
			}
		else
			{
			printTsv(rec);
			}
		}
	
	private void printAln(final SAMRecord rec)
		{
		StringBuilder L1=new StringBuilder();
		StringBuilder L2=new StringBuilder();
		StringBuilder L3=new StringBuilder();
		
		CigarIterator ci=CigarIterator.create(rec, this.indexedFastaSequenceFile);
		
		pw.println(">"+rec.getReadName()+
				"\t"+
				rec.getReferenceName()
				);
		
		while(ci.next())
			{
			int readP=ci.getReadPosition();
			int refP=ci.getReferencePosition();
			Character readBase=null;
			Character refBase=null;
			
			if(readP==-1)
				{
				L3.append('.');
				}
			else
				{
				L3.append((readBase=ci.getReadBase()));
				}
			
			if(refP==-1)
				{
				L1.append('.');
				}
			else
				{
				L1.append((refBase=ci.getReferenceBase()));
				}
			
			if(readBase==null || refBase==null)
				{
				L2.append(' ');
				}
			else if(Character.toUpperCase(readBase)==Character.toUpperCase(refBase))
				{
				L2.append('|');
				}
			else
				{
				L2.append(' ');
				}
			}
		
		pw.printf("%-30s %-8d %s %8d\n",
				rec.getReferenceName(),
				rec.getUnclippedStart(),
				L1.toString(),
				rec.getUnclippedEnd()
				);
		pw.printf("%-30s %-8s %s\n",
				"",
				"",
				L2.toString()
				);

		pw.printf("%-30s %-8d %s %8d\n",
				rec.getReadName(),
				1,
				L3.toString(),
				rec.getReadLength()
				);
		pw.println();
		}
	
	
	private void printTsv(final SAMRecord rec)
		{
		CigarIterator ci=CigarIterator.create(rec, this.indexedFastaSequenceFile);
		
			
		while(ci.next())
			{
			int readP=ci.getReadPosition();
			int refP=ci.getReferencePosition();
			Character readBase=null;
			Character refBase=null;
			
			pw.print(rec.getReadName());
			pw.print('\t');
			pw.print(rec.getFlags());
			pw.print('\t');
			if(readP==-1)
				{
				pw.print(".\t.\t.\t");
				}
			else
				{
				Integer qual=ci.getReadQual();
				pw.print(1+readP);
				pw.print('\t');
				pw.print((readBase=ci.getReadBase()));
				pw.print('\t');
				pw.print(qual==null?".":qual.toString());
				pw.print('\t');
				}
			
			pw.print(rec.getReferenceName());
			pw.print('\t');
			if(refP==-1)
				{
				pw.print(".\t.\t");
				}
			else
				{
				pw.print(refP);
				pw.print('\t');
				pw.print((refBase=ci.getReferenceBase()));
				pw.print('\t');
				}
			pw.print(ci.getCigarOperator());
			pw.print('\t');
			if(readBase==null || refBase==null)
				{
				pw.print('.');
				}
			else if(Character.toUpperCase(readBase)==Character.toUpperCase(refBase))
				{
				pw.print('=');
				}
			else
				{
				pw.print('X');
				}
			
			pw.println();
			}
		}
	
	private void scan(SAMFileReader r) 
		{
		r.setValidationStringency(super.VALIDATION_STRINGENCY);
		SAMRecordIterator iter=null;
		if(interval==null)
			{
			iter=r.iterator();
			}
		else
			{
			iter=r.queryOverlapping(interval.getSequence(),interval.getStart(), interval.getEnd());
			}
		
		while(iter.hasNext())
			{
			SAMRecord rec=iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			print(rec);
			}
		}
	
	@Override
	protected int doWork()
		{
		try {
			
			
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);
			
			if(REGION!=null)
				{
				this.interval=IntervalUtils.parseOne(
						this.indexedFastaSequenceFile.getSequenceDictionary(),
						REGION);
				if(this.interval==null)
					{
					System.err.println("Cannot parse interval "+REGION+" or chrom doesn't exists.");
					return -1;
					}
				}
			
			if(IN.isEmpty())
				{
				SAMFileReader r=new SAMFileReader(System.in);
				scan(r);
				r.close();
				}
			else
				{
				for(File f:this.IN)
					{
					SAMFileReader r=new SAMFileReader(f);
					scan(r);
					r.close();
					}
				}
			pw.flush();
			}	
		catch (FileNotFoundException e)
			{
			e.printStackTrace();
			return -1;
			}
		
		return 0;
		}

	public static void main(String[] args)
		{
		new Sam2Tsv().instanceMainWithExit(args);
		}
}
