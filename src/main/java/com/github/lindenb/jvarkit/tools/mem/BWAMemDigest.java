package com.github.lindenb.jvarkit.tools.mem;


import java.io.File;

import com.github.lindenb.jvarkit.util.picard.XPAlign;
import com.github.lindenb.jvarkit.util.picard.XPalignFactory;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.util.Log;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BWAMemDigest extends CommandLineProgram
	{
    private static final Log LOG = Log.getInstance(BWAMemDigest.class);

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="SAM or BAM input file or stdin.", optional=true)
    public File IN = null;

    private SAMFileHeader header;
    
    
    
	public BWAMemDigest()
		{
		
		}

	
	private static String bed(SAMRecord record)
		{
		return record.getReferenceName()+"\t"+
			    (record.getAlignmentStart()-1)+"\t"+
			    record.getAlignmentEnd()+"\t"+
			    (record.getReadNegativeStrandFlag()?"-":"+")
			    ;
		}
	
	private static String mate(SAMRecord record)
		{
		return record.getMateReferenceName()+"\t"+
			    (record.getMateAlignmentStart()-1)+"\t"+
			    (record.getMateNegativeStrandFlag()?"-":"+")
			    ;
		}
	
	@Override
	protected int doWork()
		{
		final float limitcigar=0.15f;
		SAMFileReader r=null;
		try {
			
			
			if(IN==null)
				{
				r=new SAMFileReader(System.in);
				}
			else
				{
				r=new SAMFileReader(IN);
				}
		r.setValidationStringency(super.VALIDATION_STRINGENCY);
		this.header=r.getFileHeader();
		XPalignFactory xPalignFactory=new XPalignFactory(this.header);
		SAMRecordIterator iter=r.iterator();	
		long readNum=0L;
		while(iter.hasNext())
			{
			SAMRecord record=iter.next();
			++readNum;
			if(!record.getReadPairedFlag()) continue;
			if(record.getProperPairFlag()) continue;
			if(record.getReadFailsVendorQualityCheckFlag()) continue;
			if(record.getDuplicateReadFlag()) continue;
			if(record.getReadUnmappedFlag()) continue;
			
			
			float countM=0f;
			float countS=0f;
			for(CigarElement c:record.getCigar().getCigarElements())
				{
				switch(c.getOperator())
					{
					case M: countM+=c.getLength(); break;
					case S: countS+=c.getLength(); break;
					default: break;
					}
				}
			
			if(countM>10 && ((countS/countM)>(1f-limitcigar) && (countS/countM)< (1f+limitcigar) ))
				{
				System.out.println(
						bed(record)+
						"\t"+readNum+
						"\tINSERT\t"+
						record.getCigarString()+"\t"+(countS/countM)
						);
				}
			
			for(XPAlign xp: xPalignFactory.getXPAligns(record))
				{
				System.out.println(
						bed(record)+
						"\t"+readNum+
						"\tXP"+
					    "\t"+record.getCigarString()+
					    "\t"+xp.getChrom()+
					    "\t"+xp.getPos()+
					    "\t"+xp.getStrand()+
					    "\t"+xp.getCigarString()
						);				
				}
			
			
			if(record.getMateUnmappedFlag())
				{
				System.out.println(
					bed(record)+
					"\t"+readNum+
				    "\tORPHAN"
					);
				continue;
				}
			
			
			System.out.println(
						bed(record)+
						"\t"+readNum+
					    "\t"+(record.getReferenceIndex()==record.getMateReferenceIndex()?"DEL":"TRANSLOC")+
					    "\t"+mate(record)
						);
			
			
			}
		}
		
		catch (Exception e)
			{
			e.printStackTrace();
			LOG.error(e);
			return -1;
			}
		finally
			{
			if(r!=null) r.close();
			}
		return 0;
		}
	
	public static void main(String[] args) {
		new BWAMemDigest().instanceMainWithExit(args);
	}
}
