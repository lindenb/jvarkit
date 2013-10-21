package com.github.lindenb.jvarkit.tools.mem;


import java.io.BufferedReader;
import java.io.File;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;
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

    @Option(shortName="BED", doc="Regions to ignore.", optional=true)
    public File IGNORE_BED = null;
 
    @Option(shortName="XBED", doc=" Extend BED Regions to ignore by 'x' bases. ", optional=true)
    public int IGNORE_EXTEND = 0;

    
    
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
		SamSequenceRecordTreeMap<Boolean> ignore=null;
		
		
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
		SAMFileHeader header=r.getFileHeader();
		
		if(IGNORE_BED!=null)
			{
			LOG.info("open "+IGNORE_BED);
			ignore=new SamSequenceRecordTreeMap<Boolean>(header.getSequenceDictionary());
			BufferedReader in=IOUtils.openFileForBufferedReading(IGNORE_BED);
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				String tokens[]=line.split("[\t]");
				if(tokens.length<3) continue;
				ignore.put(tokens[0],
						Math.max(1,Integer.parseInt(tokens[1])-(1+IGNORE_EXTEND)),
						Integer.parseInt(tokens[2])+IGNORE_EXTEND,
						Boolean.TRUE
						);
				}	
			in.close();
			}
		
		XPalignFactory xPalignFactory=new XPalignFactory(header);
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
			if(ignore!=null && ignore.containsOverlapping(
					record.getReferenceIndex(),
					record.getAlignmentStart(),
					record.getAlignmentEnd())
					)
				{
				continue;
				}
			
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
				if(ignore!=null &&
						ignore.containsOverlapping(
						xp.getChrom(),
						xp.getPos(),
						xp.getPos()
						)) continue;
				
				
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
			
			
			if(ignore!=null && ignore.containsOverlapping(
					record.getMateReferenceIndex(),
					record.getMateAlignmentStart(),
					record.getMateAlignmentStart())
					)
				{
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
