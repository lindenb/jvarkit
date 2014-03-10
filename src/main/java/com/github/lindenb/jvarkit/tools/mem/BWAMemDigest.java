package com.github.lindenb.jvarkit.tools.mem;


import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlign;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlignFactory;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.util.Log;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloserUtil;

/**
 * : ant  bwamemdigest && ~/package/bwa/bwa-0.7.4/bwa mem \
 * 		 -M  ~/tmp/DATASANGER/hg18/chr1.fa  \
 * 			  ~/tmp/DATASANGER/4368_1_1.fastq.gz  ~/tmp/DATASANGER/4368_1_2.fastq.gz 2> /dev/null | 
 * 			java -jar dist/bwamemdigest.jar BED=<(curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz" | gunzip -c | cut -f2,3,4 ) XBED=500  | tee /dev/tty | gzip --best > /tmp/jeter.mem.bed.gz

 */
public class BWAMemDigest extends AbstractCommandLineProgram
	{
    private static final Log LOG = Log.getInstance(BWAMemDigest.class);

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="SAM or BAM input file or stdin.", optional=true)
    public File IN = null;

    @Option(shortName="BED", doc="Regions to ignore.", optional=true)
    public File IGNORE_BED = null;
 
    @Option(shortName="XBED", doc=" Extend BED Regions to ignore by 'x' bases. ", optional=true)
    public int IGNORE_EXTEND = 0;
    
    
    private abstract class AbstractMemOuput
    	implements Closeable
    	{
    	
    	}
    
    
    private class DefaultMemOuput extends AbstractMemOuput
		{
    	PrintWriter out=new PrintWriter(System.out);
    	
    	public void insertion(SAMRecord record,long readNum,float countS,float countM)
	    	{
	    	System.out.println(
					bed(record)+
					"\t"+readNum+
					"\tINSERT\t"+
					record.getCigarString()+"\t"+(countS/countM)
					);
	    	}
    	
    	public void xp(SAMRecord record,long readNum,OtherCanonicalAlign xp)
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
    	public void orphan(SAMRecord record,long readNum)
	    	{
	    	System.out.println(
					bed(record)+
					"\t"+readNum+
				    "\tORPHAN"
					);
	    	}
    	public void deletion(SAMRecord record,long readNum)
	    	{
	    	System.out.println(
					bed(record)+
					"\t"+readNum+
				    "\tDEL"+
				    "\t"+mate(record)
					);
	    	}
    	public void transloc(SAMRecord record,long readNum)
	    	{
	    	System.out.println(
					bed(record)+
					"\t"+readNum+
				    "\tTRANSLOC"+
				    "\t"+mate(record)
					);
	    	}
    	@Override
    	public void close() throws IOException {
    		out.flush();
    		out.close();
    		}
    	}
    
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
		DefaultMemOuput output=new DefaultMemOuput();
		
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
				if(ignore.put(tokens[0],
						Math.max(1,Integer.parseInt(tokens[1])-(1+IGNORE_EXTEND)),
						Integer.parseInt(tokens[2])+IGNORE_EXTEND,
						Boolean.TRUE
						))
					{
					LOG.warn("BED:ignoring "+line);
					}
				}	
			in.close();
			}
		
		OtherCanonicalAlignFactory xPalignFactory=new OtherCanonicalAlignFactory(header);
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
				LOG.info("ignore "+record);
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
				output.insertion(record, readNum, countS, countM);
				}
			
			for(OtherCanonicalAlign xp: xPalignFactory.getXPAligns(record))
				{
				if(ignore!=null &&
						ignore.containsOverlapping(
						xp.getChrom(),
						xp.getPos(),
						xp.getPos()
						))
					{
					LOG.info("ignore "+record);
					continue;
					}
				
				
				output.xp(record, readNum, xp);
				}
			
			
			if(record.getMateUnmappedFlag())
				{
				output.orphan(record, readNum);
				continue;
				}
			
			
			if(ignore!=null && ignore.containsOverlapping(
					record.getMateReferenceIndex(),
					record.getMateAlignmentStart(),
					record.getMateAlignmentStart())
					)
				{
				LOG.info("ignore "+record);
				continue;
				}
			if(record.getReferenceIndex()==record.getMateReferenceIndex())
				{
				output.deletion(record, readNum);
				}
			else
				{
				output.transloc(record, readNum);
				}
			
			
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
			CloserUtil.close(output);
			}
		return 0;
		}
	
	public static void main(String[] args) {
		new BWAMemDigest().instanceMainWithExit(args);
	}
}
