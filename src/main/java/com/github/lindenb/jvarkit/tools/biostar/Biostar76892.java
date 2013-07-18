package com.github.lindenb.jvarkit.tools.biostar;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


public class Biostar76892 extends CommandLineProgram
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" fix strand of two paired reads close but on the same strand. See http://www.biostars.org/p/76892/ .";

	@Option(shortName= "osf", doc="only save pairs of reads fixed.",
    		optional=true)
	public boolean onlySaveFixed=false;
	
	@Option(shortName= "d", doc="distance beween two reads.",
    		optional=true)
	public int distance=500;
	

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM file to process.",
    		optional=false)
	public File IN=null;

    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="BAM file fixed.",
    		optional=false)
	public File OUT=null;

    
	private Log LOG=Log.getInstance(Biostar76892.class);
	
	@Override
	public String getProgramVersion() {
		return "1.0";
		}
	
	@Override
	protected int doWork()
		{
		SAMFileReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
		
			LOG.info("opening "+IN);
			sfr=new SAMFileReader(IN);
			sfr.setValidationStringency(super.VALIDATION_STRINGENCY);
			if(sfr.getFileHeader().getSortOrder()!=SortOrder.coordinate)
				{
				LOG.error("input must be sorted on coordinate ? got: "+sfr.getFileHeader().getSortOrder() );
				return -1;
				}
			
			
			SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
			sfwf.setCreateIndex(super.CREATE_INDEX);
			sfwf.setCreateMd5File(super.CREATE_MD5_FILE);
			sfwf.setMaxRecordsInRam(super.MAX_RECORDS_IN_RAM);
			if(super.TMP_DIR.isEmpty()) sfwf.setTempDirectory(super.TMP_DIR.get(0));
	
			
			LOG.info("opening "+OUT);
			SAMProgramRecord sp=new SAMProgramRecord(getClass().getSimpleName());
			sp.setProgramName(getClass().getSimpleName());
			sp.setProgramVersion(String.valueOf(getProgramVersion()));
			sp.setPreviousProgramGroupId(getClass().getSimpleName());
			sp.setCommandLine(getCommandLine().replace('\t', ' '));
			
			sfr.getFileHeader().addProgramRecord(sp);
			sfw=sfwf.makeBAMWriter(sfr.getFileHeader(), false, OUT);
			long nRecords=0;
			List<SAMRecord> buffer=new ArrayList<SAMRecord>();
			SAMRecordIterator iter=sfr.iterator();
			for(;;)
				{
				SAMRecord rec=null;
				//get next record
				if(iter.hasNext())
					{
					rec=iter.next();
					++nRecords;
					if(nRecords%1000000==0) LOG.info("records: "+nRecords);
					if(!rec.getReadPairedFlag() ||
						rec.getReadUnmappedFlag() ||
						rec.getMateUnmappedFlag() ||
						rec.getProperPairFlag() ||
						rec.getReferenceIndex()!=rec.getMateReferenceIndex() ||
						rec.getReadNegativeStrandFlag()==!rec.getMateNegativeStrandFlag()
						)
						{
						if(!onlySaveFixed) sfw.addAlignment(rec);
						continue;
						}
					}

				if(rec!=null)
					{
					int i=0;
					//cleanup buffer
					int mate_index=-1;
					while(i<buffer.size())
						{
						SAMRecord prev=buffer.get(i);
						if(prev.getReferenceIndex()!=rec.getReferenceIndex() ||
							prev.getAlignmentEnd() + distance < rec.getAlignmentStart())
							{
							if(!onlySaveFixed) sfw.addAlignment(prev);
							buffer.remove(i);
							}
						else if(prev.getReadName().equals(rec.getReadName()) && 
								(	(prev.getFirstOfPairFlag() && rec.getSecondOfPairFlag()) ||
									(rec.getFirstOfPairFlag() && prev.getSecondOfPairFlag()))
								)
							{
							mate_index=i;
							++i;
							}
						else
							{
							++i;
							}
						}
					
					if(mate_index==-1)
						{
						buffer.add(rec);
						}
					else
						{
						
						SAMRecord mate=buffer.get(mate_index);
						buffer.remove(mate_index);
						LOG.info("changing "+rec+" "+mate);
						
						if(mate.getReadNegativeStrandFlag())
							{
							mate.setReadNegativeStrandFlag(false);
							rec.setMateNegativeStrandFlag(mate.getReadNegativeStrandFlag());
							}
						else
							{
							rec.setReadNegativeStrandFlag(false);
							mate.setMateNegativeStrandFlag(rec.getReadNegativeStrandFlag());
							}
						if(!mate.getReadFailsVendorQualityCheckFlag() && !rec.getReadFailsVendorQualityCheckFlag())
							{
							mate.setProperPairFlag(true);
							rec.setProperPairFlag(true);
							}
						mate.setAttribute("rv",1);
						rec.setAttribute("rv", 1);
						sfw.addAlignment(mate);
						sfw.addAlignment(rec);
						}
					
					}
				else
					{
					for(SAMRecord r:buffer)
						{
						if(!onlySaveFixed)  sfw.addAlignment(r);
						}
					break;
					}
				}
			LOG.info("done");
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(sfw!=null)sfw.close();
			if(sfr!=null)sfr.close();
			}
		}
	
	public static void main(String[] args)throws Exception
		{
		/*
		args=new String[]{
				"/commun/data/projects/20120828.AC0KTCACXX.WHOLEGENOME1/align/CD05121/CD05121_recal.bam",
				"/commun/data/users/lindenb/jeter.bam"
				};*/
		new Biostar76892().instanceMain(args);
		}
	}
