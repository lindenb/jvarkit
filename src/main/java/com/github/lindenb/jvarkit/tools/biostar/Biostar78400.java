package com.github.lindenb.jvarkit.tools.biostar;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.transform.stream.StreamSource;


import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMTag;


public class Biostar78400 extends CommandLineProgram
	{
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class ReadGroup
		{
		@XmlAttribute(name="ID",required=true)
		public String id;
		@XmlElement(nillable=false)
		public String library;
		@XmlElement(nillable=false)
		public String platform;
		@XmlElement(nillable=false)
		public String sample;
		@XmlElement(nillable=false)
		public String platformunit;
		
		public String center;
		public String description;
		}	
	
	@XmlRootElement(name="read-groups")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class ReadGroupList
		{
		
		@XmlElement(name="flowcell")
		public List<FlowCell> flowcells=new ArrayList<FlowCell>();
		}
	
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class FlowCell
		{
		@XmlAttribute(name="name",required=true)
		public String name;
		@XmlElement(name="lane")
		public List<Lane> lanes=new ArrayList<Lane>();
		}

	@XmlAccessorType(XmlAccessType.FIELD)
	public static class Lane
		{
		@XmlAttribute(name="index")
		public int id;
		@XmlElement(name="group")
		public List<ReadGroup> readGroups=new ArrayList<ReadGroup>();
		}
	
	
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" add the read group info to the sam file on a per lane basis? .";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM file to process (or stdin).",
    		optional=true)
	public File IN=null;

    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="BAM file (or stdout).",
    		optional=true)
	public File OUT=null;
    
    @Option(shortName= "X", doc="XML desfription of the groups.",
    		optional=false)
	public File XML=null;

    
	private Log LOG=Log.getInstance(Biostar78400.class);
	
	@Override
	public String getProgramVersion() {
		return "1.0";
		}
	
	
	private Map<String, Map<Integer,String>> flowcell2lane2id=new HashMap<String, Map<Integer,String>>();
	
	@Override
	protected int doWork()
		{
		
		
		
		
		SAMFileReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			JAXBContext context = JAXBContext.newInstance(ReadGroup.class,ReadGroupList.class);
			Unmarshaller unmarshaller=context.createUnmarshaller();
			ReadGroupList rgl=unmarshaller.unmarshal(new StreamSource(XML),ReadGroupList.class).getValue();
			if(rgl.flowcells.isEmpty())
				{
				throw new PicardException("empty XML "+XML);
				}
			
			


			
			if(IN!=null)
				{
				LOG.info("opening "+IN);
				sfr=new SAMFileReader(IN);
				}
			else
				{
				LOG.info("opening stdin");
				sfr=new SAMFileReader(System.in);
				}
			sfr.setValidationStringency(super.VALIDATION_STRINGENCY);
			
			
			SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
			sfwf.setMaxRecordsInRam(super.MAX_RECORDS_IN_RAM);
			if(super.TMP_DIR.isEmpty()) sfwf.setTempDirectory(super.TMP_DIR.get(0));
	
			
			
			SAMProgramRecord sp=new SAMProgramRecord(getClass().getSimpleName());
			sp.setProgramName(getClass().getSimpleName());
			sp.setProgramVersion(String.valueOf(getProgramVersion()));
			sp.setPreviousProgramGroupId(getClass().getSimpleName());
			sp.setCommandLine(getCommandLine().replace('\t', ' '));
			sfr.getFileHeader().addProgramRecord(sp);
			
			Set<String> seenids=new HashSet<String>();
			List<SAMReadGroupRecord> samReadGroupRecords=new ArrayList<SAMReadGroupRecord>();
			for(FlowCell fc:rgl.flowcells)
				{
				Map<Integer,String> lane2id=new HashMap<Integer, String>();
				for(Lane lane:fc.lanes)
					{
					
					for(ReadGroup rg:lane.readGroups)
						{
						if(seenids.contains(rg.id))
							{
							throw new PicardException("Group id "+rg.id +" defined twice");
							}
						seenids.add(rg.id);
						 // create the read group we'll be using
				        SAMReadGroupRecord rgrec = new SAMReadGroupRecord(rg.id);
				        rgrec.setLibrary(rg.library);
				        rgrec.setPlatform(rg.platform);
				        rgrec.setSample(rg.sample);
				        rgrec.setPlatformUnit(rg.platform);
				        if (rg.center != null) rgrec.setSequencingCenter(rg.center);
				        if (rg.description != null) rgrec.setDescription(rg.description);
				        lane2id.put(lane.id,rg.id);
				        samReadGroupRecords.add(rgrec);
						}
					}
				this.flowcell2lane2id.put(fc.name,lane2id);
				}
			 sfr.getFileHeader().setReadGroups(samReadGroupRecords);
			
			boolean presorted=false;
			switch(sfr.getFileHeader().getSortOrder())
				{
				case coordinate: 
				case queryname: presorted=true;break;
				default:break;
				}
			if(OUT!=null)
				{
				LOG.info("opening "+OUT);
				sfwf.setCreateIndex(super.CREATE_INDEX);
				sfwf.setCreateMd5File(super.CREATE_MD5_FILE);
				sfw=sfwf.makeBAMWriter(sfr.getFileHeader(), presorted, OUT);
				}
			else
				{
				LOG.info("opening stdout");
				sfw=sfwf.makeSAMWriter(sfr.getFileHeader(), presorted, System.out);
				}
			
			final Pattern colon=Pattern.compile("[\\:]");
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				
				String RGID=null;
				String tokens[]=colon.split(rec.getReadName(),3);
				if(tokens.length!=3)
					{
					throw new PicardException("Cannot split "+rec.getReadName());
					}
				
				Map<Integer,String> lane2id=flowcell2lane2id.get(tokens[0]);
				if(lane2id==null) throw new PicardException("Cannot get flowcell/readgroup for "+rec.getReadName());
				try
					{
					RGID=lane2id.get(Integer.parseInt(tokens[1]));
					}
				catch (Exception e)
					{
					throw new PicardException("bad lane-Id in "+rec.getReadName());
					}
				
				if(RGID==null) 
					{
					throw new PicardException("Cannot get RGID for "+rec.getReadName());
					}
				rec.setAttribute(SAMTag.RG.name(), RGID);
				sfw.addAlignment(rec);
				}
			iter.close();
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
		new Biostar78400().instanceMain(args);
		}
	}
