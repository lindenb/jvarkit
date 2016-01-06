/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

package com.github.lindenb.jvarkit.tools.biostar;
import java.util.ArrayList;
import java.util.Collection;
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

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;


public class Biostar78400 extends AbstractBiostar78400
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(Biostar78400.class);

	
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
	
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(this.XML==null)
			{
			return wrapException("XML file missing");
			}
		final Map<String, Map<Integer,String>> flowcell2lane2id = new HashMap<String, Map<Integer,String>>();
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			JAXBContext context = JAXBContext.newInstance(ReadGroup.class,ReadGroupList.class);
			Unmarshaller unmarshaller=context.createUnmarshaller();
			ReadGroupList rgl=unmarshaller.unmarshal(new StreamSource(XML),ReadGroupList.class).getValue();
			if(rgl.flowcells.isEmpty())
				{
				throw new RuntimeException("empty XML "+XML);
				}
			
			
			sfr = openSamReader(inputName);
			
			final SAMFileHeader header = sfr.getFileHeader().clone();
			
			
			SAMProgramRecord sp=new SAMProgramRecord(getClass().getSimpleName());
			sp.setProgramName(getClass().getSimpleName());
			sp.setProgramVersion(String.valueOf(getVersion()));
			sp.setPreviousProgramGroupId(getClass().getSimpleName());
			sp.setCommandLine(getProgramCommandLine().replace('\t', ' '));
			header.addProgramRecord(sp);
			
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
							return wrapException("Group id "+rg.id +" defined twice");
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
				flowcell2lane2id.put(fc.name,lane2id);
				}
			header.setReadGroups(samReadGroupRecords);
			
		
			sfw = openSAMFileWriter(header, true);
			
			final Pattern colon=Pattern.compile("[\\:]");
			SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec= progress.watch(iter.next());
				
				String RGID=null;
				String tokens[]=colon.split(rec.getReadName(),3);
				if(tokens.length!=3)
					{
					throw new RuntimeException("Cannot split "+rec.getReadName());
					}
				
				Map<Integer,String> lane2id=flowcell2lane2id.get(tokens[0]);
				if(lane2id==null) throw new RuntimeException("Cannot get flowcell/readgroup for "+rec.getReadName());
				try
					{
					RGID=lane2id.get(Integer.parseInt(tokens[1]));
					}
				catch (Exception e)
					{
					throw new RuntimeException("bad lane-Id in "+rec.getReadName());
					}
				
				if(RGID==null) 
					{
					throw new RuntimeException("Cannot get RGID for "+rec.getReadName());
					}
				rec.setAttribute(SAMTag.RG.name(), RGID);
				sfw.addAlignment(rec);
				}
			progress.finish();
			iter.close();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(sfw);
			CloserUtil.close(sfr);
			}
		}
	
	public static void main(String[] args)throws Exception
		{
		new Biostar78400().instanceMain(args);
		}
	}
