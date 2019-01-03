/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.transform.stream.StreamSource;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;

/**

BEGIN_DOC


### Read names

Reads' name should start with the following signature:


### XML

the XML should look like this:

```

<read-groups>
<flowcell name="HS2000-1259_127">
 <lane index="1">
   <group ID="X1">
     <library>L1</library>
     <platform>P1</platform>
     <sample>S1</sample>
     <platformunit>PU1</platformunit>
     <center>C1</center>
     <description>blabla</description>
   </group>
 </lane>
</flowcell>
<flowcell name="HS2000-1259_128">
 <lane index="2">
   <group ID="x2">
     <library>L2</library>
     <platform>P2</platform>
     <sample>S2</sample>
     <platformunit>PU1</platformunit>
     <center>C1</center>
     <description>blabla</description>
   </group>
 </lane>
</flowcell>
</read-groups>

```

### Example


```

$ cat input.sam 
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
HS2000-1259_127:1:1210:15640:52255  163 ref 7   30  8M4I4M1D3M  =   37  39  
TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
HS2000-1259_128:2:1210:15640:52255  0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   
0   AAAAGATAAGGGATAAA   *

$java -jar dist/biostar78400.jar \
    -x groups.xml \
    input.sam \
   

@HD VN:1.4  SO:unsorted
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
@RG ID:X1   PL:P1   PU:P1   LB:L1   DS:blabla   SM:S1   CN:C1
@RG ID:x2   PL:P2   PU:P2   LB:L2   DS:blabla   SM:S2   CN:C1
@PG ID:Biostar78400 PN:Biostar78400 PP:Biostar78400 VN:1.0  (...)
HS2000-1259_127:1:1210:15640:52255  163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   RG:Z:X1 XX:B:S,12561,2,20,112
HS2000-1259_128:2:1210:15640:52255  0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0AAAAGATAAGGGATAAA  *   RG:Z:x2

```

END_DOC
*/


@Program(name="biostar78400",
	keywords={"sam","bam","xml","read-group"},
	biostars= {78400,302798,202358},
	description="add the read group info to the sam file on a per lane basis")
public class Biostar78400 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar78400.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-x","--xmlFile"},description="XML description of the groups.",required=true)
	private File XML = null;

	@Parameter(names={"-p","--regex"},description="Regular expression that can be used to parse read names in the incoming SAM file. Flowcell: (group 1)and the lane (group 2). Another pattern could be '[a-zA-Z0-9\\-]+:[0-9]+:([a-zA-Z0-9]+):([0-9]):[0-9]+:[0-9]+:[0-9]+.*.' (Highseq)")
	private String readNameSignatureStr = "([a-zA-Z0-9]+):([0-9]):[0-9]+:[0-9]+:[0-9]+.*";

	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
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
	public int doWork(List<String> args) {
		if(this.XML==null)
			{
			LOG.error("XML file missing");
			return -1;
			}
		final Map<String, Map<Integer,String>> flowcell2lane2id = new HashMap<String, Map<Integer,String>>();
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			final Pattern readNameSignature = Pattern.compile(this.readNameSignatureStr);
			final JAXBContext context = JAXBContext.newInstance(ReadGroup.class,ReadGroupList.class);
			final Unmarshaller unmarshaller=context.createUnmarshaller();
			final ReadGroupList rgl=unmarshaller.unmarshal(new StreamSource(XML),ReadGroupList.class).getValue();
			if(rgl.flowcells.isEmpty())
				{
				LOG.error("empty XML "+XML); return -1;
				}
			sfr = openSamReader(oneFileOrNull(args));
			
			final SAMFileHeader header = sfr.getFileHeader().clone();
			header.addComment("Processed with "+getProgramName());
			
			final Set<String> seenids=new HashSet<String>();
			final List<SAMReadGroupRecord> samReadGroupRecords=new ArrayList<SAMReadGroupRecord>();
			for(final FlowCell fc:rgl.flowcells)
				{
				final Map<Integer,String> lane2id=new HashMap<Integer, String>();
				for(final Lane lane:fc.lanes)
					{
					for(final ReadGroup rg:lane.readGroups)
						{
						if(seenids.contains(rg.id))
							{
							LOG.error("Group id "+rg.id +" defined twice"); return -1;
							}
						seenids.add(rg.id);
						 // create the read group we'll be using
						final SAMReadGroupRecord rgrec = new SAMReadGroupRecord(rg.id);
				        rgrec.setLibrary(rg.library);
				        rgrec.setPlatform(rg.platform);
				        rgrec.setSample(rg.sample);
				        rgrec.setPlatformUnit(rg.platformunit);
				        if (rg.center != null) rgrec.setSequencingCenter(rg.center);
				        if (rg.description != null) rgrec.setDescription(rg.description);
				        lane2id.put(lane.id,rg.id);
				        samReadGroupRecords.add(rgrec);
						}
					}
				if(flowcell2lane2id.containsKey(fc.name))
					{
					LOG.error("FlowCell id "+fc.name +" defined twice in XML");return -1;
					}
				flowcell2lane2id.put(fc.name,lane2id);
				}
			header.setReadGroups(samReadGroupRecords);
			
		
			sfw = this.writingBamArgs.openSAMFileWriter(this.outputFile,header, true);
			
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);
			final SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec= progress.watch(iter.next());
				final Matcher matcher = readNameSignature.matcher(rec.getReadName());
				final String flowcellStr;
				final String laneStr;
				if(matcher.matches())
					{
					flowcellStr = matcher.group(1);
					laneStr = matcher.group(2);
					}
				else
					{
					LOG.error("Read name "+rec.getReadName()+" doesn't match regular expression "+readNameSignature.pattern()+". please check options");
					return -1;
					}
				String RGID=null;
				
				final Map<Integer,String> lane2id=flowcell2lane2id.get(flowcellStr);
				if(lane2id==null) throw new RuntimeException(
						"Cannot get flowcell/readgroup for "+rec.getReadName());
				try
					{
					RGID=lane2id.get(Integer.parseInt(laneStr));
					}
				catch (final Exception e)
					{
					LOG.error("bad lane-Id in "+rec.getReadName());
					return -1;
					}
				
				if(RGID==null) 
					{
					throw new RuntimeException("Cannot get RGID for "+rec.getReadName()+" flowcell:"+flowcellStr +" lane2id:"+laneStr+ " dict:"+lane2id);
					}
				rec.setAttribute(SAMTag.RG.name(), RGID);
				sfw.addAlignment(rec);
				}
			progress.finish();
			iter.close();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfw);
			CloserUtil.close(sfr);
			}
		}
	
	public static void main(String[] args)throws Exception
		{
		new Biostar78400().instanceMainWithExit(args);
		}
	}
