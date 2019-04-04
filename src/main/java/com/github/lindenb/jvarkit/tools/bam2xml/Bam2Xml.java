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
package com.github.lindenb.jvarkit.tools.bam2xml;

import java.io.File;
import java.io.OutputStream;
import java.util.List;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.ProgressLoggerInterface;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;

/**
BEGIN_DOC



```
 samtools view samtools-1.4/examples/toy.sam |\
 	java -jar dist/bam2xml.jar  | xmllint --format -


<?xml version="1.0" encoding="UTF-8"?>
<sam>
  <header version="1.5" sort="unsorted">
    <dict size="0"/>
    <read-groups/>
    <program-records/>
  </header>
  <record id="1" flag="163" length="19" read-paired="true" proper-pair="true" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="true" first-of-pair="false" second-of-pair="true" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="7" unclipped-align-start="7" align-end="22" unclipped-align-end="22" mapq="30" mate-ref-name="ref" mate-tid="-1" mate-align-start="37" insert-size="39">
    <name>r001</name>
    <seq>TTAGATAAAGAGGATACTG</seq>
    <cigar>
      <ce op="M" length="8" read-pos="0" ref-pos="7"/>
      <ce op="I" length="4" read-pos="8"/>
      <ce op="M" length="4" read-pos="12" ref-pos="15"/>
      <ce op="D" length="1" ref-pos="19"/>
      <ce op="M" length="3" read-pos="16" ref-pos="20"/>
    </cigar>
    <attributes>
      <attribute name="XX">[S@1e81f4dc</attribute>
    </attributes>
  </record>
  <record id="2" flag="0" length="17" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="9" unclipped-align-start="8" align-end="18" unclipped-align-end="18" mapq="30">
    <name>r002</name>
    <seq>AAAAGATAAGGGATAAA</seq>
    <cigar>
      <ce op="S" length="1" read-pos="0" ref-pos="8"/>
      <ce op="I" length="2" read-pos="1"/>
      <ce op="M" length="6" read-pos="3" ref-pos="9"/>
      <ce op="P" length="1"/>
      <ce op="I" length="1" read-pos="9"/>
      <ce op="P" length="1"/>
      <ce op="I" length="1" read-pos="10"/>
      <ce op="M" length="4" read-pos="11" ref-pos="15"/>
      <ce op="I" length="2" read-pos="15"/>
    </cigar>
    <attributes/>
  </record>
  <record id="3" flag="0" length="6" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="9" unclipped-align-start="4" align-end="14" unclipped-align-end="14" mapq="30">
    <name>r003</name>
    <seq>AGCTAA</seq>
    <cigar>
      <ce op="H" length="5" ref-pos="4"/>
      <ce op="M" length="6" read-pos="0" ref-pos="9"/>
    </cigar>
    <attributes/>
  </record>
  <record id="4" flag="0" length="12" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="16" unclipped-align-start="16" align-end="40" unclipped-align-end="40" mapq="30">
    <name>r004</name>
    <seq>ATAGCTCTCAGC</seq>
    <cigar>
      <ce op="M" length="6" read-pos="0" ref-pos="16"/>
      <ce op="N" length="14" ref-pos="22"/>
      <ce op="I" length="1" read-pos="6"/>
      <ce op="M" length="5" read-pos="7" ref-pos="36"/>
    </cigar>
    <attributes/>
  </record>
  <record id="5" flag="16" length="5" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="true" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="29" unclipped-align-start="23" align-end="33" unclipped-align-end="33" mapq="30">
    <name>r003</name>
    <seq>TAGGC</seq>
    <cigar>
      <ce op="H" length="6" ref-pos="23"/>
      <ce op="M" length="5" read-pos="0" ref-pos="29"/>
    </cigar>
    <attributes/>
  </record>
  <record id="6" flag="83" length="9" read-paired="true" proper-pair="true" read-unmapped="false" mate-unmapped="false" read-reverse-strand="true" mate-reverse-strand="false" first-of-pair="true" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="37" unclipped-align-start="37" align-end="45" unclipped-align-end="45" mapq="30" mate-ref-name="ref" mate-tid="-1" mate-align-start="7" insert-size="-39">
    <name>r001</name>
    <seq>CAGCGCCAT</seq>
    <cigar>
      <ce op="M" length="9" read-pos="0" ref-pos="37"/>
    </cigar>
    <attributes/>
  </record>
  <record id="7" flag="0" length="20" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="1" unclipped-align-start="1" align-end="20" unclipped-align-end="20" mapq="30">
    <name>x1</name>
    <seq>AGGTTTTATAAAACAAATAA</seq>
    <qual>????????????????????</qual>
    <cigar>
      <ce op="M" length="20" read-pos="0" ref-pos="1"/>
    </cigar>
    <attributes/>
  </record>
  <record id="8" flag="0" length="21" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="2" unclipped-align-start="2" align-end="22" unclipped-align-end="22" mapq="30">
    <name>x2</name>
    <seq>GGTTTTATAAAACAAATAATT</seq>
    <qual>?????????????????????</qual>
    <cigar>
      <ce op="M" length="21" read-pos="0" ref-pos="2"/>
    </cigar>
    <attributes/>
  </record>
  <record id="9" flag="0" length="26" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="6" unclipped-align-start="6" align-end="27" unclipped-align-end="27" mapq="30">
    <name>x3</name>
    <seq>TTATAAAACAAATAATTAAGTCTACA</seq>
    <qual>??????????????????????????</qual>
    <cigar>
      <ce op="M" length="9" read-pos="0" ref-pos="6"/>
      <ce op="I" length="4" read-pos="9"/>
      <ce op="M" length="13" read-pos="13" ref-pos="15"/>
    </cigar>
    <attributes/>
  </record>
  <record id="10" flag="0" length="25" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="10" unclipped-align-start="10" align-end="34" unclipped-align-end="34" mapq="30">
    <name>x4</name>
    <seq>CAAATAATTAAGTCTACAGAGCAAC</seq>
    <qual>?????????????????????????</qual>
    <cigar>
      <ce op="M" length="25" read-pos="0" ref-pos="10"/>
    </cigar>
    <attributes/>
  </record>
  <record id="11" flag="0" length="24" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="12" unclipped-align-start="12" align-end="35" unclipped-align-end="35" mapq="30">
    <name>x5</name>
    <seq>AATAATTAAGTCTACAGAGCAACT</seq>
    <qual>????????????????????????</qual>
    <cigar>
      <ce op="M" length="24" read-pos="0" ref-pos="12"/>
    </cigar>
    <attributes/>
  </record>
  <record id="12" flag="0" length="23" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="14" unclipped-align-start="14" align-end="36" unclipped-align-end="36" mapq="30">
    <name>x6</name>
    <seq>TAATTAAGTCTACAGAGCAACTA</seq>
    <qual>???????????????????????</qual>
    <cigar>
      <ce op="M" length="23" read-pos="0" ref-pos="14"/>
    </cigar>
    <attributes/>
  </record>
</sam>

```

END_DOC
 */
@Program(name="bam2xml",
	keywords={"sam","bam","xml"},
	description="converts a BAM to XML")
public class Bam2Xml extends Launcher
	{
	private static final Logger LOG = Logger.build(Bam2Xml.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	public static class SAMXMLWriter implements SAMFileWriter
		{
		private long id_generator=0L;
		private ProgressLoggerInterface progress;
		private XMLStreamWriter w;
		private SAMFileHeader header;
		
		private void writeText(final XMLStreamWriter w,final String tag,Object o) throws XMLStreamException
			{
			if(o==null) return;
			final String s=String.valueOf(o);
			if(StringUtil.isBlank(s)) return;
			w.writeStartElement(tag);
			w.writeCharacters(s);
			w.writeEndElement();
			}
		
		public SAMXMLWriter(final XMLStreamWriter w,final SAMFileHeader header) throws XMLStreamException
			{
			this.w=w;
			this.header=header;
			
	        w.writeStartElement("sam");
	        
	        w.writeStartElement("header");
	        w.writeAttribute("version",header.getVersion());
	        w.writeAttribute("sort",header.getSortOrder().name());
	        
    		if(header.getCreator()!=null) w.writeAttribute("creator", header.getCreator());
	        
	        final SAMSequenceDictionary dict = header.getSequenceDictionary();
	        if(dict==null)
	        	{
	        	w.writeComment("no dictionary available");
	        	}
	        else
	        	{
	        	w.writeStartElement("dict");
	        	w.writeAttribute("size", String.valueOf(dict.size()));
	        	for(final SAMSequenceRecord ssr: dict.getSequences())
	        		{
	        		w.writeStartElement("contig");
	        		w.writeAttribute("index", String.valueOf(ssr.getSequenceIndex()));
	        		w.writeAttribute("name",ssr.getSequenceName());
	        		writeText(w,"assembly", String.valueOf(ssr.getAssembly()));
	        		writeText(w,"species", String.valueOf(ssr.getSpecies()));
	        		writeText(w,"length", String.valueOf(ssr.getSequenceLength()));
	        		writeText(w,"md5",ssr.getMd5());
	        		w.writeEndElement();
	        		}
	        	w.writeEndElement();
	        	}
	        if(!header.getReadGroups().isEmpty()) {
		        w.writeStartElement("read-groups");
		        for(final SAMReadGroupRecord g: header.getReadGroups())
			       	{
		        	w.writeStartElement("read-group");
		        	w.writeAttribute("id",g.getId());
		        	writeText(w,"sample",g.getSample());
		        	writeText(w,"description",g.getDescription());
		        	writeText(w,"flow-order",g.getFlowOrder());
		        	writeText(w,"key-sequence",g.getKeySequence());
		        	writeText(w,"platform",g.getPlatform());
		        	writeText(w,"platform-model",g.getPlatformModel());
		        	writeText(w,"platform-unit",g.getPlatformUnit());
		        	writeText(w,"program-group",g.getProgramGroup());
		        	writeText(w,"library",g.getLibrary());
		        	writeText(w,"center",g.getSequencingCenter());
		        	writeText(w,"insert-size",g.getPredictedMedianInsertSize());
		        	writeText(w,"run-date",g.getRunDate());
			       	w.writeEndElement();
			       	}
		        w.writeEndElement();
		        }
	        if(!header.getProgramRecords().isEmpty()) {
		        w.writeStartElement("program-records");
		        for(final SAMProgramRecord sp:header.getProgramRecords())
		        	{
		        	w.writeStartElement("program-record");
		        	w.writeAttribute("id",sp.getId());
		        	if(sp.getPreviousProgramGroupId()!=null) w.writeAttribute("prev-group-id",sp.getPreviousProgramGroupId());
		        	if(sp.getProgramName()!=null) w.writeAttribute("name",sp.getProgramName());
		        	if(sp.getProgramVersion()!=null) w.writeAttribute("version",sp.getProgramVersion());
		        	if(sp.getCommandLine()!=null)  w.writeCharacters(sp.getCommandLine());
		        	w.writeEndElement();
		        	}
		        w.writeEndElement();
		        }
	        
	        for(final String comm: header.getComments())
	        	{
	        	writeText(w,"comment",comm);
	        	}
	      
	        
	        w.writeEndElement();//header

			}
		
	    private void writeCharacters(final String tag,Object value) throws XMLStreamException
		    	{
	    		writeText(w, tag, value);
		    	}

		
		@Override
		public void setProgressLogger(ProgressLoggerInterface progress) {
			this.progress=progress;
			}
		@Override
		public SAMFileHeader getFileHeader() {
			return header;
			}
		@Override
		public void addAlignment(final SAMRecord rec) {
			try {
				if(progress!=null) progress.record(rec);
				++id_generator;
				w.writeStartElement("record");
				w.writeAttribute("id",String.valueOf(id_generator));
				w.writeAttribute("flag", String.valueOf(rec.getFlags()));
				w.writeAttribute("length", String.valueOf(rec.getReadLength()));
				for(final SAMFlag sf: SAMFlag.values())
					{
					w.writeAttribute(sf.name().toLowerCase().replace('_', '-'),String.valueOf(sf.isSet(rec.getFlags())));
					}
				

				if(!rec.getReadUnmappedFlag())
					{
					w.writeAttribute("ref-name",rec.getReferenceName());
					w.writeAttribute("tid",String.valueOf(rec.getReferenceIndex()));
					w.writeAttribute("align-start",String.valueOf(rec.getAlignmentStart()));
					w.writeAttribute("unclipped-align-start",String.valueOf(rec.getUnclippedStart()));
					w.writeAttribute("align-end",String.valueOf(rec.getAlignmentEnd()));
					w.writeAttribute("unclipped-align-end",String.valueOf(rec.getAlignmentEnd()));
					w.writeAttribute("mapq",String.valueOf(rec.getMappingQuality()));
					
					
					}
				if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag())
					{
					w.writeAttribute("mate-ref-name",rec.getMateReferenceName());
					w.writeAttribute("mate-tid",String.valueOf(rec.getMateReferenceIndex()));
					w.writeAttribute("mate-align-start",String.valueOf(rec.getMateAlignmentStart()));
					if(!rec.getReadUnmappedFlag())
						{
						w.writeAttribute("insert-size",String.valueOf(rec.getInferredInsertSize()));
						}
					}
				
				final String readString= rec.getReadString();
				writeCharacters("name",rec.getReadName());
				if(!SAMRecord.NULL_SEQUENCE_STRING.equals(rec.getReadString()))
					{
					writeCharacters("seq",readString);
					}
				if(!SAMRecord.NULL_QUALS_STRING.equals(rec.getBaseQualityString()))
					{
					writeCharacters("qual",rec.getBaseQualityString());
					}
				if(!rec.getReadUnmappedFlag() && rec.getCigar()!=null)
					{
					int readPos=0;
					int refpos1=rec.getUnclippedStart();
					final Cigar cigar = rec.getCigar();
					w.writeStartElement("cigar");
					for(final CigarElement ce: cigar.getCigarElements())
						{
						w.writeEmptyElement("ce");
						final CigarOperator op =ce.getOperator();
						w.writeAttribute("op",op.name());
						w.writeAttribute("length",String.valueOf(ce.getLength()));
						if(op.consumesReadBases())
							{
							w.writeAttribute("read-pos", String.valueOf(readPos));
							readPos+=ce.getLength();
							}
						
						if(op.consumesReferenceBases() || op.isClipping())
							{
							w.writeAttribute("ref-pos", String.valueOf(refpos1));
							refpos1+=ce.getLength();
							}
						}
					w.writeEndElement();
					}
				w.writeStartElement("attributes");
				for(final SAMTagAndValue attribute: rec.getAttributes())
					{
					w.writeStartElement("attribute");
	                w.writeAttribute("name", attribute.tag);
	                w.writeCharacters(String.valueOf(attribute.value));
					w.writeEndElement();
					}
				w.writeEndElement();
				
				w.writeEndElement();
			} catch (final Exception e) {
				throw new RuntimeException(e);
				}
			}
		@Override
		public void close()
			{
			try {
				w.writeEndElement();
				w.flush();
			} catch (Exception e) {
				throw new RuntimeException(e);
				}
			}
		}
	
	
	
	
		private int run(SamReader samReader)
			{    	
			OutputStream fout=null;
	        SAMRecordIterator iter=null;
	        XMLStreamWriter w= null;
	        try
		        {
		        XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
		        
		        if(this.outputFile!=null)
			        {
		        	fout= IOUtils.openFileForWriting(this.outputFile);
			        w= xmlfactory.createXMLStreamWriter(fout,"UTF-8");
			        }
		        else
		        	{
		        	w= xmlfactory.createXMLStreamWriter(stdout(),"UTF-8");
		        	}
		        w.writeStartDocument("UTF-8","1.0");
		        final SAMFileHeader header=samReader.getFileHeader();
		        final SAMXMLWriter xw =new SAMXMLWriter(w, header);
		        final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
				iter=samReader.iterator();
				while(iter.hasNext())
					{
					xw.addAlignment(progress.watch(iter.next()));
					}
				xw.close();
				w.writeEndDocument();
				if(fout!=null) fout.flush();
				} 
	    	catch (Exception e) {
	    		e.printStackTrace();
	    		LOG.error(e);
	    		return -1;
				}
	        finally
		    	{
	        	CloserUtil.close(w);
		    	CloserUtil.close(iter);
		    	CloserUtil.close(fout);
		    	}
	    	return 0;
	    	}
		
	
	@Override
	public int doWork(final List<String> args) {
			SamReader r=null;
			try
				{
				r= openSamReader(oneFileOrNull(args));
				run(r);
				return 0;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(r);
				}
			}
		
	
    public static void main(final String[] argv)
		{
	    new Bam2Xml().instanceMainWithExit(argv);
		}	


	}
