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

@Program(name="bam2xml",keywords={"sam","bam","xml"},description="converts a BAM to XML")
public class Bam2Xml extends Launcher
	{
	private static final Logger LOG = Logger.build(Bam2Xml.class).make();
	
	@Parameter(names={"-o","--out"},description="Output file or stdout")
	private File outputFile = null;

	public static class SAMXMLWriter implements SAMFileWriter
		{
		//private final SAMTagUtil tagUtil = new SAMTagUtil();
		private long id_generator=0L;
		private ProgressLoggerInterface progress;
		private XMLStreamWriter w;
		private SAMFileHeader header;
		public SAMXMLWriter(XMLStreamWriter w,final SAMFileHeader header) throws XMLStreamException
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
	        		if(ssr.getAssembly()!=null) w.writeAttribute("assembly", String.valueOf(ssr.getAssembly()));
	        		if(ssr.getSpecies()!=null) w.writeAttribute("species", String.valueOf(ssr.getSpecies()));
	        		w.writeAttribute("length", String.valueOf(ssr.getSequenceLength()));
	        		if(ssr.getMd5()!=null) w.writeAttribute("md5",ssr.getMd5());
	        		w.writeCharacters(ssr.getSequenceName());
	        		w.writeEndElement();
	        		}
	        	w.writeEndElement();
	        	}
	        
	        w.writeStartElement("read-groups");
	        for(final SAMReadGroupRecord g: header.getReadGroups())
		       	{
	        	w.writeEmptyElement("read-group");
	        	w.writeAttribute("id",g.getId());
	        	if(g.getSample()!=null) w.writeAttribute("sample",g.getSample());
	        	if(g.getDescription()!=null) w.writeAttribute("description",g.getDescription());
	        	if(g.getFlowOrder()!=null) w.writeAttribute("flow-order",g.getFlowOrder());
	        	if(g.getKeySequence()!=null) w.writeAttribute("key-sequence",g.getKeySequence());
	        	if(g.getPlatform()!=null) w.writeAttribute("platform",g.getPlatform());
	        	if(g.getPlatformModel()!=null) w.writeAttribute("platform-model",g.getPlatformModel());
	        	if(g.getPlatformUnit()!=null) w.writeAttribute("platform-unit",g.getPlatformUnit());
	        	if(g.getProgramGroup()!=null) w.writeAttribute("program-group",g.getProgramGroup());
	        	if(g.getLibrary()!=null) w.writeAttribute("library",g.getLibrary());
	        	if(g.getSequencingCenter()!=null) w.writeAttribute("center",g.getSequencingCenter());
	        	if(g.getPredictedMedianInsertSize()!=null) w.writeAttribute("insert-size",String.valueOf(g.getPredictedMedianInsertSize()));
	        	if(g.getRunDate()!=null) w.writeAttribute("run-date",String.valueOf(g.getRunDate()));
		       	}
	        w.writeEndElement();
	        
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
	        
	        for(final String comm: header.getComments())
	        	{
	        	w.writeStartElement("comment");
	        	 w.writeCharacters(comm);
	        	w.writeEndElement();
	        	}
	      
	        
	        w.writeEndElement();//header

			}
		
	    private void writeCharacters(String tag,Object value) throws XMLStreamException
		    	{
		    	if(value==null) return;
		    	w.writeStartElement(tag);
		    	w.writeCharacters(String.valueOf(value));
		    	w.writeEndElement();
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
				
				String readString= rec.getReadString();
				writeCharacters("name",rec.getReadName());
				writeCharacters("seq",readString);
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
					for(CigarElement ce: cigar.getCigarElements())
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
						
						if(op.consumesReferenceBases())
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
			} catch (Exception e) {
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
