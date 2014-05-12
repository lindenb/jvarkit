package com.github.lindenb.jvarkit.tools.bam2xml;

import java.io.File;
import java.io.FileInputStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;


import htsjdk.samtools.cmdline.CommandLineProgram;
import htsjdk.samtools.cmdline.Option;
import htsjdk.samtools.cmdline.StandardOptionDefinitions;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;


public class Bam2Xml extends CommandLineProgram
	{
		
	private static final Log log = Log.getInstance(Bam2Xml.class);
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="A BAM file to process.")
    public File INPUT=null;
    
    public boolean verboseflag=false;
    
    private void writeCharacters(XMLStreamWriter w,String tag,Object value)
    	throws XMLStreamException
    	{
    	if(value==null) return;
    	w.writeStartElement(tag);
    	w.writeCharacters(String.valueOf(value));
    	w.writeEndElement();
    	}
    
	@Override
	protected int doWork()
		{    	
        SAMFileReader samReader = null;
        SAMRecordIterator iter=null;
        try
	        {
	        XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
	        XMLStreamWriter w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
	        w.writeStartDocument("UTF-8","1.0");
	        w.writeStartElement("bam");
	        samReader=new SAMFileReader(INPUT==null?System.in:new FileInputStream(INPUT));
	        final SAMFileHeader header=samReader.getFileHeader();
	        w.writeStartElement("header");
	        
	        w.writeAttribute("version",header.getVersion());
	        
	        
	        w.writeEndElement();
	        samReader.setValidationStringency(super.VALIDATION_STRINGENCY);
			iter=samReader.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				w.writeStartElement("sam");
				w.writeAttribute("flag", String.valueOf(rec.getFlags()));
				if(verboseflag)
					{
					if(rec.getReadPairedFlag())
						{
						w.writeAttribute("paired", "true");
						}
					else
						{
						w.writeAttribute("paired", "false");
						}
					}
				writeCharacters(w,"name",rec.getReadName());
				w.writeEndElement();
				}
			
			w.writeEndElement();
			w.writeEndDocument();
			w.flush();
			w.close();
			} 
    	catch (Exception e) {
    		log.error(e);
    		return -1;
			}
        finally
	    	{
	    	if(iter!=null) iter.close();
	    	if(samReader!=null) samReader.close();
	    	}
    	return 0;
    	}
    public static void main(final String[] argv)
		{
	    new Bam2Xml().instanceMainWithExit(argv);
		}	


	}
