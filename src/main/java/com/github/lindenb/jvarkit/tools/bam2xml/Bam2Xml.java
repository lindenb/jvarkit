package com.github.lindenb.jvarkit.tools.bam2xml;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.lang.reflect.Method;
import java.util.HashSet;
import java.util.Set;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.tools.samgrep.SamGrep;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;


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
