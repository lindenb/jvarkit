package com.github.lindenb.jvarkit.tools.bam2xml;

import java.io.PrintStream;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


public class Bam2Xml extends AbstractCommandLineProgram
	{
	private boolean verboseflag=true;
    
    private void writeCharacters(XMLStreamWriter w,String tag,Object value)
    	throws XMLStreamException
    	{
    	if(value==null) return;
    	w.writeStartElement(tag);
    	w.writeCharacters(String.valueOf(value));
    	w.writeEndElement();
    	}
    
	private int run(SamReader samReader)
		{    	
        SAMRecordIterator iter=null;
        try
	        {
	        XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
	        XMLStreamWriter w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
	        w.writeStartDocument("UTF-8","1.0");
	        w.writeStartElement("bam");
	        final SAMFileHeader header=samReader.getFileHeader();
	        w.writeStartElement("header");
	        
	        w.writeAttribute("version",header.getVersion());
	        w.writeEndElement();
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
    		error(e);
    		return -1;
			}
        finally
	    	{
	    	CloserUtil.close(iter);
	    	}
    	return 0;
    	}
	
	@Override
	public String getProgramDescription() {
		return "convert bam to xml";
		}
	
	@Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"Bam2Xml";
    }
	
	@Override
	public void printOptions(PrintStream out)
		{
		super.printOptions(out);
		}

	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		SamReader r=null;
		try
			{
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			if(opt.getOptInd()==args.length)
				{
				r = srf.open(SamInputResource.of(System.in));
				}
			else if(opt.getOptInd()+1==args.length)
				{
				r = srf.open(SamInputResource.of(args[opt.getOptInd()]));
				}
			else
				{	
				error("Illegal number of arguments.");
				return -1;
				}
			run(r);
			return 0;
			}
		catch(Exception err)
			{
			error(err);
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
