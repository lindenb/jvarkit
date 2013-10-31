package com.github.lindenb.jvarkit.tools.cmpbams;

import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.tools.cmpbams.entities.BamRecord;
import com.github.lindenb.jvarkit.tools.cmpbams.entities.Comparebams;
import com.github.lindenb.jvarkit.tools.cmpbams.entities.Records;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.util.Log;


public abstract class AbstractCompareBamsHandler extends AbstractCommandLineProgram
	{
	private static final Log LOG=Log.getInstance(AbstractCompareBamsHandler.class);

	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="XML file to process. Default stdin. ",optional=true)
	public File IN=null;
	@Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="output filename. if defined, the original XML is echoed to stdout, so you can pipeline many Handlers. ",optional=true)
	public File OUT=null;
	
	private Comparebams.Header header =null;
	
	@Override
	public String getVersion() {
		return "1.0";
		}
	
	protected void finish(PrintStream out)
		{
		
		}
	
	protected void setHeader(Comparebams.Header header)
		{
		this.header=header;
		}
	
	protected Comparebams.Header getHeader() {
		return header;
		}
	
	/**
	 * 
	 * @param out where the user should print his results
	 * @return true if 'records' should be unmarshalled in the next handler
	 */
	protected abstract boolean handleRecords(
			PrintStream out,
			final Records records
			);
	
	@Override
	public int doWork()
		{
		try {
						
			PrintStream out=System.out;
			XMLInputFactory xif=XMLInputFactory.newFactory();
			XMLEventReader xmlReader=null;
			FileReader fr=null;
			JAXBContext jbc=JAXBContext.newInstance(Comparebams.class,BamRecord.class,Records.class);
			Unmarshaller unmarshaller=jbc.createUnmarshaller();
			Marshaller marshaller=jbc.createMarshaller();
			marshaller.setProperty(Marshaller.JAXB_FRAGMENT, Boolean.TRUE);

			XMLEventWriter xmlWriter=null;
			if(IN==null)
				{
				LOG.info("reading from stdin ");
				xmlReader=xif.createXMLEventReader(System.in);
				}
			else
				{
				LOG.info("opening "+IN);
				fr=new FileReader(IN);
				xmlReader=xif.createXMLEventReader(fr);
				}
			
			if(OUT!=null)
				{
				//save xml to stdout
				LOG.info("printing XML to stdout");
				XMLOutputFactory xof=XMLOutputFactory.newFactory();
				xmlWriter=xof.createXMLEventWriter(System.out);
				LOG.info("opening "+OUT);
				out=new PrintStream(OUT);
				}
			Comparebams.Header header=null;
			while(xmlReader.hasNext())
				{
				XMLEvent evt=xmlReader.peek();
				if(evt.isStartElement())
					{
					StartElement E=evt.asStartElement();
					String localName=E.getName().getLocalPart();
					
					//next record was found
					if(localName.equals("records"))
						{
						
						Records records=unmarshaller.unmarshal(xmlReader,Records.class).getValue();
						if(!handleRecords(out, records))
							{
							records=null;
							}	
						
						if(records!=null && xmlWriter!=null)
							{
							marshaller.marshal(new JAXBElement<Records>(
									E.getName(),
									Records.class, records),
									xmlWriter
									);
							}
						}
					else if(localName.equals("header"))
						{
						header=unmarshaller.unmarshal(xmlReader,Comparebams.Header.class).getValue();
						this.setHeader(header);
						if(xmlWriter!=null)
							{
							marshaller.marshal(new JAXBElement<Comparebams.Header>(
									E.getName(),
									Comparebams.Header.class, header),
									xmlWriter
									);
							}
						}
					else //consumme element
						{
						evt=xmlReader.nextEvent();
						if(xmlWriter!=null) xmlWriter.add(evt);
						}
					}
				else
					{
					evt=xmlReader.nextEvent();
					//consumme element
					if(xmlWriter!=null) xmlWriter.add(evt);
					}
				}
			
			if(xmlWriter!=null)
				{
				xmlWriter.flush();
				xmlWriter.close();
				}
			if(fr!=null) fr.close();
			if(xmlReader!=null) xmlReader.close();
			finish(out);
			out.flush();
			out.close();
			return 0;
			
		} catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		}
	

	
	}
