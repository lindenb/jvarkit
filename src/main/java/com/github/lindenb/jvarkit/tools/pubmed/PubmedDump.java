/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.pubmed;


import java.io.ByteArrayInputStream;
import java.io.PrintWriter;
import java.net.URLEncoder;
import java.util.Collection;
import java.util.List;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import htsjdk.samtools.util.CloserUtil;

/**
 * PubmedDump
 *
 */
public class PubmedDump
	extends AbstractPubmedDump
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(PubmedDump.class);

	private String tool="pubmedump";

	public PubmedDump()
		{
		
		}
	
	@Override
	public Collection<Throwable> call() throws Exception {
		PrintWriter pw=null;
		final List<String> args = getInputFiles();
		final StringBuilder query=new StringBuilder();
		if(args.isEmpty())
			{
			return wrapException("Query missing");
			}
		for(final String arg:args)
			{
			if(query.length()>0) query.append(" ");
			query.append(arg);
			}
		if(query.toString().trim().isEmpty())
			{
			return wrapException("Query is empty");
			}

		
		try
			{
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.FALSE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			xmlInputFactory.setXMLResolver(new XMLResolver() {
				@Override
				public Object resolveEntity(String publicID, String systemID, String baseURI, String namespace)
						throws XMLStreamException {
					LOG.info("ignoring "+publicID+" "+baseURI+" "+namespace);
					return new ByteArrayInputStream(new byte[0]);
				}
			});
			
			String url=
					"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="+
					URLEncoder.encode(query.toString(), "UTF-8")+	
					"&retstart=0&retmax=0&usehistory=y&retmode=xml"+
					(email==null?"":"&email="+URLEncoder.encode(email,"UTF-8"))+
					(tool==null?"":"&tool="+URLEncoder.encode(tool,"UTF-8"))
					;
			LOG.info(url);
			long count=-1;
			String WebEnv=null;
			String QueryKey=null;
			XMLEventReader r=xmlInputFactory.createXMLEventReader(new StreamSource(url));
			while(r.hasNext())
				{
				XMLEvent evt=r.nextEvent();
				if(evt.isStartElement())
					{
					String eName=evt.asStartElement().getName().getLocalPart();
					if(eName.equals("Count") && count==-1)
						{
						count=Long.parseLong(r.getElementText());
						}
					else if(eName.equals("WebEnv"))
						{
						WebEnv= r.getElementText();
						}
					else if(eName.equals("QueryKey"))
						{
						QueryKey= r.getElementText();
						}
					}
				}
			CloserUtil.close(r);
			r=null;
			
			if(count<0 || WebEnv==null || QueryKey==null)
				{
				return wrapException("Bad esearch result");
				}
			pw=super.openFileOrStdoutAsPrintWriter();
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			XMLEventWriter w=xof.createXMLEventWriter(pw);
			long nFound=0L;
			int depth=0;
			while(nFound< count)
				{
				final int ret_max=90000;
				LOG.info("nFound:"+nFound+"/"+count);
				url= "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&WebEnv="+
						URLEncoder.encode(WebEnv,"UTF-8")+
						"&query_key="+URLEncoder.encode(QueryKey,"UTF-8")+
						"&retmode=xml&retmax="+ret_max+"&retstart="+nFound+
						(email==null?"":"&email="+URLEncoder.encode(email,"UTF-8"))+
						(tool==null?"":"&tool="+URLEncoder.encode(tool,"UTF-8"))
						;
				LOG.info(url);
				int curr_count=0;
				r = xmlInputFactory.createXMLEventReader(new StreamSource(url));
				
				while(r.hasNext())
					{
					XMLEvent evt=r.nextEvent();
					
					switch(evt.getEventType())
						{
						case XMLEvent.ATTRIBUTE:
							{
							if(depth>0) w.add(evt);
							break;
							}
						case XMLEvent.START_DOCUMENT:
							{
							if(nFound==0)
								{
								w.add(evt);
								}
							break;
							}
						case XMLEvent.END_DOCUMENT:
							{
							if(nFound>= count)
								{
								w.add(evt);
								}
							break;
							}
						case XMLEvent.START_ELEMENT:
							{
							if(depth==0 && nFound==0)
								{
								w.add(evt);
								}
							else if(depth==1)
								{
								String localName= evt.asStartElement().getName().getLocalPart();
								if(!(localName.equals("PubmedArticle") || localName.equals("PubmedBookArticle")))
									{
									throw new IllegalStateException("Not PubmedArticle: "+evt);
									}
								++curr_count;
								++nFound;
								w.add(evt);
								}
							else if(depth>1)
								{
								w.add(evt);
								}
							depth++;
							break;
							}
						case XMLEvent.END_ELEMENT:
							{
							depth--;
							if(depth>0)
								{
								w.add(evt);
								}
							else if(nFound>=count)//depth ==0
								{
								w.add(evt);
								}
							
							break; 
							}
						case XMLEvent.COMMENT:break;
						case XMLEvent.PROCESSING_INSTRUCTION:break;
						case XMLEvent.DTD:
							{
							if(nFound==0) w.add(evt);
							break;	
							}
						case XMLEvent.SPACE:break;
						case XMLEvent.CHARACTERS:
							{
							if(depth>1) w.add(evt);
							break;
							}
						default:
							{
							return wrapException("XML evt no handled: "+evt);
							}
						}
					}
				
				r.close();
				if(curr_count==0)
					{
					LOG.info("Nothing found . Exiting.");
					}
				}
			w.flush();
			w.close();
			pw.flush();
			pw.close();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(pw);
			}
		}
	
	public static void main(String[] args) {
		new PubmedDump().instanceMainWithExit(args);
	}
}
