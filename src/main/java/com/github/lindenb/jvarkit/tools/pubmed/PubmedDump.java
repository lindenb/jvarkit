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


import java.io.PrintWriter;
import java.net.URLEncoder;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

/**
 * PubmedDump
 *
 */
public class PubmedDump
	extends AbstractCommandLineProgram
	{
	private String email=null;
	private String tool="pubmedump";

	private PubmedDump()
		{
		
		}
	@Override
	public String getProgramDescription() {
		return "Dump pubmed articles as XML. Can read more than ret_max=100000 articles.";
		}			
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/PubmedDump";
		}
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -e (email) user's email. Optional.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"e:"))!=-1)
			{
			switch(c)
				{
				case 'e': email=opt.getOptArg();break;
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
		StringBuilder query=new StringBuilder();
		if(opt.getOptInd()==args.length)
			{
			error("Query missing");
			return -1;
			}
		for(int i=opt.getOptInd();i< args.length;++i)
			{
			if(query.length()>0) query.append(" ");
			query.append( args[i]);
			}
		if(query.toString().trim().isEmpty())
			{
			error("Query is empty");
			return -1;
			}

		
		try
			{
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.FALSE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
						
			String url=
					"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="+
					URLEncoder.encode(query.toString(), "UTF-8")+	
					"&retstart=0&retmax=0&usehistory=y&retmode=xml"+
					(email==null?"":"&email="+URLEncoder.encode(email,"UTF-8"))+
					(tool==null?"":"&tool="+URLEncoder.encode(tool,"UTF-8"))
					;
			info(url);
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
			r.close();
			if(count<0 || WebEnv==null || QueryKey==null)
				{
				error("Bad esearch result");
				return -1;
				}
			PrintWriter pw=new PrintWriter(System.out);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			XMLEventWriter w=xof.createXMLEventWriter(pw);
			long nFound=0L;
			int depth=0;
			while(nFound< count)
				{
				final int ret_max=90000;
				info("nFound:"+nFound+"/"+count);
				url= "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&WebEnv="+
						URLEncoder.encode(WebEnv,"UTF-8")+
						"&query_key="+URLEncoder.encode(QueryKey,"UTF-8")+
						"&retmode=xml&retmax="+ret_max+"&retstart="+nFound+
						(email==null?"":"&email="+URLEncoder.encode(email,"UTF-8"))+
						(tool==null?"":"&tool="+URLEncoder.encode(tool,"UTF-8"))
						;
				info(url);
				int curr_count=0;
				r=xmlInputFactory.createXMLEventReader(new StreamSource(url));
				
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
							warning("XML evt no handled: "+evt);
							break;
							}
						}
					}
				
				r.close();
				if(curr_count==0)
					{
					info("Nothing found . Exiting.");
					}
				}
			w.flush();
			w.close();
			pw.flush();
			pw.close();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			
			}
		}
	
	public static void main(String[] args) {
		new PubmedDump().instanceMainWithExit(args);
	}
}
