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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.pubmed;


import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="pubmedgraph",
	description="Creates a Gephi-gexf graph of references-cotes for a given PMID",
	keywords={"pubmed","xml","graph"}
)
public class PubmedGraph extends Launcher
	{
	private static final Logger LOG = Logger.build(PubmedGraph.class).make();

	
	@Parameter(names="-d",description="max-depth")
	private int maxDepth=3;
	private XMLInputFactory xmlInputFactory;
	private Map<String, Article> pmid2article=new HashMap<String, PubmedGraph.Article>();
	@Parameter(names="-f",description="disable forward (cited-in)")
	private boolean disable_cited_in=false;
	@Parameter(names="-b",description="disable backward (referenced-in)")
	private boolean disable_references=false;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outFile=null;

	private class Article
		{
		String pmid;
		String title=null;
		String year=null;
		boolean visited=false;
		Set<String> cites = new HashSet<>();
		}
	
	private List<String> eLink(String pmid,String linkname) throws IOException,XMLStreamException
		{
		final List<String> L=new ArrayList<>();
		final String url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?retmode=xml&db=pubmed&dbfrom=pubmed&" +
				"id="+pmid+
				"&linkname="+linkname;
		LOG.info(url);
		StreamSource src=new StreamSource(url);
		XMLEventReader reader= this.xmlInputFactory.createXMLEventReader(src);
		int in_link=0;
		while(reader.hasNext())
			{
			XMLEvent evt=reader.nextEvent();
			if(evt.isStartElement())
				{
				String localName=evt.asStartElement().getName().getLocalPart();
				if(localName.equals("Link"))
					{
					in_link=1;
					}
				else if(localName.equals("Id") && in_link==1)
					{
					L.add(reader.getElementText().trim());
					}
				}
			else if(evt.isEndElement())
				{
				String localName=evt.asEndElement().getName().getLocalPart();
				if(localName.equals("Link"))
					{
					in_link=0;
					}
				}
			}
		reader.close();
		return L;
		}
	private void eSummary(Article a) throws IOException,XMLStreamException
		{
		final QName attName=new QName("Name");
		final QName attType=new QName("Type");
		String url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?retmode=xml&db=pubmed&" +
				"id="+a.pmid
				;
		LOG.info(url);
		StreamSource src=new StreamSource(url);
		XMLEventReader reader= this.xmlInputFactory.createXMLEventReader(src);
		int in_title=0;
		int in_date=0;
		while(reader.hasNext() && !(a.title!=null && a.year!=null))
			{
			XMLEvent evt=reader.nextEvent();
			if(evt.isStartElement())
				{
				StartElement startE=evt.asStartElement();
				String localName=startE.getName().getLocalPart();
				if(localName.equals("Item"))
					{
					Attribute name= startE.getAttributeByName(attName);
					Attribute type= startE.getAttributeByName(attType);
					if(name.getValue().equals("Title") && type.getValue().equals("String"))
						{
						in_title=1;
						}
					else if(name.getValue().equals("PubDate") && type.getValue().equals("Date"))
						{
						in_date=1;
						}
					}
				
				}
			else if(evt.isCharacters())
				{
				if(in_title==1)
					{
					a.title = evt.asCharacters().getData();
					}
				else if(in_date==1)
					{
					a.year = evt.asCharacters().getData();
					}
				in_title=0;
				in_date=0;
				}
			}
		reader.close();
		}
	
	
	
	/** Citation referenced in PubMed article. */
	private List<String> getCitesReferences(String pmid) throws IOException,XMLStreamException
		{
		if(this.disable_references) return Collections.emptyList();
		return eLink(pmid,"pubmed_pubmed_refs");
		}
	private List<String> getCitedIn(String pmid) throws IOException,XMLStreamException
		{
		if(this.disable_cited_in) return Collections.emptyList();
		return eLink(pmid,"pubmed_pubmed_citedin");
		}

	private void add_A_cites_B(Article article,String cited)
		{
		article.cites.add(cited);
		}
	
	private Article getArticleByPmid(String pmid)
		{
		Article article = this.pmid2article.get(pmid);
		if(article==null)
			{
			article=new Article();
			article.pmid=pmid;
			this.pmid2article.put(pmid,article);
			}
		return article;
		}
	
	private void run(Article article,int depth)throws IOException,XMLStreamException
		{
		if(depth>this.maxDepth) return;
		if(article.visited) return;
		article.visited=true;
		
		for(String cites: getCitesReferences(article.pmid))
			{
			add_A_cites_B(article,cites);
			run(getArticleByPmid(cites),depth+1);
			}
			
		
		for(String cited: getCitedIn(article.pmid))
			{
			Article cited_article = getArticleByPmid(cited);
			add_A_cites_B(getArticleByPmid(cited),article.pmid);
			run(cited_article,depth+1);
			}
			
		}
	
	private void gexf() throws XMLStreamException,IOException
			{
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			XMLStreamWriter w= null;
			FileWriter fw=null;
			
			if(this.outFile==null)
				{
				w=xof.createXMLStreamWriter(System.out, "UTF-8");
				}
			else
				{
				w=xof.createXMLStreamWriter((fw=new FileWriter(this.outFile)));
				}
			w.writeStartDocument("UTF-8", "1.0");
			w.writeStartElement("gexf");
			w.writeDefaultNamespace("http://www.gexf.net/1.2draft");
			w.writeAttribute("version", "1.2");
			w.writeStartElement("meta");
			  w.writeStartElement("creator");
			  w.writeCharacters("PubmedGraph  by Pierre Lindenbaum");
			  w.writeEndElement();
			 
			  w.writeStartElement("description");
			  w.writeCharacters("");
			  w.writeEndElement();
			
			w.writeEndElement();//meta
			  
			  w.writeStartElement("graph");
			  w.writeAttribute("mode", "static");
			  w.writeAttribute("defaultedgetype", "directed");
			  
			  w.writeStartElement("attributes");
			  w.writeAttribute("class", "node");
			  w.writeAttribute("mode", "static");
			  w.writeEndElement();//attributes
				
			  w.writeStartElement("attributes");                                                                                     
			  w.writeAttribute("class", "edge");
			  w.writeAttribute("mode", "static");
				  
	          w.writeEmptyElement("attribute");
				w.writeAttribute("id", "0");
				w.writeAttribute("title", "pmid");
				w.writeAttribute("type", "string");
		      w.writeEmptyElement("attribute");
				w.writeAttribute("id", "1");
				w.writeAttribute("title", "title");
				w.writeAttribute("type", "string");
			  w.writeEmptyElement("attribute");
				w.writeAttribute("id", "2");
				w.writeAttribute("title", "pubdate");
				w.writeAttribute("type", "string");
				
	          w.writeEndElement();//attributes
			  
	
			  
			  w.writeStartElement("nodes");
			  for(Article a:this.pmid2article.values())
			 	{
				w.writeStartElement("node");
				w.writeAttribute("id",a.pmid);
				if(a.title!=null && a.title.length()>30)
					{
					w.writeAttribute("label",a.title.substring(0, 27)+"...");
					}
				else
					{
					w.writeAttribute("label",String.valueOf(a.title));
					}
				w.writeStartElement("attvalues");
				
				  w.writeEmptyElement("attvalue");
					w.writeAttribute("for", "0");
					w.writeAttribute("value",a.pmid);

				  w.writeEmptyElement("attvalue");
					w.writeAttribute("for", "1");
					w.writeAttribute("value",String.valueOf(a.title));
			
				
				  w.writeEmptyElement("attvalue");
					w.writeAttribute("for", "2");
					w.writeAttribute("value",String.valueOf(a.year));

				
				w.writeEndElement();//attvalues
				w.writeEndElement();//node
			 	}
			  w.writeEndElement();//nodes
			
			  
			  w.writeStartElement("edges");
			  for(Article a:this.pmid2article.values())
			 	{
				for(String pmid: a.cites)
					{
					w.writeStartElement("edge");
					w.writeAttribute("id","E"+pmid+"_"+a.pmid);
					w.writeAttribute("type","directed");
					w.writeAttribute("source",pmid);
					w.writeAttribute("target",a.pmid);
					w.writeAttribute("weight",String.valueOf(1));
					w.writeEndElement();
					}
			 	}
			  w.writeEndElement();//edges
			  
			  w.writeEndElement();//graph
	
			
			w.writeEndElement();//gexf
			w.writeEndDocument();
			w.flush();
			if(fw!=null)
				{
				fw.flush();
				CloserUtil.close(fw);
				}
			else
				{
				System.out.flush();
				}
			}
	
	@Override
	public int doWork(final List<String> args) {
		if(args.size()!=1)
			{
			LOG.error("Illegal number of arguments");
			return -1;
			}
		try
			{
			this.xmlInputFactory=XMLInputFactory.newFactory();
			
			run(getArticleByPmid(args.get(0)),0);
			for(final Article a:this.pmid2article.values())
				{
				eSummary(a);
				}
			
			gexf();
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}

		}

	
	public static void main(String[] args)
		{
		new PubmedGraph().instanceMainWithExit(args);
		}
	}
