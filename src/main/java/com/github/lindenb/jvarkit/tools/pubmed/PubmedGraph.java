package com.github.lindenb.jvarkit.tools.pubmed;


import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;


import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

public class PubmedGraph extends AbstractCommandLineProgram
	{
	private String rootPmid=null;
	private int maxDepth=3;
	private XMLInputFactory xmlInputFactory;
	private Map<String, Article> pmid2article=new HashMap<String, PubmedGraph.Article>();
	
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
		List<String> L=new ArrayList<>();
		String url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?retmode=xml&db=pubmed&dbfrom=pubmed&" +
				"id="+pmid+
				"&linkname="+linkname;
		info(url);
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
	/** Citation referenced in PubMed article. */
	private List<String> getCitesReferences(String pmid) throws IOException,XMLStreamException
		{
		return eLink(pmid,"pubmed_pubmed_refs");
		}
	private List<String> getCitedIn(String pmid) throws IOException,XMLStreamException
		{
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
			System.err.println(this.pmid2article.size());
			}
		return article;
		}
	
	private void run(Article article,int depth)throws IOException,XMLStreamException
		{
		if(depth>this.maxDepth) return;
		if(article.visited) return;
		article.visited=true;
		System.err.println("depth="+depth);
		for(String cites: getCitesReferences(article.pmid))
			{
			add_A_cites_B(article,cites);
			run(article,depth+1);
			}
		
		for(String cited: getCitedIn(article.pmid))
			{
			Article cited_article = getArticleByPmid(cited);
			add_A_cites_B(getArticleByPmid(cited),article.pmid);
			run(cited_article,depth+1);
			}
		}
	private void graphiz(PrintStream out)
		{
		out.println("digraph G  {");
		out.println("node [ shape = point ];");
		for(String pmid : this.pmid2article.keySet())
			{
			Article article= this.pmid2article.get(pmid);
			out.println("p"+article.pmid+"[shape=point,label=\"\" ];");
			for(String k: article.cites)
				{
				out.println("p"+k+" -> p"+article.pmid+";");
				}
			}
		out.println("}");
		}
	
	@Override
	public int doWork(String[] args)
		{
		try
			{
			this.xmlInputFactory=XMLInputFactory.newFactory();
			
			run(getArticleByPmid("19505943"),0);
			graphiz(System.out);
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
	
		}
	public static void main(String[] args)
		{
		new PubmedGraph().instanceMainWithExit(args);
		}
	}
