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
package com.github.lindenb.jvarkit.tools.pubmed;


import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.net.URL;
import java.net.URLEncoder;
import java.text.Collator;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import com.sleepycat.bind.tuple.LongBinding;
import com.sleepycat.bind.tuple.StringBinding;
import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.LockMode;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.gexf.GexfConstants;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ncbi.NcbiApiKey;
import com.github.lindenb.jvarkit.util.ncbi.NcbiConstants;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

/** 
BEGIN_DOC

## About Orcid
 ORCID  (http://orcid.org) provides a persistent digital identifier that distinguishes you from every other researcher and, through integration in key research workflows such as manuscript and grant submission, supports automated linkages between you and your professional activities ensuring that your work is recognized.

You can download the papers containing some orcid Identifiers using the entrez query

http://www.ncbi.nlm.nih.gov/pubmed/?term=orcid[AUID]

I've used one of my tools pubmeddump to download the articles as XML and I wrote PubmedOrcidGraph to extract the author's orcid.
The output is a GEXF file for gephi.

## Example

using pubmed efetch output

```
java -jar dist/pubmeddump.jar --skip "MeshHeadingList ChemicalList GrantList InvestigatorList CommentsCorrectionsList ISSN DateRevised AffiliationInfo Language PublicationTypeList  ArticleDate PubmedData Abstract MedlineJournalInfo CoiStatement KeywordList Pagination ELocationID "   "orcid[AUID]" |\
java -jar dist/pubmedorcidgraph.jar -D BDB 
```

using orcid identifiers:

java -jar dist/pubmedorcidgraph.jar -D BDB --orcid 0000-0001-7751-2280 0000-0003-0677-5627 0000-0003-3628-1548 0000-0003-4530-6655 0000-0001-8007-5931 


<img src="https://pbs.twimg.com/media/Ci-h0MJWUAAvJjw.jpg"/>

## See also</h:h3>

 * http://plindenbaum.blogspot.fr/2016/05/playing-with-orcidorg-ncbipubmed-graph.html

END_DOC
*/
@Program(name="pubmedorcidgraph",
	description="Creates a graph from Pubmed and Authors' Orcid identifiers",
	keywords={"pubmed","ncbi","orcid"}
	)
public class PubmedOrcidGraph
	extends Launcher
	{
	private static final Logger LOG = Logger.build(PubmedOrcidGraph.class).make();
	private static final String NAME_NOT_FOUND="<NOT FOUND IN PUBMED>";
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-D","--berkeydb"},description="BerkeleyDB tmpDir",required=true)
	private File bdbDir = null;
	@Parameter(names={"-E","--errors"},description="Dump strange orcids (e.g: same orcird but different forename) . Default: stderr")
	private File dumpStrangeOrcidFile = null;
	private PrintStream errPrintWriter=System.err;

	
	@Parameter(names={"-d","--maxdepth"},description="Max graph depth")
	private int maxdepth = 2;
	@Parameter(names={"-orcid","--orcid"},description="Input is a set of orcids identifiers")
	private boolean input_is_orcid_id = false;
	@Parameter(names={"-links","--alllinks"},description="By default, we display only one link between two authors. Using this option will show all the links (publications)")
	private boolean all_links_between_authors = false;
	@Parameter(names={"-ea","--edgeattributes"},description="Do not show edge attributes (smaller files with less informations)")
	private boolean hide_edge_attributes = false;
	@ParametersDelegate
	private NcbiApiKey ncbiApiKey = new NcbiApiKey();

	
	private long ID_GENERATOR = 0L;
	private Database authorDatabase=null;
	private Database articleDatabase=null;
	private Database linkDatabase=null;
	private Environment environment=null;
	private Transaction txn=null;
	
	private static String _s(final String s) {
		return s==null?"":s;
	}
	
	private static class Link
		{
		String orcid1=null;
		String orcid2=null;
		Set<String> pmids = new HashSet<>();
		}
	
	private static class LinkBinding extends TupleBinding<Link>
		{
		@Override
		public Link entryToObject(TupleInput in) {
			final Link L = new Link();
			L.orcid1 = in.readString();
			L.orcid2 = in.readString();
			
			
			int c=in.readInt();
			for(int i=0;i< c;++i) 
				{
				L.pmids.add(in.readString());
				}
			return L;
			}
		@Override
		public void objectToEntry(Link L, TupleOutput out) {
			out.writeString(L.orcid1);
			out.writeString(L.orcid2);
			out.writeInt(L.pmids.size());
			for(final String o:L.pmids) {
				out.writeString(o);
				}
			}
		}
	private final LinkBinding linkBinding=new LinkBinding();
	
	private static class Article implements Comparable<Article>
		{
		String pmid=null;
		String ArticleTitle=null;
		String Year=null;
		String ISOAbbreviation=null;
		String doi  = null;
		
		@Override
		public int compareTo(Article o) {
			return pmid.compareTo(o.pmid);
			}
		

		

		}
	private static class ArticleBinding extends TupleBinding<Article>
		{
		@Override
		public Article entryToObject(TupleInput in) {
			final Article a = new Article();
			a.pmid = in.readString();
			a.ArticleTitle= in.readString();
			a.Year = in.readString();
			a.ISOAbbreviation = in.readString();
			a.doi = in.readString();
			
			return a;
			}
		@Override
		public void objectToEntry(Article a, TupleOutput w) {
			w.writeString(_s(a.pmid));
			w.writeString(_s(a.ArticleTitle));
			w.writeString(_s(a.Year));
			w.writeString(_s(a.ISOAbbreviation));
			w.writeString(_s(a.doi));
			
			}
		}
	private final ArticleBinding articleBinding = new ArticleBinding();

	

	private static class Author implements Comparable<Author>
		{
		String foreName = null;
		String lastName = null ;
		String orcid = null;
		String initials=null;
		String affiliation =null;
		boolean reviewed=false;
		int depth=-1;
		
		
		/** idea: in xml+xslt , used to presort orcid indentifiers
		 * in order to produce unique pairs of collab(orcid1,orcid2) */
		@Override
		public int compareTo(final Author o) {
			return this.orcid.compareTo(o.orcid);
			}
		
		
		void gexf(final XMLStreamWriter w) throws XMLStreamException {
			
			w.writeStartElement("node");
			w.writeAttribute("id", this.orcid);
			w.writeAttribute("label", String.join(" ",foreName,lastName));
			w.writeStartElement("attvalues");
			


			w.writeEmptyElement("attvalue");
			w.writeAttribute("for","orcid");
			w.writeAttribute("value",this.orcid);
			
			w.writeEmptyElement("attvalue");
			w.writeAttribute("for","foreName");
			w.writeAttribute("value",this.foreName);
			
			w.writeEmptyElement("attvalue");
			w.writeAttribute("for","lastName");
			w.writeAttribute("value",this.lastName);

			w.writeEmptyElement("attvalue");
			w.writeAttribute("for","initials");
			w.writeAttribute("value",this.initials);

			w.writeEmptyElement("attvalue");
			w.writeAttribute("for","affiliation");
			w.writeAttribute("value",this.affiliation);

			
			w.writeEndElement();//attvalues
			w.writeEndElement();//node
		}
		
		}	
	private static class AuthorBinding extends TupleBinding<Author>
		{
		@Override
		public Author entryToObject(final TupleInput in) {
			final Author a=new Author();
			a.foreName = in.readString();
			a.lastName = in.readString();
			a.orcid = in.readString();
			a.initials = in.readString();
			a.affiliation = in.readString();
			a.reviewed = in.readBoolean();
			a.depth = in.readInt();
			return a;
			}
		@Override
		public void objectToEntry(Author a, TupleOutput w) {
			w.writeString(_s(a.foreName));
			w.writeString(_s(a.lastName));
			w.writeString(_s(a.orcid));
			w.writeString(_s(a.initials));
			w.writeString(_s(a.affiliation));
			w.writeBoolean(a.reviewed);
			w.writeInt(a.depth);
			}
		}
	private final AuthorBinding authorBinding = new AuthorBinding();
	
	
	/** extract text from stream. Cannot use XMLEventReader.getTextContent() 
	 * when a public title contains some tag like '<sup>'
	 */
	private String textContent(final XMLEventReader r) throws XMLStreamException
		{
		final StringBuilder sb=new StringBuilder();
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isEndElement()) break;
			else if(evt.isStartElement())
				{
				sb.append(textContent(r));
				}
			else if(evt.isCharacters())
				{
				sb.append(evt.asCharacters().getData());
				}
			}
		return sb.toString();
		}
	
	private Author parseAuthor(XMLEventReader r,final String pmid,int depth)  throws XMLStreamException
		{
		final Author au = new Author();
		au.depth=depth;
		au.reviewed=false;
		while(r.hasNext()) {
			final XMLEvent evt=r.nextEvent();
			if(evt.isStartElement()) {
				final StartElement start = evt.asStartElement();
				String eltName = start.getName().getLocalPart();
				if(eltName.equals("LastName")) {
					au.lastName=r.getElementText().trim();
				} else if(eltName.equals("ForeName") || eltName.equals("FirstName")) {
					au.foreName=r.getElementText().trim();
				} else if(eltName.equals("Affiliation")) {
					au.affiliation=r.getElementText().trim();
				} else if(eltName.equals("Initials")) {
					au.initials=r.getElementText().trim();
				} else if(eltName.equals("Identifier")) {
					final Attribute source= start.getAttributeByName(new QName("Source"));
					
					if(source!=null && !source.getValue().equalsIgnoreCase("ORCID")) {
						LOG.warn("interesting, a non-orcid Identifier in pmid:"+pmid+" "+source.getValue());
						}
					else if(source!=null && source.getValue().equalsIgnoreCase("ORCID")) {
						au.orcid=r.getElementText().trim();
						final int slash = au.orcid.lastIndexOf('/');
						if(slash!=-1) au.orcid=au.orcid.substring(slash+1);
						au.orcid=au.orcid.trim();
						
						if(au.orcid.startsWith("orcid.org.")) {
							LOG.warn("Ugly orcid in "+au.orcid +" pmid:"+ pmid);
							final int dot= au.orcid.lastIndexOf('.');
							au.orcid=au.orcid.substring(dot+1);
						}
						else if(au.orcid.startsWith("orcid.org") ) {
							LOG.warn("Ugly orcid in "+au.orcid +" pmid:"+ pmid);
							au.orcid=au.orcid.substring(9);
						}
						
						
						if(au.orcid.length()==16){
							// 0000000226481522 instead of 0000-0002-2648-1522
							au.orcid= au.orcid.substring(0, 4 ) +"-"+
									  au.orcid.substring(4, 8 ) +"-"+
									  au.orcid.substring(8, 12) +"-"+
									  au.orcid.substring(12)
									  ;
						}
					}
				}
			}
			else if(evt.isEndElement() &&
				evt.asEndElement().getName().getLocalPart().equals("Author")) {
				return au;
			}
		}
		throw new IllegalStateException("should never happen");
		}
	
	private void scanArticles(InputStream in,int depth) throws XMLStreamException,IOException
		{
		XMLEventReader r=null;
		try {
			r = createXMLInputFactory().createXMLEventReader(in);
			
			while(r.hasNext()) {
				final XMLEvent evt= r.nextEvent();
				if(evt.isStartElement() )	{
					final String localName= evt.asStartElement().getName().getLocalPart();
					if(localName.equals("PubmedArticle") || localName.equals("PubmedBookArticle"))
						{
						scanArticle(localName,r,depth);
						}
				}
			}
		}
		catch(final Exception err) {
			LOG.error(err);
			}
		 finally {
			CloserUtil.close(r);
		 }
		}
	
	private boolean complete() throws IOException{
			Author authorToScan=null;
			Cursor c=null;
			try {
				final DatabaseEntry key=new DatabaseEntry();
				final DatabaseEntry data=new DatabaseEntry();
				c=this.authorDatabase.openCursor(txn, null);
				while(c.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS) {
					final Author au  = this.authorBinding.entryToObject(data);
					if(au.reviewed) continue;
					if(au.depth>=maxdepth) continue;
					if(authorToScan==null || authorToScan.depth>au.depth) {
						authorToScan=au;
						}
					}
			} catch (final Exception err) {
				LOG.error(err);
				}
			finally
				{
				CloserUtil.close(c);
				}
			if(authorToScan!=null)
				{
				LOG.info("Now scanning orcid "+authorToScan.orcid+" depth:"+authorToScan.depth);
				scanOrcid(authorToScan.orcid, authorToScan.depth);
				return true;
				}
			
		return false;
		}
	
	private void scanArticle(
			final String rootName,
			final XMLEventReader r,
			int depth
			) throws XMLStreamException,IOException {

			final Article article=new Article();

			boolean PubDate=false;
			final List<Author> authors = new ArrayList<>();
			while(r.hasNext()) {
				final XMLEvent evt=r.nextEvent();
				
				if(evt.isStartElement()) {
					final StartElement start = evt.asStartElement();
					final String eltName = start.getName().getLocalPart();
					if(article.pmid==null && eltName.equals("PMID")) {
						article.pmid=r.getElementText();
					}
					else if(article.ISOAbbreviation==null && eltName.equals("ISOAbbreviation")) {
						article.ISOAbbreviation=r.getElementText();
					}
					else if(article.ArticleTitle==null && eltName.equals("ArticleTitle")) {
						article.ArticleTitle = this.textContent(r);
					}
					else if(eltName.equals("PubDate")) {
						PubDate=true;
					}
					else if(article.doi==null && eltName.equals("ArticleId")) {
						final Attribute idType= start.getAttributeByName(new QName("IdType"));
						if(idType!=null && idType.getValue().equalsIgnoreCase("doi")) {
							article.doi = r.getElementText().trim();
						}
					}
					else if(article.Year==null && PubDate && eltName.equals("Year")) {
						article.Year=r.getElementText();
					}
					else if(eltName.equals("Author")) {
						final Author author = parseAuthor(r,article.pmid,depth+1);
						if(author!=null && author.orcid!=null) {
							authors.add(author);
							}
						}
					}
			else if(evt.isEndElement()) {
				final EndElement end = evt.asEndElement();
				final String eltName = end.getName().getLocalPart();
				if(eltName.equals("PubDate")) {
					PubDate=false;
					}
				if(eltName.equals(rootName)) break;
				}
			}//end of xml read
			
		
		
		if(authors.isEmpty()) {
			//do nothing
			return;
			}
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		Collections.sort(authors);

		
		for(int x=0;x+1< authors.size();++x)
			{
			for(int y=x+1;y< authors.size();++y)
				{
				Link L=null;
				final String orcid1=authors.get(x).orcid;
				final String orcid2 =authors.get(y).orcid;
				
				if(this.all_links_between_authors) 
					{
					LongBinding.longToEntry(++ID_GENERATOR, key);
					}
				else
					{
					StringBinding.stringToEntry(orcid1+"~"+orcid2, key);
					if(this.linkDatabase.get(txn, key, data, LockMode.DEFAULT)!=OperationStatus.NOTFOUND){
						L = this.linkBinding.entryToObject(data);
						}
					}
				if(L == null) {
					L = new Link();
					L.orcid1=orcid1;
					L.orcid2=orcid2;
					}
				L.pmids.add(article.pmid);
				this.linkBinding.objectToEntry(L, data);
				if(this.linkDatabase.put(txn, key, data)!=OperationStatus.SUCCESS) 
					{
					throw new JvarkitException.BerkeleyDbError("Cannot put in article db");
					}
				}
			}
		//for comparing names
		final Collator collator= Collator.getInstance(Locale.US);
		collator.setStrength(Collator.PRIMARY);

		StringBinding.stringToEntry(article.pmid, key);
		if(this.articleDatabase.get(txn, key, data, LockMode.DEFAULT)!=OperationStatus.NOTFOUND)
			{
			LOG.debug("Article already in database : "+article.pmid);
			}
		else
			{
			LOG.debug("inserting article "+article.pmid);
			this.articleBinding.objectToEntry(article, data);
			if(this.articleDatabase.put(txn, key, data)!=OperationStatus.SUCCESS) 
				{
				throw new JvarkitException.BerkeleyDbError("Cannot put in article db");
				}
			
			for(final Author au:authors) 
				{
				StringBinding.stringToEntry(au.orcid, key);
				
				if(this.authorDatabase.get(txn, key, data, LockMode.DEFAULT)!=OperationStatus.NOTFOUND)
					{
					LOG.debug("Author already in database : "+au.orcid);
					final Author other=this.authorBinding.entryToObject(data);
					if(!StringUtil.isBlank(other.lastName) && !StringUtil.isBlank(au.lastName) && 
							collator.compare(au.lastName,other.lastName)!=0)
						{
						this.errPrintWriter.println("Conflict\t"+au.orcid+"\t"+au.foreName+"\t"+other.foreName);
						}
					}
				else
					{
					this.authorBinding.objectToEntry(au, data);
					if(this.authorDatabase.put(txn, key, data)!=OperationStatus.SUCCESS) 
						{
						throw new JvarkitException.BerkeleyDbError("Cannot put in author db");
						}
					}
				
				}
			}
		}
	
	private Author getAuthorByOrcid(final String orcid) {
		final DatabaseEntry key=new DatabaseEntry();
		final DatabaseEntry data=new DatabaseEntry();
		StringBinding.stringToEntry(orcid, key);
		if(this.authorDatabase.get(txn, key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS) 
			{
			return this.authorBinding.entryToObject(data);
			}
		else
			{
			return null;
			}
		}
	
	private Author insertAuthor(final Author au) {
		final DatabaseEntry key=new DatabaseEntry();
		final DatabaseEntry data=new DatabaseEntry();
		StringBinding.stringToEntry(au.orcid, key);
		this.authorBinding.objectToEntry(au, data);
		if(this.authorDatabase.put(txn, key, data)!=OperationStatus.SUCCESS) {
			throw new JvarkitException.BerkeleyDbError("Cannot update author");
			}
		return au;
		}
	
	private void scanOrcid(final String orcid,int depth) throws IOException {
		if(orcid==null || orcid.trim().isEmpty())
			{
			LOG.debug("empty orcid");
			return;
			}
		if(depth>=this.maxdepth) return;
			
		//author was already scanned
		Author au = getAuthorByOrcid(orcid);
		if(au!=null) {
			if(au.depth>depth)
				{
				au.depth=depth;
				insertAuthor(au);
				}
			if(au.reviewed) return;
			}
	
		String WebEnv=null;
		String QueryKey=null;
		InputStream in=null;
		String urlstr= new StringBuilder(NcbiConstants.esearch()).append("?db=pubmed&usehistory=1&retmax=100000&term=").			 
			append(URLEncoder.encode("\""+orcid+"\"[auid]","UTF-8")).
			append(this.ncbiApiKey.getAmpParamValue()).
			toString();
			;
		
		XMLEventReader r=null;
		/* first get NCBI WebEnv */
		try
			{
			LOG.debug(urlstr);
			in = new URL(urlstr).openStream();
			r= createXMLInputFactory().createXMLEventReader(in);
			while(r.hasNext()) {
				final XMLEvent evt=r.nextEvent();
				if(evt.isStartElement() && evt.asStartElement().getName().getLocalPart().equals("WebEnv")) {
					WebEnv = r.getElementText();
					}
				else if(evt.isStartElement() && evt.asStartElement().getName().getLocalPart().equals("QueryKey")) {
					QueryKey = r.getElementText();
					}
				if(QueryKey!=null && WebEnv!=null) break;
				}
			}
		catch(final Exception err) {
			LOG.error(err);
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(in);
			r=null;in=null;
			}
		if(QueryKey==null || WebEnv==null) {
			LOG.debug("Cannot get QueryKey/WebEnv "+urlstr);
			}
		else
			{
			urlstr= new StringBuilder(NcbiConstants.efetch()).
					append("?db=pubmed&usehistory=1&retmode=xml").
					append(this.ncbiApiKey.getAmpParamValue()).
					append("&query_key=").append(QueryKey).append("&webenv=").append(WebEnv).toString();
					;
			LOG.debug(urlstr);
			
			try
				{
				in = new URL(urlstr).openStream();
				scanArticles(in,depth);
				}
			catch(Exception err) {
				LOG.error(err);
				}
			finally
				{
				CloserUtil.close(in);
				r=null;in=null;
				}
			}
		
		/* we check we found an author with this orcid */
		au = getAuthorByOrcid(orcid);
		
		if(au == null) 
			{
			au = new Author();
			au.orcid=orcid;
			au.reviewed=true;
			au.depth=depth;
			au.foreName=NAME_NOT_FOUND;
			au.lastName=NAME_NOT_FOUND;
			}
		au.reviewed=true;
		if(au.depth>depth) au.depth=depth;
		insertAuthor(au);
		}
	
	
	private XMLInputFactory createXMLInputFactory() {
		final XMLInputFactory xmlInputFactory = XMLInputFactory.newFactory();
		xmlInputFactory.setXMLResolver(new XMLResolver() {
			@Override
			public Object resolveEntity(String publicID, String systemID, String baseURI, String namespace)
					throws XMLStreamException {
				LOG.debug("Ignoring resolve Entity");
				return new ByteArrayInputStream(new byte[0]);
			}
		});
		return xmlInputFactory;
	}
	
	
	private void gexfAttDecl(
			XMLStreamWriter w,
			String key,
			String type
			)throws XMLStreamException
			{
			w.writeEmptyElement("attribute");
			w.writeAttribute("id", key);
			w.writeAttribute("title", key.replace('_', ' '));
			w.writeAttribute("type", type);
			}

	private void dumpGexf()
		{
		final XMLOutputFactory xof = XMLOutputFactory.newFactory();
		PrintWriter pw=null;
		XMLStreamWriter w=null;
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		Cursor c=null;
		try {
			pw = openFileOrStdoutAsPrintWriter(this.outputFile);
			w = xof.createXMLStreamWriter(pw);
			
			w.writeStartDocument("UTF-8","1.0");
			w.writeStartElement("gexf");
			w.writeAttribute("xmlns",GexfConstants.XMLNS);
			w.writeAttribute("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance");
			w.writeAttribute("xsi:schemaLocation",GexfConstants.XSI_SCHEMA_LOCATION);
			w.writeAttribute("version", GexfConstants.VERSION);
			
			
			/* meta */
			w.writeStartElement("meta");
				w.writeAttribute("lastmodifieddate",new SimpleDateFormat("yyyy-MM-dd").format(new Date()));
				w.writeStartElement("creator");
				  w.writeCharacters("PumedOrcidGraph");
				w.writeEndElement();
				w.writeStartElement("description");
				  w.writeCharacters("PumedOrcidGraph");
				w.writeEndElement();
			w.writeEndElement();
			
			/* graph */
			w.writeStartElement("graph");
			w.writeAttribute("mode", "static");
			w.writeAttribute("defaultedgetype", "undirected");
	
			
			/* attributes */
			w.writeStartElement("attributes");
			w.writeAttribute("class","node");
			w.writeAttribute("mode","static");
			gexfAttDecl(w,"orcid","string");
			gexfAttDecl(w,"foreName","string");
			gexfAttDecl(w,"lastName","string");
			gexfAttDecl(w,"initials","string");
			gexfAttDecl(w,"affiliation","string");
			w.writeEndElement();//attributes
			
			if(!this.hide_edge_attributes) {
			w.writeStartElement("attributes");
			w.writeAttribute("class","edge");
			w.writeAttribute("mode","static");
			gexfAttDecl(w,"pmid","string");
			gexfAttDecl(w,"title","string");
			gexfAttDecl(w,"doi","string");
			gexfAttDecl(w,"year","string");
			gexfAttDecl(w,"journal","string");
			w.writeEndElement();//attributes
			}
			
			/* nodes */
			w.writeStartElement("nodes");
			c  = this.authorDatabase.openCursor(txn, null);
			while(c.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
				{
				final Author au=this.authorBinding.entryToObject(data);
				if(NAME_NOT_FOUND.equals(au.foreName))
					{
					w.writeComment("Orcid "+au.orcid+" not found in pubmed");
					continue;
					}
				au.gexf(w);
				}
			c.close();
			w.writeEndElement();//nodes
			
			w.writeStartElement("edges");
			key=new DatabaseEntry();
			c  = this.linkDatabase.openCursor(txn, null);
			while(c.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
				{
				final Link L=this.linkBinding.entryToObject(data);
				w.writeStartElement("edge");
				w.writeAttribute("id", "E"+(++ID_GENERATOR));
				w.writeAttribute("type", "undirected");
				w.writeAttribute("weight",String.valueOf(L.pmids.size()));
				w.writeAttribute("source",L.orcid1);
				w.writeAttribute("target",L.orcid2);
				for(String pmid: L.pmids) {
					DatabaseEntry key2= new DatabaseEntry();
					DatabaseEntry data2= new DatabaseEntry();
					StringBinding.stringToEntry(pmid, key2);
					if(this.articleDatabase.get(txn, key2, data2, LockMode.DEFAULT)!=OperationStatus.SUCCESS)
						{
						throw new JvarkitException.BerkeleyDbError("cannot get article");
						}
					final Article article = this.articleBinding.entryToObject(data2);
					if(this.all_links_between_authors)
						{
						w.writeAttribute("label",String.valueOf(L.pmids.size()));
						}
					else
						{
						w.writeAttribute("label",""+article.ArticleTitle+". " +article.ISOAbbreviation+". ("+article.Year+")");
						}
					if(!this.hide_edge_attributes) {
						w.writeStartElement("attributes");
						
						w.writeEmptyElement("attribute");
							w.writeAttribute("for", "pmid");
							w.writeAttribute("value", article.pmid);
						w.writeEmptyElement("attribute");
							w.writeAttribute("for", "title");
							w.writeAttribute("value", article.ArticleTitle);
						w.writeEmptyElement("attribute");
							w.writeAttribute("for", "doi");
							w.writeAttribute("value", article.doi);
						w.writeEmptyElement("attribute");
							w.writeAttribute("for", "year");
							w.writeAttribute("value", article.Year);
						w.writeEmptyElement("attribute");
							w.writeAttribute("for", "journal");
							w.writeAttribute("value", article.ISOAbbreviation);
								
						w.writeEndElement();//attributes
						}
					break;
					}
				//
				w.writeEndElement();//edge
				}
			c.close();
			
			w.writeEndElement();//edges

			
			w.writeEndElement();
			w.writeEndDocument();
			w.flush();
			pw.flush();
			pw.close();pw=null;
			}
		catch(Exception err) {
			throw new RuntimeIOException(err);
			}
		finally 
			{
			CloserUtil.close(w);
			CloserUtil.close(pw);
			}
		}

	
	@Override
	public int doWork(final List<String> args) {
		
		
		if(!this.ncbiApiKey.isApiKeyDefined()) {
			LOG.error("NCBI API key is not defined");
			return -1;
			}
		
		InputStream in=null;
		XMLEventReader r = null;
		try {
			
			//open BDB
			final EnvironmentConfig envCfg=new EnvironmentConfig();
			envCfg.setAllowCreate(true);
			envCfg.setReadOnly(false);
			LOG.info("open BDB env...");
			this.environment= new Environment(this.bdbDir, envCfg);
			LOG.info("open BDB databases...");
			final DatabaseConfig config = new DatabaseConfig();
			config.setAllowCreate(true);
			config.setReadOnly(false);
			config.setTemporary(true);
			this.authorDatabase=this.environment.openDatabase(txn, "authors",config);
			this.articleDatabase=this.environment.openDatabase(txn, "articles",config);
			this.linkDatabase=this.environment.openDatabase(txn, "links",config);
			
			if(this.dumpStrangeOrcidFile!=null) 
				{
				this.errPrintWriter = super.openFileOrStdoutAsPrintStream(this.dumpStrangeOrcidFile);
				}
			
			if(this.input_is_orcid_id)
				{
				/* recursively scan authors */
				for(final String orcid: args) {
					scanOrcid(orcid, 0);
					}
				while(complete()) {
					// run...
					}
				}
			else{
				/* input is a efetch stream */
				String inputName=oneFileOrNull(args);
				in=(inputName==null?stdin():IOUtils.openURIForReading(inputName));
				scanArticles(in,0);
				in.close();in=null;
				}
			
			
			dumpGexf();
			this.errPrintWriter.flush();
			if(this.dumpStrangeOrcidFile!=null) 
				{
				this.errPrintWriter.close();
				}
			return 0;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(in);
			CloserUtil.close(r);
			CloserUtil.close(authorDatabase);
			CloserUtil.close(articleDatabase);
			CloserUtil.close(linkDatabase);
			CloserUtil.close(environment);
			}
		}

	
	public static void main(String[] args)
		{
		new PubmedOrcidGraph().instanceMainWithExit(args);
		}
	}
