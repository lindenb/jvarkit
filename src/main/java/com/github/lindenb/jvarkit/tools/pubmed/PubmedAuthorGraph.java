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


import java.awt.Color;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.text.Normalizer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

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
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

/** 
BEGIN_DOC

## Motivation

builds a graph with XML pubmed. Nodes are 'Article' and  'Authors'. Edges are authorships.

## Example

```
java -jar dist/pubmeddump.jar 'SCN5A redon' |\
java -jar dist/pubmedauthorgraph.jar -D BDB -singleton   > out.gexf
```

## Screenshots

https://twitter.com/yokofakun/status/1034107797439504384

![https://twitter.com/yokofakun/status/1034107797439504384](https://pbs.twimg.com/media/DlnjTj1W4AIBB6B.jpg)

https://twitter.com/yokofakun/status/1034397660189523968

![https://twitter.com/yokofakun/status/1034397660189523968](https://pbs.twimg.com/media/DlrqXqvX4AE6r27.jpg)



END_DOC
*/
@Program(name="pubmedauthorgraph",
	description="Creates a graph from Pubmed and Authors",
	keywords={"pubmed","ncbi","graph"}
	)
public class PubmedAuthorGraph
	extends Launcher
	{
	private static final Logger LOG = Logger.build(PubmedAuthorGraph.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-D","--berkeydb"},description="BerkeleyDB tmpDir",required=true)
	private File bdbDir = null;
	@Parameter(names={"-S","-singleton","--singleton"},description="Remove singleton authors (authors linked to only one paper)")
	private boolean remove_singleton_authors = false;
	@Parameter(names={"-r","--scale-articles"},description="Scale articles in function of the number of authors")
	private boolean scale_articles = false;
	@Parameter(names={"-u","--scale-authors"},description="Scale authors in function of the number of articles")
	private boolean scale_authors = false;
	@Parameter(names={"-i","--initals"},description="use author's initials to build the author-identifier. In the old pubmed record, the forename is not available.")
	private boolean use_initials_to_build_sample_id = false;	
	@Parameter(names={"-uc","--author-color"},description="viz:Color for the Authors." +ColorUtils.Converter.OPT_DESC,converter=ColorUtils.Converter.class)
	private Color authorColor = null;
	@Parameter(names={"-rc","--article-color"},description="viz:Color for the Articles." +ColorUtils.Converter.OPT_DESC,converter=ColorUtils.Converter.class)
	private Color articleColor = null;
	
	
	@ParametersDelegate
	private NcbiApiKey ncbiApiKey = new NcbiApiKey();

	
	private long ID_GENERATOR = 0L;
	private Database authorDatabase=null;
	private Database articleDatabase=null;
	private Environment environment=null;
	private Transaction txn=null;
	
	private static String _s(final String s) {
		return s==null?"":s;
	}
	

	
	private static String normalize(final String s) {
		if(StringUtil.isBlank(s)) return "";
		return Normalizer.normalize(s, Normalizer.Form.NFD).
				replace(" ", "_").replace("'", "_").
				replaceAll("\\p{InCombiningDiacriticalMarks}+", "").
				trim().
				toUpperCase();
	}
	
	private static void writeColor(final XMLStreamWriter w,final Color c) throws XMLStreamException {
		if(c==null) return;
		w.writeEmptyElement("viz", "color", GexfConstants.XMLNS_VIZ);
		w.writeAttribute("r", String.valueOf(c.getRed()));
		w.writeAttribute("g", String.valueOf(c.getGreen()));
		w.writeAttribute("b", String.valueOf(c.getBlue()));
		if(c.getAlpha()<255) {
			w.writeAttribute("a", String.valueOf(c.getAlpha()/255.0f));
			}
		}
	
	private class Article implements Comparable<Article>
		{
		String pmid=null;
		String ArticleTitle=null;
		String Year=null;
		String ISOAbbreviation=null;
		String doi  = null;
		final Set<String> authors = new HashSet<>();
		boolean singleton_flag = false;//all authors are singletons
		@Override
		public int compareTo(final Article o) {
			return pmid.compareTo(o.pmid);
			}
		
		void gexf(final XMLStreamWriter w) throws XMLStreamException {
			
			w.writeStartElement("node");
			w.writeAttribute("id","pmid:"+this.pmid);
			w.writeAttribute("label", String.valueOf(this.ArticleTitle));
			
			writeColor(w,PubmedAuthorGraph.this.articleColor);

			w.writeEmptyElement("viz", "shape", GexfConstants.XMLNS_VIZ);
			w.writeAttribute("value","square");

			if(PubmedAuthorGraph.this.scale_articles)
				{
	
				w.writeEmptyElement("viz", "size", GexfConstants.XMLNS_VIZ);
				w.writeAttribute("value",String.valueOf(1+this.authors.size()));
				}
			
			w.writeStartElement("attvalues");
			
			w.writeEmptyElement("attvalue");
				w.writeAttribute("for","nodetype");
				w.writeAttribute("value","article");

			
			w.writeEmptyElement("attvalue");
				w.writeAttribute("for", "pmid");
				w.writeAttribute("value", this.pmid);
			w.writeEmptyElement("attvalue");
				w.writeAttribute("for", "title");
				w.writeAttribute("value", this.ArticleTitle);
			w.writeEmptyElement("attvalue");
				w.writeAttribute("for", "doi");
				w.writeAttribute("value", this.doi);
			w.writeEmptyElement("attvalue");
				w.writeAttribute("for", "year");
				w.writeAttribute("value", this.Year);
			
			try {
				final int year = Integer.parseInt(this.Year);
				w.writeEmptyElement("attvalue");
				w.writeAttribute("for", "year_int");
				w.writeAttribute("value",String.valueOf(year));
				}
			catch(NumberFormatException err) 
				{
				
				}
				
				
			w.writeEmptyElement("attvalue");
				w.writeAttribute("for", "journal");
				w.writeAttribute("value", this.ISOAbbreviation);
			w.writeEmptyElement("attvalue");
				w.writeAttribute("for","count_authors");
				w.writeAttribute("value",String.valueOf(this.authors.size()));

			
			w.writeEndElement();//attvalues
			w.writeEndElement();//node
		}

		

		}
	private class ArticleBinding extends TupleBinding<Article>
		{
		@Override
		public Article entryToObject(final TupleInput in) {
			final Article a = new Article();
			a.pmid = in.readString();
			a.ArticleTitle= in.readString();
			a.Year = in.readString();
			a.ISOAbbreviation = in.readString();
			a.doi = in.readString();
			
			final int n_authors =in.readInt();
			for(int i=0;i<n_authors;++i){
				a.authors.add(in.readString());
			}
			a.singleton_flag = in.readBoolean();
			return a;
			}
		@Override
		public void objectToEntry(final Article a, final TupleOutput w) {
			w.writeString(_s(a.pmid));
			w.writeString(_s(a.ArticleTitle));
			w.writeString(_s(a.Year));
			w.writeString(_s(a.ISOAbbreviation));
			w.writeString(_s(a.doi));
			w.writeInt(a.authors.size());
			for(final String auth:a.authors)
				{
				w.writeString(auth);
				}
			w.writeBoolean(a.singleton_flag);
			}
		}
	private final ArticleBinding articleBinding = new ArticleBinding();

	

	private  class Author implements Comparable<Author>
		{
		String id= null;
		String foreName = null;
		String lastName = null ;
		String initials=null;
		String affiliation =null;
		final Set<String> pmids = new HashSet<>();
		String male = null;
		String female = null;
		String place = null;
		
		/** idea: in xml+xslt , used to presort orcid indentifiers
		 * in order to produce unique pairs of collab(orcid1,orcid2) */
		@Override
		public int compareTo(final Author o) {
			return this.id.compareTo(o.id);
			}
		
		
		void gexf(final XMLStreamWriter w) throws XMLStreamException {
			
			w.writeStartElement("node");
			w.writeAttribute("id", this.id);
			if(use_initials_to_build_sample_id) {
				w.writeAttribute("label", String.join(" ",lastName,initials));
			}
			else {
			w.writeAttribute("label", String.join(" ",lastName,foreName));
			}

			writeColor(w,PubmedAuthorGraph.this.authorColor);
			
			if(PubmedAuthorGraph.this.scale_authors)
				{
				w.writeEmptyElement("viz", "size", GexfConstants.XMLNS_VIZ);
				w.writeAttribute("value",String.valueOf(1+this.pmids.size()));
				}
			w.writeEmptyElement("viz", "shape", GexfConstants.XMLNS_VIZ);
			w.writeAttribute("value","disc");

			
			
			w.writeStartElement("attvalues");
			
			
			w.writeEmptyElement("attvalue");
			w.writeAttribute("for","nodetype");
			w.writeAttribute("value","author");
			

			if(!StringUtil.isBlank(this.foreName)) {
				w.writeEmptyElement("attvalue");
				w.writeAttribute("for","foreName");
				w.writeAttribute("value",this.foreName);
				}
			
			w.writeEmptyElement("attvalue");
			w.writeAttribute("for","lastName");
			w.writeAttribute("value",this.lastName);

			if(!StringUtil.isBlank(this.initials)) {
				w.writeEmptyElement("attvalue");
				w.writeAttribute("for","initials");
				w.writeAttribute("value",this.initials);
				}

			if(!StringUtil.isBlank(this.affiliation)) {
				w.writeEmptyElement("attvalue");
				w.writeAttribute("for","affiliation");
				w.writeAttribute("value",this.affiliation);
				}

			w.writeEmptyElement("attvalue");
			w.writeAttribute("for","count_articles");
			w.writeAttribute("value",String.valueOf(this.pmids.size()));

			if(!StringUtil.isBlank(this.male)) {
				w.writeEmptyElement("attvalue");
				w.writeAttribute("for","male");
				w.writeAttribute("value",this.male);
				}
			
			if(!StringUtil.isBlank(this.female)) {
				w.writeEmptyElement("attvalue");
				w.writeAttribute("for","female");
				w.writeAttribute("value",this.female);
				}
			
			if(!StringUtil.isBlank(this.place)) {
				w.writeEmptyElement("attvalue");
				w.writeAttribute("for","place");
				w.writeAttribute("value",this.place);
				}
			
			w.writeEndElement();//attvalues
			w.writeEndElement();//node
		}
		
		}	
	private  class AuthorBinding extends TupleBinding<Author>
		{
		@Override
		public Author entryToObject(final TupleInput in) {
			final Author a=new Author();
			a.foreName = in.readString();
			a.lastName = in.readString();
			a.id = in.readString();
			a.initials = in.readString();
			a.affiliation = in.readString();
			// pubmed gender
			a.male = in.readString();
			a.female = in.readString();
			// pubmed map
			a.place = in.readString();
			
			int n_pmid=in.readInt();
			for(int i=0;i<n_pmid;++i){
				a.pmids.add(in.readString());
			}
			return a;
			}
		@Override
		public void objectToEntry(final Author a, TupleOutput w) {
			w.writeString(_s(a.foreName));
			w.writeString(_s(a.lastName));
			w.writeString(_s(a.id));
			w.writeString(_s(a.initials));
			w.writeString(_s(a.affiliation));
			// pubmed gender
			w.writeString(_s(a.male));
			w.writeString(_s(a.female));
			// pubmed map
			w.writeString(_s(a.place));

			
			w.writeInt(a.pmids.size());
			for(final String pmid:a.pmids)
				{
				w.writeString(_s(pmid));
				}
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
	
	private Author parseAuthor(final XMLEventReader r,final StartElement root,final String pmid)  throws XMLStreamException
		{
		final Author au = new Author();
		
		Attribute att = root.getAttributeByName(PubmedGender.MALE_QNAME);
		if(att!=null) au.male=att.getValue();
		att = root.getAttributeByName(PubmedGender.FEMALE_QNAME);
		if(att!=null) au.female=att.getValue();

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
					att = root.getAttributeByName(PubmedMap.QNAME_PLACE_ATTRIBUTE);
					if(att!=null) {
						au.place = att.getValue();
						}
					au.affiliation=r.getElementText().trim();
				} else if(eltName.equals("Initials")) {
					au.initials=r.getElementText().trim();
				}
			}
			else if(evt.isEndElement() &&
				evt.asEndElement().getName().getLocalPart().equals("Author")) {
				if(StringUtil.isBlank(au.lastName)) return null;
				
				au.id= normalize(au.lastName) + "~" ;
				
				if(use_initials_to_build_sample_id)
					{
					if(StringUtil.isBlank(au.initials)) return null;
					au.id += normalize(au.initials);
					}
				else
					{
					if(StringUtil.isBlank(au.foreName)) return null;
					au.id += normalize(au.foreName);
					}
				return au;
			}
		}
		throw new IllegalStateException("should never happen");
		}
	
	private void scanArticles(final InputStream in) throws XMLStreamException,IOException
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
						scanArticle(localName,r);
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
	
	
	private void scanArticle(
			final String rootName,
			final XMLEventReader r
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
						final Author author = parseAuthor(r,start,article.pmid);
						if(author!=null) {
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
		article.authors.addAll(authors.stream().map(A->A.id).collect(Collectors.toSet()));
		
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();

		
		StringBinding.stringToEntry(article.pmid, key);
		if(this.articleDatabase.get(txn, key, data, LockMode.DEFAULT)!=OperationStatus.NOTFOUND)
			{
			LOG.debug("Article already in database : "+article.pmid);
			return;
			}

		
		Collections.sort(authors);

		// inserting authors
		for(int x=0;x< authors.size();++x)
			{
			Author au = authors.get(x);
			StringBinding.stringToEntry(au.id, key);
			// replace with old value
			if(this.authorDatabase.get(txn, key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
				{
				//LOG.debug("Author already in database : "+au.id);
				au = authorBinding.entryToObject(data);
				}
			//add pmid
			au.pmids.add(article.pmid);
			
			this.authorBinding.objectToEntry(au, data);
			if(this.authorDatabase.put(txn, key, data)!=OperationStatus.SUCCESS) {
				throw new JvarkitException.BerkeleyDbError("Cannot update author");
				}		
			}
		
		StringBinding.stringToEntry(article.pmid, key);
		this.articleBinding.objectToEntry(article, data);
		if(this.articleDatabase.put(txn, key, data)!=OperationStatus.SUCCESS) 
			{
			throw new JvarkitException.BerkeleyDbError("Cannot put in article db");
			}
		
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
		Cursor cursor = null;
		try {
			pw = openFileOrStdoutAsPrintWriter(this.outputFile);
			w = xof.createXMLStreamWriter(pw);
			
			w.writeStartDocument("UTF-8","1.0");
			w.writeStartElement("gexf");
			w.writeAttribute("xmlns",GexfConstants.XMLNS);
			w.writeAttribute("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance");
			w.writeAttribute("xmlns:viz",GexfConstants.XMLNS_VIZ);
			w.writeAttribute("xsi:schemaLocation",GexfConstants.XSI_SCHEMA_LOCATION);
			w.writeAttribute("version", GexfConstants.VERSION);
			
			
			/* meta */
			w.writeStartElement("meta");
				w.writeAttribute("lastmodifieddate",new SimpleDateFormat("yyyy-MM-dd").format(new Date()));
				w.writeStartElement("creator");
				  w.writeCharacters(PubmedAuthorGraph.class.getSimpleName());
				w.writeEndElement();
				w.writeStartElement("description");
				  w.writeCharacters(PubmedAuthorGraph.class.getSimpleName());
				w.writeEndElement();
			w.writeEndElement();
			
			/* graph */
			w.writeStartElement("graph");
			w.writeAttribute("mode", "static");
			w.writeAttribute("defaultedgetype", "directed");
	
			
			/* attributes */
			w.writeStartElement("attributes");
			w.writeAttribute("class","node");
			w.writeAttribute("mode","static");
			
			gexfAttDecl(w,"nodetype","string");
			gexfAttDecl(w,"foreName","string");
			gexfAttDecl(w,"lastName","string");
			gexfAttDecl(w,"initials","string");
			gexfAttDecl(w,"affiliation","string");
			gexfAttDecl(w,"count_articles","integer");
			//from pubmed gender
			w.writeComment("Attribute filled if used with jvarkit/pubmedgender ");
			gexfAttDecl(w,"male","integer");
			w.writeComment("Attribute filled if used with jvarkit/pubmedgender ");
			gexfAttDecl(w,"female","integer");
			//from pubmed map
			w.writeComment("Attribute filled if used with jvarkit/pubmedmap ");
			gexfAttDecl(w,"place","string");

			
			
			
			gexfAttDecl(w,"pmid","string");
			gexfAttDecl(w,"title","string");
			gexfAttDecl(w,"doi","string");
			gexfAttDecl(w,"year","string");
			gexfAttDecl(w,"year_int","integer");
			gexfAttDecl(w,"journal","string");
			gexfAttDecl(w,"count_authors","integer");

			
			w.writeEndElement();//attributes
			
			
			/* nodes */
			w.writeStartElement("nodes");
			cursor  = this.authorDatabase.openCursor(txn, null);
			while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
				{
				final Author au=this.authorBinding.entryToObject(data);
				if(this.remove_singleton_authors && au.pmids.size()<2)
					{
					continue;
					}
				au.gexf(w);
				}
			cursor.close();
			
			cursor  = this.articleDatabase.openCursor(txn, null);
			while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
				{
				final Article article =this.articleBinding.entryToObject(data);
				if(this.remove_singleton_authors )
					{
					boolean ok=false;
					for(final String author_id: article.authors) {
						DatabaseEntry key2 = new DatabaseEntry();
						DatabaseEntry data2 = new DatabaseEntry();
						StringBinding.stringToEntry(author_id, key2);
						if(this.authorDatabase.get(txn, key2, data2, LockMode.DEFAULT)!=OperationStatus.SUCCESS)
							{
							LOG.warn("cannot get "+author_id);
							continue;
							}
						final Author author = authorBinding.entryToObject(data2);
						if(author.pmids.size()<2) continue;
						ok=true;
						break;
						}
					if(!ok) {
						article.singleton_flag = true;
						articleBinding.objectToEntry(article, data);
						if(cursor.putCurrent(data)!=OperationStatus.SUCCESS)
							{
							LOG.warn("cannot update "+article.pmid);
							continue;
							}
						continue;
						}
					}
				article.gexf(w);
				}
			cursor.close();
			
			w.writeEndElement();//nodes
			
			w.writeStartElement("edges");
			key=new DatabaseEntry();
			cursor  = this.articleDatabase.openCursor(this.txn, null);
			while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
				{
				final Article article =this.articleBinding.entryToObject(data);
				if(article.singleton_flag) {
					w.writeComment("ignoring singleton article PMID: "+article.pmid);
					LOG.debug("ignoring singleton "+article.pmid);
					continue;
				}
				
				
				for(final String author_id: article.authors)
					{
					final DatabaseEntry key2 = new DatabaseEntry();
					final DatabaseEntry data2 = new DatabaseEntry();
					StringBinding.stringToEntry(author_id,key2);
					if(this.authorDatabase.get(this.txn, key2, data2, LockMode.DEFAULT)!=OperationStatus.SUCCESS) {
						LOG.warn("cannot get "+author_id);
						continue;
						}
					final Author au=this.authorBinding.entryToObject(data2);
					if(this.remove_singleton_authors && au.pmids.size()<2)
						{
						continue;
						}
					w.writeStartElement("edge");
					w.writeAttribute("id", "E"+(++ID_GENERATOR));
					w.writeAttribute("type", "directed");
					w.writeAttribute("source",author_id);
					w.writeAttribute("target","pmid:"+article.pmid);
					w.writeEndElement();//edge
					}
				//
				}
			cursor.close();
			
			w.writeEndElement();//edges

			
			w.writeEndElement();
			w.writeEndDocument();
			w.flush();
			pw.flush();
			pw.close();pw=null;
			}
		catch(final Exception err) {
			err.printStackTrace();
			LOG.error(err);
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
			
			
			/* input is a efetch stream */
			final String inputName=oneFileOrNull(args);
			in=(inputName==null?stdin():IOUtils.openURIForReading(inputName));
			scanArticles(in);
			in.close();in=null;
	
			dumpGexf();			
			return 0;
		} catch (final Throwable e) {
			e.printStackTrace();
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(in);
			CloserUtil.close(r);
			CloserUtil.close(authorDatabase);
			CloserUtil.close(articleDatabase);
			CloserUtil.close(environment);
			}
		}

	
	public static void main(final String[] args)
		{
		new PubmedAuthorGraph().instanceMainWithExit(args);
		}
	}
