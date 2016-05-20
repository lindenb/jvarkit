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
import htsjdk.samtools.util.RuntimeIOException;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.Writer;
import java.net.URL;
import java.net.URLEncoder;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Properties;
import java.util.Set;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;


public class PubmedOrcidGraph
	extends AbstractPubmedOrcidGraph
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(PubmedOrcidGraph.class);

	private XMLInputFactory xmlInputFactory=null;
	private Connection conn=null;

	private File getDerbyDirectory() {
		return new File(this.derbyFilePath);
	}

	
	private void openDerby() {
		try {
			boolean create;
			final Properties props = new Properties();
			final File derbyDir = getDerbyDirectory();
			LOG.info("open derby :" + getDerbyDirectory());
			if(derbyDir.exists()) {
				if(!derbyDir.isDirectory()) {
					throw new RuntimeIOException("derby database is not a directory : "+derbyDir);
					}
				create=false;
				}
			else
				{
				create=true;
				}
			props.setProperty("create", String.valueOf(create));
			this.conn = DriverManager.getConnection("jdbc:derby:"+derbyDir,props);
			
			if(create) {
				final Statement stmt= this.conn.createStatement();
				final String sqls[]={
						"CREATE TABLE AUTHOR(ORCID CHAR(19) NOT NULL,LASTNAME VARCHAR(50),FORENAME VARCHAR(50),AFFILIATION VARCHAR(100),PROCESSED INT,UNIQUE (ORCID))",
						"CREATE TABLE ARTICLE(PMID INTEGER NOT NULL,UNIQUE(PMID),PUBYEAR VARCHAR(10),JOURNAL VARCHAR(50),TITLE VARCHAR(200))",
						"CREATE TABLE ARTICLE2AUTHOR(ORCID CHAR(19) NOT NULL,PMID INTEGER NOT NULL,UNIQUE(ORCID,PMID))",
						"CREATE TABLE AUTHOR2AUTHOR(ORCID1 CHAR(19) NOT NULL,ORCID2 CHAR(19) NOT NULL,UNIQUE(ORCID1,ORCID2))"
						};
				for(final String sql:sqls) {
					LOG.warn(sql);
					stmt.execute(sql);
				}
				stmt.close();
			}
			this.conn.setAutoCommit(true);
		} catch (Exception e) {
			CloserUtil.close(this.conn);
			this.conn = null;
			throw new RuntimeException(e);
		}
	}
	
	private static String truncate(String s,int maxLen) {
		return s.length()>maxLen?s.substring(0,maxLen-3)+"...":s;
	}
	
	private void closeDerby() {
		CloserUtil.close(this.conn);
		this.conn = null;
		try {
			final Properties props = new Properties();
			props.setProperty("shutdown", "true");
			final File derbyDir = getDerbyDirectory();
			DriverManager.getConnection("jdbc:derby:"+derbyDir,props);
		} catch (Exception e) {
			
			}
		}

	
	/*
	private void eSummary(Article a) throws IOException,XMLStreamException
		{
		final QName attName=new QName("Name");
		final QName attType=new QName("Type");
		String url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?retmode=xml&db=pubmed&" +
				"id="+a.pmid
				;
		info(url);
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
	*/
	
	
	private Collection<Throwable> gexf() throws XMLStreamException,IOException
			{
			PreparedStatement stmt=null;
			ResultSet row=null;
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			XMLStreamWriter w= null;
			Writer fw=null;
			try {
				fw = super.openFileOrStdoutAsPrintWriter();
				w=xof.createXMLStreamWriter(fw);
				w.writeStartDocument("UTF-8", "1.0");
			
			w.writeStartElement("gexf");
			w.writeDefaultNamespace("http://www.gexf.net/1.2draft");
			w.writeAttribute("version", "1.2");
			w.writeStartElement("meta");
			  w.writeStartElement("creator");
			  w.writeCharacters(getName()+" by "+getAuthorName());
			  w.writeEndElement();
			 
			  w.writeStartElement("description");
			  w.writeCharacters(getProgramCommandLine());
			  w.writeEndElement();
			
			w.writeEndElement();//meta
			  
			  w.writeStartElement("graph");
			  w.writeAttribute("mode", "static");
			  w.writeAttribute("defaultedgetype", "undirected");
			  
			  w.writeStartElement("attributes");
			  w.writeAttribute("class", "edge");
			  w.writeAttribute("mode", "static");
			  w.writeEndElement();//attributes
				
			  w.writeStartElement("attributes");                                                                                     
			  w.writeAttribute("class", "edge");
			  w.writeAttribute("mode", "static");
				  
	          w.writeEmptyElement("attribute");
				w.writeAttribute("id", "0");
				w.writeAttribute("title", "orcid");
				w.writeAttribute("type", "string");
		      w.writeEmptyElement("attribute");
				w.writeAttribute("id", "1");
				w.writeAttribute("title", "lastName");
				w.writeAttribute("type", "string");
			  w.writeEmptyElement("attribute");
				w.writeAttribute("id", "2");
				w.writeAttribute("title", "foreName");
				w.writeAttribute("type", "string");
			  w.writeEmptyElement("attribute");
				w.writeAttribute("id", "3");
				w.writeAttribute("title", "affiliation");
				w.writeAttribute("type", "string");
	          w.writeEndElement();//attributes
			  
	
			  
			  w.writeStartElement("nodes");
			  stmt = conn.prepareStatement("SELECT ORCID,FORENAME,LASTNAME,AFFILIATION FROM AUTHOR");
			  row= stmt.executeQuery();
			  while(row.next())
			 	{
				final String orcid= row.getString(1);
				w.writeStartElement("node");
				w.writeAttribute("id",orcid);
				String forename = row.getString(2);
				forename=(forename==null?"":forename);
				String lastname = row.getString(3);
				lastname=(lastname==null?"":lastname);
				String affiliation = row.getString(4);
				affiliation=(affiliation==null?"":affiliation);
				
				w.writeAttribute("label",lastname+" "+forename);
					
				w.writeStartElement("attvalues");
				
				  w.writeEmptyElement("attvalue");
					w.writeAttribute("for", "0");
					w.writeAttribute("value",orcid);

				  w.writeEmptyElement("attvalue");
					w.writeAttribute("for", "1");
					w.writeAttribute("value",lastname);
				
				  w.writeEmptyElement("attvalue");
					w.writeAttribute("for", "2");
					w.writeAttribute("value",forename);

				  w.writeEmptyElement("attvalue");
					w.writeAttribute("for", "3");
					w.writeAttribute("value",affiliation);
	
					
					
				w.writeEndElement();//attvalues
				w.writeEndElement();//node
			 	}
			  w.writeEndElement();//nodes
			  row.close();
			  stmt.close();
			
			  stmt = conn.prepareStatement("SELECT ORCID1,ORCID2 FROM AUTHOR2AUTHOR");
			  row= stmt.executeQuery();

			  w.writeStartElement("edges");
			  long id=0;
			  while(row.next())
			 	{
				w.writeStartElement("edge");
				w.writeAttribute("id","E"+(++id));
				w.writeAttribute("type","undirected");
				w.writeAttribute("source",row.getString(1));
				w.writeAttribute("target",row.getString(2));
				w.writeEndElement();
			 	}
			  w.writeEndElement();//edges
			  row.close();
			  stmt.close();

			  
			  w.writeEndElement();//graph
	
			
			w.writeEndElement();//gexf
			
			w.writeEndDocument();
			w.flush();
			w.close();
			return RETURN_OK;
			} catch (Exception e) {
				return wrapException(e);
			} finally
				{
				CloserUtil.close(w);
				}
			}
	private static class Author
		{
		String foreName = null;
		String lastName = null ;
		String orcid = null;
		String affiliation =null;
		}
	
	private Author parseAuthor(XMLEventReader r)  throws XMLStreamException
		{
		final Author au = new Author();
		while(r.hasNext()) {
			final XMLEvent evt=r.nextEvent();
			if(evt.isStartElement()) {
				final StartElement start = evt.asStartElement();
				String eltName = start.getName().getLocalPart();
				if(eltName.equals("LastName")) {
					au.lastName=r.getElementText().trim();
				} else if(eltName.equals("ForeName")) {
					au.foreName=r.getElementText().trim();
				} else if(eltName.equals("Affiliation")) {
					au.affiliation=r.getElementText().trim();
				} else if(eltName.equals("Identifier")) {
					final Attribute source= start.getAttributeByName(new QName("Source"));
					if(source!=null && source.getValue().equalsIgnoreCase("ORCID")) {
						au.orcid=r.getElementText().trim();
						final int slash = au.orcid.lastIndexOf('/');
						if(slash!=-1) au.orcid=au.orcid.substring(slash+1);
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
	
	
	private void eSearchOrcid(final String orcid,int depth) throws XMLStreamException,IOException,SQLException {
		InputStream in=null;
		XMLEventReader r=null;
		try {	
			Set<String> pmids=new HashSet<>();
			final String url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=100&term="+
						URLEncoder.encode(orcid+"[AUID]","UTF-8");
			LOG.info("Searching:"+url);
			in = new URL(url).openStream();
			r = xmlInputFactory.createXMLEventReader(in);
			while(r.hasNext()) {
				final XMLEvent evt=r.nextEvent();
			if(evt.isStartElement()) {
				final StartElement start = evt.asStartElement();
				String eltName = start.getName().getLocalPart();
				if(eltName.equals("Id")) {
					pmids.add(r.getElementText().trim());
				}
				
				}
			}//end of xml read
			r.close();
			in.close();

			for(final String pmid : pmids) {
				fetchArticle(Long.parseLong(pmid), depth+1);
			}
		} finally {
			CloserUtil.close(in);
			CloserUtil.close(r);
		}
	}

	
	private void scanUnprocessedOrcids(int depth) throws XMLStreamException,IOException,SQLException {
		PreparedStatement stmt1=null;
		ResultSet row=null;
		try {
			for(;;){
				String orcid=null;
				stmt1 = this.conn.prepareStatement("SELECT ORCID FROM AUTHOR WHERE PROCESSED=0");
				row = stmt1.executeQuery();
				while(row.next()) {
					orcid =row.getString(1);
					break;
				}
				row.close();
				stmt1.close();
				
				if(orcid==null) break;
				LOG.info("orcid to be processed :"+orcid);
				stmt1 = this.conn.prepareStatement("UPDATE AUTHOR SET PROCESSED=1 WHERE ORCID=?");
				stmt1.setString(1, orcid);
				if(stmt1.executeUpdate()!=1) {
					throw new SQLException("Cannot update into UPPDATE");
				}
				stmt1.close();
				eSearchOrcid(orcid,depth);
			}

		} finally {
			CloserUtil.close(row);
			CloserUtil.close(stmt1);
		}
	}
	
	private void fetchArticle(long pmid,int depth) throws XMLStreamException,IOException,SQLException {
		InputStream in=null;
		PreparedStatement stmt1=null;
		ResultSet row=null;
		XMLEventReader r=null;
		try {
			int count=0;
			stmt1 = this.conn.prepareStatement("SELECT COUNT(*) FROM ARTICLE WHERE PMID=?");
			stmt1.setLong(1, pmid);
			row=stmt1.executeQuery();
			while(row.next()) {
				count = row.getInt(1);
			}
			row.close();row=null;
			stmt1.close();stmt1=null;
			if(count!=0) return;
			
			final String url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id="+pmid;
			LOG.info("Downloading pmid:"+pmid);
			in = new URL(url).openStream();
			r = xmlInputFactory.createXMLEventReader(in);
			String ArticleTitle=null;
			String Year=null;
			String ISOAbbreviation=null;
			boolean PubDate=false;
			List<Author> authors = new ArrayList<>();
			while(r.hasNext()) {
				final XMLEvent evt=r.nextEvent();
			if(evt.isStartElement()) {
				final StartElement start = evt.asStartElement();
				String eltName = start.getName().getLocalPart();
				
				if(ISOAbbreviation==null && eltName.equals("ISOAbbreviation")) {
					ISOAbbreviation=r.getElementText();
				}
				else if(ArticleTitle==null && eltName.equals("ArticleTitle")) {
					ArticleTitle=r.getElementText();
				}
				else if(eltName.equals("PubDate")) {
					PubDate=true;
				}
				else if(Year==null && PubDate && eltName.equals("Year")) {
					Year=r.getElementText();
				}
				else if(eltName.equals("Author")) {
					final Author author = parseAuthor(r);
					if(author.orcid!=null) authors.add(author);
				}
				}
			}//end of xml read
			
			LOG.info("insert article");
			stmt1 = this.conn.prepareStatement("INSERT INTO ARTICLE(PMID,PUBYEAR,JOURNAL,TITLE) VALUES(?,?,?,?)");
			stmt1.setLong(1, pmid);
			if(Year==null) {
				stmt1.setNull(2, java.sql.Types.VARCHAR);
			} else  {
				stmt1.setString(2,truncate(Year,10));
			}
			if(ISOAbbreviation==null) {
				stmt1.setNull(3, java.sql.Types.VARCHAR);
			} else  {
				stmt1.setString(3,truncate(ISOAbbreviation,50));
			}
			if(ArticleTitle==null) {
				stmt1.setNull(4, java.sql.Types.VARCHAR);
			} else  {
				stmt1.setString(4,truncate(ArticleTitle,200));
			}
			if(stmt1.executeUpdate()!=1) {
				throw new SQLException("Cannot insert into ARTICLE");
			}
			
			stmt1.close();
			
			LOG.info("insert article-author (N="+authors.size()+")");
			stmt1 = this.conn.prepareStatement("INSERT INTO ARTICLE2AUTHOR(ORCID,PMID) VALUES (?,?)");
			for(final Author au:authors) {
				stmt1.setString(1, au.orcid);
				stmt1.setLong(2, pmid);
				if(stmt1.executeUpdate()!=1) {
					throw new SQLException("Cannot insert into ARTICLE2AUTHOR");
				}
			}
			stmt1.close();
			
			
			LOG.info("insert author_author (N="+authors.size()+")");
			stmt1 = this.conn.prepareStatement("INSERT INTO AUTHOR2AUTHOR(ORCID1,ORCID2) VALUES (?,?)");
			for(int x=0;x< authors.size();++x) {
				for(int y=x+1;y< authors.size();++y) {
					if( authors.get(x).orcid.compareTo( authors.get(y).orcid)<0) {
						stmt1.setString(1, authors.get(x).orcid);
						stmt1.setString(2, authors.get(y).orcid);
					} else {
						stmt1.setString(1, authors.get(y).orcid);
						stmt1.setString(2, authors.get(x).orcid);
					}
					try {stmt1.executeUpdate(); } catch(Exception err2) {LOG.error("duplicate:",err2);}
				}
			}
			stmt1.close();
			
			LOG.info("insert authors (N="+authors.size()+")");

			int idx=0;
			while(idx<authors.size())
				{
				final Author au = authors.get(idx);
				stmt1 = this.conn.prepareStatement("SELECT COUNT(*) FROM AUTHOR WHERE ORCID=?");
				stmt1.setString(1, au.orcid);
				row=stmt1.executeQuery();
				while(row.next()) {
					count = row.getInt(1);
				}
				row.close();row=null;
				stmt1.close();stmt1=null;
				if(count!=0) {
					authors.remove(idx);
					continue;
				}
				
				stmt1 = this.conn.prepareStatement("INSERT INTO AUTHOR(ORCID,LASTNAME,FORENAME,AFFILIATION,PROCESSED) VALUES (?,?,?,?,?)");

				stmt1.setString(1, au.orcid);
				if(au.lastName!=null) {
					stmt1.setString(2, truncate(au.lastName,50));
				} else {
					stmt1.setNull(2, java.sql.Types.VARCHAR);
				}
				if(au.lastName!=null) {
					stmt1.setString(3, truncate(au.foreName,50));
				} else {
					stmt1.setNull(3, java.sql.Types.VARCHAR);
				}
				if(au.affiliation!=null) {
					stmt1.setString(4,truncate( au.affiliation,100));
				} else {
					stmt1.setNull(4, java.sql.Types.VARCHAR);
				}
				stmt1.setInt(5, 0);
				
				if(stmt1.executeUpdate()!=1) {
					throw new SQLException("Cannot insert into AUTHOR");
				}
				stmt1.close();
				++idx;
			}
			
			r.close();
			in.close();
		} finally {
			CloserUtil.close(row);
			CloserUtil.close(stmt1);
			CloserUtil.close(in);
			CloserUtil.close(r);
		}
	}
	

	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.derbyFilePath==null ||super.derbyFilePath.isEmpty())
			{
			return wrapException("Undefined Derby DB Directory: option -"+OPTION_DERBYFILEPATH);
			}
		try {
			Class.forName("org.apache.derby.jdbc.ClientDriver").newInstance();
		} catch(Exception err ){
			LOG.error("Cannot get derby driver");
			return wrapException(err);
		}
		
		try {
			this.xmlInputFactory = XMLInputFactory.newFactory();
			this.xmlInputFactory.setXMLResolver(new XMLResolver() {
				@Override
				public Object resolveEntity(String publicID, String systemID, String baseURI, String namespace)
						throws XMLStreamException {
					LOG.debug("Ignoring resolve Entity");
					return new ByteArrayInputStream(new byte[0]);
				}
			});
		} catch(Exception err) {
			return wrapException(err);
		}
		return super.initializeKnime();
	 	}

	private Collection<Throwable> doScanPmid(){
		final List<String> args = getInputFiles();
		if(args.isEmpty()) {
			LOG.info("No Pmid");
			return RETURN_OK;
		}
		
		try {
			for(final String pmidStr:args){
				final long pmid = Long.parseUnsignedLong(pmidStr);
				fetchArticle(pmid,0);
			}
			scanUnprocessedOrcids(0);
			
			return RETURN_OK;
		} catch (final Exception e) {
			return wrapException(e);
		} finally
		{
			
		}
	}
	

	@Override
	public Collection<Throwable> call() throws Exception {
		try {
			openDerby();
			final String command  = String.valueOf(super.actionStr);
			
			if(command.equals("pmid")) {
				return doScanPmid();
				} 
			else if(command.equals("gexf")) {
				return gexf();
				} 
			else {
				return wrapException("unknown command : "+command);
			}			
		} catch (Exception e) {
			return wrapException(e);
		} finally {
			closeDerby();
		}
		}

	
	public static void main(String[] args)
		{
		new PubmedOrcidGraph().instanceMainWithExit(args);
		}
	}
