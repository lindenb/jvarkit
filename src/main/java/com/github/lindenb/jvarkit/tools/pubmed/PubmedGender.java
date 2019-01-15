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

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.Collator;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;

/**
 * PubmedGender
 BEGIN_DOC
 
## Building the database

```
$ wget -O jeter.zip "https://www.ssa.gov/oact/babynames/names.zip"
$ unzip -t jeter.zip | tail
    testing: yob2009.txt              OK
    testing: yob2010.txt              OK
    testing: yob2011.txt              OK
    testing: yob2012.txt              OK
    testing: yob2013.txt              OK
    testing: yob2014.txt              OK
    testing: yob2015.txt              OK
    testing: yob1880.txt              OK
    testing: NationalReadMe.pdf       OK
No errors detected in compressed data of jeter.zip.
$ unzip -p jeter.zip yob2015.txt &gt; database.csv
```

## Example

```
$ java -jar dist/pubmeddump.jar "Lindenbaum[Author] Nantes" 2> /dev/null  | java -jar dist/pubmedgender.jar  -d jeter.csv 2> /dev/null | grep Lindenbaum -A 2 -B 1
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
--
                <Author ValidYN="Y" male="169">
                    <LastName>Lindenbaum</LastName>
                    <ForeName>Pierre</ForeName>
                    <Initials>P</Initials>
```

## Example


Get an histogram for gender ratio in pubmed:

```
$ java -jar dist/pubmeddump.jar  '("bioinformatics"[Journal]) AND ("2018-01-01"[Date - Publication] : "2018-07-06"[Date - Publication]) ' |\
	java -jar dist/pubmedgender.jar -d ./yob2017.txt  |\
	java -jar dist/xsltstream.jar -t ~/tr.xsl -n "PubmedArticle" |\
	tr ";" "\n" | sort | uniq -c |\
	java -jar dist/simpleplot.jar -su -t SIMPLE_HISTOGRAM  --title "Authors in Bioinformatics 2018"
```

with tr.xsl:

```xslt
<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet 
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:j='http://www.ibm.com/xmlns/prod/2009/jsonx'
	version='1.0'>
<xsl:output method="text"/>


<xsl:template match="/">
<xsl:apply-templates select="PubmedArticle"/>
</xsl:template>

<xsl:template match="PubmedArticle">
<xsl:apply-templates select="MedlineCitation/Article/AuthorList/Author[@female or @male]"/>
</xsl:template>

<xsl:template match="Author">
<xsl:choose>
	<xsl:when test="@female and @male and number(@female) &gt; number(@male) ">
		<xsl:text>FEMALE;</xsl:text>
	</xsl:when>
	<xsl:when test="@female and @male and number(@female) &lt; number(@male) ">
		<xsl:text>MALE;</xsl:text>
	</xsl:when>
	<xsl:when test="@female">
		<xsl:text>FEMALE;</xsl:text>
	</xsl:when>
	<xsl:when test="@male">
		<xsl:text>MALE;</xsl:text>
	</xsl:when>
</xsl:choose>
</xsl:template>

</xsl:stylesheet>
```


## See also


 * A Simple tool to get the sex ratio in pubmed :  http://plindenbaum.blogspot.fr/2010/09/simple-tool-to-get-sex-ratio-in-pubmed.html


 
 END_DOC
 */
@Program(name="pubmedgender",
	keywords={"pubmed","gender","ncbi","xml"},
	description="Add gender-related attributes in the Author tag of pubmed xml. ")
public class PubmedGender
	extends Launcher
	{
	private static final Logger LOG = Logger.build(PubmedGender.class).make();

	@Parameter(names={"-d","--database"},description="REQUIRED: A comma delimited file containing the following columns: 1) Name 2) sex (M/F) 3) Score. See http://cpansearch.perl.org/src/EDALY/Text-GenderFromName-0.33/GenderFromName.pm or https://www.ssa.gov/oact/babynames/names.zip",required=true)
	private File dataFile = null;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	/* used by pubmed graph */
	static final QName MALE_QNAME =new QName("male");
	static final QName FEMALE_QNAME =new QName("female");

	
	private static class GenderInfo {
		String male=null;
		String female=null;
	}
	private final Map<String, GenderInfo> name2gender;

	public PubmedGender()
		{
		final Collator collator= Collator.getInstance(Locale.US);
		collator.setStrength(Collator.PRIMARY);
		this.name2gender=new TreeMap<String, GenderInfo>(collator);
		}
	
	@Override
	public int doWork(final List<String> args) {
		final String inputName= oneFileOrNull(args);
		if(this.dataFile==null || !this.dataFile.exists()) {
			LOG.error("Undefined option -d");
			return -1;
		}
		OutputStream out=null;
		XMLEventReader r=null;
		InputStream in=null;
		XMLEventWriter w=null;
		BufferedReader br=null;
		try {
			final Pattern comma=Pattern.compile("[,]");
			LOG.info("load "+this.dataFile);
			this.name2gender.clear();
			br = IOUtils.openFileForBufferedReading(this.dataFile);
			String line;
			while((line=br.readLine())!=null) {
				if(line.startsWith("#")) continue;
				final String tokens[]= comma.split(line);
				if(tokens.length!=3) throw new JvarkitException.UserError("expected 3 comma-separated columns in "+line);
				tokens[0]=tokens[0].toLowerCase();
				GenderInfo gi = this.name2gender.get(tokens[0]);
				if(gi==null) {
					gi= new GenderInfo();
					this.name2gender.put(tokens[0],gi);
				}
				if(tokens[1].equals("F")) {
					gi.female=tokens[2];
				} else if(tokens[1].equals("M")) {
					gi.male=tokens[2];
				} else
				{
					throw new JvarkitException.UserError("expected 'M' or 'F' in 2nd column  in "+line);
				}
				
			}
			br.close();
			LOG.info("database contains "+this.name2gender.size());
			
			
			final XMLEventFactory xmlEventFactory = XMLEventFactory.newFactory();
			final XMLInputFactory xmlInputFactory = XMLInputFactory.newFactory();
			xmlInputFactory.setXMLResolver(new XMLResolver() {
				@Override
				public Object resolveEntity(String publicID, String systemID, String baseURI, String namespace)
						throws XMLStreamException {
					LOG.debug("Ignoring resolve Entity");
					return new ByteArrayInputStream(new byte[0]);
				}
			});
			in=(inputName==null?stdin():IOUtils.openURIForReading(inputName));
			r = xmlInputFactory.createXMLEventReader(in);
			
			out = super.openFileOrStdoutAsStream(this.outputFile);
			final XMLOutputFactory xof = XMLOutputFactory.newFactory();
			w=xof.createXMLEventWriter(out, "UTF-8");
			while(r.hasNext()) {
			final XMLEvent evt = r.nextEvent();
			if(evt.isStartElement() &&
					evt.asStartElement().getName().getLocalPart().equals("Author"))
				{
				String initials=null;
				String firstName =null;
				final List<XMLEvent> events = new ArrayList<>();
				while(r.hasNext()) {
					final XMLEvent evt2 = r.nextEvent();
					events.add(evt2);
					if(evt2.isStartElement()) {
						final String eltName = evt2.asStartElement().getName().getLocalPart();
						if((eltName.equals("ForeName") || eltName.equals("FirstName")) && r.peek().isCharacters()) {
							final XMLEvent textEvt = r.nextEvent();
							events.add(textEvt);
							firstName= textEvt.asCharacters().getData();	
							}
						else if(eltName.equals("Initials") && r.peek().isCharacters()) {
							final XMLEvent textEvt = r.nextEvent();
							events.add(textEvt);
							initials= textEvt.asCharacters().getData();	
							}
						}
					else if(evt2.isEndElement() && evt2.asEndElement().getName().getLocalPart().equals("Author")) {
						break;
					} 
				}
				
				
					
				if( firstName!=null) {
					final String tokens[]=firstName.toLowerCase().split("[ ,]+");
					firstName="";
					for(final String s:tokens)
						{
						if(s.length()> firstName.length())
							{
							firstName=s;
							}
						}
					
					if(	firstName.length()==1 ||
						(initials!=null && firstName.equals(initials.toLowerCase())))
						{
						firstName=null;
						}
					}
				
				String female=null;
				String male=null;

				if(firstName!=null) {
					final GenderInfo gi = this.name2gender.get(firstName);
					if(gi!=null) {
						male=gi.male;
						female=gi.female;
						}
					}
				
				
				final List<Attribute> attributes = new ArrayList<>();
				Iterator<?> t=evt.asStartElement().getAttributes();
				while(t.hasNext()) {
				    final Attribute att =  (Attribute)t.next();
				    if(att.getName().equals(MALE_QNAME)) continue;
				    if(att.getName().equals(FEMALE_QNAME)) continue;
					attributes.add(att);
				}
				
				if(male!=null) attributes.add(xmlEventFactory.createAttribute(MALE_QNAME, male));
				if(female!=null) attributes.add(xmlEventFactory.createAttribute(FEMALE_QNAME, female));
				
				w.add(xmlEventFactory.createStartElement(
						evt.asStartElement().getName(),
						attributes.iterator(),
						evt.asStartElement().getNamespaces()
						));
				
				
				
				for(final XMLEvent evt2:events) w.add(evt2);
				continue;
				}
			w.add(evt);
			}
			
			r.close();r=null;
			in.close();in=null;
			w.flush();w.close();w=null;
			out.flush();out.close();out=null;
			return 0;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally {
			CloserUtil.close(r);
			CloserUtil.close(in);
			CloserUtil.close(w);
			CloserUtil.close(out);
			CloserUtil.close(br);
			this.name2gender.clear();
		}

		}
		
	
	public static void main(String[] args) {
		new PubmedGender().instanceMainWithExit(args);
	}
}
