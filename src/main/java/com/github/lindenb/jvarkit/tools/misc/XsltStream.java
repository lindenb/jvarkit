package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import javax.xml.XMLConstants;
import javax.xml.namespace.QName;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.Namespace;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.Templates;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;

/**

## Example:

Dumping the Orcid from pubmed:

```
 java  -jar dist/pubmeddump.jar 'orcid[AUID]' |\
 	java -jar dist/xsltstream.jar -t pubmed2orcid.xsl -n "PubmedArticle" 
```

The XSLT stylesheet:

```xml
<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
    xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
    version='1.0'
	>

<xsl:output method="text" />


<xsl:template match="/">
<xsl:apply-templates select="PubmedArticle"/>
</xsl:template>

<xsl:template match="PubmedArticle">
<xsl:apply-templates select="MedlineCitation/Article/AuthorList/Author[Identifier/@Source='ORCID']"/>
</xsl:template>

<xsl:template match="Author">
<xsl:value-of select="LastName"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="ForeName"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="Initials"/>
<xsl:text>	</xsl:text>
<xsl:call-template name="orcid"><xsl:with-param name="s" select="Identifier[@Source='ORCID']"/></xsl:call-template>
<xsl:text>	</xsl:text>
<xsl:for-each select="Affiliation"><xsl:text> </xsl:text></xsl:for-each>
<xsl:text>	</xsl:text>
<xsl:value-of select="../../../PMID[1]"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="../../../DateCreated/Year"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="../../Journal/ISOAbbreviation"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="../../ArticleTitle"/>
<xsl:text>
</xsl:text>
</xsl:template>

<xsl:template name="orcid">
<xsl:param name="s"/>
<xsl:choose>
	<xsl:when test="starts-with($s,'http://orcid.org/')">
		<xsl:call-template name="orcid">
			<xsl:with-param name="s" select="substring($s,18)"/>
		</xsl:call-template>
	</xsl:when>
	<xsl:when test="starts-with($s,'https://orcid.org/')">
		<xsl:call-template name="orcid">
			<xsl:with-param name="s" select="substring($s,19)"/>
		</xsl:call-template>
	</xsl:when>
	<xsl:when test="starts-with($s,'https://')">
		<xsl:call-template name="orcid">
			<xsl:with-param name="s" select="substring($s,9)"/>
		</xsl:call-template>
	</xsl:when>
	<xsl:otherwise>
		<xsl:value-of select="translate($s,'-','')"/>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>

</xsl:stylesheet>
```

output:

```

```


*/
@Program(name="xsltstream",
	description="XSLT transformation for large XML files. xslt is only applied on a given subset of nodes.",
	keywords={"xml","xslt","xsl","stylesheet"}
	)
public class XsltStream extends Launcher {
	private static final Logger LOG = Logger.build(XsltStream.class).make();

	private static class QNameConverter implements IStringConverter<QName>
		{
		@Override
		public QName convert( String s) {
			String ns=null;
			String prefix =null;
			String localName=null;
			if(s.startsWith("{"))
				{
				int par=s.indexOf('}');
				if(par==-1) throw new IllegalArgumentException("closing '}' missing in '"+s+"'.");
				ns=s.substring(1, par);
				s=s.substring(par+1);
				}
			int colon = s.lastIndexOf(':');
			if(colon==-1)
				{
				localName=s;
				if(ns==null || ns.trim().isEmpty()) return new QName(localName);
				return new QName(ns, localName);
				}
			else
				{
				prefix=s.substring(0,colon);
				localName=s.substring(colon+1);
				if(ns==null || ns.trim().isEmpty()) return new QName(XMLConstants.NULL_NS_URI,localName,prefix);
				return new QName(ns, localName,prefix);
				}
			}
		}
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-t","-template"},description="XSLT template file.",required=true)
	private File templateFile = null;
	@Parameter(names={"-n","--tag","--name","-tag","-name"},description="XML node name. name has syntax '{ns}prefix:localName' or 'prefix:localName' or 'localName' or '{ns}localName' ",converter=QNameConverter.class)
	private Set<QName> targetQnames = new HashSet<>();
	@Parameter(names={"-skip","--skip"},description="Ignore those names",converter=QNameConverter.class)
	private Set<QName> skipQNames = new HashSet<>();

	
	private static String toQNAME(final QName q)
		{
		return
			(q.getPrefix()==null|| q.getPrefix().isEmpty()?"":q.getPrefix()+":")+q.getLocalPart();
		}
	
	private void skip(final StartElement startE,final XMLEventReader r)throws XMLStreamException
		{
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isStartElement())
				{
				skip(evt.asStartElement(),r);
				}
			else if(evt.isEndElement())
				{
				return;
				}
			}
			
		}
	
	private Node parseDocument(final Document owner,final StartElement startE,final XMLEventReader r) throws XMLStreamException
		{
		final QName qname=startE.getName();
		final Element root;
		
		if(qname.getNamespaceURI()==null || qname.getNamespaceURI().isEmpty())
			{
			root = owner.createElement(toQNAME(qname));
			}
		else
			{
			root = owner.createElementNS(
					qname.getNamespaceURI(),
					toQNAME(qname)
					);
			}
		Iterator<?> it= startE.getNamespaces();
		while(it.hasNext())
			{
			final Namespace att = Namespace.class.cast(it.next());
			root.setAttribute("xmlns:"+att.getPrefix(),att.getNamespaceURI());
			}
		
		it= startE.getAttributes();
		while(it.hasNext())
			{
			final Attribute att = Attribute.class.cast(it.next());
			final QName attName=att.getName();
			if(attName.getNamespaceURI()==null || attName.getNamespaceURI().isEmpty())
				{
				root.setAttribute(toQNAME(attName),att.getValue());
				}
			else
				{	
				root.setAttributeNS(
						attName.getNamespaceURI(), 
						toQNAME(attName),
						att.getValue());
				}
			}
		
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			
			if(evt.isCharacters())
				{
				root.appendChild(owner.createTextNode(evt.asCharacters().getData()));
				}
			else if(evt.isProcessingInstruction())
				{
				}
			else if(evt.isEndElement())
				{
				break;
				}
			else if(evt.isStartElement())
				{
				if(this.skipQNames.contains(evt.asStartElement().getName()))
					{
					skip(evt.asStartElement(),r);
					}
				else
					{
					root.appendChild(parseDocument(owner,evt.asStartElement(),r));
					}
				}
			else 
				{
				LOG.warn("Cannot handle "+evt);
				}
			}
		return root;
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.targetQnames.isEmpty())
			{
			LOG.error("No target name defined");
			return -1;
			}
		XMLEventReader xmlReader=null;
		OutputStream outputStream = null;
		try {
			final String inputSource = oneFileOrNull(args);
			final TransformerFactory transformerFactory = TransformerFactory.newInstance();
			final Templates template = transformerFactory.newTemplates(new StreamSource(this.templateFile));
			final Transformer transformer = template.newTransformer();
			
			outputStream = openFileOrStdoutAsStream(this.outputFile);
			
			final XMLInputFactory xif = XMLInputFactory.newFactory();
			xmlReader = (inputSource==null?
					xif.createXMLEventReader(stdin()):
					xif.createXMLEventReader(new StreamSource(new File(inputSource)))	
					);
			
			final Document dom = DocumentBuilderFactory.newInstance().
					newDocumentBuilder().
					newDocument();
			
			while(xmlReader.hasNext())
				{
				final XMLEvent evt = xmlReader.nextEvent();
				if(!evt.isStartElement())
					{
					continue;
					}
				final StartElement SE = evt.asStartElement();
				if(!this.targetQnames.contains(SE.getName()))
					{
					continue;
					}
				
				dom.appendChild(parseDocument(dom,SE,xmlReader));
				transformer.transform(
						new DOMSource(dom),
						new StreamResult(outputStream)
						);
				
				while(dom.hasChildNodes()) dom.removeChild(dom.getFirstChild());
				}
			outputStream.flush();
			outputStream.close();
			outputStream=null;
			xmlReader.close();
			xmlReader=null;
			return 0;
			}
		catch (final Exception e) {
			LOG.error(e);
			return -1;
			}
		finally {
			CloserUtil.close(xmlReader);
			CloserUtil.close(outputStream);
			}
		}
	
	public static void main(final String[] args) {
		new XsltStream().instanceMainWithExit(args);
	}
}
