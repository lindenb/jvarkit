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
package com.github.lindenb.jvarkit.tools.misc;


import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URLEncoder;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ncbi.NcbiApiKey;
import com.github.lindenb.jvarkit.util.ncbi.NcbiConstants;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC

## Annotation files

a custom annotation file can be a xml or a plain text file.

Annotation will be inserted in a `<custom-annotation>` tag (not in the official dtd) under `<Entrezgene>`

### XML files

The program will insert the first XML tag containing a XML attribute `ncbi-gene-id` corresponding to the current NCBI gene identifier.

### Plain text files

block of lines starting after `## NCBI Gene ID.` followed by the current NCBI gene identifer will be inserted.



## Examples

### Example 1

search two genes and convert to markdown using the following XSLT stylesheet.

```xslt
<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet xmlns:xsl='http://www.w3.org/1999/XSL/Transform'  version='1.0' >
<xsl:output method="text" encoding="UTF-8"/>

<xsl:template match="Entrezgene-Set">
<xsl:apply-templates select="Entrezgene"/>
</xsl:template>

<xsl:template match="Entrezgene">
<xsl:variable name="geneId" select="Entrezgene_track-info/Gene-track/Gene-track_geneid"/>
<xsl:variable name="locus" select="Entrezgene_gene/Gene-ref/Gene-ref_locus"/>
<xsl:text>**</xsl:text>
<xsl:value-of select="$locus"/>
<xsl:text>** : </xsl:text>
<xsl:value-of select="Entrezgene_summary"/>

<xsl:if test="custom-annotation">
<xsl:text>
CUSTOM-ANNOTATIONS: </xsl:text><xsl:value-of select="normalize-space(custom-annotation)"/>
</xsl:if>
<xsl:text>

</xsl:text>
</xsl:template>

</xsl:stylesheet>
```


```
$ java -jar dist/ncbigenedump.jar SCN5A NOTCH2 |\
 	xsltproc stylesheet.xsl  -
```

**NOTCH2** : This gene encodes a member of the Notch family. Members of this Type 1 transmembrane protein family share structural characteristics including an extracellular domain consisting of multiple epidermal growth factor-like (EGF) repeats, and an intracellular domain consisting of multiple, different domain types. Notch family members play a role in a variety of developmental processes by controlling cell fate decisions. The Notch signaling network is an evolutionarily conserved intercellular signaling pathway which regulates interactions between physically adjacent cells. In Drosophilia, notch interaction with its cell-bound ligands (delta, serrate) establishes an intercellular signaling pathway that plays a key role in development. Homologues of the notch-ligands have also been identified in human, but precise interactions between these ligands and the human notch homologues remain to be determined. This protein is cleaved in the trans-Golgi network, and presented on the cell surface as a heterodimer. This protein functions as a receptor for membrane bound ligands, and may play a role in vascular, renal and hepatic development. Two transcript variants encoding different isoforms have been found for this gene. [provided by RefSeq, Jan 2011]

**SCN5A** : The protein encoded by this gene is an integral membrane protein and tetrodotoxin-resistant voltage-gated sodium channel subunit. This protein is found primarily in cardiac muscle and is responsible for the initial upstroke of the action potential in an electrocardiogram. Defects in this gene are a cause of long QT syndrome type 3 (LQT3), an autosomal dominant cardiac disease. Alternative splicing results in several transcript variants encoding different isoforms. [provided by RefSeq, Jul 2008]

### Example

The xml annotation file:

```html
<html xmlns="http://www.w3.org/1999/xhtml"><head><title>My custom annots</title></head><body>
<div ncbi-gene-id="6331"><h2>SCN5A</h2><p>This is the gene, Matilde.</p></div>
<div ncbi-gene-id="4853"><h2>NOTCH2</h2><p>Hajdu-Cheney syndrome</p></div>
</body></html>
```

```
java -jar dist/ncbigenedump.jar -C annot.html NOTCH2 SCN5A | xmllint --format -

(...)
      </Xtra-Terms>
    </Entrezgene_xtra-properties>
    <custom-annotation>
      <div ncbi-gene-id="6331">
        <h2>SCN5A</h2>
        <p>This is the gene, Matilde.</p>
      </div>
    </custom-annotation>
  </Entrezgene>
</Entrezgene-Set>

```

### Example

The text annotation file `annot.md`:

```markdown
# My Annotations

last updated 2018-08-29

## NCBI Gene ID. 4853 NOTCH2

encodes a member of the Notch family

## NCBI Gene ID. 6331 SCN5A

See also SCN10A
```

```
java -jar dist/ncbigenedump.jar -C annot.md NOTCH2 SCN5A | xmllint --format -

(...)
        <Xtra-Terms_tag>PROP</Xtra-Terms_tag>
        <Xtra-Terms_value>phenotype</Xtra-Terms_value>
      </Xtra-Terms>
    </Entrezgene_xtra-properties>
    <custom-annotation>
encodes a member of the Notch family

</custom-annotation>
  </Entrezgene>
</Entrezgene-Set>
```

### Example

custom annotation and XSLT:

```
$ java -jar dist/ncbigenedump.jar -C annot.md NOTCH2 SCN5A | xsltproc tansform.xsl -
```

output:

**NOTCH2** : This gene encodes a member of the Notch family. Members of this Type 1 transmembrane protein family share structural characteristics including an extracellular domain consisting of multiple epidermal growth factor-like (EGF) repeats, and an intracellular domain consisting of multiple, different domain types. Notch family members play a role in a variety of developmental processes by controlling cell fate decisions. The Notch signaling network is an evolutionarily conserved intercellular signaling pathway which regulates interactions between physically adjacent cells. In Drosophilia, notch interaction with its cell-bound ligands (delta, serrate) establishes an intercellular signaling pathway that plays a key role in development. Homologues of the notch-ligands have also been identified in human, but precise interactions between these ligands and the human notch homologues remain to be determined. This protein is cleaved in the trans-Golgi network, and presented on the cell surface as a heterodimer. This protein functions as a receptor for membrane bound ligands, and may play a role in vascular, renal and hepatic development. Two transcript variants encoding different isoforms have been found for this gene. [provided by RefSeq, Jan 2011]

CUSTOM-ANNOTATIONS: encodes a member of the Notch family

**SCN5A** : The protein encoded by this gene is an integral membrane protein and tetrodotoxin-resistant voltage-gated sodium channel subunit. This protein is found primarily in cardiac muscle and is responsible for the initial upstroke of the action potential in an electrocardiogram. Defects in this gene are a cause of long QT syndrome type 3 (LQT3), an autosomal dominant cardiac disease. Alternative splicing results in several transcript variants encoding different isoforms. [provided by RefSeq, Jul 2008]

CUSTOM-ANNOTATIONS: See also SCN10A


## Screenshots


https://twitter.com/yokofakun/status/1034825482808819713

![https://twitter.com/yokofakun/status/1034825482808819713](https://pbs.twimg.com/media/DlxwR3fXoAIFKO2.jpg)

END_DOC
 *
 */
@Program(name="ncbigenedump",
	keywords={"ncbi","gene","xml"}, 
	description="Dump XML results from a list of gene using NCBI/Eutils"
	)
public class NcbiGeneDump
	extends Launcher
	{
	private static final Logger LOG = Logger.build(NcbiGeneDump.class).make();
	private static final QName CUSTOM_ATTRIBUTE_QNAME = new QName("ncbi-gene-id");
	private static final QName CUSTOM_ANNOT_QNAME = new QName("custom-annotation");
	private static final String CUSTOM_ANNOT_LINE_PREFIX= "## NCBI Gene ID.";
	@Parameter(names={"-e","--email"},description="optional user email")
	private String email = null;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-skip","--skip"},description="Optional set of elements names to be ignored in the output. Spaces or comma separated. .eg: 'Gene-track '")
	private String skipTagsStr = "";
	@Parameter(names={"-L","-G","--list","--genes"},description="File containing a list of genes, can be a gene name or a ncbi gene id, one per line.")
	private File userGeneFile = null;
	@Parameter(names={"-T","--taxon"},description="taxon id.")
	private String taxonId = "9606";
	@Parameter(names={"--stdin"},description="read list of genes from stdin.")
	private boolean stdinFlags = false;
	@Parameter(names={"--seconds"},description="wait 'n' seconds between each calls.")
	private int wait_seconds = 2;
	@Parameter(names={"--abort"},description="Abort with error if a gene was not found with ncbi-esearch.")
	private boolean abortOnNotFound=false;
	@Parameter(names={"-C","--custom"},description=
			"Custom annotation file. See the main documentation.")
	private File customAnnotationFile=null;

	@ParametersDelegate
	private NcbiApiKey ncbiApiKey = new NcbiApiKey();

	private static final String tool="ncbigenedump";

	public NcbiGeneDump() {
		}
	
	private void skip(final XMLEventReader r) throws XMLStreamException
		{
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			if(evt.isEndElement()) break;
			else if(evt.isStartElement())
				{
				skip(r);
				}
			}
		}
	
	private boolean isInteger(final String s) {
		try {
			new Integer(s);
			return true;
		} catch(final NumberFormatException err) {
			return false;
		}
	}
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter pw=null;
		
		
		if(!this.ncbiApiKey.isApiKeyDefined()) {
			LOG.error("NCBI API key is not defined");
			return -1;
			}
		
		
		
		final Set<String> skipTags = 
				Arrays.stream(this.skipTagsStr.split("[ \t;,]")).
				filter(S->!StringUtil.isBlank(S)).
				collect(Collectors.toSet());
		
		
		try
			{
			
			final Set<String> geneIdentifiers = new HashSet<>();
			if(this.userGeneFile!=null) {
				IOUtil.slurpLines(this.userGeneFile).stream().map(S->S.trim()).forEach(G->geneIdentifiers.add(G));
				}
			
			if(!args.isEmpty()) {
				args.stream().map(S->S.trim()).forEach(G->geneIdentifiers.add(G));
				}
			else if(this.stdinFlags)
				{
				LOG.debug("read identifiers from stdin");
				IOUtil.slurpLines(stdin()).stream().map(S->S.trim()).forEach(G->geneIdentifiers.add(G));
				}
			
			geneIdentifiers.remove("");
			
			if(geneIdentifiers.isEmpty())
				{
				LOG.warn("no gene was defined.");
				}
			
			final Set<String> geneNames = new HashSet<>(geneIdentifiers.
					stream().
					filter(G->!isInteger(G)).
					collect(Collectors.toSet())
					);
			final Set<Integer> geneIds = new HashSet<>(geneIdentifiers.
					stream().
					filter(G->isInteger(G)).
					map(G->new Integer(G)).
					collect(Collectors.toSet())
					);
			
			final XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.FALSE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			xmlInputFactory.setXMLResolver(new XMLResolver() {
				@Override
				public Object resolveEntity(String publicID, String systemID, String baseURI, String namespace)
						throws XMLStreamException {
					LOG.info("ignoring DTD : "+publicID+" "+baseURI);
					return new ByteArrayInputStream(new byte[0]);
				}
			});
			
			final int batchSize=10;
			
			
			while(!geneNames.isEmpty()) {
				LOG.debug(""+geneNames.size()+" remain for esearch...");
				final Set<String> batchNames = new HashSet<>(batchSize);
				final Iterator<String> iter = geneNames.iterator();
				while(iter.hasNext() && batchNames.size() < batchSize) {
					batchNames.add(iter.next());
					iter.remove();
					}
				
				final StringBuilder query = new StringBuilder(
					batchNames.
						stream().
						map(G->"\""+G+"\"[GENE]").
						collect(Collectors.joining(" OR " ))
					);
				
				if(!StringUtil.isBlank(taxonId)) {
					query.insert(0, "(");
					query.append(") AND \""+this.taxonId+"\"[TID]");
				}
	
				
				final String url=
						NcbiConstants.esearch()+"?db=gene&term="+
						URLEncoder.encode(query.toString(), "UTF-8")+
						ncbiApiKey.getAmpParamValue()+
						"&retstart=0&retmode=xml&retmax="+NcbiConstants.RETMAX_MAX+
						(email==null?"":"&email="+URLEncoder.encode(email,"UTF-8"))+
						(tool==null?"":"&tool="+URLEncoder.encode(tool,"UTF-8"))
						;
				LOG.info(url);
				int nFound=0;
				XMLEventReader r=xmlInputFactory.createXMLEventReader(new StreamSource(url));
				while(r.hasNext())
					{
					final XMLEvent evt=r.nextEvent();
					if(evt.isStartElement())
						{
						final  String eName=evt.asStartElement().getName().getLocalPart();
						if(eName.equals("Id"))
							{
							geneIds.add(Integer.parseInt(r.getElementText()));
							nFound++;
							}
						else if( eName.equals("QuotedPhraseNotFound"))
							{
							final String notFound = r.getElementText();
							LOG.warn("NOT FOUND :" + notFound);
							if(abortOnNotFound) 
								{
								LOG.error("The following Entrez query was not found : "+notFound);
								return -1;
								}
							}
						}
					}
				CloserUtil.close(r);
				r=null;
				
				
				if(nFound!=batchNames.size())
					{
					LOG.error("Bad esearch result. Found "+nFound+" but expected "+batchNames.size()+". was " +
							String.join(" ",batchNames));
					}
				try { if(!geneNames.isEmpty()) Thread.sleep(this.wait_seconds * 1000); } catch(final Throwable err) {}
				}
			
			pw=super.openFileOrStdoutAsPrintWriter(outputFile);
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			final XMLEventWriter w=xof.createXMLEventWriter(pw);
			final XMLEventFactory eventFactory = XMLEventFactory.newInstance();
			w.add(eventFactory.createStartDocument("UTF-8", "1.0"));
			w.add(eventFactory.createStartElement(new QName("Entrezgene-Set"), null,null));
			w.add(eventFactory.createCharacters("\n"));
			while(!geneIds.isEmpty()) {
				final Set<Integer> batchIds = new HashSet<>(batchSize);
				final Iterator<Integer> iter = geneIds.iterator();
				while(iter.hasNext() && batchIds.size() < batchSize) {
					batchIds.add(iter.next());
					iter.remove();
					}
				
				final String url= NcbiConstants.efetch()+"?"+
						"db=gene"+
						ncbiApiKey.getAmpParamValue()+
						"&retmode=xml&id="+batchIds.stream().map(I->I.toString()).collect(Collectors.joining(","))+
						(email==null?"":"&email="+URLEncoder.encode(email,"UTF-8"))+
						(tool==null?"":"&tool="+URLEncoder.encode(tool,"UTF-8"))
						;
				LOG.info(url+" remains "+geneIds.size()+" ID(s).");
				final XMLEventReader r = xmlInputFactory.createXMLEventReader(new StreamSource(url));
				boolean in_gene = false;
				Integer current_gene_id = null;
				
				while(r.hasNext())
					{
					final XMLEvent evt=r.nextEvent();
					
					switch(evt.getEventType())
						{
						case XMLEvent.ATTRIBUTE:
							{
							if(in_gene) w.add(evt);
							break;
							}
						case XMLEvent.START_DOCUMENT:
						case XMLEvent.END_DOCUMENT:
							{
							in_gene=false;
							break;
							}
						case XMLEvent.START_ELEMENT:
							{
							
							final  String localName= evt.asStartElement().getName().getLocalPart();
							if(localName.equals("Entrezgene"))
								{
								in_gene = true;
								current_gene_id = null;
								w.add(evt);
								}
							else if(in_gene && localName.equals("Gene-track_geneid"))
								{
								w.add(evt);
								final XMLEvent evt2 = r.nextEvent();
								if(!evt2.isCharacters()) throw new XMLStreamException("expected a text node for Gene-track_geneid" ,evt2.getLocation());
								current_gene_id  = new Integer(evt2.asCharacters().getData().trim());
								if(!batchIds.remove(current_gene_id)) {
									LOG.warn("found NCBI-GENE-ID "+current_gene_id+" but it is not in my list");
									}
								w.add(evt2);
								}
							else if(in_gene && skipTags.contains(localName))
								{
								skip(r);
								}
							else if(in_gene)
								{
								w.add(evt);
								}
							break;
							}
						case XMLEvent.END_ELEMENT:
							{
							if(in_gene)
								{
								final  String localName= evt.asEndElement().getName().getLocalPart();
								if(current_gene_id!=null && localName.equals("Entrezgene"))
									{
									insertCustomAnnotationsForGeneId(w,current_gene_id,eventFactory);
									}
								w.add(evt);
								if(localName.equals("Entrezgene"))
									{
									w.add(eventFactory.createCharacters("\n"));
									in_gene = false;
									current_gene_id = null;
									}
								}
							break; 
							}
						case XMLEvent.COMMENT:break;
						case XMLEvent.PROCESSING_INSTRUCTION:break;
						case XMLEvent.DTD:
							{
							break;	
							}
						case XMLEvent.SPACE:break;
						case XMLEvent.CHARACTERS:
							{
							if(in_gene) w.add(evt);
							break;
							}
						default:
							{
							throw new IllegalStateException("XML evt no handled: "+evt);
							}
						}
					}
				r.close();
				
				if(!batchIds.isEmpty())
					{
					final String msg = "The following NCBI gene identifiers were not found: "+
							batchIds.stream().
							map(I->String.valueOf(I)).
								collect(Collectors.joining(" , "))
							;
					if(abortOnNotFound) 
						{
						LOG.error(msg);
						return -1;
						}
					else
						{
						for(final Integer gid: batchIds)
							{
							w.add(eventFactory.createComment("NOT FOUND : https://www.ncbi.nlm.nih.gov/gene/ "+ gid));
							}
						LOG.warn(msg);
						}
					}
				try { if(!geneIds.isEmpty()) Thread.sleep(this.wait_seconds * 1000); } catch(final Throwable err) {}				
				}//end while ids
			w.add(eventFactory.createEndElement(new QName("Entrezgene-Set"),null));
			w.add(eventFactory.createEndDocument());
			w.flush();
			w.close();
			pw.flush();
			pw.close();
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(pw);
			}
		}
	
	private void insertCustomAnnotationsForGeneId(
			final XMLEventWriter w,
			final Integer geneid_int,
			final XMLEventFactory eventFactory
			) throws XMLStreamException,IOException
		{
		
		if(customAnnotationFile==null || geneid_int==null) return;
		final String geneid = geneid_int.toString();
		IOUtil.assertFileIsReadable(customAnnotationFile);
		
		final String suffix = IOUtil.fileSuffix(this.customAnnotationFile);
		final BufferedReader br = IOUtil.openFileForBufferedReading(this.customAnnotationFile);
		if(suffix!=null && (suffix.equals(".xml") || suffix.equals(".xml.gz") || suffix.equals(".html") || suffix.equals(".html.gz")))
			{
			final XMLEventReader r = XMLInputFactory.newInstance().createXMLEventReader(br);
			while(r.hasNext())
				{
				final XMLEvent evt = r.nextEvent();
				if(evt.isStartElement())
					{
					final StartElement R = evt.asStartElement();
					final Attribute att = R.getAttributeByName(CUSTOM_ATTRIBUTE_QNAME);
					if(att!=null && att.getValue().trim().replace(",", "").equals(geneid))
						{
						w.add(eventFactory.createStartElement(CUSTOM_ANNOT_QNAME, null, null));
						w.add(evt);
						copy(r,w);
						w.add(eventFactory.createEndElement(CUSTOM_ANNOT_QNAME, null));
						break;
						}
					}
			
				}
			r.close();
			}
		else
			{
			boolean in_gene=false;
			String line;
			while((line=br.readLine())!=null)
				{
				String line2 = line.trim().replaceAll("[ \t]+", " ");
				if(line2.startsWith(CUSTOM_ANNOT_LINE_PREFIX))
					{
					if(in_gene) break;
					line2 = line2.substring(CUSTOM_ANNOT_LINE_PREFIX.length()).trim().replace(",", "");
					line2 = line2.split("\\p{Blank}",2)[0];
					
					if(line2.equals(geneid))
						{
						in_gene = true;
						w.add(eventFactory.createStartElement(CUSTOM_ANNOT_QNAME, null, null));
						continue;
						}
					}
				if(in_gene)
					{
					w.add(eventFactory.createCharacters(line));
					w.add(eventFactory.createCharacters("\n"));
					}
				}
			if(in_gene) ; {
				w.add(eventFactory.createEndElement(CUSTOM_ANNOT_QNAME, null));
				}
			}
		br.close();
		}

	private void copy(final XMLEventReader r,final XMLEventWriter w) throws XMLStreamException
		{
		while(r.hasNext())
			{
			final XMLEvent evt = r.nextEvent();
			w.add(evt);
			switch(evt.getEventType())
				{
				case XMLEvent.END_ELEMENT: return;
				case XMLEvent.START_ELEMENT: copy(r,w);break;
				default:break;
				}
			}
		}
	
	public static void main(final String[] args) {
		new NcbiGeneDump().instanceMainWithExit(args);
	}
}
