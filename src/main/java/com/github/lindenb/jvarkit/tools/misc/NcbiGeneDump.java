/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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


import java.io.ByteArrayInputStream;
import java.io.File;
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

## Example

END_DOC
 *
 */
@Program(name="ncbigenedump",keywords={"ncbi","gene","xml"}, 
	description="Dump XML results from gene/Eutils"
	)
public class NcbiGeneDump
	extends Launcher
	{
	private static final Logger LOG = Logger.build(NcbiGeneDump.class).make();

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

	@ParametersDelegate
	private NcbiApiKey ncbiApiKey = new NcbiApiKey();

	private String tool="ncbigenedump";

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
		
		if(args.isEmpty())
			{
			LOG.error("Query missing");
			return -1;
			}
		
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
				final Set<String> batchNames = new HashSet<>(batchSize);
				final Iterator<String> iter = geneNames.iterator();
				while(iter.hasNext() && batchNames.size() < batchSize) {
					batchNames.add(iter.next());
					iter.remove();
					}
				
				final StringBuilder query = new StringBuilder(
					batchNames.
						stream().
						map(G->"\""+G+"\"[TODO]").
						collect(Collectors.joining(" OR " ))
					);
				
				if(!StringUtil.isBlank(taxonId)) {
					query.insert(0, "(");
					query.append(") AND \""+this.taxonId+"\"[TODO]");
				}
	
				
				final String url=
						NcbiConstants.esearch()+"?db=gene&term="+
						URLEncoder.encode(query.toString(), "UTF-8")+
						ncbiApiKey.getAmpParamValue()+
						"&retstart=0&retmode=xml"+
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
						}
					}
				CloserUtil.close(r);
				r=null;
				
				
				if(nFound!=batchNames.size())
					{
					LOG.error("Bad esearch result. Found "+nFound+" but expected "+batchNames.size()+". was " +
							String.join(" ",batchNames));
					}
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
				LOG.info(url);
				final XMLEventReader r = xmlInputFactory.createXMLEventReader(new StreamSource(url));
				boolean in_gene = false;
				
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
								w.add(evt);
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
								w.add(evt);
								w.add(eventFactory.createCharacters("\n"));
								final  String localName= evt.asEndElement().getName().getLocalPart();
								if(localName.equals("Entrezgene"))
									{
									in_gene = false;
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
				}//end while ids
			w.add(eventFactory.createEndElement(new QName(""),null));
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
	
	public static void main(final String[] args) {
		new NcbiGeneDump().instanceMainWithExit(args);
	}
}
