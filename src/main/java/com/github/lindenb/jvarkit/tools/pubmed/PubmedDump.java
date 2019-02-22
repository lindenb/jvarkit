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
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
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
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ncbi.NcbiApiKey;
import com.github.lindenb.jvarkit.util.ncbi.NcbiConstants;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

/**

BEGIN_DOC

## Example

```
$  java -jar dist/pubmeddump.jar "Lindenbaum P" | grep ArticleTitle
			
	<ArticleTitle>NGS library preparation may generate artifactual integration sites of AAV vectors.</ArticleTitle>
    <ArticleTitle>Mutations in FAM111B cause hereditary fibrosing poikiloderma with tendon contracture, myopathy, and pulmonary fibrosis.</ArticleTitle>
    <ArticleTitle>[The Spanish Association of Surgeon's audited teaching programme for rectal cancer. Results after six years].</ArticleTitle>
    <ArticleTitle>Common variants at SCN5A-SCN10A and HEY2 are associated with Brugada syndrome, a rare disease with high risk of sudden cardiac death.</ArticleTitle>
    <ArticleTitle>The 3rd DBCLS BioHackathon: improving life science data integration with Semantic Web technologies.</ArticleTitle>
    <ArticleTitle>Mass spectrometry-based identification of native cardiac Nav1.5 channel alpha subunit phosphorylation sites.</ArticleTitle>
    <ArticleTitle>BioStar: an online question &amp; answer resource for the bioinformatics community.</ArticleTitle>
    <ArticleTitle>Knime4Bio: a set of custom nodes for the interpretation of next-generation sequencing data with KNIME.</ArticleTitle>
    <ArticleTitle>Truncating mutations in the last exon of NOTCH2 cause a rare skeletal disorder with osteoporosis.</ArticleTitle>
    <ArticleTitle>The Gene Wiki: community intelligence applied to human gene annotation.</ArticleTitle>
    <ArticleTitle>Robust physical methods that enrich genomic regions identical by descent for linkage studies: confirmation of a locus for osteogenesis imperfecta.</ArticleTitle>
    <ArticleTitle>Association of autism with polymorphisms in the paired-like homeodomain transcription factor 1 (PITX1) on chromosome 5q31: a candidate gene analysis.</ArticleTitle>
    <ArticleTitle>Haplotypes in the gene encoding protein kinase c-beta (PRKCB1) on chromosome 16 are associated with autism.</ArticleTitle>
    <ArticleTitle>RoXaN, a novel cellular protein containing TPR, LD, and zinc finger motifs, forms a ternary complex with eukaryotic initiation factor 4G and rotavirus NSP3.</ArticleTitle>
    <ArticleTitle>CloneIt: finding cloning strategies, in-frame deletions and frameshifts.</ArticleTitle>
    <ArticleTitle>In vivo and in vitro phosphorylation of rotavirus NSP5 correlates with its localization in viroplasms.</ArticleTitle>
```

## See also

 * https://gist.github.com/lindenb/6bfb49fd8bc3dd27d99f
 * https://gist.github.com/lindenb/5d7773a93d8c2b0edbd4c01bf8834919

END_DOC

*/
@Program(name="pubmeddump",keywords={"ncbi","pubmed","xml"}, 
	description="Dump XML results from pubmed/Eutils",
	biostars= {270498,365479},
	modificationDate="20190222"
	)
public class PubmedDump
	extends Launcher
	{
	private static final Logger LOG = Logger.build(PubmedDump.class).make();

	@Parameter(names={"-e","--email"},description="optional user email")
	private String email = null;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-skip","--skip"},description="[20180302]  Optional set of elements names to be ignored in the output. Spaces or comma separated. .eg: 'AuthorList PubmedData '")
	private String skipTagsStr = "";
	@ParametersDelegate
	private NcbiApiKey ncbiApiKey = new NcbiApiKey();

	private String tool="pubmedump";

	public PubmedDump() {
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
		
		final String query=args.stream().
				collect(Collectors.joining(" ")).
				trim();
		
		if(StringUtil.isBlank(query))
			{
			LOG.error("query is empty");
			return -1;
			}
		
		final Set<String> skipTags = 
				Arrays.stream(this.skipTagsStr.split("[ \t;,]")).
				filter(S->!StringUtil.isBlank(S)).
				collect(Collectors.toSet());
		
		
		try
			{
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
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
			
			String url=
					NcbiConstants.esearch()+"?db=pubmed&term="+
					StringUtils.escapeHttp(query.toString())+
					ncbiApiKey.getAmpParamValue()+
					"&retstart=0&retmax=0&usehistory=y&retmode=xml"+
					(email==null?"":"&email="+StringUtils.escapeHttp(email))+
					(tool==null?"":"&tool="+StringUtils.escapeHttp(tool))
					;
			LOG.info(url);
			long expected_total_count=-1;
			String WebEnv=null;
			String QueryKey=null;
			XMLEventReader r=xmlInputFactory.createXMLEventReader(new StreamSource(url));
			while(r.hasNext())
				{
				final XMLEvent evt=r.nextEvent();
				if(evt.isStartElement())
					{
					final  String eName=evt.asStartElement().getName().getLocalPart();
					if(eName.equals("Count") && expected_total_count==-1)
						{
						expected_total_count=Long.parseLong(r.getElementText());
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
			CloserUtil.close(r);
			r=null;
			
			if(expected_total_count<0 || WebEnv==null || QueryKey==null)
				{
				LOG.error("Bad esearch result");
				return -1;
				}
			pw=super.openFileOrStdoutAsPrintWriter(outputFile);
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			final XMLEventWriter w=xof.createXMLEventWriter(pw);
			long total_found_so_far  = 0L;
			boolean end_document_printed=false;
			
			while(total_found_so_far< expected_total_count)
				{
				final int ret_max=90000;
				LOG.info("nFound:"+total_found_so_far+"/"+expected_total_count);
				url= NcbiConstants.efetch()+"?"+
						"db=pubmed&WebEnv="+
						StringUtils.escapeHttp(WebEnv)+
						ncbiApiKey.getAmpParamValue()+
						"&query_key="+StringUtils.escapeHttp(QueryKey)+
						"&retmode=xml&retmax="+ret_max+"&retstart="+total_found_so_far+
						(email==null?"":"&email="+StringUtils.escapeHttp(email))+
						(tool==null?"":"&tool="+StringUtils.escapeHttp(tool))
						;
				LOG.info(url);
				int curr_count=0;
				r = xmlInputFactory.createXMLEventReader(new StreamSource(url));
				int current_dom_depth =0;
				
				while(r.hasNext())
					{
					final XMLEvent evt=r.nextEvent();
					
					switch(evt.getEventType())
						{
						case XMLEvent.ATTRIBUTE:
							{
							if(current_dom_depth>0) w.add(evt);
							break;
							}
						case XMLEvent.START_DOCUMENT:
							{
							if(total_found_so_far==0)
								{
								w.add(evt);
								}
							break;
							}
						case XMLEvent.END_DOCUMENT:
							{
							if(total_found_so_far>= expected_total_count)
								{
								end_document_printed = true;
								w.add(evt);
								}
							break;
							}
						case XMLEvent.START_ELEMENT:
							{
							final  String localName= evt.asStartElement().getName().getLocalPart();

							if(current_dom_depth==0)
								{
								if(!localName.equals("PubmedArticleSet")) {
									throw new XMLStreamException("Expected <PubmedArticleSet> but got <"+localName+">",evt.getLocation());
									}
								if( total_found_so_far == 0)
									{
									w.add(evt);
									}
								current_dom_depth++;
								}
							else if(current_dom_depth==1)
								{
								if(!(localName.equals("PubmedArticle") || localName.equals("PubmedBookArticle")))
									{
									throw new IllegalStateException("Not PubmedArticle: "+evt);
									}
								++curr_count;
								++total_found_so_far;
								++current_dom_depth;
								w.add(evt);
								}
							else if(current_dom_depth>1)
								{
								if(skipTags.contains(localName))
									{
									skip(r);
									}
								else
									{
									w.add(evt);
									current_dom_depth++;
									}
								}
							else
								{
								LOG.warn("unmatched case <"+localName+"> depth:"+current_dom_depth);
								}
							break;
							}
						case XMLEvent.END_ELEMENT:
							{
							current_dom_depth--;
							if(current_dom_depth>0)
								{
								w.add(evt);
								}
							else if(total_found_so_far>=expected_total_count)//depth ==0
								{
								end_document_printed = true;
								w.add(evt);
								}
							
							break; 
							}
						case XMLEvent.COMMENT:break;
						case XMLEvent.PROCESSING_INSTRUCTION:break;
						case XMLEvent.DTD:
							{
							if(total_found_so_far==0) w.add(evt);
							break;	
							}
						case XMLEvent.SPACE:break;
						case XMLEvent.CHARACTERS:
							{
							if(current_dom_depth>1) w.add(evt);
							break;
							}
						default:
							{
							throw new IllegalStateException("XML evt no handled: "+evt);
							}
						}
					}
				
				r.close();
				if(curr_count==0)
					{
					LOG.info("Nothing found . Exiting.");
					if(!end_document_printed) {
						final XMLEventFactory xef = XMLEventFactory.newFactory();
						w.add(xef.createEndElement(new QName("PubmedArticleSet"),Collections.emptyIterator()));
						w.add(xef.createEndDocument());
						}
					break;
					}
				else
					{
					LOG.info("found "+curr_count+" total "+total_found_so_far+" expect "+expected_total_count);
					}
				}
			w.flush();
			w.close();
			pw.flush();
			pw.close();
			LOG.info("Done. found "+total_found_so_far+" / expected:" +expected_total_count+" articles.");
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
	
	public static void main(String[] args) {
		new PubmedDump().instanceMainWithExit(args);
	}
}
