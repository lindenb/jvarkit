/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.pubmed;


import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.PrintWriter;
import java.net.URLEncoder;
import java.util.List;
import java.util.stream.Collectors;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;

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

END_DOC
 *
 */
@Program(name="pubmeddump",keywords={"ncbi","pubmed","xml"}, description="Dump XML results from pubmed/Eutils")
public class PubmedDump
	extends Launcher
	{
	private static final Logger LOG = Logger.build(PubmedDump.class).make();

	@Parameter(names={"-e","--email"},description="optional user email")
	private String email = null;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	private String tool="pubmedump";

	public PubmedDump() {
		}
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter pw=null;
		
		if(args.isEmpty())
			{
			LOG.error("Query missing");
			return -1;
			}
		final String query=args.stream().collect(Collectors.joining(" ")).trim();
		
		if(query.isEmpty())
			{
			LOG.error("query is empty");
			return -1;
			}
		
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
					LOG.info("ignoring "+publicID+" "+baseURI+" "+namespace);
					return new ByteArrayInputStream(new byte[0]);
				}
			});
			
			String url=
					"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="+
					URLEncoder.encode(query.toString(), "UTF-8")+	
					"&retstart=0&retmax=0&usehistory=y&retmode=xml"+
					(email==null?"":"&email="+URLEncoder.encode(email,"UTF-8"))+
					(tool==null?"":"&tool="+URLEncoder.encode(tool,"UTF-8"))
					;
			LOG.info(url);
			long count=-1;
			String WebEnv=null;
			String QueryKey=null;
			XMLEventReader r=xmlInputFactory.createXMLEventReader(new StreamSource(url));
			while(r.hasNext())
				{
				XMLEvent evt=r.nextEvent();
				if(evt.isStartElement())
					{
					String eName=evt.asStartElement().getName().getLocalPart();
					if(eName.equals("Count") && count==-1)
						{
						count=Long.parseLong(r.getElementText());
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
			
			if(count<0 || WebEnv==null || QueryKey==null)
				{
				LOG.error("Bad esearch result");
				return -1;
				}
			pw=super.openFileOrStdoutAsPrintWriter(outputFile);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			XMLEventWriter w=xof.createXMLEventWriter(pw);
			long nFound=0L;
			int depth=0;
			while(nFound< count)
				{
				final int ret_max=90000;
				LOG.info("nFound:"+nFound+"/"+count);
				url= "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&WebEnv="+
						URLEncoder.encode(WebEnv,"UTF-8")+
						"&query_key="+URLEncoder.encode(QueryKey,"UTF-8")+
						"&retmode=xml&retmax="+ret_max+"&retstart="+nFound+
						(email==null?"":"&email="+URLEncoder.encode(email,"UTF-8"))+
						(tool==null?"":"&tool="+URLEncoder.encode(tool,"UTF-8"))
						;
				LOG.info(url);
				int curr_count=0;
				r = xmlInputFactory.createXMLEventReader(new StreamSource(url));
				
				while(r.hasNext())
					{
					XMLEvent evt=r.nextEvent();
					
					switch(evt.getEventType())
						{
						case XMLEvent.ATTRIBUTE:
							{
							if(depth>0) w.add(evt);
							break;
							}
						case XMLEvent.START_DOCUMENT:
							{
							if(nFound==0)
								{
								w.add(evt);
								}
							break;
							}
						case XMLEvent.END_DOCUMENT:
							{
							if(nFound>= count)
								{
								w.add(evt);
								}
							break;
							}
						case XMLEvent.START_ELEMENT:
							{
							if(depth==0 && nFound==0)
								{
								w.add(evt);
								}
							else if(depth==1)
								{
								String localName= evt.asStartElement().getName().getLocalPart();
								if(!(localName.equals("PubmedArticle") || localName.equals("PubmedBookArticle")))
									{
									throw new IllegalStateException("Not PubmedArticle: "+evt);
									}
								++curr_count;
								++nFound;
								w.add(evt);
								}
							else if(depth>1)
								{
								w.add(evt);
								}
							depth++;
							break;
							}
						case XMLEvent.END_ELEMENT:
							{
							depth--;
							if(depth>0)
								{
								w.add(evt);
								}
							else if(nFound>=count)//depth ==0
								{
								w.add(evt);
								}
							
							break; 
							}
						case XMLEvent.COMMENT:break;
						case XMLEvent.PROCESSING_INSTRUCTION:break;
						case XMLEvent.DTD:
							{
							if(nFound==0) w.add(evt);
							break;	
							}
						case XMLEvent.SPACE:break;
						case XMLEvent.CHARACTERS:
							{
							if(depth>1) w.add(evt);
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
					}
				}
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
	
	public static void main(String[] args) {
		new PubmedDump().instanceMainWithExit(args);
	}
}
