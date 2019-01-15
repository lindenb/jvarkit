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
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import org.apache.http.client.config.CookieSpecs;
import org.apache.http.client.config.RequestConfig;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpHead;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.message.BasicStatusLine;
import org.apache.http.StatusLine;
import org.apache.http.HttpVersion;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;

/**

## Example

```
$ java -jar dist/pubmeddump.jar 'bioinformatics 2001' 2> /dev/null |\
	java -jar dist/pubmed404.jar  2> /dev/null 

#PMID	TITLE	YEAR	URL	Status
29520589	Expression of Colocasia esculenta tuber agglutinin in Indian mustard provides resistance against Lipaphis erysimi and the expressed protein is non-allergenic.2018	http://www.fao.org/docrep/007/y0820e/y0820e00.HTM	200
29520589	Expression of Colocasia esculenta tuber agglutinin in Indian mustard provides resistance against Lipaphis erysimi and the expressed protein is non-allergenic.2018	http://www.icmr.nic.in/guide/Guidelines%20for%20Genetically%20Engineered%20Plants.pdf	-1
28482857	Horizontal gene transfer is not a hallmark of the human genome.	2017	https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0607-3	200
27899642	The UCSC Genome Browser database: 2017 update.	2017	http://genome.ucsc.edu/	200
27797935	High hospital research participation and improved colorectal cancer survival outcomes: a population-based study.	2017	http://www.bmj.com/company/products-services/rights-and-licensing/	403
25505092	NMRFAM-SPARKY: enhanced software for biomolecular NMR spectroscopy.	2015	http://pine.nmrfam.wisc.edu/download_packages.html	200
25505092	NMRFAM-SPARKY: enhanced software for biomolecular NMR spectroscopy.	2015	http://www.nmrfam.wisc.edu/nmrfam-sparky-distribution.htm	200
25428374	The UCSC Genome Browser database: 2015 update.	2015	http://genome.ucsc.edu	200
26356339	A Simple but Powerful Heuristic Method for Accelerating k-Means Clustering of Large-Scale Data in Life Science.	null	http://mlab.cb.k.u-tokyo.ac.jp/~ichikawa/boostKCP/	200
24794704	Usefulness of the Shock Index as a secondary triage tool.	2015	http://group.bmj.com/group/rights-licensing/permissions	403
24225322	Progenetix: 12 years of oncogenomic data curation.	2014	http://www.progenetix.org	200
24137000	Updates of the HbVar database of human hemoglobin variants and thalassemia mutations.	2014	http://globin.bx.psu.edu/hbvar	200
24137000	Updates of the HbVar database of human hemoglobin variants and thalassemia mutations.	2014	http://www.findbase.org	200
24137000	Updates of the HbVar database of human hemoglobin variants and thalassemia mutations.	2014	http://www.lovd.nl	200
23564938	DAMBE5: a comprehensive software package for data analysis in molecular biology and evolution.	2013	http://dambe.bio.uottawa.ca	200
22689647	SIFT web server: predicting effects of amino acid substitutions on proteins.	2012	http://sift-dna.org	200
22600740	Cyber-T web server: differential analysis of high-throughput data.	2012	http://cybert.ics.uci.edu/	200
21742331	An open source lower limb model: Hip joint validation.	2011	https://simtk.org/home/low_limb_london	200
21593132	Java bioinformatics analysis web services for multiple sequence alignment--JABAWS:MSA.	2011	http://www.compbio.dundee.ac.uk/jabaws	200
20228129	DensiTree: making sense of sets of phylogenetic trees.	2010	http://compevol.auckland.ac.nz/software/DensiTree/	404
19380317	CELLULAR OPEN RESOURCE (COR): current status and future directions.	2009	http://www.cellml.org/specifications/	200
18948284	OperonDB: a comprehensive database of predicted operons in microbial genomes.	2009	http://operondb.cbcb.umd.edu	200
18368364	Simulator for neural networks and action potentials.	2007	http://snnap.uth.tmc.edu	-1
18367465	An improved general amino acid replacement matrix.	2008	http://atgc.lirmm.fr/LG	404
18238804	Interoperability with Moby 1.0--it's better than sharing your toothbrush!	2008	http://www.biomoby.org/	200
18174178	PRALINETM: a strategy for improved multiple alignment of transmembrane proteins.	2008	http://www.ibi.vu.nl/programs/pralinewww	200
17221864	HbVar database of human hemoglobin variants and thalassemia mutations: 2007 update.	2007	http://globin.bx.psu.edu/hbvar	200
17221864	HbVar database of human hemoglobin variants and thalassemia mutations: 2007 update.	2007	http://www.goldenhelix.org/xprbase	403
(...)
```

*/
@Program(name="pubmed404",
description="Test if URL in the pubmed abstracts are reacheable.",
keywords={"pubmed","url"}
)
public class Pubmed404  extends Launcher{
	private static final Logger LOG = Logger.build(Pubmed404.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outFile=null;
	@Parameter(names={"-t","--timeout"},description="timeout in seconds")
	private int timeoutSeconds = 5;
	@Parameter(names={"-c","--collapse"},description="Only one URL per article. Print the '200/OK' first.")
	private boolean collapse = false;

	private CloseableHttpClient httpClient = null;
	private final Pattern wsSplitter = Pattern.compile("[\\[\\] ,()'\";\n\r\t]+");
	private StatusLine getStatus(final String url) {
		HttpHead httpHead = null;
		CloseableHttpResponse response = null;
		final RequestConfig requestConfig = RequestConfig.custom().
				setConnectTimeout(this.timeoutSeconds * 1000).
				setCookieSpec(CookieSpecs.STANDARD).
				build();
		try {
			httpHead = new HttpHead(url);
			httpHead.setConfig(requestConfig);
			response = this.httpClient.execute(httpHead);
			return response.getStatusLine();
			}
		catch(final java.net.UnknownHostException err) {
			return new BasicStatusLine(HttpVersion.HTTP_1_1,404,String.valueOf(err.getMessage()));
			}
		catch(final javax.net.ssl.SSLHandshakeException err) {
			return new BasicStatusLine(HttpVersion.HTTP_1_1,526,String.valueOf(err.getMessage()));
			}
		catch(final Throwable err) {
			return new BasicStatusLine(HttpVersion.HTTP_1_1,100,String.valueOf(err.getMessage()));
			}
		finally
			{
			CloserUtil.close(response);
			}
	}
	
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
	
	private void scanArticle(
			PrintWriter out,
			final String rootName,
			final XMLEventReader r
			) throws XMLStreamException,IOException {
			String article_pmid=null;
			String article_title="" ;
			String article_year=null;
			String abstractText = "";
			boolean PubDate=false;
			while(r.hasNext()) {
				final XMLEvent evt=r.nextEvent();
				
				if(evt.isStartElement()) {
					final StartElement start = evt.asStartElement();
					final String eltName = start.getName().getLocalPart();
					if(article_pmid==null && eltName.equals("PMID")) {
						article_pmid=r.getElementText();
					}
					else if(article_title.isEmpty() && eltName.equals("ArticleTitle")) {
						article_title= textContent(r);
					}
					else if(eltName.equals("PubDate")) {
						PubDate=true;
					}
					else if(article_year==null && PubDate && eltName.equals("Year")) {
						article_year=r.getElementText();
					}
					else if(eltName.equals("Abstract")) {
						abstractText = textContent(r);
						}
					}
			else if(evt.isEndElement()) {
				final EndElement end = evt.asEndElement();
				final String eltName = end.getName().getLocalPart();
				if(eltName.equals("PubDate")) {
					PubDate=false;
					}
				else if(eltName.equals(rootName)) 
					{
					break;
					}
				}
			
			}//end of xml read
		final String tokens[] = this.wsSplitter.split(abstractText);
		
		final List<String> lines = new ArrayList<>();
		for(String token: tokens)
			{
			while(!token.isEmpty()) {
				int c = token.charAt(token.length()-1);
				if(!(c=='.' || c>128)) break;
				token=token.substring(0,token.length()-1);
			}
			if(token.isEmpty()) continue;
			if(!IOUtil.isUrl(token)) {
				if(token.startsWith("http")) LOG.debug("strange url: "+token);
				continue;
			}
			final StatusLine status = getStatus(token);
			final StringBuilder sb = new StringBuilder();
			
			sb.append(article_pmid);
			sb.append('\t');
			sb.append(article_title);
			sb.append('\t');
			sb.append(article_year);
			sb.append('\t');
			sb.append(token);
			sb.append('\t');
			sb.append(status.getStatusCode());
			sb.append('\t');
			sb.append(status.getReasonPhrase());
			
			if(this.collapse && status.getStatusCode()==200)
				{
				out.println(sb.toString());
				out.flush();
				return;
				}
			else
				{
				lines.add(sb.toString());
				}
			}
		if(lines.isEmpty()) return;
		if(this.collapse)
			{
			out.println(lines.get(0));
			}
		else
			{
			for(final String line:lines) out.println(line);
			}
		out.flush();
		}	

	@Override
	public int doWork(final List<String> args) {
		final String inputName= oneFileOrNull(args);
		PrintWriter out=null;
		XMLEventReader r=null;
		InputStream in=null;
		try {
			/** create http client */
			
			
			this.httpClient = HttpClients.createSystem();//createDefault();
		
			
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
			
			out = super.openFileOrStdoutAsPrintWriter(this.outFile);
			out.println("#PMID\tTITLE\tYEAR\tURL\thttp.code\thttp.reason");
			
			while(r.hasNext()) {
				final XMLEvent evt= r.nextEvent();
				if(evt.isStartElement() )	{
					final String localName= evt.asStartElement().getName().getLocalPart();
					if(localName.equals("PubmedArticle") || localName.equals("PubmedBookArticle"))
						{
						scanArticle(out,localName,r);	
						}
					}
				}
			
			r.close();r=null;
			in.close();in=null;
			out.flush();out.close();out=null;
			return 0;
		} catch (final Exception err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(r);
			CloserUtil.close(in);
			CloserUtil.close(out);
			CloserUtil.close(this.httpClient);
		}

		}
	
public static void main(final String[] args)
	{
	new Pubmed404().instanceMainWithExit(args);
	}
}
