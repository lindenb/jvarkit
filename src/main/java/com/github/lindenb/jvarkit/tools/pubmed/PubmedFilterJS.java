/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
import java.nio.file.Path;
import java.util.List;

import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.JexlEngine;
import org.apache.commons.jexl2.MapContext;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import gov.nih.nlm.ncbi.pubmed.ObjectFactory;
import gov.nih.nlm.ncbi.pubmed.PubmedArticle;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import  jakarta.xml.bind.JAXBContext;
import  jakarta.xml.bind.JAXBElement;
import  jakarta.xml.bind.Marshaller;
import  jakarta.xml.bind.Unmarshaller;


/**

## Example

<blockquote class="twitter-tweet" lang="en"><p>Same first and last author. Did you see this before? (via <a href="https://twitter.com/joedunckley">@joedunckley</a>) <a href="http://t.co/D0CEqC8UOu">http://t.co/D0CEqC8UOu</a></p>&mdash; Nicolas Robine (@notSoJunkDNA) <a href="https://twitter.com/notSoJunkDNA/statuses/506169374953447424">August 31, 2014</a></blockquote>
<script async src="//platform.twitter.com/widgets.js" charset="utf-8"></script>

The javascript file:

```javascript
function getNamePart(author,type)
	{
	var L= author.lastNameOrForeNameOrInitialsOrSuffixOrCollectiveName;
	for(var i=0;i< L.size();++i)
		{
		if(L.get(i).getClass().getSimpleName().equals(type))
			{
			return L.get(i).value;
			}
		}
	return "";
	}

function lastName(author)
	{
	return getNamePart(author,"LastName");
	}
function foreName(author)
	{
	return getNamePart(author,"ForeName");
	}


function accept(article)
	{
	var authorList= article.medlineCitation.article.authorList;
	if(! authorList ) return false;
	var authors =authorList.author;
	if(!authors || authors.size()<3) return false;
	var first= authors.get(0);
	var last= authors.get( authors.size()-1);

	return  lastName(first).equals(lastName(last)) && 
		foreName(first).equals(foreName(last)) && 
		lastName(first).length() >0 &&
		foreName(first).length() >0 
		;
	}


accept(article);
```


```bash
java -jar dist/pubmeddump.jar '"Neuro Oncol"[TA]' |\
java -jar dist/pubmedfilterjs.jar -f authors.js |\
xmllint  --format -  | grep PMID

      <PMID Version="1">25165367</PMID>
      <PMID Version="1">25165328</PMID>
      <PMID Version="1">25165312</PMID>
      <PMID Version="1">25165305</PMID>
      <PMID Version="1">25165259</PMID>
      <PMID Version="1">25165229</PMID>
      <PMID Version="1">25165197</PMID>
```




 */
@Program(name="pubmedfilterjs",
description="Filters Pubmed XML with a javascript  (java rhino) expression. Context contain 'article' a  PubmedBookArticle or a PubmedArticle and 'index', the index in the XML file.",
keywords={"pubmed","javascript","xml","ncbi"})
public class PubmedFilterJS
	extends Launcher
	{
	private static final Logger LOG = Logger.build(PubmedFilterJS.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outFile=null;

	@Parameter(names={"-f","--scriptfile"},description="Javascript file")
	private Path javascriptFile=null;
	@Parameter(names={"-e","--expression"},description="Javascript expression")
	private String javascriptExpr=null;

	
	@SuppressWarnings("unused")
	private static ObjectFactory _fool_javac=null;
	
	public PubmedFilterJS()
		{
		
		}
	
	
	

	@Override
	public int doWork(List<String> args) {
		final String inputName= oneFileOrNull(args);
		Unmarshaller unmarshaller;
		Marshaller marshaller;
		try
			{
			if(javascriptFile==null && StringUtils.isBlank(javascriptExpr)) {
				LOG.error("-e or -f must be defined");
				return -1;
				}
			if(javascriptFile!=null && !StringUtils.isBlank(javascriptExpr)) {
				LOG.error("-e or -f are bot be defined");
				return -1;
				}
			
			final JexlEngine jexl = new JexlEngine();
	        jexl.setSilent(false); // will throw errors now for selects that don't evaluate properly
	        jexl.setLenient(false);
	        jexl.setDebug(false);
			final Expression expression= jexl.createExpression(
						this.javascriptFile!=null?IOUtils.slurpPath(this.javascriptFile):
						this.javascriptExpr
						);

			
		
						
			final JAXBContext jc = JAXBContext.newInstance("gov.nih.nlm.ncbi.pubmed");
			unmarshaller =jc.createUnmarshaller();
			marshaller =jc.createMarshaller();

			
			final XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
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
			marshaller.setProperty(Marshaller.JAXB_FRAGMENT, true);
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

			
			try(PrintWriter pw= openPathOrStdoutAsPrintWriter(this.outFile)) {
				final XMLOutputFactory xof=XMLOutputFactory.newFactory();
				final  XMLEventWriter w=xof.createXMLEventWriter(pw);
				
				
				StreamSource src=null;
				if(inputName==null)
					{
					LOG.info("Reading stdin");
					src=new StreamSource(System.in);
					}
				else 
					{
					LOG.info("Reading file");
					src=new StreamSource(new File(inputName));
					}
				
				final XMLEventReader r=xmlInputFactory.createXMLEventReader(src);
				
				final  XMLEventFactory eventFactory=XMLEventFactory.newFactory();
				
				
				long nArticles=0L;
				while(r.hasNext())
					{
					XMLEvent evt=r.peek();
					switch(evt.getEventType())
						{
						case XMLEvent.START_ELEMENT:
							{
							String localName= evt.asStartElement().getName().getLocalPart();
							Object article=null;
							JAXBElement<?> jaxbElement=null;
							if(localName.equals("PubmedArticle"))
								{
								jaxbElement= unmarshaller.unmarshal(r,PubmedArticle.class);
								article=jaxbElement.getValue();
								}
							/* no more in the latest dtd else if(localName.equals("PubmedBookArticle"))
								{
								jaxbElement= unmarshaller.unmarshal(r,PubmedBookArticle.class);
								article=jaxbElement.getValue();
								} */
							else
								{
								w.add(r.nextEvent());
								break;
								}
							
							if(article!=null)
								{
								final JexlContext context = new MapContext();
								context.set("article", article);
								context.set("index", nArticles++);
								
								if(!(boolean)expression.evaluate(context))
									{
									break;
									}
								marshaller.marshal(jaxbElement, w);
								w.add(eventFactory.createCharacters("\n"));
								}
							
							break;
							}
						default:
							{
							w.add(r.nextEvent());
							break;
							}
						}
					r.close();
					}
				w.flush();
				w.close();
				pw.flush();
				}
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
	
		}
	
	public static void main(String[] args) {
		new PubmedFilterJS().instanceMainWithExit(args);
	}
}
