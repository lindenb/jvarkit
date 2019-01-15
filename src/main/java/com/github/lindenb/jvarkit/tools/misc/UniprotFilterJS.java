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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;


import org.uniprot.*;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.SimpleBindings;
import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
BEGIN_DOC

## Example
the following script get the human (id=9606) uniprot entries having an id in ensembl:

```javascript
function accept(e)
	{
	var ok=0,i;
	// check organism is human 
	if(e.getOrganism()==null) return false;
	var L= e.getOrganism().getDbReference();
	if(L==null) return false;
	for(i=0;i<L.size();++i)
		{
		if(L.get(i).getId()=="9606") {ok=1;break;}
		}
	if(ok==0) return false;
	ok=0;
	L= e.getDbReference();
	if(L==null) return false;
	for(i=0;i<L.size();++i)
		{
		if(L.get(i).getType()=="Ensembl") {ok=1;break;}
		}
	return ok==1;
	}
accept(entry);
```


```bash
$   curl -skL "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz" | gunzip -c |\
java -jar dist/uniprotfilterjs.jar  -f filter.js > output.xml
```



END_DOC

 */


@Program(name="uniprotfilterjs",
	description= "Filters Uniprot DUMP+ XML with a javascript  (java rhino) expression. "
		+ "Context contain 'entry' an uniprot entry "
		+ "and 'index', the index in the XML file.",
		keywords={"unitprot","javascript","xjc","xml"})
public class UniprotFilterJS
	extends Launcher
	{
	private static Logger LOG=Logger.build(UniprotFilterJS.class).make();

	@SuppressWarnings("unused")
	private static ObjectFactory _fool_javac=null;
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names="-e",description=" (js expression). Optional.")
	private String scriptExpr=null;
	@Parameter(names="-f",description=" (js file). Optional.")
	private File scriptFile=null;

	
	private UniprotFilterJS()
		{
		
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		Unmarshaller unmarshaller;
		Marshaller marshaller;
		try
			{
			final ScriptEngineManager manager = new ScriptEngineManager();
			final ScriptEngine engine = manager.getEngineByName("js");
			if(engine==null)
				{
				LOG.error("not available: javascript. Use the SUN/Oracle JDK ?");
				return -1;
				}
			final CompiledScript compiledScript = super.compileJavascript(scriptExpr, scriptFile);
			
			
			
			final JAXBContext jc = JAXBContext.newInstance("org.uniprot");
			
			unmarshaller =jc.createUnmarshaller();
			marshaller =jc.createMarshaller();

			
			final XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			marshaller.setProperty(Marshaller.JAXB_FRAGMENT, true);
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
			marshaller.setProperty(Marshaller.JAXB_SCHEMA_LOCATION, "http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd");
			
			final PrintWriter pw= super.openFileOrStdoutAsPrintWriter(this.outputFile);
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			xof.setProperty(XMLOutputFactory.IS_REPAIRING_NAMESPACES, Boolean.FALSE);
			final XMLEventWriter w=xof.createXMLEventWriter(pw);
			
			
			
			final StreamSource src;
			if(args.isEmpty())
				{
				src=new StreamSource(stdin());
				}
			else if(args.size()==1)
				{
				src=new StreamSource(new File(args.get(0)));
				}
			else
				{
				LOG.error("Illegal number of arguments");
				return -1;
				}
			final XMLEventReader r=xmlInputFactory.createXMLEventReader(src);
			
			final XMLEventFactory eventFactory=XMLEventFactory.newFactory();
			
			final SimpleBindings bindings=new SimpleBindings();
			long nArticles=0L;
			while(r.hasNext())
				{
				XMLEvent evt=r.peek();
				switch(evt.getEventType())
					{
					case XMLEvent.START_ELEMENT:
						{
						StartElement sE= evt.asStartElement();
						Entry entry=null;
						JAXBElement<Entry> jaxbElement=null;
						if(sE.getName().getLocalPart().equals("entry") )
							{
							jaxbElement= unmarshaller.unmarshal(r,Entry.class);
							entry= jaxbElement.getValue();
							}
						else
							{
							w.add(r.nextEvent());
							break;
							}
						
						if(entry!=null)
							{
							
							bindings.put("entry", entry);
							bindings.put("index", nArticles++);
							final Object result=compiledScript.eval(bindings);
							
							String entryName="undefined";
							if(!entry.getAccession().isEmpty()) entryName=entry.getAccession().get(0);
							
							if(result==null) break;
							if(result instanceof Boolean)
								{
								if(Boolean.FALSE.equals(result))
									{
									w.add(eventFactory.createComment(" "+entryName+" "));
									w.add(eventFactory.createCharacters("\n"));
									break;
									}
								}
							else if(result instanceof Number)
								{
								if(((Number)result).intValue()!=1) 
									{
									w.add(eventFactory.createComment(" "+entryName+" "));
									w.add(eventFactory.createCharacters("\n"));
									break;
									}
								}
							else
								{
								LOG.warning("Script returned something that is not a boolean or a number:"+result.getClass());
								break;
								}
							
							marshaller.marshal(
									jaxbElement,
									w);
							w.add(eventFactory.createCharacters("\n"));
							}
						
						break;
						}
					case XMLEvent.SPACE:break;
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
			
			}
		}
	
	public static void main(final String[] args) {
		new UniprotFilterJS().instanceMainWithExit(args);
	}
}
