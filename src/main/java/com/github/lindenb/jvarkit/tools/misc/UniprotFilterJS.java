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
package com.github.lindenb.jvarkit.tools.misc;


import org.uniprot.*;

import com.github.lindenb.jvarkit.util.command.Command;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;

import javax.script.CompiledScript;
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

/**
 * PubmedFilterJS
 *
 */
public class UniprotFilterJS
	extends AbstractUniprotFilterJS
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(UniprotFilterJS.class);

	@SuppressWarnings("unused")
	private static ObjectFactory _fool_javac=null;
	
	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractUniprotFilterJS.AbstractUniprotFilterJSCommand
	 	{		
		
		public Collection<Throwable> call() throws Exception
			{
			final List<String> args = getInputFiles();
			Unmarshaller unmarshaller;
			Marshaller marshaller;
			try
				{
				CompiledScript compiledScript = compileJavascript();
	
				
				
				JAXBContext jc = JAXBContext.newInstance("org.uniprot");
				
				unmarshaller =jc.createUnmarshaller();
				marshaller =jc.createMarshaller();
	
				
				XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
				xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
				xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
				xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
				xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
				marshaller.setProperty(Marshaller.JAXB_FRAGMENT, true);
				marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
				marshaller.setProperty(Marshaller.JAXB_SCHEMA_LOCATION, "http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd");
				
				PrintWriter pw=new PrintWriter(System.out);
				XMLOutputFactory xof=XMLOutputFactory.newFactory();
				xof.setProperty(XMLOutputFactory.IS_REPAIRING_NAMESPACES, Boolean.FALSE);
				XMLEventWriter w=xof.createXMLEventWriter(pw);
				
				
				
				StreamSource src=null;
				if(args.isEmpty())
					{
					LOG.info("Reading stdin");
					src=new StreamSource(stdin());
					}
				else if(args.size()==1)
					{
					LOG.info("Reading file");
					src=new StreamSource(new File(args.get(0)));
					}
				else
					{
					return wrapException(getMessageBundle("illegal.number.of.arguments"));
					}
				XMLEventReader r=xmlInputFactory.createXMLEventReader(src);
				
				XMLEventFactory eventFactory=XMLEventFactory.newFactory();
				
				SimpleBindings bindings=new SimpleBindings();
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
								Object result=compiledScript.eval(bindings);
								
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
									LOG.warn("Script returned something that is not a boolean or a number:"+result.getClass());
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
				return RETURN_OK;
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				
				}
			}
	 	}
	
	public static void main(String[] args) {
		new UniprotFilterJS().instanceMainWithExit(args);
	}
}
