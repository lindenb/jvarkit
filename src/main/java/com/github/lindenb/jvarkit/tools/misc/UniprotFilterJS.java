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

import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;

import javax.script.Compilable;
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

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;



/**
 * PubmedFilterJS
 *
 */
public class UniprotFilterJS
	extends AbstractCommandLineProgram
	{
	
	@SuppressWarnings("unused")
	private static ObjectFactory _fool_javac=null;
	
	private UniprotFilterJS()
		{
		
		}
	@Override
	public String getProgramDescription() {
		return "Filters Uniprot DUMP+ XML with a javascript  (java rhino) expression. "
				+ "Context contain 'entry' an uniprot entry "
				+ "and 'index', the index in the XML file.";
		}			
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/UniprotFilterJS";
		}
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -e (js expression). Optional.");
		out.println(" -f (js file). Optional.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		

		String scriptExpr=null;
		File scriptFile=null;
		CompiledScript compiledScript=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"e:f:"))!=-1)
			{
			switch(c)
				{
				case 'e': scriptExpr=opt.getOptArg();break;
				case 'f': scriptFile=new File(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		Unmarshaller unmarshaller;
		Marshaller marshaller;
		try
			{
			ScriptEngineManager manager = new ScriptEngineManager();
			ScriptEngine engine = manager.getEngineByName("js");
			if(engine==null)
				{
				error("not available: javascript. Use the SUN/Oracle JDK ?");
				return -1;
				}
			Compilable compilingEngine = (Compilable)engine;
			if(scriptFile!=null)
				{
				info("Compiling "+scriptFile);
				FileReader r=new FileReader(scriptFile);
				compiledScript=compilingEngine.compile(r);
				r.close();
				}
			else if(scriptExpr!=null)
				{
				info("Compiling "+scriptExpr);
				compiledScript=compilingEngine.compile(scriptExpr);
				}
			else
				{
				error("Undefined script");
				return -1;
				}

			
			
			JAXBContext jc = JAXBContext.newInstance("org.uniprot");
			unmarshaller =jc.createUnmarshaller();
			marshaller =jc.createMarshaller();
			

			
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.FALSE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			marshaller.setProperty(Marshaller.JAXB_FRAGMENT, true);
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

			
			PrintWriter pw=new PrintWriter(System.out);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			xof.setProperty(XMLOutputFactory.IS_REPAIRING_NAMESPACES, Boolean.TRUE);
			XMLEventWriter w=xof.createXMLEventWriter(pw);
			w.setDefaultNamespace("http://uniprot.org/uniprot");
			
			StreamSource src=null;
			if(opt.getOptInd()==args.length)
				{
				info("Reading stdin");
				src=new StreamSource(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				info("Reading file");
				src=new StreamSource(new File(args[opt.getOptInd()]));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
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
						Object entry=null;
						JAXBElement<?> jaxbElement=null;
						if(sE.getName().getLocalPart().equals("entry") )
							{
							jaxbElement= unmarshaller.unmarshal(r,Entry.class);
							entry=jaxbElement.getValue();
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
							
							if(result==null) break;
							if(result instanceof Boolean)
								{
								if(Boolean.FALSE.equals(result)) break;
								}
							else if(result instanceof Number)
								{
								if(((Number)result).intValue()!=1) break;
								}
							else
								{
								warning("Script returned something that is not a boolean or a number:"+result.getClass());
								break;
								}
							marshaller.marshal(jaxbElement, w);
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
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			
			}
		}
	
	public static void main(String[] args) {
		new UniprotFilterJS().instanceMainWithExit(args);
	}
}
