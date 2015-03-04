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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;


import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.SimpleBindings;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;



/**
 * PubmedFilterJS
 *
 */
public class SkipXmlElements
	extends AbstractCommandLineProgram
	{
	public static interface Tag
		{
		public Tag getParent();
		/** get depth root element=1 */
		public int getDepth();
		public QName getQName();
		public String getLocalName();
		public String getNamespaceURI();
		public String getPrefix();
		public List<Attribute> getAttributes();
		}
	
	private static class TagImpl
		implements Tag
		{
		TagImpl parent=null;
		QName qname;
		List<Attribute> attributes=new ArrayList<Attribute>();
		public TagImpl(StartElement start)
			{
			this.qname=start.getName();
			Iterator<?> iter=start.getAttributes();
			while(iter.hasNext()) attributes.add(Attribute.class.cast(iter.next()));
			}
		
		@Override
		public List<Attribute> getAttributes()
			{
			return attributes;
			}
		@Override
		public QName getQName()
			{
			return qname;
			}
		@Override
		public String getLocalName()
			{
			return getQName().getLocalPart();
			}
		@Override
		public String getNamespaceURI()
			{
			return getQName().getNamespaceURI();
			}
		@Override
		public String getPrefix()
			{
			return getQName().getPrefix();
			}

		@Override
		public Tag getParent()
			{
			return this.parent;
			}
		@Override
		public int getDepth()
			{
			return 1+(getParent()==null?0:getParent().getDepth());
			}
		@Override
		public String toString()
			{
			return this.qname.toString();
			}
		}
	
	private SkipXmlElements()
		{
		
		}
	@Override
	public String getProgramDescription() {
		return "Filter XML elements with a javascript  (java rhino) expression. "
				+ "Context contain 'element' the current element. It implements" +
				"the interface Tag described in "+SkipXmlElements.class.getName();
		}			
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/SkipXmlElements";
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

			
			
			

			
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			
			PrintWriter pw=new PrintWriter(System.out);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			XMLEventWriter w=xof.createXMLEventWriter(pw);
			
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
			TagImpl curr=null;
			while(r.hasNext())
				{
				XMLEvent evt=r.peek();
				switch(evt.getEventType())
					{
					case XMLEvent.START_ELEMENT:
						{
						if(curr==null)
							{
							curr=new TagImpl(evt.asStartElement());
							}
						else
							{
							TagImpl node=new TagImpl(evt.asStartElement());
							node.parent=curr;
							curr=node;
							}
						
						boolean keep=true;
							
						bindings.put("element", curr);
						Object result=compiledScript.eval(bindings);
						
						if(result==null) 
							{
							warning("Script returned null");
							}
						else if(result instanceof Boolean)
							{
							if(Boolean.FALSE.equals(result)) keep=false;
							}
						else if(result instanceof Number)
							{
							if(((Number)result).intValue()!=1)  keep=false;
							}
						else
							{
							warning("Script returned something that is not a boolean or a number:"+result.getClass());
							}
						
						if(keep)
							{
							w.add(r.nextEvent());
							}
						else /* skip this element */
							{
							curr=curr.parent;
							int depth=0;
							while(r.hasNext())
								{
								evt=r.nextEvent();
								switch(evt.getEventType())
									{
									case XMLEvent.START_ELEMENT:
										{
										depth++;
										break;
										}
									case XMLEvent.END_ELEMENT:
										{
										depth--;
										break;
										}
									default:break;
									}
								if(depth==0) break;
								}
							}
						
						break;
						}
					case XMLEvent.END_ELEMENT :
						{
						curr=curr.parent;
						w.add(r.nextEvent());
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
		new SkipXmlElements().instanceMainWithExit(args);
	}
}
