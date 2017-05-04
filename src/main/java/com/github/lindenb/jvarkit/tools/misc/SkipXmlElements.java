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
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.script.CompiledScript;
import javax.script.SimpleBindings;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="skipxmlelements",description="Filter XML elements with a javascript  (java rhino) expression. "
		+ "Context contain 'element' the current element. It implements" +
		"the interface Tag described in  SkipXmlElements.class")

public class SkipXmlElements
	extends Launcher
	{
	private static Logger LOG=Logger.build(SkipXmlElements.class).make();
	@Parameter(names="-e",description=" (js expression). Optional.")
	private String scriptExpr=null;
	@Parameter(names="-f",description=" (js file). Optional.")
	private File scriptFile=null;

	
	private static final int KEEP_ELEMENT=1;
	@SuppressWarnings("unused")
	private static final int SKIP_ELEMENT=0;
	private static final int KEEP_ELEMENT_AND_DESCENDANTS=2;
	
	
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
		public boolean localNameIn(String...list);
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
		public boolean localNameIn(String...list)
			{
			for(String s:list)
				{	
				if(s.equals(getLocalName())) return true;
				}
			return false;
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
	public int doWork(List<String> args)
		{
		

		String scriptExpr=null;
		File scriptFile=null;
		CompiledScript compiledScript=null;
		
		try
			{
			compiledScript = super.compileJavascript(scriptExpr, scriptFile);
			
			

			
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			
			PrintWriter pw=new PrintWriter(System.out);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			XMLEventWriter w=xof.createXMLEventWriter(pw);
			
			StreamSource src=null;
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
				LOG.error("illegal.number.of.arguments");
				return -1;
				}
			XMLEventReader r=xmlInputFactory.createXMLEventReader(src);
						
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
						
						int keep=1;
						
							
						bindings.put("element", curr);
						Object result=compiledScript.eval(bindings);
						
						if(result==null) 
							{
							throw new RuntimeException("User's Script returned null");
							}
						else if((result instanceof Boolean))
							{
							keep=(Boolean.class.cast(result).booleanValue()?1:0);
							}
						else if(!(result instanceof Number))
							{
							throw new RuntimeException("User's Script didn't return a number.");
							}
						else// if(result instanceof Number)
							{
							keep = ((Number)result).intValue();
							}
						
						
						if(keep==KEEP_ELEMENT)
							{
							w.add(r.nextEvent());
							}
						else /* skip this element or keep and descendant */
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
								if(keep==KEEP_ELEMENT_AND_DESCENDANTS)
									{
									w.add(evt);
									}
								if(depth==0) break;
								}
							}
						
						break;
						}
					case XMLEvent.COMMENT:
						{
						r.nextEvent();//just consumme
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
			LOG.error(err);
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
