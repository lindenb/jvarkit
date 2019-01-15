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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.namespace.QName;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
/**
 BEGIN_DOC
 
 ## Example
 
```
 $ java -jar dist/knime2txt.jar -g out.dot  knime_3.3.2/workspace/TEST/  | xmllint --format - > out.html
 [INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/workflow.knime
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Interactive Table (#3)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/File Reader (#4)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Java Snippet Row Filter (#5)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Extract Column Header (#6)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Unpivoting (#7)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Java Snippet _simple_ (#8)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/Column Filter (#10)/settings.xml
[INFO][KnimeToText]read XML: knime_3.3.2/workspace/TEST/CSV Writer (#9)/settings.xml
```
 
```
$ xmllint out.html  | head
<?xml version="1.0"?>
<div>
  <div>
    <p>Interactive Table (#3)/settings.xml</p>
    <ul>
      <li><b>factory</b>  value:"<i>org.knime.base.node.viz.table.TableNodeFactory</i>"</li>
      <li><b>customDescription</b>  value:"<i/>"</li>
    </ul>
  </div>
  <div>
(...)
```
 
 
 
 
 END_DOC
 */
@Program(name="knime2txt",
		description="converts a Knime Workflow to a html representation.",
		keywords={"knime","workflow","convert"}
		)
public class KnimeToText extends Launcher{
	private static final String NS="http://www.knime.org/2008/09/XMLConfig";
	private static final Logger LOG = Logger.build(KnimeToText.class).make();

	
	private abstract class KnimeNode
		{
		ConfigNode parent=null;
		String key=null;
		ConfigNode asConfig() {
			return ConfigNode.class.cast(this);
			}
		EntryNode asEntry() {
			return EntryNode.class.cast(this);
			}
		boolean has(final String key)
			{
			if(!(this instanceof ConfigNode)) return false;
			return asConfig().entries.containsKey(key);
			}
		KnimeNode get(final String key)
			{
			KnimeNode n= asConfig().entries.get(key);
			if(n==null) throw new IllegalArgumentException("Cannot find "+key+" in this node. available are:"+asConfig().entries.keySet());
			return n;
			}
		EntryNode entry(final String key)
			{
			return get(key).asEntry();
			}
		ConfigNode config(final String key)
			{
			return get(key).asConfig();
			}
		List<KnimeNode> all()
			{
			return asConfig().entries.values().stream().collect(Collectors.toList());
			}
		abstract void read(StartElement root,XMLEventReader r) throws XMLStreamException;
		Document toDom()  throws Exception
			{
			Document dom=DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
			toDom(dom);
			return dom;
			}
		abstract void toDom(Node parent);
		int depth()
			{
			return (this.parent==null?0:1+this.parent.depth());
			}
		abstract void html(XMLStreamWriter w) throws XMLStreamException;
		boolean isHidden()
			{
			if(depth()<=1)
				{
				if(this.key.startsWith("node-")) return true;
				if(this.key.equals("isInactive")) return true;
				if(this.key.equals("internalObjects")) return true;
				if(this.key.equals("filestores")) return true;
				if(this.key.equals("name")) return true;
				if(this.key.equals("state")) return true;
				if(this.key.equals("internal_node_subsettings")) return true;
				if(this.key.equals("node_file")) return true;
				if(this.key.equals("ports")) return true;
				if(this.key.equals("hasContent")) return true;
				}
			if(this.parent!=null && this.parent.key.equals("nodeAnnotation")) 
				{
				if(this.key.equals("text")) return false;
				return true;
				}
			
			return false;
			}
		String getLabel()
			{
			if(!has("nodeAnnotation")) return null;
			ConfigNode nodeAnnotation = config("nodeAnnotation");
			if(!nodeAnnotation.has("text")) return null;
			return nodeAnnotation.entry("text").value;
			}
		}
	
	private class ConfigNode extends KnimeNode
		{
		final Map<String, KnimeNode> entries=new HashMap<>();
		@Override void read(StartElement root,XMLEventReader r) throws XMLStreamException
			{
			Attribute att=root.getAttributeByName(new QName("key"));
			if(att==null) throw new IllegalStateException("@key missing in entry");
			this.key = att.getValue();
			while(r.hasNext())
				{
				final XMLEvent evt=r.nextEvent();
				if(evt.isEndElement())
					{
					break;
					}
				else if(evt.isStartElement())
					{
					StartElement E = evt.asStartElement();
					String localname= E.getName().getLocalPart();
				
					if(localname.equals("config"))
						{
						final ConfigNode child = new ConfigNode();
						child.parent = this;
						child.read(E, r);
						this.entries.put(child.key, child);
						}
					else if(localname.equals("entry"))
						{
						final EntryNode child = new EntryNode();
						child.parent = this;
						child.read(E, r);
						this.entries.put(child.key, child);
						}
					else
						{
						throw new IllegalStateException("boum");
						}
					}
				}
			}
		@Override void toDom(final Node parent)
			{
			Element root=parent.getOwnerDocument().createElementNS(NS, "config");
			parent.appendChild(root);
			if(key!=null) root.setAttribute("key",key);

			for(KnimeNode kn:this.entries.values())
				{
				kn.toDom(root);
				}
			}
		
		@Override void html(XMLStreamWriter w) throws XMLStreamException
			{
			if(this.entries.isEmpty() || isHidden()) return;
			
			if(!this.entries.values().stream().filter(E->!E.isHidden()).findAny().isPresent()) return;
			
			if(this.parent==null)
				{
				w.writeStartElement("ul");
				}
			else
				{
				w.writeStartElement("li");
				w.writeStartElement("b");
				w.writeCharacters(this.key);
				w.writeEndElement();//b
				w.writeCharacters(":");
				w.writeStartElement("ul");
				}
		
			
			for(final KnimeNode kn:this.entries.values())
				{
				if(kn.isHidden()) continue;
				kn.html(w);
				}
			
		
			if(this.parent==null)
				{
				w.writeEndElement();
				}
			else
				{
				w.writeEndElement();
				w.writeEndElement();
				}
			}
		@Override
		boolean isHidden() {
			if(super.isHidden()) return true;
			return false;
			}
		}
	 class EntryNode extends KnimeNode
		{
		String type=null;
		String value=null;
		Boolean is_null=null;
		@Override void read(final StartElement root,final XMLEventReader r) throws XMLStreamException
			{
			Attribute att=root.getAttributeByName(new QName("key"));
			if(att==null) throw new IllegalStateException("@key missing in entry");
			this.key = att.getValue();
			att=root.getAttributeByName(new QName("type"));
			if(att!=null) this.type=att.getValue();
			att=root.getAttributeByName(new QName("value"));
			if(att!=null) this.value=att.getValue();

			att=root.getAttributeByName(new QName("is_null"));
			if(att!=null) this.is_null=Boolean.valueOf(att.getValue());
			while(r.hasNext())
				{
				final XMLEvent evt=r.nextEvent();
				if(evt.isEndElement())
					{
					break;
					}
				else if(evt.isStartElement())
					{
					throw new XMLStreamException("Illegal state "+evt.asStartElement().getName(), root.getLocation());
					}
				}
			}
		@Override void toDom(Node parent)
			{
			Element root=parent.getOwnerDocument().createElementNS(NS, "entry");
			parent.appendChild(root);
			if(key!=null) root.setAttribute("key",key);
			if(type!=null) root.setAttribute("type",type);
			if(is_null!=null) root.setAttribute("is_null",String.valueOf(is_null));
			}
		@Override void html(XMLStreamWriter w) throws XMLStreamException
			{
			if(this.isHidden()) return;
			w.writeStartElement("li");
			w.writeStartElement("b");
			w.writeCharacters(this.key);
			w.writeEndElement();//b
			w.writeCharacters(" ");
			
			/*
			if(type!=null) 
				{
				w.writeCharacters("  ");
				w.writeStartElement("i");
				w.writeCharacters("type");
				w.writeEndElement();
				w.writeCharacters(":\"");
				w.writeCharacters(type);
				w.writeCharacters("\"");
				}*/
			if(is_null!=null) {
				w.writeCharacters(" null:\"");
				w.writeCharacters(String.valueOf(is_null));
				w.writeCharacters("\"");
				}
			if(is_null==null || !is_null.booleanValue())
				{
				boolean escaped=false;
				String s=this.value==null?"":this.value;
				final Pattern pat=Pattern.compile("%%[0-9][0-9][0-9][0-9][0-9]");
				for(;;)
					{
					final Matcher matcher=pat.matcher(s);
					if(!matcher.find()) break;
					escaped=true;
					int ascii = Integer.parseInt(s.substring(matcher.start()+2,matcher.end()));  
					String sep =  ""+(char)((byte)ascii);
					if(ascii==0 ) sep="";
					
					s = s.substring(0,matcher.start())+ sep + s.substring(matcher.end());
					}
				if(!escaped)
					{
					w.writeCharacters(" value");
					w.writeCharacters(":\"");
					w.writeStartElement("i");
					w.writeCharacters(value);
					w.writeEndElement();
					w.writeCharacters("\"");
					}
				else
					{
					w.writeStartElement("pre");
					for(final String line: s.split("[\n]"))
						{
						w.writeCharacters(line);
						w.writeCharacters("\n");
						}
					w.writeEndElement();
					}
				
				}
			w.writeEndElement();//li
			}
		@Override
		boolean isHidden() {
			if(super.isHidden()) return true;
			return false;
			}
		}

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-g","--dot"},description="graphiz dot file.")
	private File dotFile = null;

	private ConfigNode readXmlConfig(final File settingsXml) throws Exception
		{
		LOG.info("read XML: "+settingsXml);
		ConfigNode cfg=null;
		IOUtil.assertFileIsReadable(settingsXml);
		XMLInputFactory xif = XMLInputFactory.newFactory();
		FileReader r= new FileReader(settingsXml);
		XMLEventReader xr=xif.createXMLEventReader(r);
		while(xr.hasNext())
			{
			final XMLEvent evt=xr.peek();
			if(!evt.isStartElement()) {xr.next();continue;}
			final StartElement E = xr.nextEvent().asStartElement();
			if(!E.getName().getLocalPart().equals("config"))
				{
				 throw new XMLStreamException("not <config> ?",evt.getLocation());
				}
			cfg = new ConfigNode();
			cfg.read(E,xr);
			}
		xr.close();
		r.close();
		if(cfg==null) throw new IOException("cannot get <config>");
		return cfg;
		}
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter pw=null;
		PrintWriter dotw = null;
		XMLStreamWriter w;
		try {
			final File baseDir = new File(super.oneAndOnlyOneFile(args));
			IOUtil.assertDirectoryIsReadable(baseDir);
			
			ConfigNode workflowconfig = readXmlConfig(new File(baseDir,"workflow.knime"));
			
			dotw = (this.dotFile==null?new PrintWriter(new NullOuputStream()):new PrintWriter(this.dotFile));
			pw = super.openFileOrStdoutAsPrintWriter(outputFile);
			w = XMLOutputFactory.newFactory().createXMLStreamWriter(pw);
			w.writeStartElement("div");
			dotw.println("digraph G {");

			
			for(final KnimeNode node:workflowconfig.get("nodes").all())
				{
				File node_settings_file=new File(baseDir,node.entry("node_settings_file").value);
				
				final ConfigNode nodeConfig = readXmlConfig(node_settings_file);
				w.writeStartElement("div");
				w.writeStartElement("p");
				w.writeCharacters(node.entry("node_settings_file").value);
				w.writeEndElement();
				nodeConfig.html(w);
				
				w.writeEndElement();//div
				dotw.print("n"+ node.entry("id").value +"[label=\"");
				String label= nodeConfig.getLabel();
				if(label!=null)
					{
					dotw.print(label.replace('\"', '\'') + " : ");
					}
				
				dotw.println(node.entry("node_settings_file").value.replace('\"', '\'') +
							"\"];");
				}
			for(final KnimeNode connection:workflowconfig.get("connections").all())
				{
				dotw.println("n"+ connection.entry("sourceID").value +" -> n" +  connection.entry("destID").value +";");
				}
			
			dotw.println("}");
			w.writeEndElement();
			w.flush();w.close();
			pw.flush();pw.close();pw=null;
			dotw.flush();dotw.close();

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
			CloserUtil.close(dotw);
			}
		}
	
	public static void main(String[] args) {
		new KnimeToText().instanceMainWithExit(args);
	}
	
}
