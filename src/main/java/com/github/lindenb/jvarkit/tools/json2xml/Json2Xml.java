/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.json2xml;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.file.Path;
import java.util.List;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;

/**
BEGIN_DOC

## Example

```
$ echo '[null,true,{"hello":1}]' | java -jar dist/jvarkit.jar json2xml | xmllint --format -
<?xml version="1.0" encoding="UTF-8"?>
<array xmlns="http://www.ibm.com/xmlns/prod/2009/jsonx">
  <null/>
  <boolean>true</boolean>
  <object>
    <number name="hello">1</number>
  </object>
</array>

$ echo '[null,true,{"hello":1}]' | java -jar dist/jvarkit.jar json2xml --ns "" --omit-xml-declaration 
<array><null/><boolean>true</boolean><object><number name="hello">1</number></object></array>

```



END_DOC
 */

@Program(name="json2xml",
description="convert json to xml",
keywords={"json","xml"},
creationDate="20251226",
modificationDate="20251226",
jvarkit_amalgamion = true,
generate_doc = true
)
public final class Json2Xml extends Launcher {
	private static final Logger LOG = Logger.of(Json2Xml.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--namespace","--ns","-n"},description="xml namespace (empty= no ns)")
	private String namespaceuri="http://www.ibm.com/xmlns/prod/2009/jsonx";
	@Parameter(names={"--omit-xml-declaration"},description="Omit XML declaration")
	private boolean omit_xml_decl = false;
	
	
	
	private void writeStartElement(final XMLStreamWriter w,String tag) throws XMLStreamException {
		if(StringUtils.isBlank(this.namespaceuri)) {
			w.writeStartElement(tag);
			}
		else
			{
			w.writeStartElement(this.namespaceuri, tag);
			}
		}
	
	private void writeEmptyElement(final XMLStreamWriter w,String tag) throws XMLStreamException {
		if(StringUtils.isBlank(this.namespaceuri)) {
			w.writeEmptyElement(tag);
			}
		else
			{
			w.writeEmptyElement(this.namespaceuri, tag);
			}
		}
	
	private void parseObject(final XMLStreamWriter w,final String label,final JsonReader r) throws XMLStreamException,IOException
		{
		writeStartElement(w, "object");
		if(label!=null) w.writeAttribute("name", label);
		for(;;)
			{
			if(r.peek()==JsonToken.END_OBJECT) 
				{
				w.writeEndElement();
				r.endObject();		
				break;
				}
			if(r.peek()!=JsonToken.NAME) throw new IllegalStateException(r.peek().name());
			final String s = r.nextName();
			parse(w,s,r);
			}
		}
	private void parseArray(final XMLStreamWriter w,final String label,final JsonReader r) throws XMLStreamException,IOException
		{
		
		writeStartElement(w,"array");
		if(label!=null) w.writeAttribute("name", label);
		for(;;)
			{
			if(r.peek()==JsonToken.END_ARRAY) 
				{
				w.writeEndElement();
				r.endArray();		
				break;
				}
			parse(w,null,r);
			}
		}
	
	
	private void parse(final XMLStreamWriter w,final String label,final JsonReader r) throws XMLStreamException,IOException
		{
		if(!r.hasNext()) return;
		    JsonToken token= r.peek();
		    switch(token)
		    	{
		    	case END_OBJECT://through
		    	case END_ARRAY://through
		    	case NAME: throw new IllegalStateException("unexpected "+ token);
		    	case BEGIN_OBJECT:
		    		{
		    		r.beginObject();
		    		parseObject(w,label,r);	
		    		break;
		    		}
		    	case BEGIN_ARRAY:
		    		{
		    		r.beginArray();
		    		parseArray(w,label,r);
		    		break;
		    		}
		    	case NULL:
		    		{
		    		r.nextNull();
		    		writeEmptyElement(w, "null");
		    		if(label!=null) w.writeAttribute("name", label);
		    		break;
		    		}
		    	case STRING:
		    		{
		    		writeStartElement(w, "string");
		    		if(label!=null) w.writeAttribute("name", label);
		    		w.writeCharacters(r.nextString());
		    		w.writeEndElement();
		    		break;
		    		}
		    	case NUMBER:
		    		{
		    		writeStartElement(w, "number");
		    		if(label!=null) w.writeAttribute("name", label);
		    		String s;
		    		try
		    			{
		    			s= String.valueOf(r.nextLong());
		    			}
		    		catch(Exception err)
		    			{
		    			s= String.valueOf(r.nextDouble());
		    			}

		    		
		    		w.writeCharacters(s);
		    		w.writeEndElement();
		    		break;
		    		}
		    	case BOOLEAN:
		    		{
		    		writeStartElement(w, "boolean");
		    		if(label!=null) w.writeAttribute("name", label);
		    		w.writeCharacters(String.valueOf(r.nextBoolean()));
		    		w.writeEndElement();
		    		break;
		    		}
		    	case END_DOCUMENT:
		    		{
		    		break;
		    		}
		    	default: throw new IllegalStateException(token.name());
		    	}
			
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		try {
			final String input = super.oneFileOrNull(args);
			try(Reader r=input==null || input.equals("-")?new InputStreamReader(System.in):IOUtils.openURIForBufferedReading(input)) {
			try(JsonReader jr = new JsonReader(r)) {
				jr.setLenient(true);
				final XMLOutputFactory xof = XMLOutputFactory.newFactory();
				xof.setProperty(XMLOutputFactory.IS_REPAIRING_NAMESPACES, Boolean.TRUE);
				final XMLStreamWriter w = xof.createXMLStreamWriter(System.out,"UTF-8");
				if(!StringUtils.isBlank(this.namespaceuri)) w.setDefaultNamespace(this.namespaceuri);
				if(!omit_xml_decl) w.writeStartDocument("UTF-8", "1.0");
				this.parse(w,null,jr);
				w.writeEndDocument();
				w.flush();
				w.close();
				}
			}
			return 0;
		} catch (Throwable err) {
			LOG.error(err);
			return -1;
		}

	}
	
	
	public static void main(final String[] args)
		{
		new Json2Xml().instanceMainWithExit(args);
		}
}
