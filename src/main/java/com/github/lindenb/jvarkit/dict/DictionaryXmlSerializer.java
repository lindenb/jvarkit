/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.dict;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * read/write dict from/to XML
 * @author lindenb
 *
 */
public class DictionaryXmlSerializer {
	
	
	private SAMSequenceRecord readSAMSequenceRecord(final XMLEventReader r,StartElement root)throws XMLStreamException {
		Attribute att=root.getAttributeByName(new QName("name"));
		if(att==null) 	throw new XMLStreamException("@name expected", root.getLocation());
		final String name=att.getValue();
		 att=root.getAttributeByName(new QName("length"));
		if(att==null) 	throw new XMLStreamException("@length expected", root.getLocation());
		final int length = Integer.parseInt(att.getValue());
		if(length<0) throw new XMLStreamException("positive @length expected", root.getLocation());
		
		final SAMSequenceRecord ssr=new SAMSequenceRecord(name, length);
		final Iterator<Attribute> it=root.getAttributes();
		while (it.hasNext()) {
			att = it.next();
			ssr.setAttribute(att.getName().getLocalPart(), att.getValue());
			}
		
		while(r.hasNext()) {
			final XMLEvent evt=r.nextEvent();
			if(evt.isStartElement()) {
				final StartElement se=evt.asStartElement();
				if(se.getName().getLocalPart().equals("alias")) {
					ssr.addAlternativeSequenceName(r.getElementText().trim());
					}
				else
					{
					throw new XMLStreamException("<alias> element expected", evt.getLocation());
					}
				}
			else if(evt.isEndElement()) {
				break;
				}
			else if(evt.isCharacters()) {
				if(!StringUtils.isBlank(evt.asCharacters().getData())) throw new XMLStreamException("non empty characters", evt.getLocation());
				}
			}
		return ssr;
		}
	
	
	public SAMSequenceDictionary readDictionary(final XMLEventReader r,StartElement root)throws XMLStreamException {
		if(!root.getName().getLocalPart().equals("dictionary")) throw new XMLStreamException("<dictionary> element expected", root.getLocation());
		final List<SAMSequenceRecord> L=new  ArrayList<>();
		while(r.hasNext()) {
			XMLEvent evt=r.nextEvent();
			if(evt.isStartElement()) {
				final StartElement se=evt.asStartElement();
				if(se.getName().getLocalPart().equals("sequence")) {
					L.add(readSAMSequenceRecord(r,se));
					}
				else
					{
					throw new XMLStreamException("<sequence> element expected", evt.getLocation());
					}
				}
			else if(evt.isEndElement()) {
				break;
				}
			else if(evt.isCharacters()) {
				if(!StringUtils.isBlank(evt.asCharacters().getData())) throw new XMLStreamException("non empty characters", evt.getLocation());
				}
			}
		return new SAMSequenceDictionary(L);
		}
	
	/** give a chance to add some attributes to the dictionary header */
	protected void writeOtherAttributes(final XMLStreamWriter w,final SAMSequenceDictionary dict) throws XMLStreamException {
		
		}
	
	public void writeDictionary(final XMLStreamWriter w,final SAMSequenceDictionary dict) throws XMLStreamException {
		final long genome_length = dict.getReferenceLength();
		w.writeStartElement("dictionary");
		w.writeAttribute("md5", dict.md5());
		w.writeAttribute("length", String.valueOf(genome_length));
		w.writeAttribute("count", String.valueOf(dict.size()));
		writeOtherAttributes(w,dict);
		
		final Optional<String> build=SequenceDictionaryUtils.getBuildName(dict);
		if(build.isPresent()) {
			w.writeAttribute("build",build.get());
			}
		
		long offset=0L;
		for(SAMSequenceRecord ssr:dict.getSequences()) {
			w.writeStartElement("sequence");
			w.writeAttribute("name", ssr.getSequenceName());
			w.writeAttribute("length", String.valueOf(ssr.getSequenceLength()));
			w.writeAttribute("index", String.valueOf(ssr.getSequenceIndex()));
			w.writeAttribute("offset", String.valueOf(offset));
			if(genome_length>0 ) {
				w.writeAttribute("f1", String.valueOf((offset)/(double)genome_length));
				w.writeAttribute("f2", String.valueOf((offset+ssr.getSequenceLength())/(double)genome_length));
				}
			for(Map.Entry<String, String> kv :ssr.getAttributes()) {
				if(kv.getKey().equals("name") || kv.getKey().equals(SAMSequenceRecord.SEQUENCE_NAME_TAG)) continue;
				if(kv.getKey().equals("length") || kv.getKey().equals(SAMSequenceRecord.SEQUENCE_LENGTH_TAG)) continue;
				if(kv.getKey().equals("index")) continue;
				if(kv.getKey().equals("offset")) continue;
				w.writeAttribute(kv.getKey(), kv.getValue());
				}
			final Set<String> aliases = ssr.getAlternativeSequenceNames();
			if(!aliases.isEmpty()) {
				for(String a:aliases) {
					w.writeStartElement("alias");
					w.writeCharacters(a);
					w.writeEndElement();
					}
				}
			w.writeEndElement();
			offset+=ssr.getSequenceLength();
			}
		w.writeEndElement();
		}
	}
