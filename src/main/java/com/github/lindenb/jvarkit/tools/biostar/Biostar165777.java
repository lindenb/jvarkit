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
* 2015 creation copied from http://plindenbaum.blogspot.fr/2010/11/blastxmlannotations.html

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import htsjdk.samtools.util.CloserUtil;


public class Biostar165777 extends AbstractBiostar165777
	{
	private static final String SPLIT_TOKEN="__SPLIT__";
	
	public Biostar165777() {
		}
	@Override
	public Collection<Throwable> call(String inputFile) throws Exception {
		if(getOutputFile()==null)
			{
			return wrapException("output file must be defined.");
			}
		if(!getOutputFile().getName().contains(SPLIT_TOKEN))
			{
			return wrapException("output file "+getOutputFile()+" should contains the word "+SPLIT_TOKEN);
			}
		if(getTagName()==null)
			{
			return wrapException("Tag name was not defined");
			}
		if(getCount()<1)
			{	
			return wrapException("bad count "+getCount());
			}
		XMLEventReader r = null;
		try {
			final List<FileOutputStream> fwriters = new ArrayList<>();
			final List<XMLEventWriter> writers = new ArrayList<>();
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			for(int i=0;i< getCount();++i)
				{
				File xmlout =new File(getOutputFile().getParentFile(), getOutputFile().getName().
						replaceAll(SPLIT_TOKEN, String.format("%03d", (i+1))));
				FileOutputStream fos = new FileOutputStream(xmlout);
				fwriters.add(fos);
				writers.add(xof.createXMLEventWriter(fos,"UTF-8"));
				}
			
			//READ Whole XML file
			
			final XMLInputFactory xif = XMLInputFactory.newFactory();
			if(inputFile==null)
				{
				r = xif.createXMLEventReader(new StreamSource(stdin())); 
				}
			else
				{
				r = xif.createXMLEventReader(new StreamSource(new File(inputFile))); 
				}
			int numberOfItems = 0;
			final XMLEventFactory xef = XMLEventFactory.newFactory();
			while(r.hasNext())
				{
				
				XMLEvent evt = r.nextEvent();
				if(evt.isStartElement() && evt.asStartElement().getName().getLocalPart().equals(getTagName()))
					{
					++numberOfItems;
					XMLEventWriter w = writers.get(numberOfItems%getCount());
					int depth=1;
					w.add(evt);
					while(r.hasNext())
						{
						evt=r.nextEvent();
						if(evt.isStartElement())
							{
							++depth;
							}
						else if(evt.isEndElement())
							{
							--depth;
							}
						w.add(evt);
						if(depth==0) break;
						}
					while(r.hasNext())
						{	
						evt=r.peek();
						if(!evt.isCharacters()) break;
						if(!evt.asCharacters().getData().trim().isEmpty())
							{
							break;
							}
						r.nextEvent();
						}
					continue;
					}
				if(evt.isStartDocument())
					{
					evt = xef.createStartDocument("UTF-8", "1.0");
					}
				for(XMLEventWriter w:writers) w.add(evt);
				}
			r.close();
			for(XMLEventWriter w: writers)
				{
				w.flush();
				CloserUtil.close(w);
				}
			for(FileOutputStream w:fwriters)
				{	
				CloserUtil.close(w);
				}
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally {
			
			}
		
		}
	public static void main(String[] args) {
		new Biostar165777().instanceMainWithExit(args);
	}
	}
