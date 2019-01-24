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
package com.github.lindenb.jvarkit.tools.blast;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Comparator;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import gov.nih.nlm.ncbi.blast.Iteration;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;

/**
 BEGIN_DOC
 
 ## Example
 
 ```
 $ java -jar dist/mergeblastxml.jar input1.blastn.xml  input2.blastn.xml  input2.blastn.xml > out.xml
 ``` 
 
 
 END_DOC
 */
@Program(name="mergeblastxml",biostars=246958,
description="merge XML blast results (same Iteration/Iteration_query-def in multiple xml files",
keywords={"blast","xml"}
)
public class MergeBlastXml extends Launcher {
private static final Logger LOG=Logger.build(MergeBlastXml.class).make();
private Unmarshaller unmarshaller;
private Marshaller marshaller;

@Parameter(names={"-o","--out"},description="Output SVG file or stdout")
private File outputFile=null;
@Parameter(names={"--tmpDir"},description="Tmp Directory")
private File tmpFile= IOUtils.getDefaultTmpDir();
@Parameter(names={"--maxRecordsInRam"},description="Max Records in RAM")
private int maxRecordsInRam=50000;

private class BlastIterationCodec extends AbstractDataCodec<Iteration>
	{
	@Override
	public Iteration decode(DataInputStream dis) throws IOException {
		String xmlStr;
		try
			{
			xmlStr=AbstractDataCodec.readString(dis);
			}
		catch(Exception err)
			{
			return null;
			}	
		try
			{
			return unmarshaller.unmarshal(new StreamSource(new StringReader(xmlStr)),Iteration.class).getValue();
			}
		catch(JAXBException err)
			{
			throw new RuntimeIOException(err);
			}	
		}
	@Override
	public void encode(DataOutputStream dos, Iteration hit) throws IOException {
		StringWriter sw=new StringWriter();
		try
			{
			MergeBlastXml.this.marshaller.marshal(hit,sw);
			AbstractDataCodec.writeString(dos, sw.toString());
			}
		catch(final JAXBException err)
			{
			LOG.error(err);
			throw new RuntimeIOException(err);
			}	
		}
	@Override
	public AbstractDataCodec<Iteration> clone() {
		return new BlastIterationCodec();
		}
	}

/* force javac to compile */
@SuppressWarnings("unused")
private gov.nih.nlm.ncbi.blast.ObjectFactory _ignore_for_javac=null;
	
@Override
	public int doWork(List<String> args) {
		if(args.isEmpty())
			{
			LOG.error("input xml missing");
			return -1;
			}
		XMLEventReader rx=null;
		XMLEventReader rx2=null;
		XMLEventWriter wx=null;
		SortingCollection<Iteration> sortingCollection=null;
		try {
			JAXBContext jc = JAXBContext.newInstance("gov.nih.nlm.ncbi.blast");
			this.unmarshaller=jc.createUnmarshaller();
			this.marshaller=jc.createMarshaller();
			this.marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT,true);
			this.marshaller.setProperty(Marshaller.JAXB_FRAGMENT,true);
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.FALSE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			xmlInputFactory.setXMLResolver(new XMLResolver()
				{
				@Override
				public Object resolveEntity(String arg0, String arg1, String arg2,
						String arg3) throws XMLStreamException
					{
					LOG.info("resolveEntity:" +arg0+"/"+arg1+"/"+arg2);
					return null;
					}
				});
			final Comparator<Iteration> hitComparator= (A,B) -> {
					return A.getIterationQueryDef().compareTo(B.getIterationQueryDef());
				} ;
			sortingCollection = SortingCollection.newInstance(Iteration.class, new BlastIterationCodec(),
					hitComparator, 
					this.maxRecordsInRam,
					this.tmpFile.toPath()
					);
			rx=xmlInputFactory.createXMLEventReader(new FileReader(args.get(0)));
			
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			if(this.outputFile!=null)
				{
				wx=xof.createXMLEventWriter(new StreamResult(this.outputFile));
				}
			else
				{
				wx=xof.createXMLEventWriter(new StreamResult(stdout()));
				}
			boolean in_iteration=false;
			while(rx.hasNext())
				{
				final XMLEvent evt=rx.peek();
				
				if(evt.isStartElement() && evt.asStartElement().getName().getLocalPart().equals("Iteration"))
					{
					final Iteration iteration = this.unmarshaller.unmarshal(rx, Iteration.class).getValue();
					sortingCollection.add(iteration);
					}
				else if(evt.isStartElement() && 
						evt.asStartElement().getName().getLocalPart().equals("BlastOutput_iterations"))
					{
					wx.add(rx.nextEvent());
					in_iteration=true;
					}
				else if(evt.isEndElement() && evt.asEndElement().getName().getLocalPart().equals("BlastOutput_iterations"))
					{
					for(int optind=1;optind < args.size();++optind)
						{
						LOG.info("opening "+args.get(optind));
						rx2=xmlInputFactory.createXMLEventReader(new FileReader(args.get(optind)));
						while(rx2.hasNext())
							{
							XMLEvent evt2=rx2.peek();
							if(evt2.isStartElement() && evt2.asStartElement().getName().getLocalPart().equals("Iteration"))
								{
								final Iteration iteration = this.unmarshaller.unmarshal(rx2, Iteration.class).getValue();
								sortingCollection.add(iteration);
								}
							else
								{
								rx2.nextEvent();
								}
							}
						rx2.close();
						LOG.info("close");
						}
					
					sortingCollection.doneAdding();
					sortingCollection.setDestructiveIteration(true);
					final CloseableIterator<Iteration> coliter =sortingCollection.iterator();
					final EqualRangeIterator<Iteration> eq=new EqualRangeIterator<>(coliter, hitComparator);
					while(coliter.hasNext())
						{
						final List<Iteration> L=eq.next();
						for(int i=1;i<L.size();++i)
							{
							L.get(0).getIterationHits().getHit().addAll(L.get(i).getIterationHits().getHit());
							}
						marshaller.marshal(L.get(0), wx);
						}
					eq.close();
					coliter.close();
					sortingCollection.cleanup();
					sortingCollection=null;
					wx.add(rx.nextEvent());
					in_iteration=false;
					}
				else if(in_iteration)
					{
					rx.nextEvent();//consumme text
					}
				else
					{	
					wx.add(rx.nextEvent());
					}
				}
			
			wx.flush();
			wx.close();
			return 0;
		} catch (Exception e) {
			LOG.error(e);
			if(sortingCollection!=null) 
				{
				sortingCollection.cleanup();
				}
			return -1;
		}
		finally
		{
			
		}
	}

public static void main(String[] args) {
	new MergeBlastXml().instanceMainWithExit(args);
}
}
