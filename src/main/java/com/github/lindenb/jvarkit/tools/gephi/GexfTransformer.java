/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.gephi;

import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;

/** 
BEGIN_DOC



END_DOC
*/
@Program(name="gexftr",
	description="Gexf file manipulation",
	keywords={"gexf","graph","network"},
	generate_doc=false
	)
public class GexfTransformer extends Launcher {
	private static final Logger LOG = Logger.build(GexfTransformer.class).make();
	private static int ID_GENERATOR=0;
	
	private static class Cluster{
		final int id=ID_GENERATOR++;
		final Set<String> nodes = new HashSet<>();
		Cluster() {
		}
		Cluster(final String nodeid) {
			this.nodes.add(nodeid);
			}
		Cluster(Cluster c1,Cluster c2) {
			this.nodes.addAll(c1.nodes);
			this.nodes.addAll(c2.nodes);
			}
		@Override
		public int hashCode() {
			return Integer.hashCode(id);
			}
		@Override
		public boolean equals(Object obj) {
			return obj==this || Cluster.class.cast(obj).id==id;
			}
		}
	private void slipElement(XMLEventReader r) throws XMLStreamException {
		while(r.hasNext()) {
			final XMLEvent evt = r.nextEvent();
			if(evt.isStartElement())
				{
				slipElement(r);
				}
			else if(evt.isEndElement())
				{
				return ;
				}
			}
		throw new IllegalStateException("boum");
	}
		
	@Override
	public int doWork(final List<String> args) {
		try {
			final String filename = oneFileOrNull(args);
			final File gexFile;
			final Map<String,Cluster> node2cluster = new HashMap<>();
			final Map<Integer,Cluster> clusterid2cluster = new HashMap<>();
			
			if(filename==null) {
				LOG.info("reading from stdin");
				gexFile = File.createTempFile("tmp", ".gexf");
				gexFile.deleteOnExit();
				IOUtils.copyTo(stdin(), gexFile);
				}
			else
				{
				gexFile = new File(filename);
				IOUtil.assertFileIsWritable(gexFile);
				}
			LOG.info("scanning GEXF input "+gexFile);
			// 1st pass: scan the XML
			final XMLInputFactory xif=XMLInputFactory.newFactory();
			xif.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, true);
			FileReader fr = new FileReader(gexFile);
			XMLEventReader r = xif.createXMLEventReader(fr);
			final QName id_qname = new QName("id");
			final QName source_qname = new QName("source");
			final QName target_qname = new QName("target");
			while(r.hasNext()) {
				final XMLEvent evt = r.nextEvent();
				if(evt.isStartElement()) {
					final StartElement elt =evt.asStartElement();
					final String name = elt.getName().getLocalPart();
					final Attribute node_id=elt.getAttributeByName(id_qname);
					final Attribute edge_source = elt.getAttributeByName(source_qname);
					final Attribute edge_target =elt.getAttributeByName(target_qname);
					if(name.equals("node") && 
							node_id!=null &&
							!node2cluster.containsKey(node_id.getValue())
						) {
						final Cluster cluster=new Cluster(node_id.getValue());
						node2cluster.put(node_id.getValue(), cluster);
						clusterid2cluster.put(cluster.id, cluster);
						}
					else if(name.equals("edge") && 
							edge_source!=null &&
							edge_target!=null &&
							node2cluster.containsKey(edge_source.getValue()) &&
							node2cluster.containsKey(edge_target.getValue())
							)
						{
						final Cluster cluster1 = node2cluster.get(edge_source.getValue());
						final Cluster cluster2 = node2cluster.get(edge_target.getValue());
						if(cluster1.id!=cluster2.id) 
							{
							final Cluster cluster3= new Cluster(cluster1,cluster2);
							node2cluster.put(edge_source.getValue(), cluster3);
							node2cluster.put(edge_target.getValue(), cluster3);
							clusterid2cluster.remove(cluster1.id);
							clusterid2cluster.remove(cluster2.id);
							clusterid2cluster.put(cluster3.id, cluster3);
							}
						}
				}
			}
			r.close();
			fr.close();
			
			Cluster biggest = 
					clusterid2cluster.values().
					stream().
					sorted((C1,C2)->Integer.compare(C2.nodes.size(),C1.nodes.size())).
					findFirst().orElse(new Cluster());
			
			LOG.info("number of nodes "+node2cluster.size());
			LOG.info("largest cluster size:"+biggest.nodes.size());
			
			LOG.info("saving GEXF ");

			final XMLOutputFactory xof = XMLOutputFactory.newFactory();
			final XMLEventWriter w = xof.createXMLEventWriter(stdout(),"UTF-8");
			
			fr = new FileReader(gexFile);
			r = xif.createXMLEventReader(fr);
			while(r.hasNext()) {
				final XMLEvent evt = r.nextEvent();
				
				if(evt.isStartElement()) {
					final StartElement elt =evt.asStartElement();
					final String name = elt.getName().getLocalPart();
					final Attribute node_id=elt.getAttributeByName(id_qname);
					final Attribute edge_source = elt.getAttributeByName(source_qname);
					final Attribute edge_target =elt.getAttributeByName(target_qname);
					if(name.equals("node") && 
							(node_id==null || !biggest.nodes.contains(node_id.getValue()))
						) {
						slipElement(r);
						continue;
						}
					else if(name.equals("edge") &&
							(
							edge_source==null ||
							edge_target==null ||
							!biggest.nodes.contains(edge_source.getValue()) ||
							!biggest.nodes.contains(edge_target.getValue())
							)
							)
						{
						slipElement(r);
						continue;
						}
					}
				w.add(evt);
				}
			w.flush();
			w.close();
			r.close();
			fr.close();
			
			
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
		}
	}
	
	
	public static void main(final String[] args) {
		new GexfTransformer().instanceMainWithExit(args);

	}

}
