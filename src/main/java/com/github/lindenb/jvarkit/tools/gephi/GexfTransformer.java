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
package com.github.lindenb.jvarkit.tools.gephi;

import java.io.File;
import java.io.FileReader;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gexf.GexfConstants;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;

/** 
BEGIN_DOC

## Cluster 

a **Cluster** is a set of nodes that are all connected with some edges.

## Example:

```
$ java -jar dist/pubmeddump.jar "NSP1 Rotavirus" | \
  java -jar dist/pubmedauthorgraph.jar -D BDB |\
  java -jar dist/gexftr.jar | xmllint --format - 
  
(...)
      <!--Skip edge  from HOWARD~C_R to pmid:9614866-->
      <!--Skip edge  from BRIDGER~J_C to pmid:9614866-->
      <!--Skip edge  from LÓPEZ~S to pmid:9645203-->
      <!--Skip edge  from GONZÁLEZ~R_A to pmid:9645203-->
      <!--Skip edge  from ARIAS~C_F to pmid:9645203-->
      <!--Skip edge  from TORRES-VEGA~M_A to pmid:9645203-->
    </edges>
  </graph>
  <!--Number of nodes 692-->
</gexf>


```


END_DOC
*/
@Program(name="gexftr",
	description="Gexf file manipulation",
	keywords={"gexf","graph","network"}
	)
public class GexfTransformer extends Launcher {
	private static final Logger LOG = Logger.build(GexfTransformer.class).make();
	private static int ID_GENERATOR=0;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-m","--min-nodes"},description="min number of nodes per cluster inclusive. Negative: no limit")
	private int min_num_node_per_cluster = -1;
	@Parameter(names={"-M","--max-nodes"},description="max number of nodes per cluster inclusive. Negative: no limit")
	private int max_num_node_per_cluster = -1;
	@Parameter(names={"-c","--max-cluster"},description="limit to this number of clusters. Negative: no limit")
	private int limit_count_cluster = 1;

	
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
				gexFile = File.createTempFile("tmp", ".gexf",IOUtils.getDefaultTmpDir());
				gexFile.deleteOnExit();
				IOUtils.copyTo(stdin(), gexFile);
				}
			else
				{
				LOG.info("reading from stdin");
				gexFile = new File(filename);
				IOUtil.assertFileIsWritable(gexFile);
				}
			LOG.info("scanning GEXF input "+gexFile);
			// 1st pass: scan the XML
			final Set<String> other_ns = new HashSet<>();
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
					final String ns = elt.getName().getNamespaceURI();
					if(ns==null || !(ns.equals(GexfConstants.XMLNS_VIZ) || ns.equals(GexfConstants.XMLNS)))
						{
						if(ns==null || other_ns.add(ns)) {
							LOG.warn("found namespace that is not "+GexfConstants.XMLNS+":"+ns);
							}
						slipElement(r);
						continue;
						}
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
							
							if( !cluster1.nodes.contains(edge_source.getValue()) ||
								 cluster1.nodes.contains(edge_target.getValue()) ||
								!cluster2.nodes.contains(edge_target.getValue()) ||
								 cluster2.nodes.contains(edge_source.getValue()))
								{
								throw new IllegalStateException("boum");
								}
							for(final String c1id:cluster1.nodes )
								{
								node2cluster.put(c1id, cluster3);
								}
							for(final String c2id:cluster2.nodes )
								{
								node2cluster.put(c2id, cluster3);
								}
							clusterid2cluster.remove(cluster1.id);
							clusterid2cluster.remove(cluster2.id);
							clusterid2cluster.put(cluster3.id, cluster3);
							}
						}
				}
			}
			r.close();
			fr.close();
			
			final List<Cluster> biggest_list = 
					clusterid2cluster.values().
					stream().
					filter(C->min_num_node_per_cluster<0 || C.nodes.size()>=min_num_node_per_cluster).
					filter(C->max_num_node_per_cluster<0 || C.nodes.size()<=max_num_node_per_cluster).
					sorted((C1,C2)->Integer.compare(C2.nodes.size(),C1.nodes.size())).
					limit(limit_count_cluster<0?Integer.MAX_VALUE:limit_count_cluster).
					collect(Collectors.toList())
					;
			LOG.info("number of cluster "+biggest_list.size());

			final Cluster biggest = new Cluster();
			biggest.nodes.addAll(biggest_list.stream().flatMap(C->C.nodes.stream()).collect(Collectors.toList()));
			biggest_list.clear();
			
			LOG.info("initial number of nodes "+node2cluster.size());
			LOG.info("final number of nodes:"+biggest.nodes.size());
			
			LOG.info("saving GEXF ");

			final XMLOutputFactory xof = XMLOutputFactory.newFactory();
			final OutputStream fos= super.openFileOrStdoutAsStream(this.outputFile);
			final XMLEventWriter w = xof.createXMLEventWriter(fos,"UTF-8");
			final XMLEventFactory eventFactory = XMLEventFactory.newFactory();
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
						final Cluster c1 = (node_id==null?null:node2cluster.get(node_id.getValue()));
						w.add(eventFactory.createComment("Skip node "+node_id.getValue()+
								" / cluster "+
								(c1==null?"???":"ID."+c1.id+"~size="+c1.nodes.size())));
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
						if(edge_source!=null && edge_target!=null) {
							w.add(eventFactory.createComment("Skip edge  from "+edge_source.getValue()+" to "+edge_target.getValue()));
							}
						slipElement(r);
						continue;
						}
					
					}
				
				w.add(evt);
				
				if(evt.isEndElement() && evt.asEndElement().getName().getLocalPart().equals("graph"))
					{
					w.add(eventFactory.createComment("Number of nodes "+ biggest.nodes.size()));
					}
				}
			w.flush();
			w.close();
			fos.flush();
			fos.close();
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
