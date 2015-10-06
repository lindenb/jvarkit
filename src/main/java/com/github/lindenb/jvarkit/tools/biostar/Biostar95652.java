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
package com.github.lindenb.jvarkit.tools.biostar;

import gov.nih.nlm.ncbi.gb.GBFeature;
import gov.nih.nlm.ncbi.gb.GBInterval;
import gov.nih.nlm.ncbi.gb.GBQualifier;
import gov.nih.nlm.ncbi.gb.GBSeq;
import gov.nih.nlm.ncbi.gb.GBSet;
import htsjdk.samtools.util.CloserUtil;

import java.awt.geom.Rectangle2D;
import java.io.FileWriter;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.xml.sax.InputSource;

import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.svg.SVG;



public class Biostar95652 extends AbstractBiostar95652
	{
	private static final String XLINK=com.github.lindenb.jvarkit.util.ns.XLINK.NS;
	
	@SuppressWarnings("unused")
	private gov.nih.nlm.ncbi.gb.ObjectFactory _fool_javac=null;
	private static String COLORS[]={"blue","green","orange","red","yellow","pink"};
	
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar95652.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBiostar95652.AbstractBiostar95652Command
		{    
	private Map<String,Domain> cdd2domain=new HashMap<String,Domain>();
	private Node root=new Node();
	private List<Node> leafList=new ArrayList<Node>();
	private int maxDepth=0;
	private int maxLength=0;
	private int treeWidth=500;
	private int organismWidth=300;
	private int acnWidth=200;
	private int seqWidth=1000;
	private int seqHeight=60;
	private Hershey hershey=new Hershey();
	private int maxOrganismLength=0;
	private int maxAcnLength=0;
	
	
	private class Domain
		{
		String cdd;
		String region_name;
		String color;
		}
	private class DomainRegion
		{
		Domain domain;
		int start;
		int end;
		@SuppressWarnings("unused")
		char strand;
		int length()
			{
			return end-start;
			}
		}
	private class Protein
		{
		String taxon_id=null;
		String locus;
		String definition;
		int length;
		List<DomainRegion> domains=new ArrayList<DomainRegion>();
		
		double pos2pix(int pos)
			{
			return (pos/(double)MyCommand.this.maxLength)*MyCommand.this.seqWidth;
			}
		
		void paint(XMLStreamWriter w,double x,double y) throws XMLStreamException
			{
			w.writeStartElement("g");
			double h=MyCommand.this.seqHeight*0.3;
			w.writeEmptyElement("rect");
			w.writeAttribute("class", "protein");
			w.writeAttribute("title", "Length="+this.length);
			w.writeAttribute("x", String.valueOf(x));
			w.writeAttribute("y", String.valueOf(y-h/2));
			w.writeAttribute("width", String.valueOf(pos2pix(this.length)));
			w.writeAttribute("height", String.valueOf(h));
			
			for(DomainRegion r:this.domains)
				{
				h=MyCommand.this.seqHeight*0.8;
				
				w.writeStartElement("a");
				w.writeAttribute("xlink", XLINK, "href",
						"http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="+r.domain.cdd
						);
				
				w.writeEmptyElement("rect");
				w.writeAttribute("class", "cdd"+r.domain.cdd);
				w.writeAttribute("title", r.domain.region_name+" "+r.start+"-"+r.end+" cdd"+r.domain.cdd);
				w.writeAttribute("x", String.valueOf(x+pos2pix(r.start)));
				w.writeAttribute("y", String.valueOf(y-h/2));
				w.writeAttribute("width", String.valueOf(1+pos2pix(r.end)-pos2pix(r.start)));
				w.writeAttribute("height", String.valueOf(h));
				
				w.writeEndElement();
				}
			w.writeEndElement();
			}
		
		}
	
	private class Node
		{
		double x;
		double y;
		Node parent=null;
		String name="";
		Map<String,Node> name2node=new HashMap<String,Node>();
		Protein protein=null;
		
		void insert(LinkedList<String> lineage,Protein protein)
			{
			if(lineage.isEmpty())
				{
				this.protein=protein;
				return;
				}
			String organism=lineage.removeFirst();
			
			Node b=this.name2node.get(organism);
			if(b==null)
				{
				b=new Node();
				b.parent=this;
				b.name=organism;
				this.name2node.put(organism,b);
				}
			b.insert(lineage, protein);
			}
		public int getDepth()
			{
			int c=0;
			Node n=this;
			while(n!=null)
				{
				n=n.parent;
				c++;
				}
			return c;
			}
		void simplify()
			{
			for(Node c:new ArrayList<Node>(name2node.values()))
				{
				c.simplify();
				}
			if(protein!=null) return;
			if(this.name2node.size()!=1) return;
			if(parent==null) return;
			Node next=this.name2node.values().iterator().next();
			if(next.protein!=null) return;
			parent.name2node.remove(this.name);
			next.name=this.name+"; "+next.name;
			parent.name2node.put(next.name, next);
			next.parent=parent;
			}
		void compile()
			{
			if(protein!=null)
				{
				MyCommand.this.leafList.add(this);
				MyCommand.this.maxDepth=Math.max(MyCommand.this.maxDepth, getDepth());
				MyCommand.this.maxLength=Math.max(MyCommand.this.maxLength, protein.length);
				
				
				MyCommand.this.maxOrganismLength=Math.max(this.name.length(),MyCommand.this.maxOrganismLength);
				MyCommand.this.maxAcnLength=Math.max(this.protein.locus.length(),MyCommand.this.maxAcnLength);
				}
			for(Node c:this.name2node.values()) c.compile();
			}
		int getWeitht()
			{
			if(protein!=null) return 1;
			int n=0;
			for(Node c:this.name2node.values())
				{
				n+=c.getWeitht();
				}
			return n;
			}
		void compileXY(double y0,double y1)
			{
			double weigh0=this.getWeitht();
			double y=y0;
			for(Node c:this.name2node.values())
				{
				double h=((y1-y0)/weigh0)*c.getWeitht();
				c.y=y+h/2;
				c.compileXY(y, y+h);
				y+=h;
				}
			
			this.x=(MyCommand.this.treeWidth/(double)MyCommand.this.maxDepth)*this.getDepth();
			}
		void paint(XMLStreamWriter w) throws XMLStreamException
			{
			if(parent!=null)
				{
				w.writeEmptyElement("polyline");
				w.writeAttribute("class", "tree");
				w.writeAttribute("title", this.name);
				w.writeAttribute("points", 
						this.x+","+this.y+" "+
						((this.x+this.parent.x)/2)+","+this.y+" "+
						((this.x+this.parent.x)/2)+","+((this.y+this.parent.y)/2)+" "+
						(this.parent.x)+","+((this.y+this.parent.y)/2)+" "+
						this.parent.x+","+this.parent.y
						);
				}
			if(protein!=null)
				{
				Rectangle2D.Double r=new Rectangle2D.Double();
				r.height= MyCommand.this.seqHeight*0.5;
				r.x=this.x;
				r.width=(this.name.length()/(double)MyCommand.this.maxOrganismLength)*(MyCommand.this.organismWidth-10);
				r.y=this.y-r.height/2;
				if(protein.taxon_id!=null)
					{
					w.writeStartElement("a");
					w.writeAttribute("xlink", XLINK, "href",
							"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="+protein.taxon_id
							);
					}
				w.writeEmptyElement("path");
				w.writeAttribute("class", "organism");
				w.writeAttribute("title", this.name);
				w.writeAttribute("d",MyCommand.this.hershey.svgPath(this.name, r) );
				if(protein.taxon_id!=null)
					{
					w.writeEndElement();
					}
				//
				
				w.writeStartElement("a");
				w.writeAttribute("xlink", XLINK, "href",
						"http://www.ncbi.nlm.nih.gov/protein/"+protein.locus
						);
				r.x+= MyCommand.this.organismWidth;
				r.width=MyCommand.this.acnWidth-10;
				w.writeEmptyElement("path");
				w.writeAttribute("class", "acn");
				w.writeAttribute("title", this.protein.definition);
				w.writeAttribute("d",MyCommand.this.hershey.svgPath(this.protein.locus, r) );
				w.writeEndElement();
				
				
				this.protein.paint(w,
						MyCommand.this.treeWidth+
						MyCommand.this.organismWidth+
						MyCommand.this.acnWidth
						,this.y);
				}
			for(Node c:this.name2node.values()) c.paint(w);
			}
		}
	
	
	@Override
		public Collection<Throwable> call() throws Exception {
			final List<String> args=getInputFiles();
			XMLStreamWriter w=null;
			FileWriter fw=null;
			try
				{
				if(args.isEmpty())
					{
					return wrapException("protein ID missing");
					}
				JAXBContext context = JAXBContext.newInstance("gov.nih.nlm.ncbi.gb");
				Unmarshaller unmarshaller=context.createUnmarshaller();
				for(String arg:args)
					{
					String uri="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=gb&retmode=xml&id="+
							URLEncoder.encode(arg,"UTF-8");
					
					LOG.info("Reading from "+uri);
					GBSet gbset=(GBSet)unmarshaller.unmarshal(new InputSource(uri));
					if(gbset.getGBSeq().isEmpty())
						{
						LOG.info("Nothing in "+uri);
						continue;
						}
					GBSeq gbseq=gbset.getGBSeq().get(0);
					Protein protein=new Protein();
					protein.length=Integer.parseInt(gbseq.getGBSeqLength());
					protein.locus=gbseq.getGBSeqLocus();
					protein.definition=gbseq.getGBSeqDefinition();
					for(GBFeature feat:gbseq.getGBSeqFeatureTable().getGBFeature())
						{
						if(feat.getGBFeatureIntervals().getGBInterval().isEmpty()) continue;
						
						String cdd=null;
						String region_name=null;
						for(GBQualifier qual:feat.getGBFeatureQuals().getGBQualifier())
							{
							if(qual.getGBQualifierName().equals("db_xref") &&
									qual.getGBQualifierValue().startsWith("CDD:"))
								{
								cdd=qual.getGBQualifierValue().substring(4);
								}
							else if(qual.getGBQualifierName().equals("db_xref") &&
									qual.getGBQualifierValue().startsWith("taxon:"))
								{
								protein.taxon_id=qual.getGBQualifierValue().substring(6);
								}
							else if(qual.getGBQualifierName().equals("region_name"))
								{
								region_name=qual.getGBQualifierValue();
								}
							}
						if(cdd==null || region_name==null)
							{
							continue;
							}
						Domain domain= cdd2domain.get(cdd);
						if(domain==null)
							{
							domain=new Domain();
							domain.cdd=cdd;
							domain.region_name=region_name;
							domain.color=COLORS[cdd2domain.size()%COLORS.length];
							cdd2domain.put(domain.cdd, domain);
							}
						for(GBInterval interval:feat.getGBFeatureIntervals().getGBInterval())
							{
							if(interval.getGBIntervalFrom()==null || interval.getGBIntervalTo()==null) continue;
							DomainRegion region=new DomainRegion();
							region.domain=domain;
							int start=Integer.parseInt(interval.getGBIntervalFrom());
							int end=Integer.parseInt(interval.getGBIntervalTo());
							if(start<end)
								{
								region.start=start;
								region.end=end;
								region.strand='+';
								}
							else
								{
								region.start=end;
								region.end=start;
								region.strand='-';
								}
							protein.domains.add(region);
							}
						LinkedList<String> lineage=new LinkedList<String>(
								Arrays.asList(gbseq.getGBSeqTaxonomy().split("[;][ ]*")));
						lineage.add(gbseq.getGBSeqOrganism());
						Collections.sort(protein.domains,new Comparator<DomainRegion>()
							{
							@Override
							public int compare(DomainRegion o1, DomainRegion o2)
								{
								return o2.length()-o1.length();
								}
							});
						this.root.insert(lineage,protein);
						}
					
					}
				root.simplify();
				root.compile();
				root.x=0;
				root.y= (this.leafList.size()*seqHeight)/2.0;
				root.compileXY(0,this.leafList.size()*seqHeight);
				
				XMLOutputFactory xof=XMLOutputFactory.newFactory();
				xof.setProperty(XMLOutputFactory.IS_REPAIRING_NAMESPACES, Boolean.TRUE);
				
				if(getOutputFile()==null)
					{
					w=xof.createXMLStreamWriter(stdout(), "UTF-8");
					}
				else
					{
					fw=new FileWriter(getOutputFile());
					w=xof.createXMLStreamWriter(fw);
					}
				w.writeStartDocument("UTF-8", "1.0");
				w.writeStartElement("svg");
				w.writeDefaultNamespace(SVG.NS);
				w.writeNamespace("xlink", XLINK);
				w.writeAttribute("version", "1.1");
				w.writeAttribute("width",String.valueOf(2+this.treeWidth+this.organismWidth+this.acnWidth+this.seqWidth));
				w.writeAttribute("height",String.valueOf(2+this.leafList.size()*seqHeight));
				w.writeComment(this.getProgramCommandLine());
				w.writeComment("Version:"+getVersion());
				w.writeComment("Author:"+getFactory().getAuthorName());
				w.writeComment("URL:"+getFactory().getOnlineDocUrl());
				
				
				w.writeStartElement("defs");
				
				w.writeStartElement("linearGradient");
					w.writeAttribute("id","grad01");
					w.writeAttribute("x1","50%");
					w.writeAttribute("x2","50%");
					w.writeAttribute("y1","0%");
					w.writeAttribute("y2","100%");
					w.writeEmptyElement("stop");
						w.writeAttribute("offset","0%");
						w.writeAttribute("style","stop-color:black;stop-opacity:1;");
					w.writeEmptyElement("stop");
						w.writeAttribute("offset","50%");
						w.writeAttribute("style","stop-color:white;stop-opacity:1;");
					w.writeEmptyElement("stop");
						w.writeAttribute("offset","100%");
						w.writeAttribute("style","stop-color:black;stop-opacity:1;");
				w.writeEndElement();
	
				for(Domain cdd:this.cdd2domain.values())
					{
					w.writeStartElement("linearGradient");
						w.writeAttribute("id","grad"+cdd.cdd);
						w.writeAttribute("x1","50%");
						w.writeAttribute("x2","50%");
						w.writeAttribute("y1","0%");
						w.writeAttribute("y2","100%");
						w.writeEmptyElement("stop");
							w.writeAttribute("offset","0%");
							w.writeAttribute("style","stop-color:"+cdd.color+";stop-opacity:1;");
						w.writeEmptyElement("stop");
							w.writeAttribute("offset","50%");
							w.writeAttribute("style","stop-color:white;stop-opacity:1;");
						w.writeEmptyElement("stop");
							w.writeAttribute("offset","100%");
							w.writeAttribute("style","stop-color:"+cdd.color+";stop-opacity:1;");
					w.writeEndElement();
					}
				
				
				w.writeEndElement();//defs
				
				w.writeStartElement("style");
				w.writeCharacters(
						"svg {fill:none; stroke:black;}\n"+
						".protein { stroke:red;}\n"+
						".tree { stroke:black;fill:none;stroke-width:2}\n"+
						".organism { stroke:black;fill:none;stroke-width:2}\n"+
						".acn { stroke:blue;fill:none;stroke-width:2}\n"+
						".protein {fill:url(#grad01);stroke:black;}\n"
						);
				for(Domain cdd:this.cdd2domain.values())
					{
					w.writeCharacters(".cdd"+cdd.cdd+" {fill:url(#grad"+cdd.cdd+");stroke:orange;stroke-width:3;fill-opacity:0.8;}\n");
					}
				w.writeEndElement();//style
				
				w.writeStartElement("g");
				this.root.paint(w);
				w.writeEndElement();//g
				w.writeEndElement();//svg
				w.writeEndDocument();
				w.flush();
				w.close();
				LOG.info("Done");
				return Collections.emptyList();
				}
			catch(Exception err)
				{
				LOG.error(err);
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(fw);
				CloserUtil.close(w);
				}
			}
		}
	public static void main(String[] args) {
		new Biostar95652().instanceMainWithExit(args);
	}
}
