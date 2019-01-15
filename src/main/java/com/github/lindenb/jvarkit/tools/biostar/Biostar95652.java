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
package com.github.lindenb.jvarkit.tools.biostar;

import gov.nih.nlm.ncbi.gb.GBFeature;
import gov.nih.nlm.ncbi.gb.GBInterval;
import gov.nih.nlm.ncbi.gb.GBQualifier;
import gov.nih.nlm.ncbi.gb.GBSeq;
import gov.nih.nlm.ncbi.gb.GBSet;

import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.xml.XMLConstants;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.parsers.SAXParserFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.transform.sax.SAXSource;

import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ncbi.NcbiApiKey;
import com.github.lindenb.jvarkit.util.ncbi.NcbiConstants;
import com.github.lindenb.jvarkit.util.svg.SVG;

/**
BEGIN_DOC

## Example

```bash
$ java -jar dist/biostar95652.jar \
   NP_077719.2 \
   XP_513697.3 XP_001114248.1 \
   XP_540266.3 XP_002686160.2 \
   NP_035058.2 NP_077334.1 \
   NP_001238962.1 NP_001108566.1 > result.svg
```

Result:

![Hosted by imgur.com](http://i.imgur.com/SYn6IAal.png)

END_DOC

 */

@Program(name="biostar95652",
	keywords={"genbank","svg","tree","evolution"},
	description="Drawing a schematic genomic context tree.",
	biostars=95652)
public class Biostar95652 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar95652.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private NcbiApiKey ncbiApiKey = new NcbiApiKey();

	private static final String XLINK="http://www.w3.org/1999/xlink";

	@SuppressWarnings("unused")
	private gov.nih.nlm.ncbi.gb.ObjectFactory _fool_javac=null;
	private static String COLORS[]={"blue","green","orange","red","yellow","pink"};
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
		final List<DomainRegion> domains=new ArrayList<DomainRegion>();
		
		double pos2pix(int pos)
			{
			return (pos/(double)Biostar95652.this.maxLength)*Biostar95652.this.seqWidth;
			}
		
		void paint(XMLStreamWriter w,double x,double y) throws XMLStreamException
			{
			w.writeStartElement("g");
			double h=Biostar95652.this.seqHeight*0.3;
			w.writeEmptyElement("rect");
			w.writeAttribute("class", "protein");
			w.writeAttribute("title", "Length="+this.length);
			w.writeAttribute("x", String.valueOf(x));
			w.writeAttribute("y", String.valueOf(y-h/2));
			w.writeAttribute("width", String.valueOf(pos2pix(this.length)));
			w.writeAttribute("height", String.valueOf(h));
			
			for(final DomainRegion r:this.domains)
				{
				h=Biostar95652.this.seqHeight*0.8;
				
				w.writeStartElement("a");
				w.writeAttribute("xlink", XLINK, "href",
						"https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="+r.domain.cdd
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
				Biostar95652.this.leafList.add(this);
				Biostar95652.this.maxDepth=Math.max(Biostar95652.this.maxDepth, getDepth());
				Biostar95652.this.maxLength=Math.max(Biostar95652.this.maxLength, protein.length);
				
				
				Biostar95652.this.maxOrganismLength=Math.max(this.name.length(),Biostar95652.this.maxOrganismLength);
				Biostar95652.this.maxAcnLength=Math.max(this.protein.locus.length(),Biostar95652.this.maxAcnLength);
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
			
			this.x=(Biostar95652.this.treeWidth/(double)Biostar95652.this.maxDepth)*this.getDepth();
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
				r.height= Biostar95652.this.seqHeight*0.5;
				r.x=this.x;
				r.width=(this.name.length()/(double)Biostar95652.this.maxOrganismLength)*(Biostar95652.this.organismWidth-10);
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
				w.writeAttribute("d",Biostar95652.this.hershey.svgPath(this.name, r) );
				if(protein.taxon_id!=null)
					{
					w.writeEndElement();
					}
				//
				
				w.writeStartElement("a");
				w.writeAttribute("xlink", XLINK, "href",
						"https://www.ncbi.nlm.nih.gov/protein/"+protein.locus
						);
				r.x+= Biostar95652.this.organismWidth;
				r.width=Biostar95652.this.acnWidth-10;
				w.writeEmptyElement("path");
				w.writeAttribute("class", "acn");
				w.writeAttribute("title", this.protein.definition);
				w.writeAttribute("d",Biostar95652.this.hershey.svgPath(this.protein.locus, r) );
				w.writeEndElement();
				
				
				this.protein.paint(w,
						Biostar95652.this.treeWidth+
						Biostar95652.this.organismWidth+
						Biostar95652.this.acnWidth
						,this.y);
				}
			for(Node c:this.name2node.values()) c.paint(w);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
	
		try
			{
			if(args.isEmpty())
				{
				LOG.error("protein ID missing");
				return -1;
				}
			if(!this.ncbiApiKey.isApiKeyDefined()) {
				LOG.error("NCBI API key is not defined");
				return -1;
				}
			
			JAXBContext context = JAXBContext.newInstance("gov.nih.nlm.ncbi.gb");
			final Unmarshaller unmarshaller=context.createUnmarshaller();
			//  https://stackoverflow.com/questions/31293624
			try {
				unmarshaller.setProperty(XMLConstants.ACCESS_EXTERNAL_DTD, "all");
				}
			catch(final Throwable err)
				{
				LOG.warn("Ignoring :" +err.getMessage());
				}
			for(final String arg:args)
				{
				String uri=
						NcbiConstants.efetch()+
						"?db=protein&rettype=gb&retmode=xml&id="+
						StringUtils.escapeHttp(arg)+
						this.ncbiApiKey.getAmpParamValue();
				
				LOG.info("Reading from "+uri);
				
				// https://stackoverflow.com/questions/24460892/
				SAXParserFactory spf = SAXParserFactory.newInstance();
			    spf.setValidating(false); // Not required for JAXB/XInclude
			    final XMLReader xr = spf.newSAXParser().getXMLReader();
			    final SAXSource source = new SAXSource(xr, new InputSource(uri));
				GBSet gbset=(GBSet)unmarshaller.unmarshal(source);
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
			PrintStream ps = super.openFileOrStdoutAsPrintStream(outputFile);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			xof.setProperty(XMLOutputFactory.IS_REPAIRING_NAMESPACES, Boolean.TRUE);
			XMLStreamWriter w=xof.createXMLStreamWriter(ps, "UTF-8");
			w.writeStartDocument("UTF-8", "1.0");
			w.writeStartElement("svg");
			w.writeDefaultNamespace(SVG.NS);
			w.writeNamespace("xlink", XLINK);
			w.writeAttribute("version", "1.1");
			w.writeAttribute("width",String.valueOf(2+this.treeWidth+this.organismWidth+this.acnWidth+this.seqWidth));
			w.writeAttribute("height",String.valueOf(2+this.leafList.size()*seqHeight));
			w.writeComment(this.getProgramCommandLine());
			w.writeComment("Version:"+getVersion());
			w.writeComment("Author: Pierre lindenbaum Phd");
			
			
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
			ps.close();
			ps=null;
			LOG.info("Done");
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
		new Biostar95652().instanceMainWithExit(args);
	}
}
