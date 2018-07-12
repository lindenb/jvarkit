package com.github.lindenb.jvarkit.tools.treepack;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.io.InputStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.XMLConstants;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;


public class TreePackApp
	extends Launcher
	{
	private static final Logger LOG = Logger.build(TreePackApp.class).make();

	@Parameter(names="-x",description="{width}x{height} dimension of the output rectangle")
	protected Rectangle viewRect=new Rectangle(1000,1000);
	protected final String SVG="http://www.w3.org/2000/svg";
	protected final String XLINK="http://www.w3.org/1999/xlink";
	protected final String JVARKIT_NS="https://github.com/lindenb/jvarkit";
	protected Node root=null;
	protected Hershey hershey=new Hershey();   

	
	private class Node
		implements TreePack
		{
		private Node parent=null;
		private Rectangle2D bounds=new Rectangle2D.Double(-10,-10,10,10);
		private String label="";
		private Map<String,Node> children=null;
		private Double weight=null;
		
		protected Node(Node parent)
			{
			this.parent=parent;
			}
		
		
		public final boolean isLeaf()
			{
			return children==null || children.isEmpty();
			}
		
		protected String convertWeightToString()
			{
			return TreePackApp.this.convertWeightToString(getWeight());
			}


		
		public String getPath()
			{
			StringBuilder b=new StringBuilder();
			Node curr=this;
			while(curr!=null)
				{
				if(b.length()!=0) b.insert(0,"/");
				b.insert(0, curr.getLabel());
				curr=curr.parent;
				}
			return b.toString();
			}
		
		public String getLabel()
			{
			return this.label;
			}
		
		
		@Override
		public void setBounds(Rectangle2D bounds)
			{
			this.bounds=bounds;
			}
		
		public Rectangle2D getInsetBounds()
			{
			final double ratio=0.9;
			final Rectangle2D r0=getBounds();
			Rectangle2D.Double r1=new Rectangle2D.Double();
			r1.width=r0.getWidth()*ratio;
			r1.x=r0.getX()+(r0.getWidth()-r1.width)/2;
			r1.height=r0.getHeight()*ratio;
			r1.y=r0.getY()+(r0.getHeight()-r1.height)/2;
			return r1;
			}
	
		@Override
		public Rectangle2D getBounds()
			{
			return this.bounds;
			}
	
		public double getWeight()
			{
			if(isLeaf())
				{
				return this.weight==null?1.0:this.weight;
				}
			else
				{
				double N=0;
				for(Node c:this.children.values())
					{
					if(c.getWeight()<=0) throw new IllegalStateException(getPath());
					N+=c.getWeight();
					}
				return N;
				}
			}

		protected Rectangle2D getTitleFrame()
			{
			Rectangle2D r=getBounds(); 
			Rectangle2D.Double frame= new Rectangle2D.Double();
			frame.height=r.getHeight()*0.1;
			frame.y=r.getY();
			frame.width=r.getWidth();
			frame.x=r.getX();
			return frame;
			}
		
		protected Rectangle2D getChildrenFrame()
			{
			Rectangle2D r=getBounds(); 
			Rectangle2D.Double frame= new Rectangle2D.Double();
			frame.height=r.getHeight()*0.9;//yes again after
			frame.y=r.getMaxY()-frame.height;
			frame.height=r.getHeight()*0.85;//yes again 
			frame.width=r.getWidth()*0.95;
			frame.x=r.getX()+(r.getWidth()-frame.width)/2.0;
			return frame;
			}
		
		public void layout(TreePacker packer)
			{
			if(getWeight()<=0 || isLeaf()) return ;
			List<TreePack> L=new ArrayList<TreePack>(children.values());
			packer.layout(L, getChildrenFrame() );
			for(Node c:children.values())
				{
				c.layout(packer);
				}
			}
		
		public int getDepth()
			{
			int d=0;
			Node p=this;
			while(p.parent!=null)
				{
				p=p.parent;
				++d;
				}
			return d;
			}
		
		
		public void svg(XMLStreamWriter w) throws XMLStreamException
		   {
		   if(getWeight()<=0)
			   {
			   LOG.info("weight < 0 for "+getPath());
			   return ;
			   }
		   final Rectangle2D bounds=this.getBounds();
		   Rectangle2D insets=this.getInsetBounds();
		   Rectangle2D frameUsed=bounds;
		   
		   
		   if(bounds.getWidth()<=1 || bounds.getHeight()<=1)
			   {
			   return;
			   }
		   
		   w.writeStartElement("svg","g",SVG);
		   w.writeAttribute("title",getPath()+"="+convertWeightToString());
		   w.writeAttribute("jvarkit",JVARKIT_NS,"path",getPath());
		   w.writeAttribute("jvarkit",JVARKIT_NS,"weight",convertWeightToString());
		  	if(!isLeaf())
		   		{
		   		w.writeAttribute("jvarkit",JVARKIT_NS,"count",""+ this.children.size());
		   		}
		   
		   w.writeEmptyElement("svg", "rect", SVG);
		   w.writeAttribute("class", "r"+(getDepth()%2));
		   w.writeAttribute("x",String.valueOf(frameUsed.getX()));
		   w.writeAttribute("y",String.valueOf(frameUsed.getY()));
		   w.writeAttribute("width",String.valueOf(frameUsed.getWidth()));
		   w.writeAttribute("height",String.valueOf(frameUsed.getHeight()));
			   

		   if(!isLeaf())
			   {
			   
			   w.writeEmptyElement("svg", "path", SVG);
			   w.writeAttribute("d",
						hershey.svgPath(getLabel()+"="+convertWeightToString(),
								this.getTitleFrame())
						);
			   }
		   

		   
		   
		   
		   if(isLeaf())
			   {
			   	Rectangle2D f_up=new Rectangle2D.Double(
			   			insets.getX(),insets.getY(),
			   			insets.getWidth(),insets.getHeight()/2.0
			   			);
			   	w.writeEmptyElement("svg", "path", SVG);
				w.writeAttribute("d",
						hershey.svgPath(getLabel(),f_up)
						);	
				w.writeAttribute("class","lbla"+(getDepth()%2));

				Rectangle2D f_down=new Rectangle2D.Double(
						insets.getX(),insets.getCenterY(),
						insets.getWidth(),insets.getHeight()/2.0
			   			);
				w.writeEmptyElement("svg", "path", SVG);
				w.writeAttribute("d",
						hershey.svgPath(convertWeightToString(),f_down)
						);
				w.writeAttribute("class","lblb"+(getDepth()%2));
			   }
		   else
			   {
			   for(Node child:children.values())
				   {
				   child.svg(w);
				   }
			   }
		   
		   w.writeEndElement();//g
		   }
		
		}
	

	
	protected String getCascadingStylesheet()
		{
		return "svg {fill:none;stroke:black;stroke-width:0.5px;}\n"+
				".r0 {fill:none;stroke:black;stroke-width:0.5px;}\n"+
				".r1 {fill:none;stroke:black;stroke-width:0.5px;}\n"+
				".lbla0 {stroke:black;stroke-width:0.3px;}\n"+
				".lblb0 {stroke:red;stroke-width:0.3px;}\n"+
				".lbla1 {stroke:gray;stroke-width:0.3px;}\n"+
				".lblb1 {stroke:red;stroke-width:0.3px;}\n"+
				""
				;
		}
	
	protected void svg() throws XMLStreamException
	  	{
		LOG.info("writing svg");
		XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
		XMLStreamWriter w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
		w.writeStartDocument("UTF-8","1.0");
		w.writeStartElement("svg","svg",SVG);
		
		
		
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "svg", SVG);
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "xlink", XLINK);
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "jvarkit", JVARKIT_NS);
		w.writeAttribute("width",String.valueOf(this.viewRect.getWidth()));
		w.writeAttribute("height",String.valueOf(this.viewRect.getHeight()));
		
		
		w.writeStartElement("svg","defs",SVG);
		w.writeStartElement("svg","style",SVG);
		w.writeAttribute("type","text/css");
		w.writeCharacters(getCascadingStylesheet());
		w.writeEndElement();//svg:style
		w.writeEndElement();//svg:def

		
		w.writeStartElement("svg","title",SVG);
		w.writeCharacters(getProgramName());
		w.writeEndElement();
		w.writeComment("Cmd-Line:"+getProgramCommandLine());
		w.writeComment("Version "+getVersion());
		root.svg(w);
		w.writeEndElement();//svg
		w.writeEndDocument();
		w.flush();
	  	}
	
	protected TreePackApp()
		{
		
		}
	
	protected String convertWeightToString(double weight)
		{
		NumberFormat format=NumberFormat.getInstance();
		return format.format((long)weight);
		}
	
	protected String intervalToString(int value,int step)
		{
		value=((int)(value/(double)step));
		return String.valueOf((int)(value)*step+"-"+((value+1)*step));
		}

	private void layout()
		{
		LOG.info("layout");
		this.root.setBounds(new Rectangle2D.Double(0, 0, this.viewRect.getWidth(), this.viewRect.getHeight()));
		TreePacker packer=new TreePacker();
		this.root.layout(packer);
		}
	
	private void skip(XMLEventReader r) throws XMLStreamException
		{
		while(r.hasNext())
			{
			XMLEvent evt=r.nextEvent();
			if(evt.isStartElement())
				{
				skip(r);
				}
			else if(evt.isEndElement())
				{
				return;
				}
			}
		throw new XMLStreamException();
		}

	private Node parseNode(XMLEventReader r,StartElement E,Node parent) throws XMLStreamException
		{
		Node node=new Node(parent);
		Attribute att=E.getAttributeByName(new QName("label"));
		node.label=(att!=null?att.getValue():"");
		while(r.hasNext())
			{
			XMLEvent evt=r.nextEvent();
			if(evt.isStartElement())
				{
				QName qName=evt.asStartElement().getName();
				if(qName.getLocalPart().equals("node"))
					{
					Node child=parseNode(r,evt.asStartElement(),node);
					if(node.children==null) node.children=new  HashMap<String,Node>();
					node.children.put(node.label,child);
					}
				else
					{
					skip(r);
					}
				}
			else if(evt.isEndElement())
				{
				if(node.children==null || node.children.isEmpty())
					{
					att=E.getAttributeByName(new QName("weight"));
					if(att==null)
						{
						node.weight=1.0;
						}
					else
						{
						node.weight=Double.parseDouble(att.getValue());
						if(node.weight<=0) throw new XMLStreamException("bad @weight:"+node.weight, E.getLocation());
						}
					}
				return node;
				}
			}
		throw new IllegalStateException();
		}
	
	
	private Node scan(InputStream in) throws XMLStreamException
		{
		Node n=null;
		XMLInputFactory xif=XMLInputFactory.newFactory();
		XMLEventReader r=xif.createXMLEventReader(in);
		while(r.hasNext())
			{
			XMLEvent evt=r.nextEvent();
			if(evt.isStartElement())
				{
				QName qName=evt.asStartElement().getName();
				if(qName.getLocalPart().equals("node"))
					{
					n=parseNode(r,evt.asStartElement(),null);
					}
				}
			}
		r.close();
		return n;
		}
	
	
	
	@Override
	public int doWork(List<String> args) {
		if(viewRect.x<=0 ||this.viewRect.y<=0 )
			{
			LOG.error("bad viewRect "+viewRect);
			return -1;
			}
		
		try
			{
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				this.root=scan(stdin());
				}
			else if(args.size()==1)
				{
				String filename=args.get(0);
				LOG.info("Reading from "+filename);
				InputStream in=IOUtils.openURIForReading(filename);
				this.root=scan(in);
				CloserUtil.close(in);
				}
			else
				{
				LOG.error("Illegal number of arguments.");
				return -1;
				}
			layout();
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
	}
