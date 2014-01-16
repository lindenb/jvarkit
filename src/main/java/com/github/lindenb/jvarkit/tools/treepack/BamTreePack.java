package com.github.lindenb.jvarkit.tools.treepack;

import java.awt.Font;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.XMLConstants;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.sun.org.apache.xerces.internal.impl.xs.XMLSchemaException;

public class BamTreePack extends AbstractCommandLineProgram
	{
	private Rectangle viewRect=new Rectangle(1000,1000);
	private final String SVG="http://www.w3.org/2000/svg";
	private final String XLINK="http://www.w3.org/1999/xlink";
	private List<NodeFactory> nodeFactoryChain=new ArrayList<NodeFactory>();
	private interface NodeFactory
		{
		public Node create(Node parent);
		}
	
	private class SimpleCountNodeFactory implements NodeFactory
		{
		private class MyNode extends BamTreePack.Node
			{
			MyNode(Node parent)
				{
				super(parent);
				}
			@Override
			public void watch(SAMFileHeader header,SAMRecord seq)
				{
				++count;
				}
			}
		public Node create(Node parent)
			{
			return new MyNode(parent);
			}
		}
	
	private class SampleNodeFactory implements NodeFactory
		{
		private class MyNode extends BamTreePack.Node
			{
			private Map<String,Node> chrom2node=new HashMap<String,Node>();
			private Counter<String> chrom2count=new Counter<String>();
			private MyNode(Node parent)
				{
				super(parent);
				}
			
			@Override
			public void watch(SAMFileHeader header,SAMRecord seq)
				{
				++count;
				SAMReadGroupRecord g=seq.getReadGroup();
				String sample=(g==null?"*":g.getSample());
				if(sample==null || sample.isEmpty()) sample="*";
				Node c=chrom2node.get(sample);
				if(c==null)
					{
					c=getFactoryByDepth(this.getDepth()).create(this);
					chrom2node.put(sample, c);
					}
				chrom2count.incr(sample);
				c.watch(header, seq);
				}
			}
		public Node create(Node parent)
			{
			return new MyNode(parent);
			}
		}
	
	private abstract class Node
		implements TreePack
		{
		Node parent=null;
		long count=0L;
		Rectangle2D bounds=new Rectangle2D.Double();		
		
		public abstract void watch(SAMFileHeader header,SAMRecord seq);
		
		protected Node(Node parent)
			{
			this.parent=parent;
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
		
		@Override
		public double getWeight()
			{
			return 0;
			}
		@Override
		public void setBounds(Rectangle2D bounds)
			{
			
			}
		@Override
		public Rectangle2D getBounds()
			{
			return null;
			}
		}
	
	private BamTreePack()
		{
		
		}
	
	 private NodeFactory getFactoryByDepth(int d)
		 {
		 return (d<this.nodeFactoryChain.size()?
				 this.nodeFactoryChain.get(d):
				 new SimpleCountNodeFactory()
		 		);
		 }
	 
	 private void svg(XMLStreamWriter w,Node f)throws Exception
	   {
	   String textPath=null;
	   String selector=null;
	   Font font=new Font(f.getFontFamily(),Font.BOLD,1+(int)Math.max(f.bounds.width,f.bounds.height));
	   w.writeStartElement("svg","g",SVG);
	   selector=f.getStyle("stroke",null);
	   if(selector!=null)  w.writeAttribute("stroke",selector);
	   selector=f.getStyle("fill",null);
	   if(selector!=null)  w.writeAttribute("fill",selector);
	   selector=f.getStyle("stroke-width",String.valueOf(Math.max(0.2,2/(f.getDepth()+1.0))));
	   w.writeAttribute("stroke-width",selector);
	   
	   
	   w.writeEmptyElement("svg", "rect", SVG);
	   if(f.getDescription()!=null) w.writeAttribute("title",f.getDescription());
	  
	   w.writeAttribute("x",String.valueOf(f.bounds.getX()));
	   w.writeAttribute("y",String.valueOf(f.bounds.getY()));
	   w.writeAttribute("width",String.valueOf(f.bounds.getWidth()));
	   w.writeAttribute("height",String.valueOf(f.bounds.getHeight()));
	 
	   if(!f.isLeaf())
		   {
		   Branch branch=(Branch)f;
		   for(Frame c:branch.children)
			   {
			   svg(w,c);
			   }
		   if(f.getLabel()!=null)
			   {
			   Rectangle2D.Double bbox=f.getInnerRect();
			   if(bbox.getMaxY()+6 < f.bounds.getMaxY())
				   {
				   bbox=new Rectangle2D.Double(
						   f.bounds.x,bbox.getMaxY()+2,
						   f.bounds.width,
						   (f.bounds.getMaxY()-bbox.getMaxY())-4
						   );
				   textPath=fitText(f.getLabel(), font, bbox);
				   }
			   else if(bbox.getMaxX()+6 < f.bounds.getMaxX())
				   {
				   bbox=new Rectangle2D.Double(bbox.getMaxX(),f.bounds.y,(f.bounds.getMaxX()-bbox.getMaxX()),f.bounds.height);
				   textPath=fitText(f.getLabel(), font, bbox);
				   }
			   }
		   
		   }
	   else
		   {
		   if(f.getLabel()!=null)
			   {
			   textPath=fitText(f.getLabel(), font, f.getInnerRect());
			   }
		   }
	   
	   if(textPath!=null)
		   {
		   w.writeEmptyElement("svg", "path", SVG);
		   w.writeAttribute("d",textPath);
		   selector=f.getStyle("font-stroke",null);
		   if(selector!=null) w.writeAttribute("font-stroke","none");
		   selector=f.getStyle("font-fill",null);
		   if(selector!=null) w.writeAttribute("font-fill",selector);
		   }
	   
	   w.writeEndElement();//g
	   }
 
  private void svg(Frame f) throws XMLSchemaException
  	{
	XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
	XMLStreamWriter w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
	w.writeStartDocument("UTF-8","1.0");
	w.writeStartElement("svg","svg",SVG);
	w.writeAttribute("style", "fill:none;stroke:black;stroke-width:0.5px;");
	w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "svg", SVG);
	w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "xlink", XLINK);
	w.writeAttribute("width",String.valueOf(this.viewRect.getWidth()));
	w.writeAttribute("height",String.valueOf(this.viewRect.getHeight()));
	w.writeStartElement("svg","title",SVG);
	w.writeCharacters("made with TreeMapMaker (c) Pierre Lindenbaum");
	w.writeEndElement();
	svg(w,f);
	w.writeEndElement();//svg
	w.writeEndDocument();
	w.flush();
  	}

	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -v invert.");
		out.println(" -S (name) add sample.");
		out.println(" -f (file) read file containing sample names. Optional.");
		out.println(" -r remove variant if there is not any called genotype on the line. Optional.");
		out.println(" -E unknown user sample is not a fatal error . Optional.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "ErvS:f:"))!=-1)
			{
			switch(c)
				{
				default: 
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		return -1;
		}
	
	public static void main(String[] args)
		{
		new BamTreePack().instanceMainWithExit(args);

		}

	}
