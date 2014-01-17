package com.github.lindenb.jvarkit.tools.treepack;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.XMLConstants;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;


import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Hershey;


public abstract class AbstractTreePackCommandLine<WATCH>
	extends AbstractCommandLineProgram
	{
	protected Rectangle viewRect=new Rectangle(1000,1000);
	protected final String SVG="http://www.w3.org/2000/svg";
	protected final String XLINK="http://www.w3.org/1999/xlink";
	protected List<NodeFactory> nodeFactoryChain=new ArrayList<NodeFactory>();
	protected AbstractNode root=null;
	protected Hershey hershey=new Hershey();   

	
	public abstract class NodeFactory
		{
		public  abstract String getName();
		public  abstract BranchNode createBranch(BranchNode parent);
		}
	
	protected abstract class AbstractNode
		implements TreePack
		{
		private AbstractNode parent=null;
		private Rectangle2D bounds=new Rectangle2D.Double();
		
		protected AbstractNode(AbstractNode parent)
			{
			this.parent=parent;
			}
		
		public AbstractNode getParent()
			{
			return parent;
			}
		
		
		public int getDepth()
			{
			int d=0;
			AbstractNode p=this;
			while(p.parent!=null)
				{
				p=p.parent;
				++d;
				}
			return d;
			}
		
		@Override
		public void setBounds(Rectangle2D bounds)
			{
			this.bounds=bounds;
			}
		
		public Rectangle2D getInsetBounds()
			{
			final double ratio=0.95;
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
	
		@Override
		public abstract double getWeight();
		
		public abstract void watch(WATCH object);
		public abstract boolean isLeaf();
		
		
		public abstract void layout(TreePacker packer);
		
		public void svg(XMLStreamWriter w,String label)throws XMLStreamException
		   {
		   if(this.isLeaf()) return;
		   final Rectangle2D bounds=this.getBounds();
		   w.writeStartElement("svg","g",SVG);
		   if(label!=null) w.writeAttribute("title",label);
		 
		  
		   BranchNode me=(BranchNode)this;
		   for(Object key:me.children.keySet())
			   {
			   AbstractNode c=me.children.get(key);
			   if(!c.isLeaf())
				   {
				   c.svg(w,String.valueOf(key));
				   }
			   else
			   	{
				w.writeEmptyElement("svg", "path", SVG);
				w.writeAttribute("d", hershey.svgPath(String.valueOf(getWeight()),c.getInsetBounds()));
			   	}
			   }
		   
			  

		   w.writeEmptyElement("svg", "rect", SVG);
		  
		   w.writeAttribute("x",String.valueOf(bounds.getX()));
		   w.writeAttribute("y",String.valueOf(bounds.getY()));
		   w.writeAttribute("width",String.valueOf(bounds.getWidth()));
		   w.writeAttribute("height",String.valueOf(bounds.getHeight()));

		   
		   w.writeEndElement();//g
		   }
		}
	
	
	
	protected class LeafNode
		extends AbstractNode
		{
		long count;
		
		LeafNode(AbstractNode parent)
			{
			super(parent);
			}
		@Override
		public void watch(WATCH object) {
			++count;
			}
		@Override
		public boolean isLeaf() {
			return true;
			}
		@Override
		public double getWeight() {
			return count;
			}
		@Override
		public void layout(TreePacker packer) {
			//nothing
			}
		}
	
	
	
	protected abstract class BranchNode
		extends AbstractNode
			{
			long count;
			private Map<Object,AbstractNode> children=new HashMap<Object,AbstractNode>();
			protected NodeFactory factory=null;
			
			
			
			protected BranchNode(NodeFactory factory,AbstractNode parent)
				{
				super(parent);
				this.factory=factory;
				}
			
			public List<AbstractNode> getChildren()
				{
				return new ArrayList<AbstractNode>(this.children.values());
				}
			
			@Override
			public double getWeight()
				{
				double N=0;
				for(TreePack c:this.getChildren()) N+=c.getWeight();
				return N;
				}
			protected AbstractNode createChildren(Object key)
				{
				if(getDepth()+1< AbstractTreePackCommandLine.this.nodeFactoryChain.size())
					 {
					 return nodeFactoryChain.get(getDepth()+1).createBranch(this);
					 }
				 else
					 {
					 return  new LeafNode(this);
					 }
				}
			
			protected AbstractNode get(Object key)
				{
				AbstractNode c=this.children.get(key);
				if(c==null)
					{
					c =createChildren(key);
					this.children.put(key, c);
					}
				return c;
				}

			@Override
			public final boolean isLeaf()
				{
				return false;
				}
			@Override
			public void layout(TreePacker packer)
				{
				Rectangle2D frame=getParent()==null ?
						new Rectangle2D.Double(0,0,viewRect.width,viewRect.height):
						getParent().getBounds();
				List<TreePack> L=new ArrayList<TreePack>(children.values());
				packer.layout(L, frame 	);
				for(AbstractNode c:getChildren())
					{
					c.layout(packer);
					}
				}
			}
	
	protected void svg() throws XMLStreamException
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
		w.writeCharacters(getProgramName());
		w.writeEndElement();
		w.writeComment(getProgramCommandLine());
		w.writeComment(getVersion());
		root.svg(w,null);
		w.writeEndElement();//svg
		w.writeEndDocument();
		w.flush();
	  	}
	
	protected AbstractTreePackCommandLine()
		{
		
		}
	protected void layout()
		{
		this.root.setBounds(new Rectangle2D.Double(0, 0, this.viewRect.getWidth(), this.viewRect.getHeight()));
		this.root.layout(new TreePacker());
		}
	protected abstract List<NodeFactory> getAllAvailableFactories();
	
	protected int buildFactoryChain(String path)
		{
		Map<String,NodeFactory> factorymap=new HashMap<String,NodeFactory>();
		for(NodeFactory f:getAllAvailableFactories())
			{
			factorymap.put(f.getName(), f);
			}
		
		if(path==null)
			{
			warning("No path defined!");
			this.root=new LeafNode(null);
			}
		else
			{
			for(String p: path.split("[ \t/,;]+"))
				{
				NodeFactory nf=factorymap.get(p.toLowerCase());
				if(nf==null)
					{
					error("Cannot get type \'"+p+"\' in "+factorymap.keySet());
					return -1;
					}
				if(this.nodeFactoryChain.isEmpty())
					{
					this.root=nf.createBranch(null);
					}
				this.nodeFactoryChain.add(nf);
				}
			if(this.root==null)
				{
				error("Wrong path "+path);
				return -1;
				}
			}
		return 0;
		}
	
	
	}
