package com.github.lindenb.jvarkit.tools.treepack;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.io.PrintStream;
import java.text.NumberFormat;
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
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.svg.SVG;


public abstract class AbstractTreePackCommandLine<WATCH>
	extends AbstractCommandLineProgram
	{
	protected Rectangle viewRect=new Rectangle(1000,1000);
	protected final String XLINK="http://www.w3.org/1999/xlink";
	protected final String JVARKIT_NS="https://github.com/lindenb/jvarkit";
	protected List<NodeFactory> nodeFactoryChain=new ArrayList<NodeFactory>();
	protected AbstractNode root=null;
	protected Hershey hershey=new Hershey();   

	
	public abstract class NodeFactory
		{
		public  abstract String getName();
		public abstract String getDescription();
		public  abstract BranchNode createBranch(BranchNode parent,String label);
		@Override
		public String toString()
			{
			return getName();
			}
		}
	
	protected abstract class AbstractNode
		implements TreePack
		{
		private AbstractNode parent=null;
		private Rectangle2D bounds=new Rectangle2D.Double(-10,-10,10,10);
		private String label="";
		
		
		protected AbstractNode(AbstractNode parent,String label)
			{
			this.parent=parent;
			this.label=label;
			}
		
		public boolean isRoot()
			{
			return getParent()==null;
			}
		
		protected String convertWeightToString()
			{
			return AbstractTreePackCommandLine.this.convertWeightToString(getWeight());
			}


		
		public String getPath()
			{
			StringBuilder b=new StringBuilder();
			AbstractNode curr=this;
			while(curr!=null)
				{
				if(b.length()!=0) b.insert(0,"/");
				b.insert(0, curr.getLabel());
				curr=curr.getParent();
				}
			return b.toString();
			}
		
		public String getLabel()
			{
			return this.label;
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
	
		@Override
		public abstract double getWeight();
		
		public abstract void watch(WATCH object);
		public abstract boolean isLeaf();
		
		
		public abstract void layout(TreePacker packer);		
		
		public void svg(XMLStreamWriter w) throws XMLStreamException
		   {
		   if(getWeight()<=0)
			   {
			   info("weight < 0 for "+getPath());
			   return ;
			   }
		   final Rectangle2D bounds=this.getBounds();
		   Rectangle2D insets=this.getInsetBounds();
		   Rectangle2D frameUsed=bounds;
		   
		   
		   if(bounds.getWidth()<=1 || bounds.getHeight()<=1)
			   {
			   return;
			   }
		   
		   w.writeStartElement("svg","g",SVG.NS);
		   w.writeAttribute("title",getPath()+"="+convertWeightToString());
		   w.writeAttribute("jvarkit",JVARKIT_NS,"path",getPath());
		   w.writeAttribute("jvarkit",JVARKIT_NS,"weight",convertWeightToString());
		   	if(!isLeaf())
		   		{
		   		w.writeAttribute("jvarkit",JVARKIT_NS,"category",((BranchNode)this).factory.getName());
		   		w.writeAttribute("jvarkit",JVARKIT_NS,"count",""+((BranchNode)this).children.size());
		   		}
		   
		   w.writeEmptyElement("svg", "rect", SVG.NS);
		   w.writeAttribute("class", "r"+(getDepth()%2));
		   w.writeAttribute("x",String.valueOf(frameUsed.getX()));
		   w.writeAttribute("y",String.valueOf(frameUsed.getY()));
		   w.writeAttribute("width",String.valueOf(frameUsed.getWidth()));
		   w.writeAttribute("height",String.valueOf(frameUsed.getHeight()));
			   

		   if(!isLeaf())
			   {
			   BranchNode me=((BranchNode)this);
			   
			   w.writeEmptyElement("svg", "path", SVG.NS);
			   w.writeAttribute("d",
						hershey.svgPath(getLabel()+"="+convertWeightToString(),
								me.getTitleFrame())
						);
			   }
		   

		   
		   
		   
		   if(isLeaf())
			   {
			   	Rectangle2D f_up=new Rectangle2D.Double(
			   			insets.getX(),insets.getY(),
			   			insets.getWidth(),insets.getHeight()/2.0
			   			);
			   	w.writeEmptyElement("svg", "path", SVG.NS);
				w.writeAttribute("d",
						hershey.svgPath(getLabel(),f_up)
						);	
				w.writeAttribute("class","lbla"+(getDepth()%2));

				Rectangle2D f_down=new Rectangle2D.Double(
						insets.getX(),insets.getCenterY(),
						insets.getWidth(),insets.getHeight()/2.0
			   			);
				w.writeEmptyElement("svg", "path", SVG.NS);
				w.writeAttribute("d",
						hershey.svgPath(convertWeightToString(),f_down)
						);
				w.writeAttribute("class","lblb"+(getDepth()%2));
			   }
		   else
			   {
			   BranchNode me=(BranchNode)this;
			   for(String key:me.children.keySet())
				   {
				   me.children.get(key).svg(w);
				   }
			   }
		   
		   w.writeEndElement();//g
		   }
		
		}
	
	
	
	protected abstract class LeafNode
		extends AbstractNode
		{
		LeafNode(AbstractNode parent,String label)
			{
			super(parent,label);
			}
		
		@Override
		public final boolean isLeaf() {
			return true;
			}
		
		@Override
		public void layout(TreePacker packer) {
			//nothing
			}
		}
	protected class DefaultLeafNode
		extends LeafNode
			{
			long count;
			
			DefaultLeafNode(AbstractNode parent,String label)
				{
				super(parent,label);
				}
			@Override
			public void watch(WATCH object) {
				++count;
				}
			
			@Override
			public double getWeight() {
				return count;
				}
			}
	
	
	protected abstract class BranchNode
		extends AbstractNode
			{
			long count;
			private Map<String,AbstractNode> children=new HashMap<String,AbstractNode>();
			protected NodeFactory factory=null;
			
			
			
			
			protected BranchNode(NodeFactory factory,AbstractNode parent,String label)
				{
				super(parent,label);
				this.factory=factory;
				}
			
			protected String findKeyByValue(AbstractNode child)
				{
				for(String key:children.keySet())
					{
					if(children.get(key)==child) return key;
					}
				return null;
				}
			
			public List<AbstractNode> getChildren()
				{
				return new ArrayList<AbstractNode>(this.children.values());
				}
			
			@Override
			public double getWeight()
				{
				double N=0;
				for(TreePack c:this.getChildren())
					{
					if(c.getWeight()<=0) throw new IllegalStateException(getPath());
					N+=c.getWeight();
					}
				return N;
				}
			
			protected LeafNode createLeafNode(String label)
				{
				return  AbstractTreePackCommandLine.this.createLeafNode(this,label);
				}
			
			protected AbstractNode createChildren(String key)
				{
				if(getDepth()+1< AbstractTreePackCommandLine.this.nodeFactoryChain.size())
					 {
					 return nodeFactoryChain.get(getDepth()+1).createBranch(this,key);
					 }
				 else
					 {
					 return  this.createLeafNode(key);
					 }
				}
			
			protected AbstractNode get(String key)
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
			
			@Override
			public void layout(TreePacker packer)
				{
				if(getWeight()<=0) return ;
				List<TreePack> L=new ArrayList<TreePack>(children.values());
				packer.layout(L, getChildrenFrame() );
				for(AbstractNode c:getChildren())
					{
					c.layout(packer);
					}
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
	  	info("writing svg");
		XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
		XMLStreamWriter w= xmlfactory.createXMLStreamWriter(System.out,"UTF-8");
		w.writeStartDocument("UTF-8","1.0");
		w.writeStartElement("svg","svg",SVG.NS);
		
		
		
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "svg", SVG.NS);
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "xlink", XLINK);
		w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "jvarkit", JVARKIT_NS);
		w.writeAttribute("width",String.valueOf(this.viewRect.getWidth()));
		w.writeAttribute("height",String.valueOf(this.viewRect.getHeight()));
		
		
		w.writeStartElement("svg","defs",SVG.NS);
		w.writeStartElement("svg","style",SVG.NS);
		w.writeAttribute("type","text/css");
		w.writeCharacters(getCascadingStylesheet());
		w.writeEndElement();//svg:style
		w.writeEndElement();//svg:def

		
		w.writeStartElement("svg","title",SVG.NS);
		w.writeCharacters(getProgramName());
		w.writeEndElement();
		w.writeComment("Cmd-Line:"+getProgramCommandLine());
		w.writeComment("Version "+getVersion());
		root.svg(w);
		w.writeEndElement();//svg
		w.writeEndDocument();
		w.flush();
	  	}
	
	protected AbstractTreePackCommandLine()
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

	protected void layout()
		{
		info("layout");
		this.root.setBounds(new Rectangle2D.Double(0, 0, this.viewRect.getWidth(), this.viewRect.getHeight()));
		TreePacker packer=new TreePacker();
		this.root.layout(packer);
		}
	protected abstract List<NodeFactory> getAllAvailableFactories();
	
	protected NodeFactory findFactoryByName(String s)
		{
		for(NodeFactory f:getAllAvailableFactories())
			{
			if(f.getName().equalsIgnoreCase(s)) return f;
			}
		return null;
		}
	protected String getGetOptDefault()
		{
		return super.getGetOptDefault()+"e:x:";
		}
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -e 'path1 path2 path3 path4 ... A list of observed components.");
		out.println("    The current list of available components is:");
		out.println("    #name\tdescription");
		for(NodeFactory f:getAllAvailableFactories())
			{
			out.println("    "+f.getName()+"\t\""+f.getDescription()+"\"");
			}
		out.println(" -x {width}x{height} dimension of the output rectangle. default:"+this.viewRect.width+"x"+this.viewRect.height);
		super.printOptions(out);
		}
	@Override
	protected GetOptStatus handleOtherOptions(int c, GetOpt opt)
		{
		switch(c)
			{
			case 'e':
				if(buildFactoryChain(opt.getOptArg())!=0) return GetOptStatus.EXIT_FAILURE;
				return GetOptStatus.OK;
			case 'x':
					{
					String s=opt.getOptArg();
					int x=s.indexOf('x');
					if(x==-1)
						{
						error("'x' missing in "+s);
						return GetOptStatus.EXIT_FAILURE;
						}
					this.viewRect.x=Integer.parseInt(s.substring(0, x));
					this.viewRect.y=Integer.parseInt(s.substring(x+1));
					if(viewRect.x<=0 ||this.viewRect.y<=0 )
						{
						error("bad viewRect "+viewRect);
						return GetOptStatus.EXIT_FAILURE;
						}
					return GetOptStatus.OK;
					}
			default:return super.handleOtherOptions(c, opt);
			}
		
		}
	
	protected LeafNode createLeafNode(BranchNode parent,String label)
		{
		return new DefaultLeafNode(null, label);
		}
	
	protected int buildFactoryChain(String path)
		{
		
		if(path==null)
			{
			warning("No path defined!");
			this.root=createLeafNode(null,"");
			}
		else
			{
			for(String p: path.split("[ \t/,;]+"))
				{
				NodeFactory nf=findFactoryByName(p);
				if(nf==null)
					{
					error("Cannot get type \'"+p+"\' in "+getAllAvailableFactories());
					return -1;
					}
				if(this.nodeFactoryChain.isEmpty())
					{
					this.root=nf.createBranch(null,"");
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
