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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.treepack;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicReference;

import javax.script.CompiledScript;
import javax.script.ScriptException;
import javax.script.SimpleBindings;
import javax.xml.XMLConstants;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.svg.SVG;



public abstract class AbstractTreePackCommandLine
	extends Launcher
	{
	private static final Logger LOG = Logger.build(AbstractTreePackCommandLine.class).make();

	
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	protected File outputFile = null;

	
	protected Rectangle viewRect=new Rectangle(1000,1000);
	private final String XLINK="http://www.w3.org/1999/xlink";
	private final String JVARKIT_NS="https://github.com/lindenb/jvarkit";
	protected RootNodeFactory nodeFactoryChain=new RootNodeFactory();
	protected final TreeNode rootNode =new TreeNode("(root)",nodeFactoryChain,null,null);
	private final Hershey hershey=new Hershey();   
	protected final SimpleBindings bindings = new SimpleBindings();
	private final AtomicReference<String> styleReference = new AtomicReference<String>(null);

	/** Node Factory */
	protected abstract class NodeFactory
		{
		public abstract String getName();
		public NodeFactory next=null;
		public NodeFactory prev=null;
		@Override
		public String toString()
			{
			return getName();
			}
		public NodeFactory append(final NodeFactory child) {
			if(next!=null) return next.append(child);
			this.next=child;
			child.prev=this;
			return child;
		}
		public TreeNode increment(final TreeNode parentNode,String key) {
			final TreeNode child;
			if(parentNode.factory!=this.prev) throw new IllegalStateException();
			if(!parentNode.children.containsKey(key)) {
				 final String style = styleReference.get();
				 styleReference.set(null);
				 child =new TreeNode(key, this,parentNode,style);
				 parentNode.children.put(key,child);
			} else {
				child = parentNode.children.get(key);
			}
			child.count++;
			return child;
			}
		public abstract void watch(final TreeNode parentNode,final Object o);
		}
	
	/** Root Node Factory */
	protected final class RootNodeFactory extends NodeFactory
		{
		public String getName() { return "ALL";}
		public void watch(final TreeNode parentNode,final Object o) {
			if(parentNode!=rootNode) throw new IllegalStateException();
			if(next==null) return;
			next.watch(parentNode, o);
			}
		}
	/** Root Node Factory */
	protected class JsNodeFactory extends NodeFactory
		{
		protected final String name;
		protected final CompiledScript script;
		public JsNodeFactory(final String name,CompiledScript script) {
			this.name=name;
			this.script=script;
		}
		public String getName() { return name;}
		public void watch(final TreeNode parentNode,final Object o) {
			if(parentNode.factory!=this.prev) throw new IllegalStateException();
			try {
				final Object jso=script.eval(bindings);
				if(jso==null) return;
				final String key = String.valueOf(jso);
				if(key==null || key.trim().isEmpty()) return;
				final TreeNode child = increment(parentNode, key);
				if(next!=null) next.watch(child, o);
			} catch (ScriptException e) {
				throw new RuntimeException(e);
			}
			
			}
		}
	
	protected final class TreeNode
		implements TreePack
		{
		private Rectangle2D bounds=new Rectangle2D.Double(-10,-10,10,10);
		private final String label;
		/** css style */
		private final String style;
		protected final TreeNode parent;
		protected final NodeFactory factory;
		private Map<String,TreeNode> children=new HashMap<String,TreeNode>();
		private long count=0L;
		
		protected TreeNode(
				final String label,
				final NodeFactory factory,
				final TreeNode parent,
				final String style
				)
			{
			this.label=label;
			this.factory = factory; 
			this.parent = parent;
			this.style = style;
			if(this.factory==null) throw new IllegalStateException("Factory is null");
			if(this.label==null || this.label.trim().isEmpty()) throw new IllegalStateException("label is null");
			}
		
		public final boolean isRoot() {
			return this.parent == null;
		}
		
		protected String convertWeightToString()
			{
			return AbstractTreePackCommandLine.this.convertWeightToString(getWeight());
			}


		
		public String getPath()
			{
			final StringBuilder b=new StringBuilder();
			TreeNode curr=this;
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
		
		public TreeNode getParent()
			{
			return parent;
			}
		
		public int getDepth()
			{
			int d=0;
			TreeNode p=this;
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
			final Rectangle2D.Double r1=new Rectangle2D.Double();
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
		public double getWeight() {
			if(isLeaf()) {
				return this.count;
			}
			else
				{
				double N=0;
				for(final String k:this.children.keySet())
					{
					final TreePack c = this.children.get(k);
					if(c.getWeight()<=0) throw new IllegalStateException(getPath());
					N+=c.getWeight();
					}
				return N;
				}
		}
		
		public boolean isLeaf(){
			return this.factory.next==null;
		}
		
		private Rectangle2D getTitleFrame()
			{
			final Rectangle2D r=getBounds(); 
			final Rectangle2D.Double frame= new Rectangle2D.Double();
			frame.height=r.getHeight()*0.1;
			frame.y=r.getY();
			frame.width=r.getWidth();
			frame.x=r.getX();
			return frame;
			}
		private Rectangle2D getChildrenFrame()
			{
			final Rectangle2D r=getBounds(); 
			if(isRoot()) return r;
			final Rectangle2D.Double frame= new Rectangle2D.Double();
			frame.height=r.getHeight()*0.9;//yes again after
			frame.y=r.getMaxY()-frame.height;
			frame.height=r.getHeight()*0.85;//yes again 
			frame.width=r.getWidth()*0.95;
			frame.x=r.getX()+(r.getWidth()-frame.width)/2.0;
			return frame;
			}
	
	public void layout(final TreePacker packer)
		{
		if(isLeaf()) return;
		if(getWeight()<=0) return ;
		final List<TreePack> L=new ArrayList<TreePack>(children.values());
		packer.layout(L, getChildrenFrame() );
		for(final TreeNode c:this.children.values())
			{
			c.layout(packer);
			}
		}			
	
		
		public void svg(final XMLStreamWriter w) throws XMLStreamException
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
		   
		   w.writeStartElement("svg","g",SVG.NS);
		   w.writeAttribute("title",getPath()+"="+convertWeightToString());
		   
		   
		   w.writeAttribute("jvarkit",JVARKIT_NS,"path",getPath());
		   w.writeAttribute("jvarkit",JVARKIT_NS,"weight",convertWeightToString());
		   if(!isLeaf())
		   		{
		   		w.writeAttribute("jvarkit",JVARKIT_NS,"category",factory.getName());
		   		w.writeAttribute("jvarkit",JVARKIT_NS,"count",""+children.size());
		   		}
		   
		   w.writeEmptyElement("svg", "rect", SVG.NS);
		   
		   w.writeAttribute("class", "r"+(getDepth()%2));
		   if(!(this.style==null || this.style.trim().isEmpty())) {
			   w.writeAttribute("style",this.style);
		   }
		   w.writeAttribute("x",String.valueOf(frameUsed.getX()));
		   w.writeAttribute("y",String.valueOf(frameUsed.getY()));
		   w.writeAttribute("width",String.valueOf(frameUsed.getWidth()));
		   w.writeAttribute("height",String.valueOf(frameUsed.getHeight()));
			   
		   /* write banner */
		   if(!isLeaf() && !isRoot())
			   {			   
			   w.writeEmptyElement("svg", "path", SVG.NS);
			   w.writeAttribute("d",
						hershey.svgPath(getLabel()+"="+convertWeightToString(),
								this.getTitleFrame())
						);
			   }
		   

		   
		   if(isLeaf())
			   {
			   final Rectangle2D f_up=new Rectangle2D.Double(
			   			insets.getX(),insets.getY(),
			   			insets.getWidth(),insets.getHeight()/2.0
			   			);
			   	w.writeEmptyElement("svg", "path", SVG.NS);
				w.writeAttribute("d",
						hershey.svgPath(getLabel(),f_up)
						);	
				w.writeAttribute("class","lbla"+(getDepth()%2));

				final Rectangle2D f_down=new Rectangle2D.Double(
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
			   for(final String key:this.children.keySet())
				   {
				   this.children.get(key).svg(w);
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
	
	
	protected void svg(final File out) throws XMLStreamException,IOException
	  	{
	  	LOG.info("writing svg");
		XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
		XMLStreamWriter w= xmlfactory.createXMLStreamWriter(
				super.openFileOrStdoutAsStream(out),"UTF-8");
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
		w.writeComment("Author Pierre Lindenbaum");
		rootNode.svg(w);
		w.writeEndElement();//svg
		w.writeEndDocument();
		w.flush();
	  	}
	
	protected AbstractTreePackCommandLine()
		{
		this.bindings.put("style",this.styleReference);
		}
	
	protected String convertWeightToString(double weight)
		{
		NumberFormat format=NumberFormat.getInstance();
		return format.format((long)weight);
		}
	protected String intervalToString(int value,int step)
		{
		value=((int)(value/(double)step));
		return String.valueOf((value)*step+"-"+((value+1)*step));
		}

	protected void layout()
		{
		LOG.info("layout");
		this.rootNode.setBounds(new Rectangle2D.Double(0, 0, this.viewRect.getWidth(), this.viewRect.getHeight()));
		final TreePacker packer=new TreePacker();
		this.rootNode.layout(packer);
		}

	protected void setDimension(final String s){
		final int x=s.indexOf('x');
		if(x==-1)
			{
			throw new IllegalArgumentException("'x' missing in "+s);
			}
		this.viewRect.x=Integer.parseInt(s.substring(0, x));
		this.viewRect.y=Integer.parseInt(s.substring(x+1));
		if(viewRect.x<=0 ||this.viewRect.y<=0 )
			{
			throw new IllegalArgumentException("bad viewRect "+s);
			}
		}
	
	
	
	
	}
