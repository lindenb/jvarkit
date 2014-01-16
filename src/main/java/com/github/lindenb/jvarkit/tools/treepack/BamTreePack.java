package com.github.lindenb.jvarkit.tools.treepack;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.XMLConstants;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

public class BamTreePack extends AbstractCommandLineProgram
	{
	private Rectangle viewRect=new Rectangle(1000,1000);
	private final String SVG="http://www.w3.org/2000/svg";
	private final String XLINK="http://www.w3.org/1999/xlink";
	private List<NodeFactory> nodeFactoryChain=new ArrayList<NodeFactory>();
	private enum FactoryType {sample,chrom,mapq,properlypaired};
	private Node root=null;
	private interface NodeFactory
		{
		public FactoryType getType();
		public Node create(Branch parent);
		}
	
	
	
	/** Group by Sample */
	private class SampleNodeFactory implements NodeFactory
		{
		private class MyNode extends Branch
			{
			private MyNode(Branch parent)
				{
				super(parent);
				}
			
			@Override
			public void watch(SAMFileHeader header,SAMRecord seq)
				{
				SAMReadGroupRecord g=seq.getReadGroup();
				String sample=(g==null?"*":g.getSample());
				if(sample==null || sample.isEmpty()) sample="*";
				this.get(sample).watch(header, seq);
				}
			}
		@Override
		public FactoryType getType()	{return FactoryType.sample;}
		public Node create(Branch parent)
			{
			return new MyNode(parent);
			}
		}
	
	/** Group by Chromosome */
	private class ChromosomeNodeFactory implements NodeFactory
		{
		private class MyNode extends Branch
			{
			private MyNode(Branch parent)
				{
				super(parent);
				}
			
			@Override
			public void watch(SAMFileHeader header,SAMRecord seq)
				{
				String chr=seq.getReferenceName();
				if(chr==null) chr=SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
				this.get(chr).watch(header, seq);
				}
			}
		@Override
		public FactoryType getType()	{return FactoryType.chrom;}
		public Node create(Branch parent)
			{
			return new MyNode(parent);
			}
		}
	
	/** Group by MAPQ */
	private class MappingQualityFactory implements NodeFactory
		{
		private class MyNode extends Branch
			{
			private MyNode(Branch parent)
				{
				super(parent);
				}
			
			@Override
			public void watch(SAMFileHeader header,SAMRecord seq)
				{
				String qualStr="*";
				if(!seq.getReadUnmappedFlag())
					{
					int qual=seq.getMappingQuality();
					if(qual==255)
						{
						qualStr="?";
						}
					else
						{
						qual=((int)(qual/10.0))*10;
						qualStr=String.valueOf(qual);
						}
					}
				this.get(qualStr).watch(header, seq);
				}
			}
		@Override
		public FactoryType getType()	{return FactoryType.mapq;}

		public Node create(Branch parent)
			{
			return new MyNode(parent);
			}
		}

	/** Group by MAPQ */
	private class PropertyPairedFactory implements NodeFactory
		{
		private class MyNode extends Branch
			{
			private MyNode(Branch parent)
				{
				super(parent);
				}
			
			@Override
			public void watch(SAMFileHeader header,SAMRecord seq)
				{
				String s="";
				if(seq.getReadPairedFlag())
					{
					s=(seq.getReadPairedFlag()?"properly-paired":"not-property-paired");
					}
				this.get(s).watch(header, seq);
				}
			}
		@Override
		public FactoryType getType()	{return FactoryType.properlypaired;}

		public Node create(Branch parent)
			{
			return new MyNode(parent);
			}
		}
	
	
	
	private abstract class Node
		extends AbstractTreePack
		{
		protected Node(Node parent)
			{
			super(parent);
			}
		public abstract void watch(SAMFileHeader header,SAMRecord seq);
		public abstract boolean isLeaf();
		public void layout()
			{
			
			}
		
		void svg(XMLStreamWriter w,String label)throws XMLStreamException
		   {
		   if(this.isLeaf()) return;
		   final Rectangle2D bounds=this.getBounds();
		   w.writeStartElement("svg","g",SVG);
		   	
		   w.writeEmptyElement("svg", "rect", SVG);
		   if(label!=null) w.writeAttribute("title",label);
		  
		   w.writeAttribute("x",String.valueOf(bounds.getX()));
		   w.writeAttribute("y",String.valueOf(bounds.getY()));
		   w.writeAttribute("width",String.valueOf(bounds.getWidth()));
		   w.writeAttribute("height",String.valueOf(bounds.getHeight()));
		 
		  
		   Branch me=(Branch)this;
		   for(String key:me.children.keySet())
			   {
			   Node c=me.children.get(key);
			   if(!c.isLeaf())
				   {
				   c.svg(w,key);
				   }
			   }
		   if(label!=null)
			   {
			   //todo hershey
			   }
			   
			   
		 
		   
		   w.writeEndElement();//g
		   }
		}
	
	
	
	private abstract class Branch extends Node
		{
		private Map<String,Node> children=new HashMap<String,Node>();
	
		public Branch(BamTreePack.Branch parent)
			{
			super(parent);
			}
		public Collection<Node> getChildren()
			{
			return children.values();
			}
		
		public void layout()
			{
			Rectangle2D frame=getParent()==null ?
					new Rectangle2D.Double(0,0,viewRect.width,viewRect.height):
					getParent().getBounds();
			List<TreePack> L=new ArrayList<TreePack>(children.values());
			new  TreePacker().layout(L, frame 	);
			}
		
		@Override
		public double getWeight()
			{
			double N=0;
			for(TreePack c:this.getChildren()) N+=c.getWeight();
			return N;
			}
		
		protected Node get(String key)
			{
			Node c=this.children.get(key);
			if(c==null)
				{
				 if(getDepth()< BamTreePack.this.nodeFactoryChain.size())
					 {
					 c=nodeFactoryChain.get(getDepth()).create(this);
					 }
				 else
					 {
					 c=new Leaf(this);
					 }
				this.children.put(key, c);
				}
			return c;
			}

		@Override
		public final boolean isLeaf()
			{
			return false;
			}

		}
	

	private class Leaf extends Node
		{
		protected long count=0L;
		public Leaf(Branch parent)
			{
			super(parent);
			}
		@Override
		public void watch(SAMFileHeader header, SAMRecord seq)
			{
			++count;
			}
		@Override
		public double getWeight()
			{
			return count;
			}
		
		@Override
		public final boolean isLeaf()
			{
			return true;
			}
		}

	
	
	private BamTreePack()
		{
		
		}
	
	
	 
 
	  private void svg(Node root) throws XMLStreamException
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
		root.svg(w,null);
		w.writeEndElement();//svg
		w.writeEndDocument();
		w.flush();
	  	}

	 private void scan(SAMFileReader sfr)
		 {
		 sfr.setValidationStringency(ValidationStringency.LENIENT);
		 SAMFileHeader header=sfr.getFileHeader();
		 SAMRecordIterator iter=sfr.iterator();
		 while(iter.hasNext())
			 {
			 root.watch(header,iter.next());
			 }
		 }
	  
	@Override
	public void printOptions(PrintStream out)
		{
		super.printOptions(out);
		}
	
	private static void register(Map<String,NodeFactory> map,NodeFactory nf)
		{
		map.put(nf.getType().name(), nf);
		}
	
	@Override
	public int doWork(String[] args)
		{
		Map<String,NodeFactory> factorymap=new HashMap<String,NodeFactory>();
		register(factorymap,new ChromosomeNodeFactory());
		register(factorymap,new SampleNodeFactory());
		register(factorymap,new PropertyPairedFactory());
		register(factorymap,new MappingQualityFactory());
		String path=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "e:"))!=-1)
			{
			switch(c)
				{
				case 'e': path=opt.getOptArg();break;
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
		if(path==null)
			{
			warning("No path defined!");
			root=new Leaf(null);
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
					root=nf.create(null);
					}
				this.nodeFactoryChain.add(nf);
				}
			if(root==null)
				{
				error("Wrong path");
				return -1;
				}
			}
		
		SAMFileReader sfr=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading stdin");
				sfr=new SAMFileReader(System.in);
				scan(sfr);
				CloserUtil.close(sfr);
				}
			else
				{
				for(int optind =opt.getOptInd(); optind < args.length;++optind)
					{
					File f=new File(args[optind]);
					info("Reading stdin");
					sfr=new SAMFileReader(f);
					scan(sfr);
					CloserUtil.close(sfr);
					}
				}
			root.setBounds(new Rectangle2D.Double(0, 0, this.viewRect.getWidth(), this.viewRect.getHeight()));
			root.layout();
			svg(root);
			return 0;
			}
		catch (Exception e)
			{
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
			}
		}
	
	public static void main(String[] args)
		{
		new BamTreePack().instanceMainWithExit(args);

		}

	}
