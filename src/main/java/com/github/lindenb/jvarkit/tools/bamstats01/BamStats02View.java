
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
package com.github.lindenb.jvarkit.tools.bamstats01;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.regex.Pattern;

import javax.swing.AbstractAction;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.treepack.TreePack;
import com.github.lindenb.jvarkit.tools.treepack.TreePacker;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFlag;

/**

BEGIN_DOC


Explore those statistics with a graphical user interface. Left-click shows a popup menu used to include/exclude categories.



### Synopsis


```
$java -jar dist/bamstats02view.jar output_of_bamstats02.tsv
```



### Example


```
$java -jar dist/bamstats02view.jar output.tsv
```


![img](http://i.imgur.com/xC6grmV.jpg)

![img](http://i.imgur.com/cjItYY0.jp)



### See also

BamStats02


END_DOC
*/

@Program(name="bamstats02view",description="Statistics about the flags and reads in a BAM. Visualize date from BamStats02")
public class BamStats02View
	extends Launcher
	{
	private static final Logger LOG = Logger.build(BamStats02View.class).make();
    
   
	private static boolean isEmpty(String s)
		{
		return s==null || s.trim().isEmpty() || s.equals(".");
		}
	
    private interface NodeFactory
		{
		public String toString();
		String label(final CategoryAndCount cac);
		}
    
    private static class StringPropNodeFactory
		implements NodeFactory
		{
    	final BamStats02.STRING_PROPS prop;
    	StringPropNodeFactory(final BamStats02.STRING_PROPS prop)
    		{
    		this.prop=prop;
    		}
    	@Override
    	public String label(CategoryAndCount cac) {
    		String s= cac._strings[this.prop.ordinal()];
    		return isEmpty(s)?null:s;
    		}
    	@Override
    	public String toString() {
    		return this.prop.toString();
    		}
		}
    private static class IntPropNodeFactory
		implements NodeFactory
		{
    	final BamStats02.INT_PROPS prop;
    	IntPropNodeFactory(final BamStats02.INT_PROPS prop)
    		{
    		this.prop=prop;
    		}
    	@Override
    	public String label(CategoryAndCount cac) {
    		int s= cac._ints[this.prop.ordinal()];
    		return s<0?null:String.valueOf(s);
    		}
    	@Override
    	public String toString() {
    		return this.prop.toString();
    		}
    	
		}

    
    private static class SamFlagNodeFactory
    	implements NodeFactory
		{
    	SAMFlag flg;
		SamFlagNodeFactory(SAMFlag flg)
			{
			this.flg=flg;
			}
		@Override
		public String toString()
			{
			return "SAM Flag: "+flg.getLabel();
			}
		@Override
		public String label(CategoryAndCount cac)
			{
			boolean is_set=this.flg.isSet(cac.flag);
			//no chromosome, we should not consider this
			if(isEmpty(cac._strings[BamStats02.STRING_PROPS.chromosome.ordinal()]))
				{
				switch(this.flg)
					{
					case DUPLICATE_READ: 
					case READ_REVERSE_STRAND: 
					case SUPPLEMENTARY_ALIGNMENT:
					case SECONDARY_ALIGNMENT:
						return null;
					default:break;
					}
				}
			
			switch(this.flg)
				{
				case FIRST_OF_PAIR: return (is_set?"R1":"R2");
				case DUPLICATE_READ: return (is_set?"DUPLICATE":"UNSET");
				case READ_PAIRED: return (is_set?"PAIRED":"NORMAL");
				case READ_UNMAPPED: return (is_set?"READ UNMAPPED":"READ MAPPED");
				case MATE_UNMAPPED: return (is_set?"MATE UNMAPPED":"MATE MAPPED");
				case READ_REVERSE_STRAND: return (is_set?"READ MINUS":"READ PLUS");
				case MATE_REVERSE_STRAND: return (is_set?"MATE MINUS":"MATE PLUS");
				default: return (this.flg.isSet(cac.flag)?this.flg.getLabel():"UNSET");
				}
			
			}
		}
    
    
    
    /** used in GUI, Category+final count */
    private static class CategoryAndCount extends BamStats02.Category
    	{
    	long count=0L;
    	}
    
    /** GUI window */
   
	@SuppressWarnings("serial")
	private static class ViewDialog extends JFrame
    	{
    	private List<JComboBox<NodeFactory>> comboBoxes=new ArrayList<JComboBox<NodeFactory>>();
    	private JPanel drawingArea=null;
    	private List<CategoryAndCount> _categoryAndCounts=null;
    	private List<CategoryAndCount> _filteredCategoryAndCounts=null;
    	private AbstractTreePack root=null;
    	private final NumberFormat fmt = new DecimalFormat("#,###");
    	/** http://www.perbang.dk/rgbgradient/ */
    	private Color papers[]=new Color[]{
    	      new Color(135, 152, 222),  
    	      new Color(135, 152, 222),
    	      new Color(165, 189, 232),
    	      new Color(180, 207, 237),
    	      new Color(195, 226, 242),
    	      new Color(211, 245, 248)
    			}; 
    	
    	private class NodeFilter
    		{
    		NodeFactory nodeFactory;
    		String label;
    		NodeFilter _next=null;
    		boolean accept(CategoryAndCount cac)
    			{
    			if(cac==null) return false;
    			String s= this.nodeFactory.label(cac);
    			if(s==null || !s.equals(this.label)) return false;
    			return _next!=null?_next.accept(cac):true;
    			}
    		}
    	
    	/** abstract implementation of TreePack */
        private abstract class AbstractTreePack implements TreePack
    		{
        	
    		String label=null;
    		NodeFactory nodeFactory=null;
    		AbstractTreePack parent=null;
    		List<DefaultTreePack> childNodes=new ArrayList<>();
    		
    		
    		AbstractTreePack findTreePackAt(int x,int y)
    			{	
    			for(DefaultTreePack c:this.childNodes)
    				{
    				if(!c.getBounds().contains(x, y)) continue;
    				
    				return c.findTreePackAt(x,y);
    				}
    			return this;
    			}
    		
    		public int getDepth()
    			{
    			return parent==null?0:1+parent.getDepth();
    			}
    		
    		public NodeFilter getNodeFilter()
				{
				if(this.nodeFactory==null) return null;//e.g: root node
				NodeFilter filter=new NodeFilter();
				filter.label=this.label;
				filter.nodeFactory=this.nodeFactory;
				filter._next=null;
				return filter;
				}
    		
    		public NodeFilter getNodeFilters()
    			{
    			NodeFilter filter=null;
    			AbstractTreePack curr=this;
    			while(curr!=null)
    				{
    				NodeFilter filter2= curr.getNodeFilter();//e.g: root node
    				if(filter2!=null)
    					{
    					if(filter==null)
    						{
    						filter=filter2;
    						}
    					else
    						{
    						filter2._next=filter;
    						filter=filter2;
    						}
    					}
    				curr=curr.parent;
    				}
    			
    			return filter;
    			}
    		
    		
    		
    		public Rectangle2D getInsetArea()
    			{
    			Rectangle2D area=getBounds();
    			if(parent==null) return area;
    			double dw=area.getWidth()*0.05;
    			double dh=area.getHeight()*0.05;
    			return new Rectangle2D.Double(
    					area.getX()+dw,
    					area.getY()+dh,
    					area.getWidth()-dw*2.0,
    					area.getHeight()-dh*2.0
    					);
    			}
    		
    		public Rectangle2D getLabelArea()
    			{
    			Rectangle2D area=getInsetArea();
    			double labelH=area.getHeight()*0.1;
    			return new Rectangle2D.Double(
    					area.getX(),
    					area.getY()+labelH*0.05,
    					area.getWidth(),
    					labelH-(labelH*0.05)*2
    					);
    			}
    		
    		public Rectangle2D getChildrenArea()
    			{
    			Rectangle2D area=getInsetArea();
    			if(this.label==null || this.label.trim().isEmpty())
    				{
    				return area;
    				}
    			double labelH=area.getHeight()*0.1;
    			return new Rectangle2D.Double(
    					area.getX(),
    					area.getY()+labelH,
    					area.getWidth(),
    					area.getHeight()-labelH
    					);
    			}
    		}

    	
    	private class RootTreePack
    		extends AbstractTreePack
    		{
    		@Override
    		public Rectangle2D getBounds()
    			{
    			return new Rectangle2D.Double(0,0,drawingArea.getWidth()+1,drawingArea.getHeight()+1);
    			}
    		
    		@Override
    		public void setBounds(Rectangle2D bounds)
    			{
    			throw new IllegalStateException();
    			}
    		@Override
    		public double getWeight()
    			{
    			long weight=0L;
    			for(CategoryAndCount cac: ViewDialog.this._filteredCategoryAndCounts)
    				{
    				weight+=cac.count;
    				}
    			return weight;
    			}
    		
    		}
    	
    	 private class DefaultTreePack
    	 	extends AbstractTreePack
	     	{
    		 long weight=0L;
	     	Rectangle2D rec=new Rectangle2D.Double();
	     		     	
	     	@Override
	     	public Rectangle2D getBounds() {
	     		return rec;
	     		}
	     	@Override
	     	public void setBounds(Rectangle2D bounds) {
	     		this.rec=bounds;
	     		}
	     	@Override
	     	public double getWeight() {
	     		return weight;
	     		}
	     	}
    	
    	ViewDialog(List<CategoryAndCount> categories )
    		{
    		super("BamStats02");
    		this._categoryAndCounts= Collections.unmodifiableList(categories);
    		this._filteredCategoryAndCounts= this._categoryAndCounts;
    		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    		
    		JPanel topPane=new JPanel(new GridLayout(1,0,1,1));
    		
    	
    		JPanel contentPane=new JPanel(new BorderLayout(5, 5));
    		
    		contentPane.add(topPane,BorderLayout.NORTH);
    		
    		this.drawingArea=new JPanel(null,true)
    			{
    			protected void paintComponent(java.awt.Graphics c)
					{
					paintDrawingArea(Graphics2D.class.cast(c));
					}
    			@Override
    			public String getToolTipText(MouseEvent event)
    				{
    				if(ViewDialog.this.root==null) return null;
    				AbstractTreePack p = ViewDialog.this.root.findTreePackAt(event.getX(),event.getY());
    				if(p==null) return null;
    				return p.label+" ("+(long)p.getWeight()+")";
    				}
    			};
    		this.drawingArea.setToolTipText("");
    		this.drawingArea.addComponentListener(new ComponentAdapter()
    				{
	    			@Override
	    			public void componentResized(ComponentEvent e) {
	    				ViewDialog.this.root=null;
	    				ViewDialog.this.drawingArea.repaint();
	    				}
    				});
    		this.drawingArea.addMouseListener(new MouseAdapter() 
    			{
    			@Override
    			public void mousePressed(MouseEvent event) {
    				if(!event.isPopupTrigger()) return;
    				if(ViewDialog.this.root==null) return;
    				final AbstractTreePack focused = ViewDialog.this.root.findTreePackAt(event.getX(),event.getY());
    				if(focused==null) return;
					final NodeFilter nodeFilters=focused.getNodeFilters();
					final NodeFilter thisNodeFilter =focused.getNodeFilter();

    				JPopupMenu popup=new JPopupMenu(focused.label+" ("+(long)focused.getWeight()+")");
    				AbstractAction action=new AbstractAction("Focus on this category")
						{
						@Override
						public void actionPerformed(ActionEvent e)
							{
							
							if(nodeFilters==null) return;
							List<CategoryAndCount> L=new Vector<>();
							//create a new list of variant, set as the main filtered list
							for(CategoryAndCount cac: ViewDialog.this._categoryAndCounts)
								{
								if(nodeFilters.accept(cac))
									{
									L.add(cac);
									}
								}
							ViewDialog.this._filteredCategoryAndCounts=L;
							ViewDialog.this.root=null;
							ViewDialog.this.drawingArea.repaint();
							}
						};
    				action.setEnabled(nodeFilters!=null);
    				popup.add(new JMenuItem(action));
    				
    				
    				action=new AbstractAction("Hide this category")
						{
						@Override
						public void actionPerformed(ActionEvent e)
							{
							
							if(thisNodeFilter==null) return;
							List<CategoryAndCount> L=new Vector<>();
							//filter the existing a new list of variant, set as the main filtered list
							for(CategoryAndCount cac: ViewDialog.this._filteredCategoryAndCounts)
								{
								if(!thisNodeFilter.accept(cac))
									{
									L.add(cac);
									}
								}
							ViewDialog.this._filteredCategoryAndCounts=L;
							ViewDialog.this.root=null;
							ViewDialog.this.drawingArea.repaint();
							}
						};
					action.setEnabled(thisNodeFilter!=null);
	 				popup.add(new JMenuItem(action));
    				
    				action=new AbstractAction("Show All")
						{
						@Override
						public void actionPerformed(ActionEvent e)
							{
							//restore original list of variants
							ViewDialog.this.root=null;
							ViewDialog.this._filteredCategoryAndCounts=ViewDialog.this._categoryAndCounts;
							ViewDialog.this.drawingArea.repaint();
							}
						};
    				
    				popup.add(new JMenuItem(action));
    				popup.show(drawingArea, event.getX(),event.getY());
    				}
    			});
        	this.drawingArea.setOpaque(true);
        	this.drawingArea.setBackground(Color.WHITE);
        	contentPane.add(this.drawingArea, BorderLayout.CENTER);
    		
    		setContentPane(contentPane);
    		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
    		
    		List<NodeFactory> nodeFactories=new ArrayList<NodeFactory>();
    		nodeFactories.add(null);//mean: remove this factory
    		for(final BamStats02.STRING_PROPS p: BamStats02.STRING_PROPS.values())
    			{
    			nodeFactories.add(new StringPropNodeFactory(p));
    			}
    		for(final BamStats02.INT_PROPS p: BamStats02.INT_PROPS.values())
				{
    			if(p== BamStats02.INT_PROPS.inTarget) continue;//special, see below
				nodeFactories.add(new IntPropNodeFactory(p));
				}
    		
    		nodeFactories.add(new NodeFactory()
				{
				@Override
				public String toString()
					{
					return "IN/OFF target";
					}
				@Override
				public String label(CategoryAndCount cac)
					{
					if(cac.getChrom()==null) return null;
					int v= cac._ints[BamStats02.INT_PROPS.inTarget.ordinal()];
					if(v==1) return "IN_TARGET";
					if(v==0) return "OFF_TARGET";
					return null;
					}
				});
    		
    		nodeFactories.add(new NodeFactory()
				{
				@Override
				public String toString()
					{
					return "Disjoint Chromosomes";
					}
				@Override
				public String label(final CategoryAndCount cac)
					{
					String c1= cac._strings[BamStats02.STRING_PROPS.chromosome.ordinal()];
					String c2= cac._strings[BamStats02.STRING_PROPS.mate_chromosome.ordinal()];
					if(isEmpty(c1)) return null;
					if(isEmpty(c2)) return null;
					 if(c1.equals(c2))
						{
						return null;
						}
					else if(c1.compareTo(c2)<=0)
						{
						return c1+":"+c2;
						}
					else
						{
						return c2+":"+c1;
						}
					}
				});
    		
    		
    		for(SAMFlag sf:SAMFlag.values())
    			{
    			nodeFactories.add(new SamFlagNodeFactory(sf));
    			}
    		
    		ActionListener repainAction=new ActionListener()
				{
				@Override
				public void actionPerformed(ActionEvent arg0)
					{
					ViewDialog.this.root=null;
					drawingArea.repaint();
					}
				};
    		for(int i=0;i<6;++i)
    			{
    			JComboBox<NodeFactory> combox=new JComboBox<NodeFactory>(
    					nodeFactories.toArray(new NodeFactory[nodeFactories.size()])
    					);
    			combox.setFont(new Font("Dialog",Font.PLAIN,9));
    			combox.addActionListener(repainAction);
    			this.comboBoxes.add(combox);
    			topPane.add(combox);
    			}
    		
    		}
    	
    	
    	private List<NodeFactory> getSelectedNodeFactories()
    		{
    		List<NodeFactory> nodeFactories=new ArrayList<NodeFactory>();
    		for(JComboBox<NodeFactory> combox:this.comboBoxes)
    			{
    			int selIndex= combox.getSelectedIndex();
    			if(selIndex==-1) continue;
    			NodeFactory nodeFactory = combox.getItemAt(selIndex);
    			if(nodeFactory==null) continue;
    			nodeFactories.add(nodeFactory);
    			}
    		return nodeFactories;
    		}
    	
    	private void paintDrawingArea(Graphics2D g)
    		{
    		
    		if(this.root==null)
    			{
    			this.root=new RootTreePack();
    			List<NodeFactory> nodeFactories=getSelectedNodeFactories();
        		Rectangle2D.Double area=new Rectangle2D.Double(0, 0, this.drawingArea.getWidth(), this.drawingArea.getHeight());
        		
        		layoutPack(this.root,area,this._filteredCategoryAndCounts,nodeFactories,0);
    			}
    		
    		
    		paintPack(this.root,g);
    		}
    	
    	/** recursively paint each tree pack */
    	private void paintPack(AbstractTreePack root,Graphics2D g)
    		{
    		Rectangle2D area=(Rectangle2D)root.getBounds().clone();
			Hershey hershey=new Hershey();
			
			/* paint background */
			g.setColor(Color.LIGHT_GRAY);
			g.fill(area);
			g.setColor(papers[root.getDepth()%papers.length]);
			g.fill(new Rectangle2D.Double(
					area.getX(),area.getY(),
					(area.getWidth()>5?area.getWidth()-2:area.getWidth()),
					(area.getHeight()>5?area.getHeight()-2:area.getHeight())
					));
			g.setColor(Color.BLACK);
			g.draw(area);
			
			/* paint label */
			if(root.label!=null && !root.label.trim().isEmpty())
				{
				Rectangle2D lblArea=root.getLabelArea();
				double w= lblArea.getHeight()*root.label.length();
				if(lblArea.getWidth()> w)
					{
					
					lblArea=new Rectangle2D.Double(
							lblArea.getCenterX()-w/2.0,
							lblArea.getY(),
							w,
							lblArea.getHeight()
							);
					}
				g.setColor(root.getDepth()<=1?Color.WHITE:Color.BLACK);
				hershey.paint(g,
						root.label,
						lblArea
						);
				}
			
			/* paint children */
			area = root.getChildrenArea();
			if(!root.childNodes.isEmpty())
				{
				for(DefaultTreePack tp:root.childNodes)
					{
					paintPack(DefaultTreePack.class.cast(tp),g);
					}
				}
			else
				{
				g.setColor(root.getDepth()<=1?Color.WHITE:Color.BLACK);
				hershey.paint(g, fmt.format((long)root.getWeight()), area);
				}
    		}
    	
    	/** layout tree pack */
    	private void layoutPack(
    			AbstractTreePack root,
    			Rectangle2D area,
    			List<CategoryAndCount> categoryAndCounts,
    			List<NodeFactory> factories,
    			int factoryIndex
    			)
	    	{
    		if(root==null || factoryIndex>=factories.size()) return;
    		if(categoryAndCounts.isEmpty()) return;
    		NodeFactory nodeFactory=factories.get(factoryIndex);
			Map<String,TreePack> label2node=new HashMap<>();
			Map<String,List<CategoryAndCount>> label2treepackList=new HashMap<>();
			for(CategoryAndCount cac:categoryAndCounts)
				{
				String label=nodeFactory.label(cac);
				if(label==null || label.isEmpty()) continue;
				DefaultTreePack dtp= (DefaultTreePack)label2node.get(label);
				
				List<CategoryAndCount> cacL= label2treepackList.get(label);
				if(dtp==null) 
					{
					dtp =new DefaultTreePack();
					root.childNodes.add(dtp);
					dtp.label=label;
					dtp.parent=root;
					dtp.nodeFactory=nodeFactory;
					label2node.put(label,dtp);
					
					cacL = new ArrayList<>();
					label2treepackList.put(label, cacL);
					}
				cacL.add(cac);
				dtp.weight+=cac.count;
				}
			List<TreePack> children= new ArrayList<>(label2node.values());
			TreePacker packer=new TreePacker();
			
			area= root.getChildrenArea();
			packer.layout(children, area);
			
			for(final DefaultTreePack tp:root.childNodes)
				{
				layoutPack(
						DefaultTreePack.class.cast(tp),
						tp.getBounds(),
						label2treepackList.get(tp.label),
						factories,
						factoryIndex+1
						);
				label2treepackList.remove(tp.label);
				}
	    	}
    	
    	}
    
	private BamStats02View()
		{
		
		}
	
	
	
	@Override
	public int doWork(final List<String> args) {
		
			
		
		
			{
			List<CategoryAndCount> categories = new ArrayList<>();
			BufferedReader in=null;
			try {
				final Pattern tab=Pattern.compile("[\t]");
				final String inputFilename = oneAndOnlyOneFile(args);
				System.err.println("[LOG] Reading "+inputFilename);
				in=IOUtils.openURIForBufferedReading(inputFilename);
				String line;
				while((line=in.readLine())!=null)
					{
					if(line.startsWith("#")) continue;
					if(line.isEmpty()) continue;
					if(categories.size()%100000==0) System.err.println("[LOG] Reading N="+categories.size());
					CategoryAndCount cat=new CategoryAndCount();
					String tokens[]=tab.split(line);
					int col=0;
					for(BamStats02.STRING_PROPS p:BamStats02.STRING_PROPS.values())
						{
						if(col>=tokens.length) throw new IOException("Not enough columns in "+line);
						cat.set(p, tokens[col++]);
						}
					for(BamStats02.INT_PROPS p:BamStats02.INT_PROPS.values())
						{
						if(col>=tokens.length) throw new IOException("Not enough columns in "+line);
						String s= tokens[col++];
						cat.set(p,isEmpty(s)?-1:Integer.parseInt(s));
						}
					
					cat.flag=0;
					for(final SAMFlag sf:SAMFlag.values())
						{
						if(col>=tokens.length) throw new IOException("Not enough columns in "+line);
						final String ok=tokens[col++];
						if(ok.equals("0"))
							{
							//nothing
							}
						else if(ok.equals("1"))
							{
							cat.flag+=sf.intValue();
							}
						else
							{
							System.err.println("Illegal Flag Column for "+sf+" in "+line);
							System.exit(-1);
							}
						}
					if(col>=tokens.length) throw new IOException("Not enough columns in "+line);
					cat.count=Long.parseLong(tokens[col++]);
					categories.add(cat);
					}
				LOG.info("[LOG] Down reading "+inputFilename+" N="+categories.size());
				final ViewDialog f=new ViewDialog(categories);
				final Dimension screen=Toolkit.getDefaultToolkit().getScreenSize();
				f.setBounds(50, 50, screen.width-100, screen.height-100);
				SwingUtilities.invokeAndWait(new Runnable()
						{
						@Override
						public void run() {
							f.setVisible(true);
						}
					});
				return RETURN_OK;
				}
			catch (Exception e)
				{
				LOG.error(e);
				return -1;
				}
			}
		}
	
	public static void main(String[] args) {
		new BamStats02View().instanceMain(args);//no exit
	}

	}
