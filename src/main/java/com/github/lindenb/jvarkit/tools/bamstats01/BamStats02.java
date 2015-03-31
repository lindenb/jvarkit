
/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.tools.treepack.TreePack;
import com.github.lindenb.jvarkit.tools.treepack.TreePacker;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFlag;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

public class BamStats02
	extends AbstractKnimeApplication
	{
	private File bedFile=null;
	private IntervalTreeMap<Boolean> intervals=null;
    private Counter<Category> counter=new Counter<>();
   
	
    private interface NodeFactory
		{
		public String toString();
		String label(CategoryAndCount cac);
		}
    private static class SamFlagNodeFactory
    	implements NodeFactory
		{
		SamFlag flg;
		SamFlagNodeFactory(SamFlag flg)
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
			return (this.flg.isSet(cac.flag)?this.flg.getLabel():"UNSET");
			}
		}
    
    private static abstract class AbstractTreePack implements TreePack
		{
		String label=null;
		AbstractTreePack parent=null;
		abstract List<CategoryAndCount> getCategoryAndCountList();
		boolean containsSelection=false;
		
		void paint(Graphics2D g,Rectangle2D area,CategoryAndCount selected,List<NodeFactory> factories,int factoryIndex)
			{
			System.err.println("Repainting");
			
			Hershey hershey=new Hershey();
			g.setColor(this.containsSelection?Color.ORANGE:Color.WHITE);
			g.fill(area);
			g.setColor(Color.BLACK);
			g.draw(area);
			
			double dw=area.getWidth()*0.05;
			double dh=area.getHeight()*0.05;
			//shrink
			area=new Rectangle2D.Double(
					area.getX()+dw,
					area.getY()+dh,
					area.getWidth()-dw*2.0,
					area.getHeight()-dh*2.0
					);
			
			if(label!=null && !label.trim().isEmpty())
				{
				
				double labelH=area.getHeight()*0.1;
				hershey.paint(g, this.label,
						new Rectangle2D.Double(
								area.getX(),
								area.getY()+labelH*0.05,
								area.getWidth(),
								labelH-(labelH*0.05)*2
								) 
						);
				
				area =new Rectangle2D.Double(
						area.getX(),
						area.getY()+labelH,
						area.getWidth(),
						area.getHeight()-labelH*2
						);
				}
			
			if(factoryIndex<factories.size() && !this.getCategoryAndCountList().isEmpty())
				{
				NodeFactory f=factories.get(factoryIndex);
				Map<String,TreePack> label2node=new HashMap<>();
				for(CategoryAndCount cac:this.getCategoryAndCountList())
					{
					String label=f.label(cac);
					if(label==null || label.isEmpty()) continue;
					DefaultTreePack dtp= (DefaultTreePack)label2node.get(label);
					if(dtp==null) 
						{
						dtp =new DefaultTreePack();
						dtp.label=label;
						dtp.categoriesAndCount=new ArrayList<>();
						dtp.parent=this;
						label2node.put(label,dtp);
						}
					if(cac==selected) dtp.containsSelection=true;
					dtp.categoriesAndCount.add(cac);
					dtp.weight+=cac.count;
					}
				List<TreePack> children= new ArrayList<>(label2node.values());
				TreePacker packer=new TreePacker();
				packer.layout(children, area);
				for(TreePack tp:children)
					{
					DefaultTreePack.class.cast(tp).paint(g, tp.getBounds(),selected, factories,factoryIndex+1);
					}
				}
			else
				{
				hershey.paint(g, String.valueOf(getWeight()), area);
				}
			}
		}
    
    private static class DefaultTreePack extends AbstractTreePack
    	{
    	List<CategoryAndCount> categoriesAndCount;
    	Rectangle2D rec=new Rectangle2D.Double();
    	double weight=0.0;
    	@Override
    	public Rectangle2D getBounds() {
    		return rec;
    		}
    	@Override
    	public double getWeight() {
    		return weight;
    		}
    	@Override
    	public void setBounds(Rectangle2D bounds) {
    		this.rec=bounds;
    		}
    	@Override
    	List<CategoryAndCount> getCategoryAndCountList()
    		{
    		return categoriesAndCount;
    		}
    	}
    
    
    private abstract class SampleNameCategoryFilter
    	implements NodeFactory
		{
    	List<DefaultTreePack> filter(List<CategoryAndCount> cats)
	    	{
	    	Set<String> labels=new HashSet<>();
	    	for(CategoryAndCount c:cats) labels.add(c.sampleName);
	    	List<DefaultTreePack> array=new ArrayList<>(labels.size());
	    	for(String label:labels)
	    		{
	    		DefaultTreePack dtp=new DefaultTreePack();
	    		dtp.label=label;
	    		for(CategoryAndCount c:cats)
	    			{
	    			if(c.sampleName.equals(label));
	    			}
	    		}
	    	return array;
	    	}
		}

    
    private static class Category
    	{
    	String filename=".";
    	String sampleName=".";
    	String chromosome=".";
    	int mapq;
    	int inTarget;
    	int flag;
    	
    	private void print(PrintWriter pw)
    		{
    		pw.print(filename);
    		pw.print("\t");
    		pw.print(sampleName);
    		pw.print("\t");
    		pw.print(chromosome);
    		pw.print("\t");
    		pw.print(mapq);
    		pw.print("\t");
    		pw.print(inTarget);
    		for(SamFlag flg:SamFlag.values())
    			{
    			pw.print("\t");
    			pw.print(flg.isSet(this.flag)?1:0);
    			}
    		}
    	
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((chromosome == null) ? 0 : chromosome.hashCode());
			result = prime * result
					+ ((filename == null) ? 0 : filename.hashCode());
			result = prime * result + flag;
			result = prime * result + inTarget;
			result = prime * result + mapq;
			result = prime * result
					+ ((sampleName == null) ? 0 : sampleName.hashCode());
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Category other = (Category) obj;
			if (chromosome == null) {
				if (other.chromosome != null)
					return false;
			} else if (!chromosome.equals(other.chromosome))
				return false;
			if (filename == null) {
				if (other.filename != null)
					return false;
			} else if (!filename.equals(other.filename))
				return false;
			if (flag != other.flag)
				return false;
			if (inTarget != other.inTarget)
				return false;
			if (mapq != other.mapq)
				return false;
			if (sampleName == null) {
				if (other.sampleName != null)
					return false;
			} else if (!sampleName.equals(other.sampleName))
				return false;
			return true;
		}
    	
    	}

    
    private static class CategoryAndCount extends Category
    	{
    	long count=0;
    	}
    
    private static class ViewDialog extends JFrame
    	{
    	private CatTableModel catTableModel=null;
    	private List<JComboBox<NodeFactory>> comboBoxes=new ArrayList<JComboBox<NodeFactory>>();
    	private JPanel drawingArea=null;
    	private JTable catTable;
    	
    	private class RootTreePack
    		extends AbstractTreePack
    		{
    		@Override
    		List<CategoryAndCount> getCategoryAndCountList()
    			{
    			return catTableModel.getRows();
    			}
    		@Override
    		public Rectangle2D getBounds()
    			{
    			return new Rectangle2D.Double(0,0,drawingArea.getWidth()+1,drawingArea.getHeight()+1);
    			}
    		@Override
    		public double getWeight()
    			{
    			double n=0;
    			for(CategoryAndCount cac:catTableModel.getRows())
    				{
    				n+=cac.count;
    				}
    			return n;
    			}
    		@Override
    		public void setBounds(Rectangle2D bounds)
    			{
    			throw new IllegalStateException();
    			}
    		
    		}
    	
    	private class CatTableModel extends AbstractGenericTable<CategoryAndCount>
    		{
    		CatTableModel(List<CategoryAndCount> rows)
    			{
    			super(rows);
    			}
    		
    		@Override
    		public int getColumnCount()
    			{
    			return 6+SamFlag.values().length;
    			}
    		
    		@Override
    		public Class<?> getColumnClass(int columnIndex) {

    			switch(columnIndex)
    				{
    				case 0: return String.class;
    				case 1: return String.class;
    				case 2: return String.class;
    				case 3: return Integer.class;
    				case 4: return Boolean.class;
    				default:
    					columnIndex-=5;
    					if(columnIndex< SamFlag.values().length)
    						{
    						return Boolean.class;
    						}
    					return Long.class;
    				}
    			  }
    		
    		@Override
    		public Object getValueOf(CategoryAndCount o, int columnIndex)
    			{
    			switch(columnIndex)
    				{
    				case 0: return o.filename;
    				case 1: return o.sampleName;
    				case 2: return o.chromosome;
    				case 3: return o.mapq;
    				case 4: return (o.inTarget==1?Boolean.TRUE:(int)o.inTarget==0?Boolean.FALSE:null);
    				default:
    					columnIndex-=5;
    					if(columnIndex< SamFlag.values().length)
    						{
    						return SamFlag.values()[columnIndex].isSet(o.flag);
    						}
    					return o.count;
    				}
    			}
    		
    		public String getColumnName(int column)
    			{
    			switch(column)
					{
					case 0: return "File";
					case 1: return "Sample";
					case 2: return "Chrom";
					case 3: return "Mapq";
					case 4: return "In/Off target";
					default:
						column-=5;
						if(column< SamFlag.values().length)
							{
							return SamFlag.values()[column].name();
							}
						return "Count";
					}
    			}
    		
    		}
    	
    	ViewDialog(List<CategoryAndCount> categories )
    		{
    		super("BamStats02");
    		this.catTableModel=new CatTableModel(categories);
    		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    		this.catTable=new JTable(this.catTableModel);
    		this.catTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    		JScrollPane scroll=new JScrollPane(this.catTable);
    		
    		JPanel topPane=new JPanel(new GridLayout(1,0,1,1));
    		
    	
    		JPanel contentPane=new JPanel(new BorderLayout(5, 5));
    		
    		contentPane.add(topPane,BorderLayout.NORTH);
    		
    		this.drawingArea=new JPanel(null,true)
    			{
    			protected void paintComponent(java.awt.Graphics c)
					{
					paintDrawingArea(Graphics2D.class.cast(c));
					}
    			};
        	this.drawingArea.setOpaque(true);
        	this.drawingArea.setBackground(Color.WHITE);
        	JSplitPane split=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,scroll, this.drawingArea);
        	contentPane.add(split, BorderLayout.CENTER);
    		
    		setContentPane(contentPane);
    		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
    		
    		List<NodeFactory> nodeFactories=new ArrayList<NodeFactory>();
    		nodeFactories.add(null);
    		nodeFactories.add(new NodeFactory()
				{
				@Override
				public String toString()
					{
					return "Filename";
					}
				@Override
				public String label(CategoryAndCount cac)
					{
					return cac.filename;
					}
				});
    		
    		nodeFactories.add(new NodeFactory()
				{
				@Override
				public String toString()
					{
					return "Sample Name";
					}
				@Override
				public String label(CategoryAndCount cac)
					{
					return cac.sampleName;
					}
				});
   		
    		nodeFactories.add(new NodeFactory()
				{
				@Override
				public String toString()
					{
					return "Chromosome";
					}
				@Override
				public String label(CategoryAndCount cac)
					{
					if(cac.chromosome==null || cac.chromosome.isEmpty() || cac.chromosome.equals(".")) return null;
					return cac.chromosome;
					}
				});
    		
    		nodeFactories.add(new NodeFactory()
				{
				@Override
				public String toString()
					{
					return "MAPQ";
					}
				@Override
				public String label(CategoryAndCount cac)
					{
					if(cac.chromosome==null || cac.chromosome.isEmpty() || cac.chromosome.equals(".")) return null;
					return String.valueOf(cac.mapq);
					}
				});
    		
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
					if(cac.chromosome==null || cac.chromosome.isEmpty() || cac.chromosome.equals(".")) return null;
					if(cac.inTarget==1) return "IN_TARGET";
					if(cac.inTarget==0) return "OFF_TARGET";
					return null;
					}
				});
    		
    		for(SamFlag sf:SamFlag.values())
    			{
    			nodeFactories.add(new SamFlagNodeFactory(sf));
    			}
    		
    		ActionListener repainAction=new ActionListener()
				{
				@Override
				public void actionPerformed(ActionEvent arg0)
					{
					drawingArea.repaint();
					}
				};
    		for(int i=0;i<5;++i)
    			{
    			JComboBox<NodeFactory> combox=new JComboBox<NodeFactory>(
    					nodeFactories.toArray(new NodeFactory[nodeFactories.size()])
    					);
    			combox.setFont(new Font("Dialog",Font.PLAIN,9));
    			combox.addActionListener(repainAction);
    			this.comboBoxes.add(combox);
    			topPane.add(combox);
    			}
    		this.catTable.getSelectionModel().addListSelectionListener(new ListSelectionListener()
				{
				@Override
				public void valueChanged(ListSelectionEvent arg0)
					{
					drawingArea.repaint();
					}
				});
    		}
    	
    	void paintDrawingArea(Graphics2D g)
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
    		RootTreePack rt=new RootTreePack();
    		Rectangle2D.Double area=new Rectangle2D.Double(0, 0, this.drawingArea.getWidth(), this.drawingArea.getHeight());
    		
    		int r=this.catTable.getSelectedRow();
    		if(r!=-1) r=this.catTable.convertRowIndexToModel(r);    		
    		rt.paint(g, area, (r==-1?null:this.catTableModel.getElementAt(r)),nodeFactories,0);
    		}
    	}
    
	public BamStats02()
		{
		
		}
	
	public void setBedFile(File bedFile) {
		this.bedFile = bedFile;
	}
	
	@Override
	public String getProgramDescription()
		{
		return "Statistics about the flags and reads in a BAM.";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"BamStats02";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-o (file) output. Default: stdout");
		super.printOptions(out);
		}
	
	
	
	private void run(String filename,SamReader r)
		{
		SAMRecordIterator iter=null;
		try
			{
			iter=r.iterator();
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(r.getFileHeader().getSequenceDictionary());
			while(iter.hasNext())
				{
				SAMRecord record=progress.watch(iter.next());
				Category cat=new Category();
				cat.filename=filename;
				cat.flag = record.getFlags();
				cat.inTarget = -1;
				cat.mapq=0;
				SAMReadGroupRecord g=record.getReadGroup();
				String sampleName=(g!=null?g.getSample():null);
				cat.sampleName=(sampleName==null?".":sampleName);
				if(!record.getReadUnmappedFlag())
					{
					cat.mapq=(int)(Math.ceil(record.getMappingQuality()/10.0)*10);
					cat.chromosome=record.getReferenceName();
					if(this.intervals!=null)
						{
						if(this.intervals.containsOverlapping(
								new Interval(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd())
								))
								{
								cat.inTarget=1;
								}
						else
								{
								cat.inTarget=0;
								}
						}
					
					}
				this.counter.incr(cat);
				}
			progress.finish();
			}
		finally
			{
			CloserUtil.close(iter);
			}
		
		}
	
	@SuppressWarnings("resource")
	@Override
	public int executeKnime(List<String> args)
		{
		SamReader samFileReader=null;
		PrintWriter out=null;
		try
			{
			if(bedFile!=null)
				{
				this.intervals=new IntervalTreeMap<Boolean>();
				Pattern tab=Pattern.compile("[\t]");
				String line;
				BufferedReader bedIn=IOUtils.openFileForBufferedReading(bedFile);
				while((line=bedIn.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					String tokens[]=tab.split(line,5);
					if(tokens.length<3) throw new IOException("bad bed line in "+line+" "+this.bedFile);
					String crhrom= tokens[0];
					int chromStart1= 1+Integer.parseInt(tokens[1]);//add one
					int chromEnd1= Integer.parseInt(tokens[2]);
					intervals.put(new Interval(crhrom, chromStart1, chromEnd1),Boolean.TRUE);
					}
				bedIn.close();
				}
			

			
			SamReaderFactory srf=SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.SILENT);
			
			if(args.isEmpty())
				{
				info("Reading from stdin");
				samFileReader= srf.open(SamInputResource.of(System.in));
				run("stdin",samFileReader);
				samFileReader.close();
				samFileReader=null;
				}
			else
				{
				for(String filename:IOUtils.unrollFiles(args))
					{
					info("Reading from "+filename);
					samFileReader=srf.open(new File(filename));
					run(filename,samFileReader);
					samFileReader.close();
					samFileReader=null;
					}
				}
			
			
			out = 	(
					getOutputFile()==null ?
					new PrintWriter(System.out) :
					new PrintWriter(getOutputFile())
					);
			out.print("#filename");
    		out.print("\t");
    		out.print("sampleName");
    		out.print("\t");
    		out.print("chromosome");
    		out.print("\t");
    		out.print("mapq");
    		out.print("\t");
    		out.print("inTarget");
    		for(SamFlag flg:SamFlag.values())
    			{
    			out.print("\t");
    			out.print(flg.name());
    			}
    		out.print("\t");
    		out.print("count");
			out.println();
			
			for(Category cat:this.counter.keySetDecreasing())
				{
				cat.print(out);
	    		out.print("\t");
	    		out.print(this.counter.count(cat));
				out.println();
				}
			
			out.flush();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(samFileReader);
			CloserUtil.close(out);
			}
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:B:"))!=-1)
			{
			switch(c)
				{
				case 'o': setOutputFile(opt.getOptArg());break;
				case 'B': this.setBedFile(new File(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
		}
	
	public static void main(String[] args)
		{
		if(!(args.length==2 && args[0].equals("gui")))
			{
			new BamStats02().instanceMainWithExit(args);
			}
		else
			{
			List<CategoryAndCount> categories = new ArrayList<>();
			BufferedReader in=null;
			try {
				Pattern tab=Pattern.compile("[\t]");
				in=IOUtils.openURIForBufferedReading(args[1]);
				String line;
				while((line=in.readLine())!=null)
					{
					if(line.startsWith("#")) continue;
					if(line.isEmpty()) continue;
					CategoryAndCount cat=new CategoryAndCount();
					String tokens[]=tab.split(line);
					int col=0;
					cat.filename  =tokens[col++];
					cat.sampleName=tokens[col++];
					cat.chromosome=tokens[col++];
					cat.mapq=Integer.parseInt(tokens[col++]);
					cat.inTarget=Integer.parseInt(tokens[col++]);
					cat.flag=0;
					for(SamFlag sf:SamFlag.values())
						{
						String ok=tokens[col++];
						if(ok.equals("0"))
							{
							//nothing
							}
						else if(ok.equals("1"))
							{
							cat.flag+=sf.getFlag();
							}
						else
							{
							System.err.println("Illegal Flag Column for "+sf+" in "+line);
							System.exit(-1);
							}
						}
					cat.count=Long.parseLong(tokens[col++]);
					categories.add(cat);
					}
				final ViewDialog f=new ViewDialog(categories);
				Dimension screen=Toolkit.getDefaultToolkit().getScreenSize();
				f.setBounds(50, 50, screen.width-100, screen.height-100);
				SwingUtilities.invokeAndWait(new Runnable()
						{
						@Override
						public void run() {
							f.setVisible(true);
						}
					});
				}
			catch (Exception e)
				{
				e.printStackTrace();
				System.exit(-1);
				}
			}
		}
	}
