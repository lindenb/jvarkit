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
package com.github.lindenb.jvarkit.tools.sigframe;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.LinearGradientPaint;
import java.awt.Paint;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Pattern;

import javax.swing.AbstractAction;
import javax.swing.ActionMap;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JDesktopPane;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import javax.swing.filechooser.FileFilter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;
import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;



/**
 
BEGIN_DOC

##Input file:

tab delimited text file, compressed with **tabix/bgzip** and indexed with **tabix/tabix.**
Columns are:
* chrom (string)
* start (int)
* end (int)
* value (double)
* color (optional: 3 comma-separated integers [0-255] for R,G,B )

```tsv
chrM	0	10	0.005661	120,120,120
chrM	0	10	0.114468	
chrM	0	10	-0.877466
chrM	0	10	1.456670	120,220,120
chrM	0	10	-1.720711
chrM	1	11	-1.427848
chrM	1	11	-1.433984
chrM	1	11	1.891983
chrM	1	11	2.255700
```

## Synopsis

```
$ java -jar dist/sigframe.jar
```

or start from the command-line

```
$ java -jar dist/sigframe.jar -R ref.fasta tabix1.gz tabix2.gz  tabix3.gz 
```

END_DOC
 */
@SuppressWarnings("serial")

public class SigFrame
	extends JFrame
	{
	private static final Logger LOG=Logger.build(SigFrame.class).make();
	/** desktop */
	private JDesktopPane desktop;
	/** map string to action */
	private ActionMap actionMap=new ActionMap();
	/** reference genome */
	private SAMSequenceDictionary genome=null;
	/** about message */
	private String aboutMessage="<h2 align='center'>Pierre Lindenbaum PhD @yokofakun</h2>";
	
	
	/** tool */
	private enum Tool
		{
		ZOOM,
		SHOW_DATA_FOR_FILE,
		CHANGE_HEIGHT
		};
	
	/** a chrom/start/end/value */
	static private class SigData
		{
		String name;
		double value=0;
		String chrom;
		int start;
		int end;
		String color_string;
		private Color _color;
		
		public int getStart() { return start;}
		public int getEnd() { return end;}
		public String getChromosome() { return chrom;}
		public String getName() {return name;}
		public double getValue() { return value;}
		public Color getColor()
			{
			if(_color==null)
				{
				_color=Color.BLACK;
				if(color_string!=null)
					{
					String tokens[]=color_string.split("[,]");
					_color=new Color(
							Integer.parseInt(tokens[0]),
							Integer.parseInt(tokens[1]),
							Integer.parseInt(tokens[2])
							);
					
					}
				}
			return _color;
			}
		
		}
	
	/** a vector of signalfiles */
	class SigFiles extends Vector<SigDataFileTabixReader>
		{
		
		}
	
	/** A wrapper around a file indexed with tabix */
	class SigDataFileTabixReader extends AbstractTabixObjectReader<SigData>
		{
		private File file;
		SigDataFileTabixReader(final File file ) throws IOException
			{
			super(file.getPath());
			this.file=file;
			}
		
		public File getFile()
			{
			return file;
			}
		public boolean hasDataForChrom(final SAMSequenceRecord rec)
			{
			return iterator(rec.getSequenceName(), 1, 1+rec.getSequenceLength()).hasNext();
			}
		
		@Override
		protected Iterator<SigData> iterator(Iterator<String> delegate)
			{
			return new MyIterator(delegate);
			}
		
		/** convert line to  SigData */
	    private class MyIterator
    	extends AbstractMyIterator
	    	{
	    	private final Pattern tab=Pattern.compile("[\t]");
	    	MyIterator(final Iterator<String> delegate)
	    		{
	    		super(delegate);
	    		}
	    	
	    	@Override
	    	public SigData next()
	    		{
	    		if(!hasNext()) throw new IllegalStateException();
	    		final String tokens[]=this.tab.split(super.delegate.next());
	    		final SigData d=new SigData();
	    		d.chrom=tokens[0];
	    		d.start=Integer.parseInt(tokens[1]);
	    		d.end=Integer.parseInt(tokens[2]);
	    		d.value=Double.parseDouble(tokens[3]);
	    		if(tokens.length>4) d.color_string=tokens[4];
	    		return d;
	    		}
	    	}

		
		}
	
	
	/**
	 * LoadAtStartup
	 *
	 */
	private static class LoadAtStartup
		extends WindowAdapter
		{
		private Set<File> sources = new HashSet<File>();
		private SigFrame frame;
		LoadAtStartup(SigFrame frame)
			{
			this.frame=frame;
			}
		@Override
		public void windowOpened(WindowEvent e)
			{
			frame.removeWindowListener(this);
			if(sources.isEmpty()) return;
			frame.openTabixFiles(sources);	
			}
		}
	
	/**
	 * IFrame
	 *
	 */
	private class IFrame extends JInternalFrame
		{
		protected IFrame()
			{
			super("iframe",true,true,true,true);
			JPanel content=new JPanel(new BorderLayout(2,2));
			content.setBorder(new EmptyBorder(2, 2, 2, 2));
			setContentPane(content);
			this.addInternalFrameListener(new InternalFrameAdapter()
				{
				@Override
				public void internalFrameOpened(InternalFrameEvent e)
					{
					IFrame.this.removeInternalFrameListener(this);
					onOpen();
					}
				});
			}
		
		void onOpen()
			{
			final float ratio=0.9f;
			final Dimension dim= SigFrame.this.desktop.getSize();
			final Dimension d2=new Dimension(dim);
			d2.width=(int)(ratio*dim.width);
			d2.height=(int)(ratio*dim.height);
			this.setBounds(
					(int)(Math.random()*(dim.width-d2.width)),
					(int)(Math.random()*(dim.height-d2.height)),
					d2.width, d2.height);
			}
		}
	
	/**
	 * Browser: a browser is an iframe containing 
	 * multiple chromView, each chromView contains a file(tabix) view
	 *
	 */
	private class Browser extends IFrame
		implements MouseListener,MouseMotionListener
		{
		private final int TOP_DISTANCE=100;
		private final int LEFT_DISTANCE=50;
		private final int BORDER_SIZE=5;
		private final int DEFAULT_FILE_VIEW_HEIGHT=150;
		private JScrollBar hScrollBar;
		private JScrollBar vScrollBar;
		private JPanel drawingArea;
		private BufferedImage offscreen=null;
		private Point mouseStart=null;
		private Point mousePrev=null;
		private int max_chrom_length=0;
		private int view_start=0;
		private int view_length=0;
		private int position_in_history=0;
		private List<History> history=new ArrayList<History>();
		private List<ChromView> chromViews=new ArrayList<ChromView>();
		private Tool currentTool=Tool.ZOOM;
		private ChromView selectedChromView=null;
		private ChromView.FileView selectedFileFiew=null;
		private SAMSequenceDictionary samSequenceDictionary;
		/** all tabix files associated to this reader */
		private SigFiles sigilesList;
		
		/** select tool */
		private class SelectToolAction
			extends AbstractAction
			{
			SelectToolAction(Tool tool)
				{
				super(tool.name());
				this.putValue("sigframe.tool", tool);
				}
			@Override
			public void actionPerformed(ActionEvent e)
				{
				currentTool=(Tool)this.getValue("sigframe.tool");
				}
			}
		
		/** browser history */
		private class History
			{
			int view_start=0;
			int view_length=0;
			History(int start,int length)
				{
				this.view_start=start;
				this.view_length=length;
				}
			}
		
		/** Graphics Context */
		private class GC
			{
			int step=0;
			Graphics2D g;
			int top_y;
			int width;
			int height;
			
			GC()
				{
				top_y=vScrollBar.getValue();
				width=drawingArea.getWidth();
				height=drawingArea.getHeight();
				}
			}
		
		/** abstract view: offsets Y for the view */
		private abstract class View
			{
			int y_pixels=0;
			int height_pixels=0;
			public int getMaxY()
				{
				return getTopY()+getHeight();
				}
			public int getTopY()
				{
				return this.y_pixels;
				}
			public int getHeight()
				{
				return this.height_pixels;
				}
			
			public int getTopY(GC gc)
				{
				return this.getTopY()+topMargin()-gc.top_y;
				}
			}
		
		/**
		 * ChromView
		 */
		private class ChromView
			extends View
			{
			/** background color for this view */
			private Color background=Color.WHITE;
			/** chromosome for this view */
			private SAMSequenceRecord chrom;
			/** one file view for each tabix */
			private List<FileView> fileViews=new ArrayList<FileView>();
			
			
			
				
			/** inside a chromosome view: one view per BED file */
			class FileView
					extends View
				{
				/** associated tabix file */
				SigDataFileTabixReader tabixFile;
				/** min value Y */
				double min_value;
				/** max value Y */
				double max_value;
				
				
				public FileView(SigDataFileTabixReader sigFile)
					{
					this.tabixFile=sigFile;
					this.height_pixels=DEFAULT_FILE_VIEW_HEIGHT;
					this.min_value=-1.0;//sigChrom.getMinValue();
					this.max_value=1.0;//sigChrom.getMaxValue();
					}
				
				/** return owner */
				public ChromView getChromView()
					{
					return ChromView.this;
					}
				
				
				/** return owner */
				public Browser getBrowser()
					{
					return getChromView().getBrowser();
					}
				
				/** return associated chromosome */
				public SAMSequenceRecord getChromosome()
					{
					return getChromView().getChromosome();
					}
				
				/** returns wether that fileview contains x,y */
				public boolean contains(GC gc,int x,int y)
					{
					int left= leftMargin();
					if(y<topMargin()) return false;
					if(getTopY(gc)> gc.height) return false;
					if(x <left) return false;
					return new Rectangle(left, getTopY(gc), gc.width-left, getHeight()).contains(x,y);
					}
				
				/** fills a popup menu invoked in that menu */
				public void fillPopup(JPopupMenu popup)
					{
					popup.add(new JSeparator());
					popup.add(new AbstractAction("Adjust File Height....")
						{
						@Override
						public void actionPerformed(ActionEvent e)
							{
							Integer h=getBrowser().askHeight();
							if(h==null) return;
							FileView.this.height_pixels=h;
							getBrowser().recalcSize();
							}});
					popup.add(new AbstractAction("File Min/Max value....")
						{
						@Override
						public void actionPerformed(ActionEvent e)
							{
							double values[]=askMinValue(min_value,max_value);
							if(values==null) return;
							FileView.this.min_value=values[0];
							FileView.this.max_value=values[1];
							getBrowser().offscreen=null;
							getBrowser().drawingArea.repaint();
							}});
						}
				
				/** paint that file view */
				private void paint(GC gc)
					{
					if(this.getMaxY()< gc.top_y) return;
					if(this.getTopY()>= (gc.top_y+(gc.height-topMargin()))) return;
					Graphics2D g=gc.g;
					Hershey hershey=new Hershey();
					Shape clip=g.getClip();
					//g.setClip(new Rectangle(leftMargin()/2,getTopY()-topMargin(),leftMargin()/2,getHeight()));
					g.setFont(new Font("Courier",Font.PLAIN,9));
					g.setColor(Color.BLACK);
	        		AffineTransform old=g.getTransform();
	        		AffineTransform tr= AffineTransform.getTranslateInstance(leftMargin()/2, getTopY()-gc.top_y+topMargin());
	        		tr.rotate(Math.PI/2);
	        		g.setTransform(tr);
	        		hershey.paint(g, this.tabixFile.getFile().getName(), 0,0,this.getHeight(), leftMargin()/2);
	        		//g.drawString(this.tabixFile.getFile().getName(),0,0);
	        		g.setTransform(old);
					g.setClip(clip);
					
					
					g.setColor(Color.WHITE);
					g.fillRect(
							leftMargin(),
							getTopY(),
							gc.width-leftMargin(),
							getHeight());
					
					/* paint gradient */
					for(int side=0;side<2;++side)
						{
						float y_0=value2pixel(0f);
						float y_value2=value2pixel(2.0*(side==0?1.0:-1.0));
						float fractions[]=new float[]{0.0f,1.0f};
						if(y_value2==y_0) continue;
						Color colors[]=new Color[]{Color.WHITE,new Color(255,210,210)};
						LinearGradientPaint grad=new LinearGradientPaint(0,y_0,0,y_value2,fractions,colors);
						Paint oldp= g.getPaint();
						g.setPaint(grad);
						Rectangle r=new Rectangle();
						r.x=leftMargin();
						r.width=gc.width-leftMargin();
						
						if(side==0)
							{
							r.y=getTopY();
							r.height=(int)(y_0-getTopY());
							if(r.height<=0) continue;
							}
						else
							{
							r.y=(int)y_0;
							r.height=(int)(getMaxY()-y_0);
							if(r.y>= getMaxY() || r.height<=0) continue;
							}
						
						g.fill(r);
						g.setPaint(oldp);
						}
						
					/* paint y axis */
					if(this.min_value<=0 && this.max_value>=0)
						{
						g.setColor(Color.RED);
						int y_axis=value2pixel(0);
						g.drawLine(
								leftMargin(),
								y_axis,
								gc.width,
								y_axis);
						}
					/* paint data */
					Iterator<SigData> iter=this.tabixFile.iterator(
							ChromView.this.getChromosome().getSequenceName(),
							(Browser.this.view_start+1),
							1+Math.min(getChromosome().getSequenceLength(),Browser.this.view_start+Browser.this.view_length)
							);
					/* loop over data */
					while(iter.hasNext())
						{
						SigData data=iter.next();
						if(data.getEnd()< Browser.this.view_start) continue;
						if(data.getStart()>(Browser.this.view_start+Browser.this.view_length)) continue;//break ?
						if(data.getValue()< this.min_value || data.getValue()> this.max_value) continue;
						g.setColor(data.getColor());
						int x0= getBrowser().base2pixel(data.getStart(), gc);
						int x1= getBrowser().base2pixel(data.getEnd(), gc);
						if(x1==x0) x1++;
						g.fillRect(
								x0,
								this.value2pixel(data.getValue()),
								(x1-x0),
								5);
						
						}
					//iter.close();
					/* draw frame */
					g.setColor(Color.GRAY);
					g.drawRect(
							leftMargin(),
							getTopY(),
							gc.width-leftMargin(),
							getHeight());
					}
				
				/** convert value to Y -pixel*/
				private int value2pixel(double value)
					{
					return this.getTopY()+
							(int)(((this.max_value-value)/(this.max_value-this.min_value))*this.getHeight());
					}
				}			

			/** constructor for this chromosome */
			ChromView(SAMSequenceRecord chrom )
				{
				this.chrom=chrom;
				}
			
			/** get associated browser */
			public Browser getBrowser()
				{
				return Browser.this;
				}
			
			/** get chromosome */
			public SAMSequenceRecord getChromosome()
				{
				return chrom;
				}
			
			/** return true wether this chromosome contains (x,y) */
			public boolean contains(GC gc,int x,int y)
				{
				if(y< topMargin()) return false;
				if(getTopY(gc)> gc.height) return false;
				return new Rectangle(0, getTopY(gc), gc.width, getHeight()).contains(x,y);
				}
			
			/** fills a popup for this chromosome */
			public void fillPopup(JPopupMenu popup)
				{
				popup.add(new JSeparator());
				popup.add(new AbstractAction("Adjust Chromosome Height....")
					{
					@Override
					public void actionPerformed(ActionEvent e)
						{
						Integer h=getBrowser().askHeight();
						if(h==null) return;
						for(FileView fv:fileViews)
							{
							fv.height_pixels=h;
							}
						getBrowser().recalcSize();
						}});
				
				}
			
			/** paint this chrom view */
			private void paint(GC gc)
				{
				if(this.getMaxY()< gc.top_y) return;
				if(this.getTopY()>= (gc.top_y+(gc.height-topMargin()))) return;
				
				Graphics2D g=gc.g;
				g.translate(0, -gc.top_y);
				g.setColor(background);
				g.fillRect(0, getTopY(), gc.width, this.getHeight());
				Hershey hershey=new Hershey();
				Shape clip=g.getClip();
				//g.setClip(new Rectangle(0,getTopY()-topMargin(),leftMargin()/2,getHeight()));
				//g.setFont(new Font("Courier",Font.BOLD,14));
				g.setColor(Color.BLACK);
        		AffineTransform old=g.getTransform();
        		AffineTransform tr= AffineTransform.getTranslateInstance(
        				leftMargin(),
        				getTopY()-gc.top_y+topMargin()
        				);
        		tr.rotate(Math.PI/2);
        		g.setTransform(tr);
        		//hershey.paint(g, getChromosome().getSequenceName(),0,0,this.getHeight()+100,100+leftMargin()/2);
        		//g.drawString(getChromosome().getSequenceName(),0,0);
        		hershey.paint(g, getChromosome().getSequenceName(),0,0,this.getHeight(),leftMargin()/2);
        		g.setTransform(old);
				g.setClip(clip);
				
				// paint all file views
				for(FileView fv:this.fileViews)
					{
					fv.paint(gc);
					}
				g.translate(0, gc.top_y);
				}
			}
		

		/** browser constuctor dictinary and set of tabix Files */
		Browser(SAMSequenceDictionary dict,SigFiles sigilesList)
			{
			super();
			this.samSequenceDictionary=dict;
			this.sigilesList=sigilesList;
			this.setTitle("Browser :"+this.sigilesList.size()+" file(s)");
			this.setMinimumSize(new Dimension(leftMargin()+10,topMargin()+10));
			JToolBar toolbar=new JToolBar();
			this.getContentPane().add(toolbar,BorderLayout.NORTH);
			JPanel pane=new JPanel(new BorderLayout());
			this.getContentPane().add(pane,BorderLayout.CENTER);
			
			this.addInternalFrameListener(new InternalFrameAdapter()
				{
				@Override
				public void internalFrameClosed(InternalFrameEvent e)
					{
					for(SigDataFileTabixReader f: Browser.this.sigilesList) f.close();
					Browser.this.doMenuClose();
					}
				});
			
			this.addComponentListener(new ComponentAdapter()
				{
				@Override
				public void componentResized(ComponentEvent e)
					{
					adjustScrollBars();
					drawingArea.repaint();
					}
				});
			
			this.vScrollBar=new JScrollBar(JScrollBar.VERTICAL);
			pane.add(this.vScrollBar,BorderLayout.EAST);
			this.hScrollBar=new JScrollBar(JScrollBar.HORIZONTAL);
			pane.add(this.hScrollBar,BorderLayout.SOUTH);
			
			/* create drawing area */
			this.drawingArea=new JPanel(null,true)
				{
				/* create a tooltip text */
				@Override
				public String getToolTipText(MouseEvent event)
					{
					/* create genome-position */
					GC gc=new GC();
					StringBuilder b=new StringBuilder();
					int b1= pixel2base(event.getX(), gc);
					if(b1!=-1)
						{
						b.append("position:").append(niceNumber(b1));
						int b2= pixel2base(event.getX()+1, gc);
						if(b2!=-1 && b2!=b1)
							{
							b.append("-").append(niceNumber(b2));
							}
						}
					/* find chromview at this location */
					ChromView cv=findChromView(event.getX(), event.getY());
					if(cv!=null)
						{
						b.insert(0,cv.getChromosome().getSequenceName());
						ChromView.FileView fv=findFileView(event.getX(), event.getY());
						if( fv!=null)
							{
							b.insert(0,fv.tabixFile.getURI()+" ");
							}
						}
					return b.length()==0?null:b.toString();
					}
				/* paint component */
				@Override
				protected void paintComponent(Graphics g1)
					{
					if(offscreen==null ||
							offscreen.getWidth()!=drawingArea.getWidth() ||
							offscreen.getHeight()!=drawingArea.getHeight()
							)
							{
							offscreen=new BufferedImage(
								drawingArea.getWidth(),
								drawingArea.getHeight(),
								BufferedImage.TYPE_INT_RGB
								);
							GC gc= new GC();
							gc.top_y=vScrollBar.getValue();
							gc.width=drawingArea.getWidth();
							gc.height=drawingArea.getHeight();
							gc.g=offscreen.createGraphics();
							paintDrawingArea(gc);
							gc.g.dispose();
							}
					Graphics2D.class.cast(g1).drawImage(offscreen, 0, 0, drawingArea);
					}
				};
			this.drawingArea.setToolTipText("");
			this.drawingArea.setOpaque(true);
			this.drawingArea.setBackground(Color.WHITE);
			pane.add(this.drawingArea,BorderLayout.CENTER);
			
			/* create chrom-views */
			for(SAMSequenceRecord chrom: this.samSequenceDictionary.getSequences())
				{
				ChromView cv=new ChromView(chrom);
				for(SigDataFileTabixReader sf: sigilesList)
					{
					if(!sf.hasDataForChrom(chrom)) continue;
					this.max_chrom_length = Math.max(this.max_chrom_length,chrom.getSequenceLength()+1);
				
					cv.fileViews.add(cv.new FileView(sf));
					}
				cv.background=(chromViews.size()%2==0?Color.WHITE:Color.LIGHT_GRAY);
				chromViews.add(cv);
				}
			this.view_length=this.max_chrom_length;
			this.recalcSize();
			
			AdjustmentListener al=new AdjustmentListener()
				{
				@Override
				public void adjustmentValueChanged(AdjustmentEvent e)
					{
					if(e.getAdjustable()==hScrollBar)
						{
						view_start=hScrollBar.getValue();
						}
					if(e.getValueIsAdjusting())
						{
						//handle?
						}
					offscreen=null;
					drawingArea.repaint();
					}
				};
			vScrollBar.addAdjustmentListener(al);
			hScrollBar.addAdjustmentListener(al);
			drawingArea.addMouseListener(this);
			drawingArea.addMouseMotionListener(this);
			
			this.history.add(new History(this.view_start,this.view_length));
			this.position_in_history=0;
			
			/* create history buttons */
			AbstractAction action=new AbstractAction("Prev")
				{	
				@Override
				public void actionPerformed(ActionEvent e)
					{
					goHistory(position_in_history-1);
					}
				};
			action.setEnabled(false);
			this.getActionMap().put("ACTION.PREV.HISTORY", action);
			toolbar.add(new JButton(action));
				
			action=new AbstractAction("Next")
				{	
				@Override
				public void actionPerformed(ActionEvent e)
					{
					goHistory(1);
					}
				};
			action.setEnabled(false);
			this.getActionMap().put("ACTION.NEXT.HISTORY", action);
			toolbar.add(new JButton(action));
			
			/* create tools buttons */
			ButtonGroup butGroup=new ButtonGroup();
			for(Tool tool: Tool.values())
				{
				JToggleButton toggle=new JToggleButton(new SelectToolAction(tool));
				butGroup.add(toggle);
				toolbar.add(toggle);
				}
			JMenuBar bar=new JMenuBar();
			this.setJMenuBar(bar);
			JMenu menu=new JMenu("Options");
			bar.add(menu);
			menu.add(new AbstractAction("Set Height")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					Integer h=Browser.this.askHeight();
					if(h==null) return;
					for(ChromView cv:chromViews)
						{
						for(ChromView.FileView fv:cv.fileViews)
							{
							fv.height_pixels=h;
							}
						}
					recalcSize();
					}
				});
			menu.add(new AbstractAction("Set Min/Max")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					double h[]=Browser.this.askMinValue(-100, 100);
					if(h==null) return;
					for(ChromView cv:chromViews)
						{
						for(ChromView.FileView fv:cv.fileViews)
							{
							fv.min_value=h[0];
							fv.max_value=h[1];
							}
						}
					offscreen=null;
					drawingArea.repaint();
					}
				});
			bar.add(new JSeparator());
			bar.add(menu);
			menu.add(new AbstractAction("Close")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					Browser.this.doMenuClose();
					}
				});
			}
		
		void doMenuClose()
			{
			//this.setVisible(false);
			//this.dispose();
			}
		
		int leftMargin()
			{
			return LEFT_DISTANCE;
			}
		
		int topMargin()
			{
			return TOP_DISTANCE;
			}
		
		private void goHistory(int new_position_history)
			{
			if(new_position_history<0 || new_position_history>=history.size() || history.isEmpty()) return;
			position_in_history=new_position_history;
			Browser.this.view_start=history.get(position_in_history).view_start;
			Browser.this.view_length=history.get(position_in_history).view_length;
			offscreen=null;
			drawingArea.repaint();
			this.getActionMap().get("ACTION.PREV.HISTORY").setEnabled(position_in_history>0);
			this.getActionMap().get("ACTION.NEXT.HISTORY").setEnabled(position_in_history+1< history.size());
			}	
		
		
		private void paintDrawingArea(GC gc)
			{
			Graphics2D g=gc.g;
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, gc.width, gc.height);
			for(gc.step=0;gc.step<1;++gc.step)
				{
				paintAxis(gc);
				Shape oldClip=g.getClip();
				g.setClip(new Rectangle(0,topMargin(),gc.width,gc.height-topMargin()));
				g.translate(0, topMargin());
				for(ChromView cv: this.chromViews)
					{
					cv.paint(gc);
					}
				g.translate(0, -topMargin());
				g.setClip(oldClip);
				}
			}
	
		private void recalcSize()
			{
			int y=0;
			for(ChromView cv:this.chromViews)
				{
				cv.y_pixels=y;
				for(ChromView.FileView fv:cv.fileViews)
					{
					y+=BORDER_SIZE;
					fv.y_pixels=y;
					y+= fv.height_pixels;
					}
				y+=BORDER_SIZE;
				cv.height_pixels=(y-cv.y_pixels);
				}
			adjustScrollBars();
			}
		
		
		private void adjustScrollBars()
			{
			//adjust horizontal
			this.hScrollBar.setValues(
				this.view_start, 
				this.view_length,
				0,
				this.max_chrom_length-this.view_length
				);
			ChromView last=this.chromViews.get(this.chromViews.size()-1);
			int top_y= this.vScrollBar.getValue();
			int max_y= last.getMaxY();
			int extend= drawingArea.getHeight()-topMargin();
			if(extend<=0) extend=1;
			if(max_y<extend) max_y=extend;
			
			//adjust vertical
			this.vScrollBar.setValues(
				top_y,
				extend,
				0,
				max_y
				);
			}
		
		private void paintAxis(GC gc)
			{
			AffineTransform old;
			Graphics2D g=gc.g;
		    final int num_steps=10;
		    int step=(int)(Math.pow(10,Math.ceil(Math.log10(this.view_length)))/num_steps);
		    if(step>=10 && step*2>=this.view_length)
		            {
		            step/=10;
		            }
		    
		    if(step>0)
		            {
		            //create a scale with vertical bars
		            int ticks = view_start - view_start%step;
		            
		            while(ticks<= view_start+view_length)
		                    {
		                    float x= (float)base2pixel(ticks,gc);
		                    if(x>=leftMargin() && x< gc.width)
		                            {
		                    		g.setColor(Color.BLACK);
		                    		g.fillRect((int)x, 0, 3, gc.height);
		                    		g.setColor(Color.GRAY); 
		                    		g.fillRect((int)x, 0, 2, gc.height);
		                           
		                    		g.setColor(Color.BLACK);
		                    		old=g.getTransform();
		                    		AffineTransform tr= AffineTransform.getTranslateInstance(x+2, 10);
		                    		tr.rotate(Math.PI/2);
		                    		g.setTransform(tr);
		                    		g.drawString(niceNumber(ticks),0,0);
		                    		g.setTransform(old);
		                           
		                            }
		                    
		                    if(step>=10)
		                            {
		                            for(int t2=ticks;t2<ticks+step;t2+=step/10)
		                                    {
		                                    float x2= base2pixel(t2,gc);
		                                    if(!(x2>=leftMargin() && x2< gc.width)) continue;
		                                    g.setColor(Color.LIGHT_GRAY);
		                                    g.drawLine((int)x2, 0, (int)x2,  gc.height);
		                                    
		                                    /*
		                                    g.setColor(Color.BLACK);
				                    		old=g.getTransform();
				                    		AffineTransform tr= AffineTransform.getTranslateInstance(x2+2, 10);
				                    		tr.rotate(Math.PI/2);
				                    		g.setTransform(tr);
				                    		g.drawString(niceNumber(t2),0,0);
				                    		g.setTransform(old);
		                                    */
		                                    }
		                            }
		                    
		                    ticks+=step;
		                    }
		            	}
		            }
			
		  private String niceNumber(final int i)
	          {
	          String s=String.valueOf(i);
	          StringBuilder b= new StringBuilder(s.length());
	          for(int j=0;j< s.length();j++)
	                  {
	                  if(j!=0 && j%3==0) b.insert(0, ' ');
	                  b.insert(0,s.charAt(s.length()-1-j));
	                  }
	          return b.toString();
	          }

		
		
		    private int base2pixel(final int genome_pos,GC gc)
		    	{
		    	int left=leftMargin();
		    	return left+(int)(((genome_pos- this.view_start)/(double)this.view_length)*(gc.width-left));
		    	}
		    
		    private int pixel2base(final int x1,GC gc)
		    	{
		    	int left=leftMargin();
		    	if(x1< left) return -1;
		    	return this.view_start+(int)(((x1-left)/(double)(gc.width-left))*this.view_length);
		    	}
		    
			
			@Override
			public void mouseClicked(MouseEvent e)
				{
				
				}
	
			@Override
			public void mouseEntered(MouseEvent e)
				{
				
				}
	
			@Override
			public void mouseExited(MouseEvent e)
				{
				}
	
			@Override
			public void mousePressed(MouseEvent e)
				{
				this.mouseStart=null;
				this.mousePrev=null;
				this.selectedChromView=findChromView(e.getX(),e.getY());
				this.selectedFileFiew=findFileView(e.getX(),e.getY());
				if(e.isPopupTrigger() || e.isControlDown())
					{
					JPopupMenu menu=new JPopupMenu();
					if(this.selectedChromView!=null)
						{
						this.selectedChromView.fillPopup(menu);
						if(this.selectedFileFiew!=null)
							{
							this.selectedFileFiew.fillPopup(menu);
							}
						}
					menu.show(drawingArea, e.getX(), e.getY());
					return;
					}
				
				this.mouseStart=new Point(e.getX(),e.getY());
				this.mousePrev=null;
				switch(currentTool)
					{
					case CHANGE_HEIGHT:
					
					case SHOW_DATA_FOR_FILE:
						if(selectedFileFiew==null)
							{
							mouseStart=null;
							}
						break;
					case ZOOM: break;
					default:break;
					}
				}
			
			@Override
			public void mouseDragged(MouseEvent e)
				{
				if(mouseStart!=null && e.getX()> (leftMargin()) && e.getX()<this.drawingArea.getWidth())
					{
					int topY= topMargin();
					int heightY= drawingArea.getHeight()-topMargin();
					switch(currentTool)
						{
						case CHANGE_HEIGHT:
						case SHOW_DATA_FOR_FILE:
							{
							GC gc=new GC();
							topY=selectedFileFiew.getTopY(gc);
							heightY= selectedFileFiew.getHeight();
							break;
							}
						default:break;
						}
					Graphics2D g=(Graphics2D)drawingArea.getGraphics();
					g.setXORMode(drawingArea.getBackground());
					if(currentTool==Tool.CHANGE_HEIGHT)
						{
						int mid=topY+heightY/2;
						int radius;
						if(mousePrev!=null)
							{
							radius= (int) Point2D.distance(mouseStart.x, mid, mousePrev.x, mousePrev.y);
							g.drawOval(mouseStart.x-radius, mid-radius, radius*2, radius*2);
							}
						mousePrev=new Point(e.getX(),e.getY());
						radius= (int) Point2D.distance(mouseStart.x, mid, mousePrev.x, mousePrev.y);
						g.drawOval(mouseStart.x-radius, mid-radius, radius*2, radius*2);
						}
					else
						{
						
						if(mousePrev!=null)
							{
							g.fillRect(
								Math.min(mousePrev.x, mouseStart.x),
								topY,
								Math.abs(mousePrev.x - mouseStart.x),
								heightY
								);
							}
						mousePrev=new Point(e.getX(),e.getY());
						g.fillRect(
								Math.min(mousePrev.x, mouseStart.x),
								topY,
								Math.abs(mousePrev.x - mouseStart.x),
								heightY
								);
						
						}
					g.setPaintMode();
					}
				}
			
			@Override
			public void mouseReleased(MouseEvent e)
				{
				if(mouseStart!=null &&
					mousePrev!=null &&
					!mouseStart.equals(mousePrev))
					{
					switch(currentTool)
						{
						case ZOOM:
							{
							GC gc=new GC();
							int chromStart=pixel2base(Math.max(leftMargin(),Math.min(mouseStart.x, mousePrev.x)), gc);
							int chromEnd=pixel2base(Math.min(drawingArea.getWidth(),Math.max(mouseStart.x, mousePrev.x)), gc);
							if(chromStart< chromEnd)
								{
								History h=new History(chromStart, chromEnd-chromStart);
								
								history.add(h);
								position_in_history++;
								while(position_in_history+1>history.size() )
									{
									history.remove(history.size()-1);
									}
								this.view_start=h.view_start;
								this.view_length=h.view_length;
								this.getActionMap().get("ACTION.PREV.HISTORY").setEnabled(position_in_history>0);
								this.getActionMap().get("ACTION.NEXT.HISTORY").setEnabled(false);
								adjustScrollBars();
								}
							break;
							}
						case CHANGE_HEIGHT:
							{
							GC gc=new GC();
							int mid=selectedFileFiew.getTopY(gc)+selectedFileFiew.getHeight()/2;
							int radius= (int) Point2D.distance(mouseStart.x, mid, mousePrev.x, mousePrev.y);
							radius*=2;
							if(radius< 20 || radius>1000) break;
							selectedFileFiew.height_pixels=radius;
							recalcSize();
							break;
							}
						case SHOW_DATA_FOR_FILE:
							{
							GC gc=new GC();
							int chromStart=pixel2base(Math.max(leftMargin(),Math.min(mouseStart.x, mousePrev.x)), gc);
							int chromEnd=pixel2base(Math.min(drawingArea.getWidth(),Math.max(mouseStart.x, mousePrev.x)), gc);
							Iterator<SigData> iter=selectedFileFiew.tabixFile.iterator(
									selectedFileFiew.getChromosome().getSequenceName(),
									chromStart,
									chromEnd
									);
							Table table=new Table(iter);
							CloserUtil.close(iter);
							SigFrame.this.desktop.add(table);
							table.setVisible(true);
							break;
							}
						}
					}
				mouseStart=null;
				mousePrev=null;
				offscreen=null;
				selectedChromView=null;
				selectedFileFiew=null;
				drawingArea.repaint();
				}
	
	
			@Override
			public void mouseMoved(MouseEvent e)
				{
				
				}
			
			/** find chrom view at given location */
			private ChromView findChromView(int x,int y)
				{
				GC gc=new GC();
				for(ChromView cv:this.chromViews)
					{
					if(cv.contains(gc,x,y)) return cv;
					}
				return null;
				}
			
			/** find file-view at given location */
			private ChromView.FileView findFileView(int x,int y)
				{
				GC gc=new GC();
				for(ChromView cv:this.chromViews)
					{
					for( ChromView.FileView fv:cv.fileViews)
						{
						if(fv.contains(gc,x,y)) return fv;
						}
					}
				return null;
				}
			
			private Integer askHeight()
				{
				JSpinner spin=new JSpinner(new SpinnerNumberModel(DEFAULT_FILE_VIEW_HEIGHT, 50, 1000, 1));
				JPanel pane=new JPanel(new FlowLayout(FlowLayout.LEADING));
				JLabel label=new JLabel("new Height:",JLabel.RIGHT);
				label.setLabelFor(spin);
				pane.add(label);
				pane.add(spin);
				if(JOptionPane.showConfirmDialog(
						Browser.this, pane,"Height",
						JOptionPane.OK_CANCEL_OPTION
						)!=JOptionPane.OK_OPTION)
					{
					return null;
					}
				
				return SpinnerNumberModel.class.cast(spin.getModel()).getNumber().intValue();
				}
			private double[] askMinValue(double min,double max)
				{
				JSpinner spinMin=new JSpinner(new SpinnerNumberModel(min,Math.min(min, -20.0),Math.max(max, 20.0),0.001));
				JPanel pane=new JPanel(new FlowLayout(FlowLayout.LEADING));
				JLabel label=new JLabel("new Min Value:",JLabel.RIGHT);
				label.setLabelFor(spinMin);
				pane.add(label);
				pane.add(spinMin);
				JSpinner spinMax=new JSpinner(new SpinnerNumberModel(max,Math.min(min, -20.0),Math.max(max, 20.0),0.001));
				label=new JLabel("new Max Value:",JLabel.RIGHT);
				label.setLabelFor(spinMax);
				pane.add(label);
				pane.add(spinMax);
				do
					{
					if(JOptionPane.showConfirmDialog(
							Browser.this, pane,"Height",
							JOptionPane.OK_CANCEL_OPTION
							)!=JOptionPane.OK_OPTION)
						{
						return null;
						}
					
				    min= SpinnerNumberModel.class.cast(spinMin.getModel()).getNumber().doubleValue();
				    max= SpinnerNumberModel.class.cast(spinMax.getModel()).getNumber().doubleValue();
					} while(min>max);
				return new double[]{min,max};
				}
			}
	
	/**
	 * Table
	 */
	private class Table extends IFrame
		{
		/**
		 * Model for this table
		 *
		 */
		class SigTableModel
			extends AbstractGenericTable<SigData>
			{
			public SigTableModel(Iterator<SigData> iter)
				{
				while(iter.hasNext() && this.getRowCount()<10000)
					{
					this.getRows().add(iter.next());
					}
				}
			
			@Override
			public int getColumnCount()
				{
				return 6;
				}
		
			
			@Override
			public Class<?> getColumnClass(int columnIndex) {
				switch(columnIndex)
					{
					case 0: return String.class;
					case 1: return Integer.class;
					case 2: return Integer.class;
					case 3: return String.class;
					case 4: return Double.class;
					case 5: return String.class;
					}
				return null;
				}
			@Override
			public Object getValueOf(SigData d, int columnIndex)
				{
				switch(columnIndex)
					{
					case 0: return d.getChromosome();
					case 1: return d.getStart();
					case 2: return d.getEnd();
					case 3: return d.getName();
					case 4: return d.getValue();
					case 5: return d.color_string;
					}
				return null;
				}
			
			@Override
			public String getColumnName(int columnIndex)
				{
				switch(columnIndex)
					{
					case 0: return "Chrom";
					case 1: return "Start";
					case 2: return "End";
					case 3: return "Name";
					case 4: return "Value";
					case 5: return "Color";
					}
				return null;
				}
			
			
			}
		private JTable table;
		Table(Iterator<SigData> iter)
			{
			SigTableModel model=new SigTableModel(iter);
			this.table=new JTable(model);
			table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			getContentPane().add(new JScrollPane(this.table),BorderLayout.CENTER);
			if(iter.hasNext())
				{
				getContentPane().add(new JLabel("Warning: maximum number of rows reached"),
						BorderLayout.NORTH);
				}
			}
		}
	
	
	private SigFrame()
		{
		super("SigFrame");
		setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
		addWindowListener(new WindowAdapter()
			{
			@Override
			public void windowClosing(WindowEvent e)
				{
				doMenuQuit();
				}
			});
		
		JPanel contentPane=new JPanel(new BorderLayout());
		contentPane.setBorder(new EmptyBorder(4,4, 4, 4));
		setContentPane(contentPane);
		this.desktop=new JDesktopPane();
		contentPane.add(desktop);
		
		
		
		
		AbstractAction action=new AbstractAction("Quit")
			{
			@Override
			public void actionPerformed(ActionEvent e)
				{
				doMenuQuit();
				}
			};
		this.actionMap.put("ACTION_QUIT", action);
		
		
		action=new AbstractAction("Open Reference...")
			{
			@Override
			public void actionPerformed(ActionEvent e)
				{
				doMenuOpenGenome();
				}
			};
		action.putValue(AbstractAction.SHORT_DESCRIPTION, "Load reference");
		this.actionMap.put("ACTION_OPEN_FAIDX", action);
		
		action=new AbstractAction("Open File...")
			{
			@Override
			public void actionPerformed(ActionEvent e)
				{
				doMenuOpen();
				}
			};
		action.putValue(AbstractAction.SHORT_DESCRIPTION, "Load Files....");
		action.setEnabled(false);
		this.actionMap.put("ACTION_OPEN", action);
		
		
		/*
		action=new AbstractAction("Open URL...")
			{
			@Override
			public void actionPerformed(ActionEvent e)
				{
				doMenuLoadURL();
				}
			};
		action.putValue(AbstractAction.SHORT_DESCRIPTION, "Load from remote URL");
		this.actionMap.put("ACTION_LOADURL", action);
		*/
		
		
		JMenuBar bar=new JMenuBar();
		setJMenuBar(bar);
		JMenu menu=new JMenu("File");
		bar.add(menu);
		menu.add(new AbstractAction("About...")
			{
			@Override
			public void actionPerformed(ActionEvent ae)
				{
				JOptionPane.showMessageDialog(
				 SigFrame.this,
                  "<html><body><h1 align='center'>"+
                  "SigFrame"+
                  "</h1>"+aboutMessage+"</body></html>",
                  "About...",JOptionPane.PLAIN_MESSAGE);
				}
			});
		  

		menu.add(new JSeparator());
		menu.add(this.actionMap.get("ACTION_OPEN_FAIDX"));
		menu.add(this.actionMap.get("ACTION_OPEN"));
		menu.add(new JSeparator());
		menu.add(this.actionMap.get("ACTION_QUIT"));
		}
	
	private void doMenuQuit()
		{
		for(JInternalFrame ji:this.desktop.getAllFrames())
			{
			ji.setVisible(false);
			ji.dispose();
			}
		this.setVisible(false);
		this.dispose();
		}
	
	public SAMSequenceDictionary getSamSequenceDictionary()	
		{
		return genome;
		}
	private void doMenuOpenGenome()
		{
		JFileChooser chooser=new JFileChooser();
		chooser.setFileFilter(new FileFilter()
			{
			@Override
			public String getDescription()
				{
				return "faidx file";
				}
			
			@Override
			public boolean accept(File f)
				{
				return f.isDirectory() || (f.isFile() && f.getName().endsWith(".dict") || f.getName().endsWith(".fai"));
				}
			});
		if(chooser.showOpenDialog(this)!=JFileChooser.APPROVE_OPTION) return;
		try {
			doMenuOpenGenome(chooser.getSelectedFile());
			} 
		catch (Exception e)
			{
			this.actionMap.get("ACTION_OPEN").setEnabled(false);
			showError(e);
			}
		}
	private void doMenuOpenGenome(File f) throws IOException
		{	
		LOG.info("open faidx "+f);
		this.actionMap.get("ACTION_OPEN").setEnabled(false);
		this.genome= SAMSequenceDictionaryExtractor.extractDictionary(f);
		this.actionMap.get("ACTION_OPEN").setEnabled(true);
		}
	
	
	private void doMenuOpen()
		{
		if(this.genome==null) return ;
		JFileChooser chooser=new JFileChooser();
		chooser.setMultiSelectionEnabled(true);
		chooser.setFileFilter(new FileFilter()
			{
			@Override
			public String getDescription()
				{
				return "VCF file indexed with tabix";
				}
			
			@Override
			public boolean accept(File f)
				{
				return f.isDirectory() || TabixFileReader.isValidTabixFile(f);
				}
			});
		if(chooser.showOpenDialog(this)!=JFileChooser.APPROVE_OPTION) return;
		File files[]=chooser.getSelectedFiles();
		if(files==null || files.length==0) return;
		
		
		openTabixFiles(Arrays.asList(files));
			
		}
	
	
	private void openTabixFiles(Collection<File> inputs)
		{
		SigFiles storage=new SigFiles();
		try
			{
			for(File f:inputs) 
				{
				LOG.info("opening "+f);
				storage.add(new SigDataFileTabixReader(f));
				}
			Browser browser=new Browser(this.genome,storage);
			SigFrame.this.desktop.add(browser);
			browser.setVisible(true);
			}
		catch(Exception err)
			{
			err.printStackTrace();
			showError(err);
			return;
			}
		}
	
	private void showError(Object o)
		{
		
		}
	
	@Program(name="sigframe",
			description="SigFrame displays CGH/ position+values in a GUI",
			keywords={"cgh","gui","visualization"}
			)
	public static class Main extends Launcher
		{
		@Parameter(names="-R",description="Reference indexed fasta file",required=true)
		private File referenceFile;
		private SigFrame app=new SigFrame();
		private Main()
			{
			}
		@Override
		public String getProgramName()
			{
			return "SigFrame";
			}
		

		@Override
		public int doWork(List<String> args) {
			if(this.referenceFile==null)
			{
				LOG.error("ref missing");
				return -1;
			}
			try
				{
				app.doMenuOpenGenome(this.referenceFile);
				app.aboutMessage="Author: Pierre Lindenbaum "+
						"<br>Version:"+getVersion();
				
				JFrame.setDefaultLookAndFeelDecorated(true);
				JDialog.setDefaultLookAndFeelDecorated(true);
				
			   /** loop over tabix file */
	           if(!args.isEmpty())
	                {
	        	   	if(app.getSamSequenceDictionary()==null)
	        	   		{
	        	   		LOG.error("No indexed genome defined");
	        	   		return -1;
	        	   		}
	            	LoadAtStartup startup= new LoadAtStartup(app);
	                for(final String fname:args)
	                        {
	                		final File f=new File(fname);
	                		if(!TabixFileReader.isValidTabixFile(f))
	                			{
	                			LOG.error("Not a valid tabix file:"+f);
	                			return -1;
	                			}
	                		startup.sources.add(f);
	                        }
	                app.addWindowListener(startup);
	                }
	           
	        	Dimension screen= Toolkit.getDefaultToolkit().getScreenSize();
	        	app.setBounds(
	        		50,50,
	        		screen.width-100,
	        		screen.height-100
	        		);
				SwingUtilities.invokeAndWait(new Runnable()
					{	
					@Override
					public void run()
						{
						app.setVisible(true);
						}
					});
				return 0;
				}
			catch (Exception e)
				{
				e.printStackTrace();
				return -1;
				}
			}
		private void startInstance(final String[] args)
			{
			instanceMain(args);
			}
		}
	public static void main(String[] args)
		{
		new Main().startInstance(args);
		}
	}