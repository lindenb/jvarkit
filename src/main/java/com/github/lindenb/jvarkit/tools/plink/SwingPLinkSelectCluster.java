
/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.plink;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.Vector;
import java.util.function.Function;
import java.util.regex.Pattern;

import javax.swing.AbstractAction;
import javax.swing.AbstractButton;
import javax.swing.ButtonGroup;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JToggleButton;
import javax.swing.SwingUtilities;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;

import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.swing.PreferredDirectory;
import com.github.lindenb.jvarkit.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

GUI selecting the samples of a MDS file generated with plink.

## Example:

END_DOC
 */
@Program(name="swingplinkselectcluster",
description="Swing-based Plink/MDS sample selector",
keywords={"plink","sample","swing"},
creationDate = "20231123",
modificationDate="20231123",
jvarkit_amalgamion =  true
)
public class SwingPLinkSelectCluster extends Launcher {
private static final Logger LOG = Logger.build(SwingPLinkSelectCluster.class).make();
private static final String ACTION_XY_KEY= "plink.xy";
private static final String ACTION_TOOL= "select.tool";
private  enum ToolType {AND,NOT,XOR}

private static class Column {
	final int index;
	final String label;
	double min = 0.0;
	double max = 0.0;
	Column(int index,String label) {
		this.index  = index;
		this.label = label;
		}
	}


private static class ColXY {
	final Column columnX;
	final Column columnY;
	ColXY(final Column columnX,final Column columnY) {
		this.columnX = columnX;
		this.columnY = columnY;
		}
	Column get(int i) {
		return i==0?columnX:columnY;
		}
	String getName() {
		return columnX.label+" / "+columnY.label;
	}
}



private static class Sample {
	String fid;
	String iid;
	double[] C;
	//
	boolean selectedInTable = false;
	boolean selectForOutput = true;
	}

@SuppressWarnings("serial")
private static class SampleTableModel extends AbstractTableModel{
	final List<Sample> rows;
	final int countCcols;
	 SampleTableModel(final List<Sample> rows,final int countCcols) {
		this.rows= rows;
		this.countCcols = countCcols;
	 	}
	 
	@Override
	public Class<?> getColumnClass(int columnIndex) {
		 switch(columnIndex) {
			case 0:
			case 1: return String.class;
			default: return Double.class;
			}
	 	}
	@Override
	public int getRowCount() {
		return this.rows.size();
		}
	@Override
	public int getColumnCount() {
		return this.countCcols+2;
		}
	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		return false;
		}
	@Override
	public String getColumnName(int columnIndex) {
		switch(columnIndex) {
			case 0: return "FID";
			case 1: return "IID";
			default: return "C"+((columnIndex-2)+1);
			}
		}
	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		final Sample o = this.rows.get(rowIndex);
		switch(columnIndex) {
			case 0: return o.fid;
			case 1: return o.iid;
			default: return o.C[columnIndex-2];
			}
		}
	}

@SuppressWarnings("serial")
private static class XFrame extends JFrame {
	final List<Sample> samples;
	final ButtonGroup buttonGroupXY;
	final ButtonGroup buttonGroupTool;
	final JPanel drawingArea ;
	final JTable XYTable;
	private boolean dirtyFlag=false;
	private File saveAsFile = null;
	XFrame(final List<Column> columns,
			final List<Sample> samples
			)
		{
		super(SwingPLinkSelectCluster.class.getSimpleName());
		
		
		final List<ColXY> columnsXY = new Vector<>();
		for(int x=0;x< columns.size();++x) {
			for(int y=x+1;y< columns.size();++y) {
				columnsXY.add( new ColXY(columns.get(x),columns.get(y)));
			}
		}
		this.samples= samples;
		super.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		this.addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				if(dirtyFlag) {
					saveAs();
					}
				XFrame.this.setVisible(false);
				XFrame.this.dispose();
				}
			});
		
		
		final JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		final JMenu menu=new JMenu("File");
		menuBar.add(menu);
		menu.add(new JMenuItem(new AbstractAction("Save...") {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				saveAs(saveAsFile);
				}
			}));
		menu.add(new JMenuItem(new AbstractAction("Save As") {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				saveAs();
				}
			}));
		menu.add(new JSeparator());
		menu.add(new JMenuItem(new AbstractAction("Quit") {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				if(dirtyFlag) {
					saveAs();
					}
				XFrame.this.setVisible(false);
				XFrame.this.dispose();
				}
			}));
		final JPanel mainPanel=new JPanel(new BorderLayout());
		this.setContentPane(mainPanel);
		
		final JPanel topPane = new JPanel(new FlowLayout()); 
		mainPanel.add(topPane,BorderLayout.NORTH);
		
		this.buttonGroupXY = new ButtonGroup();
		this.buttonGroupTool = new ButtonGroup();
		
		this.drawingArea = new JPanel(null) {
			@Override
			public String getToolTipText(MouseEvent event) {
				final ColXY columnXY = getCurrentColXY();
				if(columnXY==null) return "";
				StringBuilder sb = null;
				for(int i=0;i< samples.size();i++) {
					final Sample sn = samples.get(i);
					final Point2D pt = sampleToPixel(columnXY,sn);
					
					if(pt.distance(event.getX(), event.getY()) < 3) {
						if(sb==null) sb=new StringBuilder();
						else if(sb.length()>0) sb.append(" ");
						sb.append(sn.iid);
						}
					}
				return sb==null?null:sb.toString();
				}
			@Override
			protected void paintComponent(Graphics g) {
				paintDrawingArea(Graphics2D.class.cast(g));
				}
			};
		drawingArea.setToolTipText("");
		drawingArea.setPreferredSize(new Dimension(100,100));
		drawingArea.setOpaque(true);
		final MouseAdapter mouse = new MouseAdapter() {
			Point prev = null;
			GeneralPath gpath;
			@Override
			public void mousePressed(MouseEvent e) {
				prev = new Point(e.getX(),e.getY());
				gpath =new GeneralPath();
				gpath.moveTo(prev.getX(), prev.getY());
				}
			@Override
			public void mouseDragged(MouseEvent e) {
				Graphics2D g = (Graphics2D)drawingArea.getGraphics();
				g.setColor(Color.MAGENTA);
				g.draw(new Line2D.Double(prev.getX(),prev.getY(),e.getX(),e.getY()));
				prev = new Point(e.getX(),e.getY());
				gpath.lineTo(e.getX(),e.getY());
				}
			@Override
			public void mouseReleased(MouseEvent e) {
				gpath.closePath();
				final ToolType tt = getCurrentTooltype();
				final ColXY cxy = getCurrentColXY();
				for(int i=0;cxy!=null && i< samples.size();i++) {
					final Sample sn = samples.get(i);
					final Point2D pt = sampleToPixel(cxy, sn);
					if(gpath.contains(pt.getX(),pt.getY())) {
						switch(tt) {
							case AND: sn.selectForOutput=true; break;
							case NOT: sn.selectForOutput=false; break;
							case XOR: sn.selectForOutput=!sn.selectForOutput; break;
							}
						}
					}
				drawingArea.repaint();
				dirtyFlag=true;
				gpath= null;
				prev=null;
				}
			};
		this.drawingArea.addMouseListener(mouse);
		this.drawingArea.addMouseMotionListener(mouse);
		
		
		final JPanel paneTool = new JPanel(new FlowLayout());
		paneTool.add(new JLabel("Tool:"));
		topPane.add(paneTool);
		
		for(ToolType toolType: ToolType.values() ) {
			final AbstractAction action = new AbstractAction(toolType.name()) {
					@Override
					public void actionPerformed(ActionEvent e) {
					}
				};
			action.putValue(ACTION_TOOL, toolType);
			final JToggleButton button=new  JToggleButton(action);
			
			this.buttonGroupTool.add(button);
			this.buttonGroupTool.setSelected(button.getModel(), true);
			paneTool.add(button);
			}
		
		final JPanel paneFlowXY = new JPanel(new FlowLayout());
		paneFlowXY.add(new JLabel("XY:"));
		topPane.add(paneFlowXY);
		for(ColXY cxy: columnsXY ) {
			final AbstractAction action = new AbstractAction(cxy.getName()) {
					@Override
					public void actionPerformed(ActionEvent e) {
						drawingArea.repaint();
					}
				};
			action.putValue(ACTION_XY_KEY, cxy);
			final JToggleButton button=new  JToggleButton(action);
			this.buttonGroupXY.add(button);
			this.buttonGroupXY.setSelected(button.getModel(), true);
			paneFlowXY.add(button);
			}
		
		
		
		
		
			JPanel pane2 = new JPanel(new BorderLayout());
			
			this.XYTable = new JTable(new SampleTableModel(samples,columns.size()));
			this.XYTable.getSelectionModel().addListSelectionListener(E->{
				if(E.getValueIsAdjusting()) return;
				for(Sample sample: samples) {
					sample.selectedInTable= false;
					}
				final  int[] idxy  = this.XYTable.getSelectedRows();
				for(int i1:idxy) {
					int i2 = this.XYTable.convertRowIndexToModel(i1);
					if(i2<0) continue;
					samples.get(i2).selectedInTable=true;
					}
				drawingArea.repaint();
				});
			final DefaultTableCellRenderer renderer = new DefaultTableCellRenderer() {
				public java.awt.Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
					java.awt.Component c = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
					if(!isSelected) {
						c.setBackground(row%2==0?Color.WHITE:Color.LIGHT_GRAY);
						if(column==0 || column==1) {
							row = table.convertRowIndexToModel(row);
							if(row==-1) return c;
							final SampleTableModel tm = (SampleTableModel)XYTable.getModel();
							c.setBackground(tm.rows.get(row).selectForOutput?Color.GREEN:Color.RED);
							}
						}
					return c;
					}
				};
				this.XYTable.setDefaultRenderer(Object.class, renderer);
				this.XYTable.setDefaultRenderer(String.class, renderer);
				this.XYTable.setDefaultRenderer(Double.class, renderer);
			pane2.add(new JScrollPane(this.XYTable));
			
		
		
		mainPanel.add(
				new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,this.drawingArea,pane2),
				BorderLayout.CENTER
				);
		}
	
	
	private boolean saveAs(final File fname) {
		if(fname==null) {
			return saveAs();
			}
		try( PrintWriter pw = IOUtils.openFileForPrintWriter(fname)) {
			for(int i=0;i< this.samples.size();i++) {
				final Sample sn = this.samples.get(i);
				pw.print(sn.fid);
				pw.print("\t");
				pw.print(sn.iid);
				pw.print("\t");
				pw.print(sn.selectForOutput?".":"EXCLUDED");
				pw.println();
				}
			pw.flush();
			}
		catch(final Throwable err) {
			ThrowablePane.show(this, err);
			return false;
			}
		dirtyFlag=false;
		saveAsFile=fname;
		return true;
		}
	
	
	private boolean saveAs() {
		final JFileChooser chooser = new JFileChooser(PreferredDirectory.get(SwingPLinkSelectCluster.class));
		int ret=chooser.showSaveDialog(this);
		if(ret==JFileChooser.CANCEL_OPTION) return false;
		
		final File f = chooser.getSelectedFile();
		if(f.exists() && JOptionPane.showConfirmDialog(this, "File "+f.getName()+" exists. Overwite ?", "Overwite ?", JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE, null)!=JOptionPane.OK_OPTION)
			{
			return false;
			}
		if(saveAs(f)) {
			PreferredDirectory.update(f);
			return true;
			}
		return false;
		}
	
	private Point2D sampleToPixel(ColXY columnXY,Sample sn) {
		final double width = drawingArea.getWidth();
		final double height = drawingArea.getHeight();
		Column col = columnXY.get(0);
		double x = sn.C[col.index];
		x = ((x-col.min)/(col.max-col.min))*width;
		col =columnXY.get(1);
		double y = sn.C[col.index];
		y = ((y-col.min)/(col.max-col.min))*height;
		return new Point2D.Double(x,y);
		}
	
	private ColXY getCurrentColXY() {
		 for (Enumeration<AbstractButton> iter = this.buttonGroupXY.getElements();
				iter.hasMoreElements();) {
	            final JToggleButton button = (JToggleButton)iter.nextElement();

	            if (button.isSelected()) {
	                return ColXY.class.cast(AbstractAction.class.cast(button.getAction()).getValue(ACTION_XY_KEY));
	            }
	        }
		 return null;
		}

	private ToolType getCurrentTooltype() {
		 for (Enumeration<AbstractButton> iter = this.buttonGroupTool.getElements();
				iter.hasMoreElements();) {
	            final JToggleButton button = (JToggleButton)iter.nextElement();

	            if (button.isSelected()) {
	                return ToolType.class.cast(AbstractAction.class.cast(button.getAction()).getValue(ACTION_TOOL));
	            }
	        }
		 return ToolType.XOR;
		}
	

	private void paintDrawingArea(Graphics2D g) {
		final int width = this.drawingArea.getWidth();
		final int height = this.drawingArea.getHeight();
		g.setColor(Color.LIGHT_GRAY);
		g.fillRect(0, 0,width,height);
		g.setColor(Color.DARK_GRAY);
		g.drawRect(0, 0,width-1,height-1);
		final ColXY columnXY = getCurrentColXY();
		if(columnXY==null) return;
		final int radius =5;
		
		for(int side=0;side<2;++side) {
			for(int i=0;i< this.samples.size();i++) {
				final Sample sn = this.samples.get(i);
				final Point2D pt = sampleToPixel(columnXY, sn);
				
				
				if(side==1) {
					if(sn.selectedInTable) {
						g.setColor(Color.BLUE);
						g.draw(new Ellipse2D.Double(
								pt.getX()-(radius+1),
								pt.getY()-(radius+2),
								(radius+1)*2,
								(radius+1)*2)
								);
						}
					}
				else
					{
					
					if(!sn.selectForOutput) {
						g.setColor(Color.RED);
						g.draw(new Line2D.Double(pt.getX()-radius, pt.getY(), pt.getX()+radius, pt.getY()));
						g.draw(new Line2D.Double(pt.getX(), pt.getY()-radius, pt.getX(), pt.getY()+radius));
						}
					else
						{
						g.setColor(Color.GREEN);
						g.fill(new Ellipse2D.Double(pt.getX()-radius, pt.getY()-radius, radius*2,radius*2));
						g.setColor(Color.GRAY);
						g.draw(new Ellipse2D.Double(pt.getX()-radius, pt.getY()-radius, radius*2,radius*2));
						}
					}
				}
			if(this.XYTable.getSelectionModel().isSelectionEmpty()) break;
			}
		}

}


@Override
public int doWork(List<String> args) {
	try {
		final List<Column> columns = new Vector<>();
		final List<Sample> rows = new Vector<>();
		final String input = super.oneAndOnlyOneFile(args);
		try(BufferedReader br=IOUtils.openPathForBufferedReading(Paths.get(input))) {
			final String hdr = br.readLine();
			if(hdr==null) throw new IOException("Cannot read header from "+input);
			final Pattern ws = Pattern.compile("\\s+");
			final Function<String, List<String>> splitter = S->Arrays.asList(ws.split(S.trim()));;
			final FileHeader header= new FileHeader(hdr,splitter);
			for(int i=0;i< header.size();i++) {
				final String col = header.get(i);
				if(col.matches("C[0-9]+")) {
					columns.add(new Column(columns.size(),col));
					}
				}
			if(columns.size()<2) {
				LOG.error("Zero or only one column matching C[0-9]+ was found in "+header);
				return -1;
				}
			for(;;) {
				final String line = br.readLine();
				if(line==null) break;
				final FileHeader.RowMap row= header.toMap(line);
				final Sample sample = new Sample();
				sample.fid = row.get("FID");
				sample.iid = row.get("IID");
				sample.C = new double[columns.size()];
				for(Column column: columns) {
					sample.C[column.index] = Double.parseDouble(row.get(column.label));
					}
				rows.add(sample);
				}
			}
		
		if(rows.isEmpty()) {
			LOG.error("no row in  "+input);
			return -1;
			}
		
		for(Column col:columns) {
			col.min = rows.stream().mapToDouble(R->R.C[col.index]).min().getAsDouble();
			col.max = rows.stream().mapToDouble(R->R.C[col.index]).max().getAsDouble();
			double diff = col.max - col.min;
			diff= diff*0.05;
			col.min-=diff;
			col.max+=diff;
		}

		
		final XFrame frame = new XFrame(columns,rows);
		final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
		frame.setBounds(50, 50, screen.width-100, screen.height-100);

		
		SwingUtilities.invokeAndWait(()->{
			frame.setVisible(true);
			});

		return 0;
		}
	catch(Throwable err) {
		LOG.error(err);
		return -1;
		}
	}

public static void main(String[] args) {
	new SwingPLinkSelectCluster().instanceMain(args);
}
}
