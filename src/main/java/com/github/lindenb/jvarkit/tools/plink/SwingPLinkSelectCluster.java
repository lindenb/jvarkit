
/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
import java.awt.Shape;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Vector;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.table.AbstractTableModel;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.swing.SAMRecordPanel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.swing.AbstractGenericTableModel;
import com.github.lindenb.jvarkit.swing.PropertyChangeObserver;
import com.github.lindenb.jvarkit.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.iterator.LineIterators;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFGenotypesTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFInfoTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVariantsTableModel;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;

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
private static class Column {
	final int index;
	final String label;
	final int index_in_file;
	double min = 0.0;
	double max = 0.0;
	Column(int index,String label, int index_in_file) {
		this.index  = index;
		this.label = label;
		this.index_in_file= index_in_file;
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
	boolean selected=false;
	int state = 0;
	}

private static class SampleTableModel extends AbstractTableModel{
	final ColXY columnXY;
	final List<Sample> rows;
	 SampleTableModel(final ColXY columnXY,final List<Sample> rows) {
		this.rows= rows;
		this.columnXY = columnXY;
	 	}
	 
	@Override
	public Class<?> getColumnClass(int columnIndex) {
		 switch(columnIndex) {
			case 0:
			case 1: return String.class;
			case 2:
			case 3: return Double.class;
			default: return Object.class;
			}
	 	}
	@Override
	public int getRowCount() {
		return this.rows.size();
		}
	@Override
	public int getColumnCount() {
		return 4;
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
			case 2: return this.columnXY.get(0).label;
			case 3: return this.columnXY.get(1).label;
			default: return null;
			}
		}
	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		final Sample o = this.rows.get(rowIndex);
		switch(columnIndex) {
			case 0: return o.fid;
			case 1: return o.iid;
			case 2: return o.C[this.columnXY.get(0).index];
			case 3: return o.C[this.columnXY.get(1).index];
			default: return null;
			}
		}
	}

@SuppressWarnings("serial")
private static class XFrame extends JFrame {
	final List<Sample> samples;
	final List<ColXY> columnsXY;
	final JPanel drawingArea ;
	final List<JTable> XYTables = new Vector<>();
	final JTabbedPane jTabbedPane ;
	XFrame(final List<ColXY> columnsXY,
			final List<Sample> samples
			)
		{
		super(SwingPLinkSelectCluster.class.getSimpleName());
		this.columnsXY=columnsXY;
		this.samples= samples;
		super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		
		final JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		JMenu menu=new JMenu("File");
		menuBar.add(menu);
		menu.add(new JSeparator());
		menu.add(new JMenuItem(new AbstractAction("Quit") {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				XFrame.this.setVisible(false);
				XFrame.this.dispose();
				}
			}));
		final JPanel mainPanel=new JPanel(new BorderLayout());
		this.setContentPane(mainPanel);
		
		this.drawingArea = new JPanel(null) {
			@Override
			public String getToolTipText(MouseEvent event) {
				final int width = this.getWidth();
				final int height = this.getHeight();
				int idx = jTabbedPane.getSelectedIndex();
				if(idx<0) return "";
				StringBuilder sb = new StringBuilder();
				SampleTableModel tm = (SampleTableModel)XYTables.get(idx).getModel();
				for(int i=0;i< samples.size();i++) {
					final Sample sn = samples.get(i);
					Column col = tm.columnXY.get(0);
					double x = sn.C[col.index];
					x = ((x-col.min)/(col.max-col.min))*width;
					col = tm.columnXY.get(1);
					double y = sn.C[col.index];
					y = ((y-col.min)/(col.max-col.min))*height;
					if(Point.distance(event.getX(), event.getY(), x, y) < 3) {
						if(sb.length()>0) sb.append(" ");
						sb.append(sn.iid);
						}
					}
				return sb.toString();
				}
			@Override
			protected void paintComponent(Graphics g) {
				paintDrawingArea(Graphics2D.class.cast(g));
				}
			};
		drawingArea.setToolTipText("x");
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
				int idx = jTabbedPane.getSelectedIndex();
				SampleTableModel tm = (SampleTableModel)XYTables.get(idx).getModel();
				for(int i=0;i< samples.size();i++) {
					final Sample sn = samples.get(i);
					Column col = tm.columnXY.get(0);
					double x = sn.C[col.index];
					x = ((x-col.min)/(col.max-col.min))*drawingArea.getWidth();
					col = tm.columnXY.get(1);
					double y = sn.C[col.index];
					y = ((y-col.min)/(col.max-col.min))*drawingArea.getHeight();
					if(gpath.contains(x, y)) {
						sn.state = 1;
						}
					}
				drawingArea.repaint();
				gpath= null;
				prev=null;
				}
			};
		this.drawingArea.addMouseListener(mouse);
		this.drawingArea.addMouseMotionListener(mouse);
		this.jTabbedPane = new JTabbedPane();
		for(ColXY cxy: columnsXY ) {
			JPanel pane = new JPanel(new BorderLayout());
			this.jTabbedPane.addTab(cxy.getName(), pane);
			
			final JTable table = new JTable(new SampleTableModel(cxy,samples));
			table.getSelectionModel().addListSelectionListener(E->{
				if(E.getValueIsAdjusting()) return;
				resetSelectedSamples();
				int[] idxy  = table.getSelectedRows();
				for(int i1:idxy) {
					int i2 = table.convertRowIndexToModel(i1);
					if(i2<0) continue;
					samples.get(i2).selected=true;
					}
				drawingArea.repaint();
				});
			this.XYTables.add(table);
			pane.add(new JScrollPane(table));
			}
		
		this.jTabbedPane.addChangeListener((AE)->drawingArea.repaint());
		
		mainPanel.add(
				new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,this.drawingArea,jTabbedPane),
				BorderLayout.CENTER
				);
		}
	
	private void resetSelectedSamples() {
		for(Sample sample: samples) sample.selected= false;
	}
	
	private void paintDrawingArea(Graphics2D g) {
		final int width = this.drawingArea.getWidth();
		final int height = this.drawingArea.getHeight();
		g.setColor(Color.LIGHT_GRAY);
		g.fillRect(0, 0,width,height);
		g.setColor(Color.DARK_GRAY);
		g.drawRect(0, 0,width-1,height-1);
		int idx = this.jTabbedPane.getSelectedIndex();
		if(idx<0) return;
		SampleTableModel tm = (SampleTableModel)this.XYTables.get(idx).getModel();
		for(int i=0;i< this.samples.size();i++) {
			final Sample sn = this.samples.get(i);
			Column col = tm.columnXY.get(0);
			double x = sn.C[col.index];
			x = ((x-col.min)/(col.max-col.min))*width;
			col = tm.columnXY.get(1);
			double y = sn.C[col.index];
			y = ((y-col.min)/(col.max-col.min))*height;
			g.setColor(sn.selected?Color.RED:Color.BLUE);
			if(sn.state==0) {
				g.draw(new Line2D.Double(x-5, y, x+5, y));
				g.draw(new Line2D.Double(x, y-5, x, y+5));
				}
			else
				{
				g.fill(new Ellipse2D.Double(x-5, y-5, 10, 10));
				}
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
					columns.add(new Column(columns.size(),col,i));
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
		final List<ColXY> colsXY = new Vector<>();
		for(int x=0;x< columns.size();++x) {
			for(int y=x+1;y< columns.size();++y) {
				colsXY.add( new ColXY(columns.get(x),columns.get(y)));
			}
		}
		
		final XFrame frame = new XFrame(colsXY,rows);
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
