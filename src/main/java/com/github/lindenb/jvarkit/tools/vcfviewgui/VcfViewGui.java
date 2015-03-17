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
package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;
import java.util.regex.Pattern;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JCheckBox;
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
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableModel;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;

import com.github.lindenb.jvarkit.util.igv.IgvSocket;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.util.picard.SamFlag;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;


/** base class for VCF or BAM */
abstract class AbstractFileRef
	{
	File filePath;
	boolean fileisIndexed=false;
	abstract boolean isVCF();
	abstract boolean isBAM();
	boolean isIndexed()
		{
		return this.fileisIndexed;
		}
	}
/**
 * VCF header, codec and header
 */
class VCFFileRef extends AbstractFileRef
	{
	VCFUtils.CodecAndHeader codhead;
	@Override
	final boolean isVCF() { return true;}
	@Override
	final boolean isBAM() { return false;}
	}

/**
 * BAm reference
 */
class BamFileRef extends AbstractFileRef
	{
	SAMFileHeader header;
	@Override
	final boolean isVCF() { return false;}
	@Override
	final boolean isBAM() { return true;}
	}




/**
 * AbstractRefInternalFrame
 */
abstract class AbstractRefInternalFrame extends JInternalFrame
	{
	private static final long serialVersionUID = 1L;
	private AbstractFileRef refFile;
	AbstractRefInternalFrame(AbstractFileRef refFile)
		{
		super(refFile.filePath.getName()+" Indexed:"+refFile.isIndexed(),true,false,true,true);
		this.refFile=refFile;
		}
	public AbstractFileRef getRefFile()
		{
		return refFile;
		}
	}

/**
 * Main FRAME
 */
class VCFInternalFrame extends AbstractRefInternalFrame
	{
	private static final long serialVersionUID = 1L;

	JTable jTable;
	VCFTableModel tableModel;
	InfoTreeModel infoTreeModel;
	private JTree infoTree;
	GenotypeTableModel genotypeTableModel;
	private ListSelectionListener selList;
	VCFInternalFrame(VCFFileRef ref)
		{
		super(ref);
		JPanel mainPane=new JPanel(new BorderLayout(5,5));
		setContentPane(mainPane);
		JTabbedPane tabbedPane=new JTabbedPane();
		mainPane.add(tabbedPane,BorderLayout.CENTER);

		
		
		
		JPanel pane=new JPanel(new BorderLayout(5,5));
		tabbedPane.addTab("VCF", pane);
		
		this.tableModel=new VCFTableModel(ref);
		this.jTable=new JTable(tableModel);
		this.jTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		this.jTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		
		JScrollPane scroll1=new JScrollPane(this.jTable);
		
		this.infoTreeModel=new InfoTreeModel();
		this.infoTree=new JTree(this.infoTreeModel);
		
		
		this.genotypeTableModel=new GenotypeTableModel();
		JTable tGen=new JTable(this.genotypeTableModel);
		
		
		JSplitPane splitH=new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
				new JScrollPane(this.infoTree),
				new JScrollPane(tGen)
				);
	
		JSplitPane splitVert=new JSplitPane(JSplitPane.VERTICAL_SPLIT,scroll1,splitH);
		
		this.jTable.getSelectionModel().addListSelectionListener(new ListSelectionListener()
			{
			@Override
			public void valueChanged(ListSelectionEvent e)
				{
				if(e.getValueIsAdjusting()) return;
		        int row = jTable.getSelectedRow();
		        VariantContext ctx;
				if(row==-1 || (ctx=tableModel.getVariantContext(row))==null)
					{
					infoTreeModel.setContext(null,VCFInternalFrame.this.getVCFFileRef().codhead.header);
					genotypeTableModel.setContext(null);
					}
				else
					{
					infoTreeModel.setContext(ctx,VCFInternalFrame.this.getVCFFileRef().codhead.header);
					genotypeTableModel.setContext(ctx);
					}
				
				}
			});
		
		pane.add(splitVert);
		
		pane=new JPanel(new BorderLayout(5,5));
		
		ByteArrayOutputStream baos=new ByteArrayOutputStream();
		VariantContextWriter w=VCFUtils.createVariantContextWriterToOutputStream(baos);
		w.writeHeader(ref.codhead.header);
		
		tabbedPane.addTab("Header", pane);
		JTextArea area=new JTextArea(new String(baos.toByteArray()));
		area.setCaretPosition(0);
		area.setEditable(false);
		pane.add(new JScrollPane(area),BorderLayout.CENTER);
		
		this.selList=new ListSelectionListener()
			{
			@Override
			public void valueChanged(ListSelectionEvent e)
				{
				if(e.getValueIsAdjusting()) return;
				listSelectionChanged();
				}
			};
		
		this.addInternalFrameListener(new InternalFrameAdapter()
			{
			@Override
			public void internalFrameActivated(InternalFrameEvent e)
				{
				jTable.getSelectionModel().addListSelectionListener(selList);
				}
			@Override
			public void internalFrameDeactivated(InternalFrameEvent e) {
				jTable.getSelectionModel().removeListSelectionListener(selList);
				}
			});
		}
	
	VCFFileRef getVCFFileRef()
		{
		return VCFFileRef.class.cast(super.getRefFile());
		}
	
	private void listSelectionChanged()
		{
		int row=jTable.getSelectedRow();
		if(row==-1 || this.getDesktopPane()==null) return;
		VariantContext ctx=this.tableModel.getVariantContext(row);
		
		if(ctx==null) return;
		
		for(JInternalFrame jif:this.getDesktopPane().getAllFrames())
			{
			if(jif==this) continue;
			if(jif.getClass()!=this.getClass() ) continue;
			VCFInternalFrame other=VCFInternalFrame.class.cast(jif);
			int row2=-1;
			for(int i=0;i< other.tableModel.getRowCount();++i)
				{
				VariantContext ctx2=other.tableModel.getVariantContext(i);
				if(ctx.getChr().equals(ctx2.getChr()) &&
						ctx.getStart()==ctx2.getStart()
						)
					{
					row2=i;
					break;
					}
				}
			if(row2==-1)
				{
				other.jTable.getSelectionModel().clearSelection();
				}
			else
				{
				other.jTable.getSelectionModel().setSelectionInterval(row2, row2);
				other.jTable.scrollRectToVisible((other.jTable.getCellRect(row2,0, true)));
				}
			}
		}
	
	}

class InfoTreeModel
	extends DefaultTreeModel
	{
	private static Log LOG=Log.getInstance(VcfViewGui.class);
	private Pattern pipeVep=Pattern.compile("[\\|]");
	private Pattern pipeSnpEff=Pattern.compile("[\\|\\(\\)]");

	private static final long serialVersionUID = 1L;
	public InfoTreeModel()
		{
		super(new DefaultMutableTreeNode("INFO",true));
		}
	
	private DefaultMutableTreeNode getTreeNodeRoot()
		{
		return DefaultMutableTreeNode.class.cast(getRoot());
		}	
	
	public void setContext(VariantContext ctx,VCFHeader header)
		{
		getTreeNodeRoot().removeAllChildren();
		if(ctx!=null)
			{
			for(String key:ctx.getAttributes().keySet())
				{
				Object v=ctx.getAttribute(key);
				Object o[];
				if(v==null)
					{
					o=new Object[]{null};
					}
				else if(v instanceof java.util.Collection)
					{
					o=((java.util.Collection<?>)v).toArray();
					}
				else if(v.getClass().isArray())
					{
					o=(Object[])v;
					}
				else
					{
					o=new Object[]{v};
					}
				if(o.length==0)
					{
					
					}
				else if(o.length==1)
					{
					if(key.equals("CSQ") || key.equals("EFF"))
						{
						List<String> columns= key.equals("CSQ")?getCSQCols(header):getEFFCols(header);
						DefaultMutableTreeNode n1=new DefaultMutableTreeNode(
								"<html><b>"+key+"</html>"
								,true);
						getTreeNodeRoot().add(n1);
						String tokens[]=(key.equals("CSQ")?pipeVep:pipeSnpEff).split(String.valueOf(o));
						for(int i=0;i< tokens.length && i< columns.size();++i)
							{
							DefaultMutableTreeNode n2=new DefaultMutableTreeNode(
									"<html><b>"+columns.get(i)+"</b>:"+tokens[i]+"</html>"
									,false);
							n1.add(n2);
							}
						}
					else
						{
						DefaultMutableTreeNode n=new DefaultMutableTreeNode(
								"<html><b>"+key+"</b>:"+o[0]+"</html>"
								,false);
						getTreeNodeRoot().add(n);
						}
					}
				else 
					{
					if(key.equals("CSQ") || key.equals("EFF"))
						{
					
						List<String> columns= key.equals("CSQ")?getCSQCols(header):getEFFCols(header);
						DefaultMutableTreeNode n1=new DefaultMutableTreeNode(
								"<html><b>"+key+"</html>"
								,true);
						getTreeNodeRoot().add(n1);
						int index=0;
						for(Object v2:o)
							{
							DefaultMutableTreeNode n2=new DefaultMutableTreeNode(
									String.valueOf(++index),
									true
									);
							n1.add(n2);
							String tokens[]= (key.equals("CSQ")?pipeVep:pipeSnpEff).split(String.valueOf(v2));
							for(int i=0;i< tokens.length && i< columns.size();++i)
								{
								DefaultMutableTreeNode n3=new DefaultMutableTreeNode(
										"<html><b>"+columns.get(i)+"</b>:"+tokens[i]+"</html>"
										,false);
								n2.add(n3);
								}
							}
						}
					else
						{
						DefaultMutableTreeNode n1=new DefaultMutableTreeNode(
								"<html><b>"+key+"</b></html>"
								,true);
						getTreeNodeRoot().add(n1);
						for(Object v2:o)
							{
							DefaultMutableTreeNode n2=new DefaultMutableTreeNode(
									String.valueOf(v2),
									false
									);
	
							n1.add(n2);
							}
						}
					}
				
				}
			}
		this.fireTreeStructureChanged();
		}
	private List<String> getCSQCols(VCFHeader header)
			{
			VCFInfoHeaderLine ihl=header.getInfoHeaderLine("CSQ");
			if(ihl==null) return Collections.emptyList();
			String description=ihl.getDescription();
			String chunck=" Format:";
			int i=description.indexOf(chunck);
			if(i==-1)
				{
				LOG.warn("Cannot find "+chunck+ " in "+description);
				return Collections.emptyList();
				}
			description=description.substring(i+chunck.length()).replaceAll("[ \'\\.\\(\\)]+","").trim();
			String tokens[]=pipeVep.split(description);
			ArrayList<String> L=new ArrayList<String>(tokens.length);
			for(String s:tokens)
				{
				if(s.trim().isEmpty()) continue;
				L.add(s);
				}
			return L;
			}
	private List<String> getEFFCols(VCFHeader header)
		{
		VCFInfoHeaderLine ihl=header.getInfoHeaderLine("EFF");
		if(ihl==null) return Collections.emptyList();
		String description=ihl.getDescription();
		String chunck="Format:";
		int i=description.indexOf(chunck);
		if(i==-1)
			{
			LOG.warn("Cannot find "+chunck+ " in "+description);
			return Collections.emptyList();
			}
		description=description.substring(i+chunck.length()).replace('(','|').replaceAll("[ \'\\.)\\[\\]]+","").trim();
		String tokens[]=pipeSnpEff.split(description);
		ArrayList<String> L=new ArrayList<String>(tokens.length);
		for(String s:tokens)
			{
			if(s.trim().isEmpty()) continue;
			L.add(s);
			}
		return L;
		}
	
	
	public void fireTreeStructureChanged()
		{
		fireTreeStructureChanged(this, new Object[]{getRoot()}, null, null);
		}	
	}

class GenotypeTableModel
extends AbstractTableModel
	{
	private static final long serialVersionUID = 1L;
	private List<Genotype> genotypes=new Vector<Genotype>();
	private static final String COLS[]=new String[]{"SAMPLE","TYPE","DP","GQ"};
	public GenotypeTableModel()
		{
		
		}
	
	public void setContext(VariantContext ctx)
		{
		if(ctx==null)
			{
			this.genotypes.clear();
			}
		else
			{
			this.genotypes=new Vector<Genotype>(ctx.getGenotypes());
			}
		super.fireTableDataChanged();
		}
	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		return false;
		}
	@Override
	public Class<?> getColumnClass(int columnIndex)
		{
		switch(columnIndex)
			{
			case 0: return String.class;
			case 1: return String.class;
			case 2: return Integer.class;
			case 3: return Integer.class;
			default: return Object.class;
			}
		}
	@Override
	public String getColumnName(int column)
		{
		return COLS[column];
		}
	@Override
	public Object getValueAt(int rowIndex, int columnIndex)
		{
		Genotype g=this.genotypes.get(rowIndex);
		if(g==null ) return null;
		switch(columnIndex)
			{
			case 0: return g.getSampleName();
			case 1: return g.getType()==null?null:g.getType().name();
			case 2: return g.hasDP()?g.getDP():null;
			case 3: return g.hasGQ()?g.getGQ():null;
			default: return null;
			}
		}
	@Override
	public int getRowCount() {
		return genotypes.size();
		}
	@Override
	public int getColumnCount() {
		return  COLS.length;
		}
}



class VCFTableModel
	extends AbstractTableModel
	{
	private static final long serialVersionUID = 1L;

	private VCFFileRef ref;
	private List<String> rows=new Vector<String>();
	private Pattern tab=Pattern.compile("[\t]");
	
	
	VCFTableModel(VCFFileRef ref)
		{
		this.ref=ref;
		}
	
	@Override
	public int getRowCount()
		{
		return rows.size();
		}
	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex)
		{
		return false;
		}
	
	VariantContext getVariantContext(int rowIndex)
		{
		String s=this.rows.get(rowIndex);
		return this.ref.codhead.codec.decode(s);
		}
	
	@Override
	public Object getValueAt(int rowIndex, int columnIndex)
		{
		String s=this.rows.get(rowIndex);
		String tokens[]=tab.split(s,columnIndex+2);
		return columnIndex<tokens.length?tokens[columnIndex]:null;
		}
	
	@Override
	public String getColumnName(int column)
		{
		VCFHeader.HEADER_FIELDS headers[]= VCFHeader.HEADER_FIELDS.values();
		if(column< headers.length)
			{
			return headers[column].name();
			}
		if(column==headers.length)
			{
			return "FORMAT";
			}
		column-=(headers.length+1);
		return ref.codhead.header.getSampleNamesInOrder().get(column);
		}
	
	@Override
	public int getColumnCount()
		{
		int n=VCFHeader.HEADER_FIELDS.values().length;
		if(this.ref.codhead.header.getSampleNamesInOrder().isEmpty())
			{
			return n;
			}
		return n+1+ref.codhead.header.getSampleNamesInOrder().size();
		}
	
	synchronized void updateRow(List<String> L)
		{
		this.rows.clear();
		this.rows.addAll(L);
		this.fireTableDataChanged();
		}
	}
class VCFPos
	{
	String chrom;
	int start;
	int end=-1;
	@Override
	public String toString() {
		return chrom+":"+start+(end==-1?"":"-"+end);
		}
	}




@SuppressWarnings("serial")
class BamInternalFrame extends AbstractRefInternalFrame
	{
	private static final long serialVersionUID = 1L;

	JTable jTable;
	BamTableModel tableModel;
	BamFileRef ref;
	FlagTableModel infoTableModel;
	ReadGroupTableModel groupTableModel;
	SAMTagAndValueModel genotypeTableModel;
	private ListSelectionListener selList;
	BamInternalFrame(BamFileRef ref)
		{
		super(ref);
		this.ref=ref;
		JPanel mainPane=new JPanel(new BorderLayout(5,5));
		setContentPane(mainPane);
		JTabbedPane tabbedPane=new JTabbedPane();
		mainPane.add(tabbedPane,BorderLayout.CENTER);

		
		
		
		JPanel pane=new JPanel(new BorderLayout(5,5));
		tabbedPane.addTab("BAM", pane);
		
		this.tableModel=new BamTableModel();
		this.jTable=createTable(tableModel);
		this.jTable.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		this.jTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		
		JScrollPane scroll1=new JScrollPane(this.jTable);
		
		this.infoTableModel=new FlagTableModel();
		JTable tInfo=createTable(this.infoTableModel);
		
		
		this.genotypeTableModel=new SAMTagAndValueModel();
		JTable tGen=createTable(this.genotypeTableModel);
		
		
		this.groupTableModel=new ReadGroupTableModel();
		JTable tGrp=createTable(this.groupTableModel);
	
		
		JPanel splitH=new JPanel(new GridLayout(1, 0,5,5));
		splitH.add(new JScrollPane(tInfo));
		splitH.add(new JScrollPane(tGen));
		splitH.add(new JScrollPane(tGrp));
			
	
		JSplitPane splitVert=new JSplitPane(JSplitPane.VERTICAL_SPLIT,scroll1,splitH);
		
		this.jTable.getSelectionModel().addListSelectionListener(new ListSelectionListener()
			{
			@Override
			public void valueChanged(ListSelectionEvent e)
				{
				if(e.getValueIsAdjusting()) return;
		        int row = jTable.getSelectedRow();
		        SAMRecord ctx;
				if(row==-1 || (ctx=tableModel.getElementAt(row))==null)
					{
					infoTableModel.setContext(null);
					genotypeTableModel.setContext(null);
					groupTableModel.setContext(null);
					}
				else
					{
					infoTableModel.setContext(ctx);
					genotypeTableModel.setContext(ctx);
					groupTableModel.setContext(ctx);
					}
				
				}
			});
		
		pane.add(splitVert);
		
		//header as text
		pane=new JPanel(new BorderLayout(5,5));
		tabbedPane.addTab("Header", pane);
		JTextArea area=new JTextArea(String.valueOf(ref.header.getTextHeader()));
		area.setCaretPosition(0);
		area.setEditable(false);
		pane.add(new JScrollPane(area),BorderLayout.CENTER);
		
		//dict
		pane=new JPanel(new BorderLayout(5,5));
		tabbedPane.addTab("Reference", pane);
		JTable dictTable=createTable(new SAMSequenceDictionaryTableModel(ref.header.getSequenceDictionary()));
		dictTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		pane.add(new JScrollPane(dictTable),BorderLayout.CENTER);
		
		
		this.selList=new ListSelectionListener()
			{
			@Override
			public void valueChanged(ListSelectionEvent e)
				{
				if(e.getValueIsAdjusting()) return;
				listSelectionChanged();
				}
			};
		
		this.addInternalFrameListener(new InternalFrameAdapter()
			{
			@Override
			public void internalFrameActivated(InternalFrameEvent e)
				{
				jTable.getSelectionModel().addListSelectionListener(selList);
				}
			@Override
			public void internalFrameDeactivated(InternalFrameEvent e) {
				jTable.getSelectionModel().removeListSelectionListener(selList);
				}
			});
		}
	private static final Color color1=Color.WHITE;
	private static final Color color2=new Color(255,255,230);
	
	private static JTable createTable(TableModel m)
		{
		JTable t=new JTable(m)
			{
			@Override
			public String getToolTipText(MouseEvent event)
				{
				JTable t=(JTable)event.getSource();
				int x= t.rowAtPoint(event.getPoint()); if(x==-1) return null;
				int y= t.columnAtPoint(event.getPoint()); if(y==-1) return null;
				Object o= t.getValueAt(x, y);
				if(o==null) return null;
				return String.valueOf(o);
				}
			
			
			};
		t.setToolTipText("");
		t.setShowVerticalLines(false);
		DefaultTableCellRenderer render=new DefaultTableCellRenderer()
			{
			@Override
			public Component getTableCellRendererComponent(JTable table,
					Object value, boolean isSelected, boolean hasFocus,
					int row, int column) {
				Component c= super.getTableCellRendererComponent(table, value, isSelected, hasFocus,
						row, column);
				if(!isSelected && !hasFocus) this.setBackground(row%2==0?color1:color2);
				if(value !=null && value instanceof Boolean)
					{
					if(Boolean.TRUE.equals(value)) this.setText("\u2612");
					else if(Boolean.FALSE.equals(value)) this.setText("");
					}
				return c;
				}
			};
		render.setOpaque(true);
		for(int i=0;i< t.getColumnModel().getColumnCount();++i)
			{
			t.getColumnModel().getColumn(i).setCellRenderer(render);
			}
		return t;
		}
	
	private void listSelectionChanged()
		{
		int row=jTable.getSelectedRow();
		if(row==-1 || this.getDesktopPane()==null) return;
		SAMRecord ctx=this.tableModel.getElementAt(row);
		
		if(ctx==null) return;
		
		for(JInternalFrame jif:this.getDesktopPane().getAllFrames())
			{
			if(jif==this) continue;
			if(jif.getClass()!=this.getClass() ) continue;
			BamInternalFrame other=BamInternalFrame.class.cast(jif);
			int row2=-1;
			for(int i=0; !ctx.getReadUnmappedFlag() &&
					i< other.tableModel.getRowCount();
					++i)
				{
				SAMRecord ctx2=other.tableModel.getElementAt(i);
				if(ctx2.getReadUnmappedFlag()) continue;
				if(ctx.getReferenceName().equals(ctx2.getReferenceName()) &&
						ctx.getAlignmentStart()<=ctx2.getAlignmentStart()
						)
					{
					row2=i;
					break;
					}
				}
			if(row2==-1)
				{
				other.jTable.getSelectionModel().clearSelection();
				}
			else
				{
				other.jTable.getSelectionModel().setSelectionInterval(row2, row2);
				other.jTable.scrollRectToVisible((other.jTable.getCellRect(row2,0, true)));
				}
			}
		}
	
	}

class FlagTableModel
	extends AbstractTableModel
	{
	private static final long serialVersionUID = 1L;
	private int flag=0;
	public FlagTableModel()
		{
		}
	
	public void setContext(SAMRecord rec)
		{
		setContext(rec==null ?0:rec.getFlags());
		}
	
	public void setContext(int flag)
		{
		this.flag=flag;
		this.fireTableDataChanged();
		}
	
	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		return false;
		}
	@Override
	public Class<?> getColumnClass(int columnIndex)
		{
		switch(columnIndex)
			{
			case 0: return String.class;
			case 1: return Boolean.class;
			}
		return Object.class;
		}
	@Override
	public String getColumnName(int column)
		{
		return column==0?"Flag":"ON/OFF";
		}
	@Override
	public Object getValueAt(int rowIndex, int columnIndex)
		{
		SamFlag f=SamFlag.values()[rowIndex];
		return (columnIndex==0?
				f.name():
				f.isSet(this.flag)
				);
		}
	@Override
	public int getRowCount() {
		return SamFlag.values().length;
		}
	@Override
	public int getColumnCount() {
		return 2;
		}
	}

@SuppressWarnings("serial")
class SAMTagAndValueModel
extends AbstractGenericTable<SAMTagAndValue>
	{
	private static final long serialVersionUID = 1L;
	private static final String COLS[]=new String[]{"KEY","DESC","VALUE","TYPE"};
	public SAMTagAndValueModel()
		{
		
		}
	
	public void setContext(SAMRecord ctx)
		{
		if(ctx==null)
			{
			super.clear();
			}
		else
			{
			super.setRows(ctx.getAttributes());
			
			}
		}
	
	@Override
	public Object getValueOf(SAMTagAndValue o, int columnIndex)
		{
		if(o==null) return null;
		switch(columnIndex)
			{
			case 0: return o.tag;
			case 1: return getDescription(o.tag);
			case 2: return o.value;
			case 3: return o.value==null?null:o.value.getClass().getSimpleName();
			
			}
		return null;
		}
	
	@Override
	public Class<?> getColumnClass(int columnIndex) {
		switch(columnIndex)
			{
			case 0: 
			case 1:
			case 3: return String.class;
			default: return Object.class;
			}
		}
	
	@Override
	public String getColumnName(int column)
		{
		return COLS[column];
		}
	
	@Override
	public int getColumnCount() {
		return  COLS.length;
		}
	
	private static final Map<String,String> defs=new HashMap<String,String>()
			{{{
			put("AM","The smallest template-independent mapping quality of segments in the rest ");
		  put("AS","Alignment score generated by aligner ");
		  put("BC","Barcode sequence, with any quality scores stored in the QT tag. ");
		  put("CC","Reference name of the next hit; \"=\" for the same chromosome ");
		  put("CM","Edit distance between the color sequence and the color reference (see also NM)");
		  put("CO","Free-text comments ");
		  put("CQ","Color read quality on the original strand of the read. Same encoding as QUAL; same length as CS.");
		  put("CS","Color read sequence on the original strand of the read. The primer base must be included.");
		  put("CT","Complete read annotation tag, used for consensus annotation dummy features.");
		  put("E2","The 2nd most likely base calls. Same encoding and same length as QUAL.");
		  put("FI","The index of segment in the template.");
		  put("FS","Segment suffix.");
		  put("FZ","Flow signal intensities on the original strand of the read. ");
		  put("LB","Library. Value to be consistent with the header   RG-LB tag if @RG is present.");
		  put("H0","Number of perfect hits");
		  put("H1","Number of 1-difference hits (see also NM)");
		  put("H2","Number of 2-difference hits ");
		  put("HI","Query hit index, indicating the alignment record is the i-th one stored in SAM");
		  put("IH","Number of stored alignments in SAM that contains the query in the current record");
		  put("MD","String for mismatching positions.");
		  put("MQ","Mapping quality of the mate/next segment ");
		  put("NH","Number of reported alignments that contains the query in the current record");
		  put("NM","Edit distance to the reference, including ambiguous bases but excluding clipping");
		  put("OQ","Original base quality (usually before recalibration).");
		  put("OP","Original mapping position (usually before realignment) ");
		  put("OC","Original CIGAR (usually before realignment) ");
		  put("PG","Program. Value matches the header  PG-ID tag if @PG is present. ");
		  put("PQ","Phred likelihood of the template, conditional on both the mapping being correct ");
		  put("PT","Read annotations for parts of the padded read sequenc");
		  put("PU","Platform unit. Value to be consistent with the header RG-PU tag if @RG is present.");
		  put("QT","Phred quality of the barcode sequence in the  BC or RT tag. Same encoding as  QUAL. ");
		  put("Q2","Phred quality of the mate/next segment sequence in the  R2 tag. Same encoding as QUAL.");
		  put("R2","Sequence of the mate/next segment in the template. ");
		  put("RG","Read group. Value matches the header RG-ID tag if   @RG is present in the header. ");
		  put("RT","Deprecated alternative to  BC tag originally used at Sanger. ");
		  put("SM","Template-independent mapping quality ");
		  put("TC","The number of segments in the template.");
		  put("U2","Phred probility of the 2nd call being wrong conditional on the best being wrong. The same encoding as QUAL. ");
		  put("UQ","Phred likelihood of the segment, conditional on the mapping being correct ");
			}}};
	/* curl -s "https://raw.github.com/samtools/hts-specs/master/SAMv1.tex" | grep -E '\\$' | grep '{\\tt' | grep ' & ' | cut -d '&' -f 1,3 | sed -e  's/{\\tt /if(key.equals("/' -e 's/} \& /")) return "/' | sed 's/\\\\$/";/' */
	private static String getDescription(String key)
		{
		if(key==null) return null;
		 if(key.startsWith("X")) return "Reserved fields for end users (together with Y? and Z?) ";
		 
		 return defs.get(key);
		}
	}

@SuppressWarnings("serial")
class ReadGroupTableModel
	extends  AbstractTableModel
	{
	private List<Object> rows=new Vector<Object>();
	ReadGroupTableModel()
		{
		
		}
	
	private void add(String key,Object value)
		{
		if(value==null) return;
		rows.add(key);
		rows.add(value);
		}
	
	public void setContext(SAMRecord rec)
		{
		rows.clear();
		SAMReadGroupRecord g=(rec==null?null:rec.getReadGroup());
		if(g!=null)
			{
			add("Library",g.getLibrary());
			add("Sample",g.getSample());
			add("KeySequence",g.getKeySequence());
			add("Description",g.getDescription());
			add("FlowOrder",g.getFlowOrder());
			add("Platform",g.getPlatform());
			add("PlatformUnit",g.getPlatformUnit());
			add("PredictedMedianInsertSize",g.getPredictedMedianInsertSize());
			add("ReadGroupId",g.getReadGroupId());
			add("Sample",g.getSample());
			add("SequencingCenter",g.getSequencingCenter());
			add("RunDate",g.getRunDate());
			}
		fireTableDataChanged();
		}
	
	@Override
	public int getRowCount()
		{
		return rows.size()/2;
		}
	
	@Override
	public int getColumnCount() {
		return 2;
		}
	
	@Override
	public String getColumnName(int column) {
		return (column==0?"Key":"Value");
		}
	@Override
	public Class<?> getColumnClass(int columnIndex) {
		return (columnIndex==0?String.class:Object.class);
		}
	
	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		return rows.get(rowIndex*2+columnIndex);
		}
	
	}

class BamTableModel
	extends  AbstractGenericTable<SAMRecord>
	{
	private enum COLS{
		QNAME,FLAG,RNAME,POS,
		MAPQ,CIGAR,RNEXT,PNEXT,
		TLEN,SEQ,QUAL};
	private static final long serialVersionUID = 1L;
	//private BamFileRef ref;
	
	
	BamTableModel()
		{
		}
	
	
	@Override
	public Object getValueOf(SAMRecord o, int columnIndex)
		{
		if(o==null) return null;
		switch(COLS.values()[columnIndex])
			{
			case QNAME: return o.getReadName();
			case FLAG: return o.getFlags();
			case RNAME: return (o.getReadUnmappedFlag()?null:o.getReferenceName());
			case POS: return (o.getReadUnmappedFlag()?null:o.getAlignmentStart());
			case MAPQ: return (o.getReadUnmappedFlag()?null:o.getMappingQuality());
			case CIGAR: return (o.getReadUnmappedFlag()?null:o.getCigarString());
			case RNEXT: return (!o.getReadPairedFlag() || o.getMateUnmappedFlag()?null:o.getMateReferenceName());
			case PNEXT: return (!o.getReadPairedFlag() || o.getMateUnmappedFlag()?null:o.getMateAlignmentStart());
			case TLEN:
				{
				if(!o.getReadPairedFlag()) return null;
				if(o.getReadUnmappedFlag()) return null;
				if(o.getMateUnmappedFlag()) return null;
				if(o.getReferenceIndex()!=o.getMateReferenceIndex()) return null;
				return o.getInferredInsertSize();
				}
			case SEQ: return o.getReadString();
			case QUAL: return o.getBaseQualityString();
			}
		return null;
		}
	
	@Override
	public String getColumnName(int column)
		{
		return COLS.values()[column].name();
		}
	
	@Override
	public int getColumnCount()
		{
		return COLS.values().length;
		}
	
	synchronized void updateRow(List<SAMRecord> L)
		{
		super.rows.clear();
		super.rows.addAll(L);
		this.fireTableDataChanged();
		}
	}


class WorkbenchFrame extends JDialog
	{
	private static Log LOG=Log.getInstance(VcfViewGui.class);

	private static final long serialVersionUID = 1L;
	private JDesktopPane desktopPane;
	private List<AbstractRefInternalFrame> allInternalFrames=new Vector<AbstractRefInternalFrame>();
	private ReloadVCF dataLoaderVCFThread;
	private ReloadBam dataLoaderBamThread;
	private SpinnerNumberModel numFetchModel;
	private SpinnerNumberModel numSecondsModel;
	private JTextField selectRgnField;
	private List<JCheckBox> requiredFlags=new Vector<JCheckBox>();
	private List<JCheckBox> filteringFlags=new Vector<JCheckBox>();

	
	//private JTextField jexlField;
	class ReloadBam extends Thread
		{
		List<SAMRecord> rows=new Vector<SAMRecord>();
		int maxRows;
		int maxTimeSeconds=10;
		Interval reg=null;
		private BamInternalFrame bamintf;
		private Set<SamFlag> g_flag_off=new HashSet<SamFlag>();
		private int g_flag_on=0;

		
		private void updateRows() throws InterruptedException,InvocationTargetException
			{
			if(dataLoaderBamThread!=this) return;
			LOG.info("updating "+rows.size());
			SwingUtilities.invokeAndWait(new Runnable()
				{
				@Override
				public void run() {
					bamintf.tableModel.setRows(rows);
					}
				});
			
			}
		
		
		@Override
		public void run()
			{
			SamReader samReader=null;
			SAMRecordIterator iter=null;
			try
				{
				SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
				
				for(int i=0;i< allInternalFrames.size() && dataLoaderBamThread==this;++i)
					{
					if(!allInternalFrames.get(i).getRefFile().isBAM()) continue;
					this.bamintf = (BamInternalFrame)allInternalFrames.get(i);
					this.rows.clear();
					try
						{
						LOG.info("reading "+bamintf.ref.filePath.getPath()+" region \""+reg+"\"");
						samReader=srf.open(bamintf.ref.filePath);
						if(bamintf.ref.isIndexed() && reg!=null)
							{
							iter=samReader.queryOverlapping(
									reg.getSequence(),
									reg.getStart(),
									reg.getEnd()
									);
							}
						else
							{
							iter=samReader.iterator();
							}
						long start=System.currentTimeMillis();
						while(   dataLoaderBamThread==this &&
								iter!=null && iter.hasNext() &&
								rows.size()<maxRows &&
								(System.currentTimeMillis()-start)< maxTimeSeconds*1000
								)
							{
							SAMRecord rec=iter.next();
							for(SamFlag f:g_flag_off)
								{
								if(f.isSet(rec.getFlags()))
									{
									rec=null;break;
									}
								}
							
							if(rec==null) continue;
							if( (g_flag_on & rec.getFlags())!=g_flag_on) continue;//OK checked
							rows.add(rec);
							}
						}
					catch(Exception err)
						{
						err.printStackTrace();
						this.rows.clear();
						}
					finally
						{
						CloserUtil.close(iter);
						CloserUtil.close(samReader);
						samReader=null;
						iter=null;
						}
					
					updateRows();
					}
					
				}
			catch(Exception err)
				{
				err.printStackTrace();
				}
			finally
				{
				CloserUtil.close(samReader);
				rows.clear();
				}
			}
		}

	
	
	//private JTextField jexlField;
	class ReloadVCF extends Thread
		{
		List<String> rows=new Vector<String>();
		int maxRows;
		VCFPos reg;
		private VCFInternalFrame vcfi;
		
		private void updateRows() throws InterruptedException,InvocationTargetException
			{
			if(dataLoaderVCFThread!=this) return;
			LOG.info("updating "+rows.size());
			SwingUtilities.invokeAndWait(new Runnable()
				{
				@Override
				public void run() {
					vcfi.tableModel.updateRow(rows);
					}
				});
			
			}
		
		
		@Override
		public void run()
			{
			TabixReader tabixReader=null;
			BufferedReader br=null;
			try
				{
				
				for(int i=0;i< allInternalFrames.size() &&
						dataLoaderVCFThread==this;++i)
					{
					if(!allInternalFrames.get(i).getRefFile().isVCF()) continue;
					this.vcfi = (VCFInternalFrame)allInternalFrames.get(i);
					String line;
					this.rows.clear();
					try
						{
						LOG.info("reading "+vcfi.getVCFFileRef().filePath.getPath()+" region:"+reg);
						tabixReader=new TabixReader(vcfi.getVCFFileRef().filePath.getPath());
						if(!vcfi.getVCFFileRef().isIndexed())
							{
							br = IOUtils.openFileForBufferedReading(vcfi.getVCFFileRef().filePath);
							while(   dataLoaderVCFThread==this &&
									((line=br.readLine())!=null) &&
									rows.size()<maxRows)
								{
								if(line.startsWith(VCFHeader.METADATA_INDICATOR)) continue;
								if(line.startsWith(VCFHeader.HEADER_INDICATOR)) continue;
								rows.add(line);
								}
							CloserUtil.close(br);
							br=null;
							}
						else if(reg!=null)
							{
							TabixReader.Iterator iter=tabixReader.query(reg.toString());
							while(   dataLoaderVCFThread==this &&
									iter!=null && (line=iter.next())!=null &&
									rows.size()<maxRows)
								{
								rows.add(line);
								}
							tabixReader=null;
							}
						else
							{

							while(dataLoaderVCFThread==this &&
									(line=tabixReader.readLine())!=null &&
									rows.size()<maxRows)
								{
								if(line.startsWith(VCFHeader.METADATA_INDICATOR)) continue;
								if(line.startsWith(VCFHeader.HEADER_INDICATOR)) continue;
								
								this.rows.add(line);
								}
							}
						}
					catch(Exception err)
						{
						err.printStackTrace();
						this.rows.clear();
						}
					finally
						{
						if(tabixReader!=null) tabixReader.close();
						tabixReader=null;
						CloserUtil.close(br);
						}
					
					updateRows();
					}
					
				}
			catch(Exception err)
				{
				err.printStackTrace();
				}
			finally
				{
				if(tabixReader!=null) tabixReader.close();
				rows.clear();
				}
			}
		}
	
	
	@SuppressWarnings("serial")
	WorkbenchFrame(final List<AbstractFileRef> allFileRefs)
		{
		super((JFrame)null,"WorkbenchFrame",ModalityType.APPLICATION_MODAL);
		this.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		
		
		 
		
		addWindowListener(new WindowAdapter()
			{
			@Override
			public void windowOpened(WindowEvent e)
				{
				removeWindowListener(this);
				for(AbstractFileRef vfr:allFileRefs)
					{
					addRefFile(
							Toolkit.getDefaultToolkit().getScreenSize(),
							vfr);
					}
				reloadFrameContent();
				}
			});
		
		addWindowListener(new WindowAdapter()
			{
			@Override
			public void windowClosing(WindowEvent e)
				{
				doMenuClose();
				}
			});
		JMenuBar bar=new JMenuBar();
		setJMenuBar(bar);
		
		JPanel contentPane=new JPanel(new BorderLayout(5,5));
		setContentPane(contentPane);
		this.desktopPane=new JDesktopPane();
		contentPane.add(this.desktopPane, BorderLayout.CENTER);
		
		JPanel top=new JPanel(new FlowLayout(FlowLayout.LEADING));
		contentPane.add(top, BorderLayout.NORTH);
		
		JLabel lbl=new JLabel("Max Rows:",JLabel.LEADING);
		JSpinner spinner=new JSpinner(this.numFetchModel=new SpinnerNumberModel(1000, 1, 1000000, 10));
		lbl.setLabelFor(spinner);
		top.add(lbl);
		top.add(spinner);
		
		lbl=new JLabel("Timeout (secs):",JLabel.LEADING);
		spinner=new JSpinner(this.numSecondsModel=new SpinnerNumberModel(2, 1, 10000, 1));
		lbl.setLabelFor(spinner);
		top.add(lbl);
		top.add(spinner);

		
		
		lbl=new JLabel("JEXL:",JLabel.LEADING);
		/*jexlField=new JTextField(20);
		lbl.setLabelFor(jexlField);
		
		top.add(lbl);
		top.add(jexlField);*/
		
		lbl=new JLabel("Region:",JLabel.LEADING);
		selectRgnField=new JTextField(20);
		lbl.setLabelFor(selectRgnField);
		AbstractAction action=new AbstractAction("Select")
			{
			private static final long serialVersionUID = 1L;

			@Override
			public void actionPerformed(ActionEvent a)
				{
				if(
					/* (jexlField.getText().trim().isEmpty() || parseJex(jexlField.getText().trim())!=null) && */
				   (selectRgnField.getText().trim().isEmpty() ||
				   parseOne(selectRgnField.getText())!=null))
					{
					reloadFrameContent();
					}
				else
					{
					LOG.info("Bad input "+selectRgnField.getText());
					}
				}
			};
		selectRgnField.addActionListener(action);
		//jexlField.addActionListener(action);

		
		top.add(lbl);
		top.add(selectRgnField);
		top.add(new JButton(action));
		
		JMenu menu=new JMenu("File");
		getJMenuBar().add(menu);
		action=new AbstractAction("Load VCF")
			{
			
			@Override
			public void actionPerformed(ActionEvent e) {
				doMenuLoadVCF();
				}
			};
		//top.add(new JButton(""));
		menu.add(action);
		
		action=new AbstractAction("Quit")
			{
			
			@Override
			public void actionPerformed(ActionEvent e) {
				setVisible(false);
				dispose();
				}
			};
		menu.add(action);
		
		
		
		menu=new JMenu("Flags");
		bar.add(menu);
			for(SamFlag flag:SamFlag.values())
			{
			JCheckBox cbox=new JCheckBox("Require "+flag);
			requiredFlags.add(cbox);
			menu.add(cbox);
			}
		menu.add(new JSeparator());
		for(SamFlag flag:SamFlag.values())
			{
			JCheckBox cbox=new JCheckBox("Filter out "+flag);
			filteringFlags.add(cbox);
			menu.add(cbox);
			}
		}
	private void addRefFile(final Dimension d,AbstractFileRef ref)
		{
		LOG.info("Reading "+ref.filePath);
		int w=(int)(d.width*0.8);
		int h=(int)(d.height*0.8);
		AbstractRefInternalFrame iFrame;
		if(ref.isVCF())
			{
			VCFInternalFrame vif= new VCFInternalFrame((VCFFileRef) ref);
			iFrame=vif;
			vif.jTable.addMouseListener(new MouseAdapter(){
				@Override
				public void mouseClicked(MouseEvent e) {
				      if (e.getClickCount() == 2) {
				         JTable t = (JTable)e.getSource();
				         int row = t.getSelectedRow();
				         if(row==-1) return;
				         VCFTableModel tm=(VCFTableModel)t.getModel();
				         showIgv(tm.getValueAt(row, 0),tm.getValueAt(row, 1));
						}
					}
			});
			
			}
		else if(ref.isBAM())
			{
			BamInternalFrame bif=new BamInternalFrame((BamFileRef) ref);
			iFrame=bif;
			bif.jTable.addMouseListener(new MouseAdapter(){
				@Override
				public void mouseClicked(MouseEvent e) {
				      if (e.getClickCount() == 2) {
				         JTable t = (JTable)e.getSource();
				         int row = t.getSelectedRow();
				         if(row==-1) return;
				         BamTableModel tm=(BamTableModel)t.getModel();
				         showIgv(tm.getValueAt(row, 0),tm.getValueAt(row, 1));
						}
					}
				});
			}
		else
			{
			return;
			}
		this.allInternalFrames.add(iFrame);
		iFrame.setBounds(
				Math.max((int)((d.width-w)*Math.random()),0),
				Math.max((int)((d.height-h)*Math.random()),0),
				w, h);
		desktopPane.add(iFrame);
		iFrame.setVisible(true);
		}
	
	private void doMenuLoadVCF()
		{
		Dimension d=this.desktopPane.getSize();
		File vcf[]=WorkbenchFrame.selectNgsFiles(this);
		if(vcf==null || vcf.length==0) return;
		try
			{
			for(File f:vcf)
				{
				AbstractFileRef ref=create(f);
				addRefFile(d,ref);
				}
			}
		catch(Exception err)
			{
			err.printStackTrace();
			JOptionPane.showMessageDialog(this, "Error:"+err.getMessage());
			}
		reloadFrameContent();
		}
	/** create a new VCF file ref */
    static VCFFileRef createVCFRef(File vcfFile) throws IOException
		{
    	IOUtil.assertFileIsReadable(vcfFile);
		VCFFileRef vfr=new VCFFileRef();
		vfr.filePath=vcfFile;
		LineIterator r=IOUtils.openFileForLineIterator(vcfFile);
		vfr.codhead=VCFUtils.parseHeader(r);
		CloserUtil.close(r);
		if(VCFUtils.isValidTabixFile(vcfFile))
			{
			vfr.fileisIndexed=true;
			}
		return vfr;
		}

    static BamFileRef createBamRef(File bamFile) throws IOException
		{
    	IOUtil.assertFileIsReadable(bamFile);
    	BamFileRef bamref=new BamFileRef();
    	bamref.filePath=bamFile;
    	SamReader r=null;
    	try
    		{
    		SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
    		r=srf.open(bamFile);
    		bamref.header=r.getFileHeader();
    		bamref.fileisIndexed=r.hasIndex();
    		}
    	catch(Exception err)
    		{
    		bamref.header=new SAMFileHeader();
    		}
    	finally
    		{	
    		CloserUtil.close(r);
    		}
    	return bamref;
    	}
	
    static AbstractFileRef create(File vcfFile) throws IOException
		{
		if(vcfFile.getName().endsWith(".vcf.gz") || vcfFile.getName().endsWith(".vcf"))
			{
			return createVCFRef(vcfFile);
			}
		else if(vcfFile.getName().endsWith(".bam"))
			{
			return createBamRef(vcfFile);
			}
		else
			{
			throw new IOException("Unkown File Type");
			}
		}

    
	
	static File[] selectNgsFiles(Component owner)
		{
		Preferences prefs=Preferences.userNodeForPackage(VcfViewGui.class);
		String dirStr=prefs.get("last.directory", null);
		File lastDir=null;
		if(dirStr!=null) lastDir=new File(dirStr);
		JFileChooser chooser=new JFileChooser(lastDir);
		chooser.setFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "VCF indexed with tabix, Indexed BAM";
			}
			
			@Override
			public boolean accept(File f) {
				if(f.isDirectory()) return true;
				if(VCFUtils.isTabixVcfFile(f)) return true;
				
				String name=f.getName();
				if(!name.endsWith(".bam")) return false;
				if(new File(f.getParentFile(),name+".bai").exists()) return true;
				return new File(
						f.getParentFile(),
						name.substring(0,name.length()-4)+".bai"
						).exists();
				};
			});
		chooser.setMultiSelectionEnabled(true);
		if(chooser.showOpenDialog(owner)!=JFileChooser.APPROVE_OPTION)
			{
			LOG.info("user pressed cancel");
			return new File[0];
			}
		File fs[]=chooser.getSelectedFiles();
		if(fs!=null)
			{
			if(fs.length>0 && fs[0].getParentFile()!=null)
				{
				prefs.put("last.directory", fs[0].getParentFile().getPath());
				try { prefs.sync();}catch(BackingStoreException err){}
				}
			return fs;
			}
		return new File[0];
		}
	
	private synchronized void reloadFrameContent()
		{
		reloadVCFFrameContent();
		reloadBamFrameContent();
		}
	
	private void reloadVCFFrameContent()
		{
		if(dataLoaderVCFThread!=null)
			{
			try { dataLoaderVCFThread.interrupt();}
			catch(Exception err) {}
			}
		dataLoaderVCFThread=new ReloadVCF();
		dataLoaderVCFThread.maxRows=numFetchModel.getNumber().intValue();
		dataLoaderVCFThread.reg=parseOne(this.selectRgnField.getText());
		dataLoaderVCFThread.start();
		}
	
	private synchronized void reloadBamFrameContent()
		{
		if(dataLoaderBamThread!=null)
			{
			try { dataLoaderBamThread.interrupt();}
			catch(Exception err) {}
			}
		dataLoaderBamThread=new ReloadBam();
		
		for(SamFlag flag:SamFlag.values())
			{
			if(this.requiredFlags.get(flag.ordinal()).isSelected())
				{
				dataLoaderBamThread.g_flag_on |=flag.getFlag();
				}
			if(this.filteringFlags.get(flag.ordinal()).isSelected())
				{
				dataLoaderBamThread.g_flag_off.add(flag);
				}
			}
		
		dataLoaderBamThread.maxRows=numFetchModel.getNumber().intValue();
		dataLoaderBamThread.maxTimeSeconds=numSecondsModel.getNumber().intValue();
		VCFPos pos=parseOne(this.selectRgnField.getText());
		dataLoaderBamThread.reg=null;
		if(pos!=null)
			{
			dataLoaderBamThread.reg=new Interval(pos.chrom, pos.start, pos.start+1000);
			}
		dataLoaderBamThread.start();
		}

	private void doMenuClose()
		{
		this.setVisible(false);
		this.dispose();
		}
	
	IgvSocket igvSocket=null;
	
	//http://plindenbaum.blogspot.fr/2011/07/controlling-igv-through-port-my.html
	private void showIgv(final Object chrom,final Object pos)
		{
		if(igvSocket==null || igvSocket.getPort()<0) return;
		Thread thread=new Thread()
			{	
			@Override
			public void run()
				{
				PrintWriter out=null;
				try
					{
					out = igvSocket.getWriter();
					igvSocket.getReader();
					out.println("goto "+chrom+":"+pos);
					try{ Thread.sleep(5*1000);}
					catch(InterruptedException err2) {}
					 }
				catch(Exception err)
					{
					LOG.error(err);
					}
				finally
					{
					CloserUtil.close(igvSocket);
					}
				}
			};
		thread.start();
		}
	
	public static VCFPos parseOne(String s)
		{
		s=s.trim();
		int colon=s.indexOf(':');
		if(colon<1 || colon+1==s.length()) return null;
		int hyphen=s.indexOf('-',colon+1);
		try {
			VCFPos pos=new VCFPos();
			pos.chrom=s.substring(0,colon);
			
			
			if(hyphen==-1)
				{
				pos.start=Integer.parseInt(s.substring(colon+1));
				}
			else
				{
				pos.start=Integer.parseInt(s.substring(colon+1),hyphen);
				pos.end=Integer.parseInt(s.substring(hyphen+1));
				if(pos.start<1 || pos.end<pos.start) return null;
				}
			return pos;
		} catch (Exception e)
			{
			return null;
			}
		
		}
	}

/**
 * 
 * VcfViewGui
 *
 */
public class VcfViewGui
	extends AbstractCommandLineProgram
	{

    
    @Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"VcfViewGui";
    	}
    
    @Override
	public String getProgramDescription() {
		return "Simple java-Swing-based VCF+Bam viewer.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -O IGV Host. default:"+IgvSocket.DEFAULT_HOST);
		out.println(" -P IGV Port. default:"+IgvSocket.DEFAULT_PORT+" negative number to disable.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		@SuppressWarnings("resource")
		IgvSocket igvSocket=new IgvSocket();
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"O:P:"))!=-1)
			{
			switch(c)
				{
				case 'O': igvSocket.setHost(opt.getOptArg());break;
				case 'P': igvSocket.setPort(Integer.parseInt(opt.getOptArg()));break;
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
		
		
		
		try
			{
			List<File> IN=new ArrayList<File>();
			JFrame.setDefaultLookAndFeelDecorated(true);
			JDialog.setDefaultLookAndFeelDecorated(true);
			List<AbstractFileRef> vfrs=new ArrayList<AbstractFileRef>();
			
			
			for(int i=opt.getOptInd();i< args.length;++i)
				{
				String filename=args[i];
				IN.add(new File(filename));
				}
				
			
			
			if(IN.isEmpty())
				{
				info("NO VCF provided; Opening dialog");
				
				File fs[]=WorkbenchFrame.selectNgsFiles(null);
				if(fs!=null)
					{
					IN.addAll(Arrays.asList(fs));
					}
				
				}
			
			
			for(File in:IN)
				{
				vfrs.add(WorkbenchFrame.create(in));
				}
			
			info("showing Workbench frame");
			Dimension screen=Toolkit.getDefaultToolkit().getScreenSize();
			final WorkbenchFrame f=new WorkbenchFrame(vfrs);
			f.igvSocket=igvSocket;
			f.setBounds(50, 50, screen.width-100, screen.height-100);
			SwingUtilities.invokeAndWait(new Runnable()
				{
				
				@Override
				public void run() {
					f.setVisible(true);
					}
				});
			
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfViewGui().instanceMainWithExit(args);
		}

	}
