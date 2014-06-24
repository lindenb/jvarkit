package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Vector;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;
import java.util.regex.Pattern;

import javax.swing.AbstractAction;
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
import javax.swing.JScrollPane;
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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import com.github.lindenb.jvarkit.util.igv.IgvSocket;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

/**
 * VCF header, codec and header
 */
class VCFFileRef
	{
	File vcfFile;
	VCFUtils.CodecAndHeader codhead;
	}

/**
 * Main FRAME
 */
class VCFInternalFrame extends JInternalFrame
	{
	private static final long serialVersionUID = 1L;

	JTable jTable;
	VCFTableModel tableModel;
	VCFFileRef ref;
	InfoTreeModel infoTreeModel;
	private JTree infoTree;
	GenotypeTableModel genotypeTableModel;
	private ListSelectionListener selList;
	VCFInternalFrame(VCFFileRef ref)
		{
		super(ref.vcfFile.getName(),true,false,true,true);
		this.ref=ref;
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
					infoTreeModel.setContext(null,VCFInternalFrame.this.ref.codhead.header);
					genotypeTableModel.setContext(null);
					}
				else
					{
					infoTreeModel.setContext(ctx,VCFInternalFrame.this.ref.codhead.header);
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
								System.err.println(""+i+")  "+columns.get(i)+"\t"+v2+"  ## "+tokens[i]);
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
		System.err.println("EFF "+description);
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

class VCFFrame extends JDialog
	{
	private static Log LOG=Log.getInstance(VcfViewGui.class);

	private static final long serialVersionUID = 1L;
	private JDesktopPane desktopPane;
	private List<VCFFileRef> vcfFileRefs;
	private List<VCFInternalFrame> vcfInternalFrames=new Vector<VCFInternalFrame>();
	private Reload dataLoaderThread;
	private SpinnerNumberModel numFetchModel;
	private JTextField selectRgnField;
	//private JTextField jexlField;
	class Reload extends Thread
		{
		List<String> rows=new Vector<String>();
		int maxRows;
		VCFPos reg;
		private VCFInternalFrame vcfi;
		
		private void updateRows() throws InterruptedException,InvocationTargetException
			{
			if(dataLoaderThread!=this) return;
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
			
			try
				{
				
				for(int i=0;i< vcfInternalFrames.size() && dataLoaderThread==this;++i)
					{
					this.vcfi = vcfInternalFrames.get(i);
					String line;
					this.rows.clear();
					try
						{
						LOG.info("reading "+vcfi.ref.vcfFile.getPath()+" region:"+reg);
						tabixReader=new TabixReader(vcfi.ref.vcfFile.getPath());
						if(reg!=null)
							{
							
							TabixReader.Iterator iter=tabixReader.query(reg.toString());
							while(   dataLoaderThread==this &&
									iter!=null && (line=iter.next())!=null &&
									rows.size()<maxRows)
								{
								rows.add(line);
								}
							tabixReader=null;
							}
						else
							{
							LOG.info(""+(dataLoaderThread==this));

							while(dataLoaderThread==this &&
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
	VCFFrame(List<VCFFileRef> vcfFileRefs)
		{
		super((JFrame)null,"VCF View ("+vcfFileRefs.size()+" files)",ModalityType.APPLICATION_MODAL);
		this.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		this.vcfFileRefs=vcfFileRefs;
		
		
		 
		
		addWindowListener(new WindowAdapter()
			{
			@Override
			public void windowOpened(WindowEvent e)
				{
				removeWindowListener(this);
				for(VCFFileRef vfr:VCFFrame.this.vcfFileRefs)
					{
					addVCFFile(
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
		}
	
	private void addVCFFile(Dimension d,VCFFileRef vfr)
		{
		LOG.info("Reading "+vfr.vcfFile);
		int w=(int)(d.width*0.8);
		int h=(int)(d.height*0.8);
		VCFInternalFrame iFrame=new VCFInternalFrame(vfr);
		iFrame.setBounds(
				Math.max((int)((d.width-w)*Math.random()),0),
				Math.max((int)((d.height-h)*Math.random()),0),
				w, h);
		desktopPane.add(iFrame);
		vcfInternalFrames.add(iFrame);
		iFrame.setVisible(true);
		iFrame.jTable.addMouseListener(new MouseAdapter(){
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
		LOG.info(iFrame.getBounds());
		}
	
	private void doMenuLoadVCF()
		{
		Dimension d=this.desktopPane.getSize();
		File vcf[]=VCFFrame.selectVcfFiles(this);
		if(vcf==null || vcf.length==0) return;
		try
			{
			for(File f:vcf)
				{
				VCFFileRef ref=create(f);
				this.vcfFileRefs.add(ref);
				addVCFFile(d,ref);
				}
			}
		catch(Exception err)
			{
			err.printStackTrace();
			JOptionPane.showMessageDialog(this, "Error:"+err.getMessage());
			}
		reloadFrameContent();
		}
	
    static VCFFileRef create(File vcfFile) throws IOException
		{
		VCFFileRef vfr=new VCFFileRef();
		vfr.vcfFile=vcfFile;
		LineIterator r=IOUtils.openFileForLineIterator(vcfFile);
		vfr.codhead=VCFUtils.parseHeader(r);
		CloserUtil.close(r);
		return vfr;
		}

	
	
	static File[] selectVcfFiles(Component owner)
		{
		Preferences prefs=Preferences.userNodeForPackage(VcfViewGui.class);
		String dirStr=prefs.get("last.directory", null);
		File lastDir=null;
		if(dirStr!=null) lastDir=new File(dirStr);
		JFileChooser chooser=new JFileChooser(lastDir);
		chooser.setFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "VCF indexed with tabix";
			}
			
			@Override
			public boolean accept(File f) {
				if(f.isDirectory()) return true;
				return VCFUtils.isTabixVcfFile(f);
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
		if(dataLoaderThread!=null)
			{
			try { dataLoaderThread.interrupt();}
			catch(Exception err) {}
			}
		dataLoaderThread=new Reload();
		dataLoaderThread.maxRows=numFetchModel.getNumber().intValue();
		dataLoaderThread.reg=parseOne(this.selectRgnField.getText());
		dataLoaderThread.start();
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
    	return "https://github.com/lindenb/jvarkit/wiki/VcfViewGui";
    	}
    
    @Override
	public String getProgramDescription() {
		return "Simple java-Swing-based VCF viewer.";
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
			List<VCFFileRef> vfrs=new ArrayList<VCFFileRef>();
			
			
			for(int i=opt.getOptInd();i< args.length;++i)
				{
				String filename=args[i];
				IN.add(new File(filename));
				}
				
			
			
			if(IN.isEmpty())
				{
				info("NO VCF provided; Opening dialog");
				
				File fs[]=VCFFrame.selectVcfFiles(null);
				if(fs!=null)
					{
					IN.addAll(Arrays.asList(fs));
					}
				
				}
			
			
			for(File in:IN)
				{
				vfrs.add(VCFFrame.create(in));
				}
			
			info("showing VCF frame");
			Dimension screen=Toolkit.getDefaultToolkit().getScreenSize();
			final VCFFrame f=new VCFFrame(vfrs);
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
