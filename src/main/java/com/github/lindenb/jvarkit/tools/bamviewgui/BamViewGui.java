package com.github.lindenb.jvarkit.tools.bamviewgui;

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
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.net.Socket;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Logger;

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
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.util.picard.SamFlag;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;


class BamFileRef
	{
	File bamFile;
	SAMFileHeader header;
	}

@SuppressWarnings("serial")
class BamInternalFrame extends JInternalFrame
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
		super(ref.bamFile.getName(),true,false,true,true);
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

@SuppressWarnings("serial")
class BamFrame extends JDialog
	{
	private static Logger LOG=Logger.getLogger("jvarkit");

	private static final long serialVersionUID = 1L;
	private JDesktopPane desktopPane;
	private List<BamFileRef> BamFileRefs;
	private List<BamInternalFrame> BamInternalFrames=new Vector<BamInternalFrame>();
	private Reload dataLoaderThread;
	private SpinnerNumberModel numFetchModel;
	private SpinnerNumberModel numSecondsModel;
	private JTextField selectRgnField;
	private List<JCheckBox> requiredFlags=new Vector<JCheckBox>();
	private List<JCheckBox> filteringFlags=new Vector<JCheckBox>();

	
	//private JTextField jexlField;
	class Reload extends Thread
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
			if(dataLoaderThread!=this) return;
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
			SAMFileReader tabixReader=null;
			SAMRecordIterator iter=null;
			try
				{
				
				for(int i=0;i< BamInternalFrames.size() && dataLoaderThread==this;++i)
					{
					this.bamintf = BamInternalFrames.get(i);
					this.rows.clear();
					try
						{
						LOG.info("reading "+bamintf.ref.bamFile.getPath()+" region \""+reg+"\"");
						tabixReader=new SAMFileReader(bamintf.ref.bamFile);
						if(reg!=null)
							{
							iter=tabixReader.queryOverlapping(
									reg.getSequence(),
									reg.getStart(),
									reg.getEnd()
									);
							}
						else
							{
							iter=tabixReader.iterator();
							}
						long start=System.currentTimeMillis();
						while(   dataLoaderThread==this &&
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
						CloserUtil.close(tabixReader);
						tabixReader=null;
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
				if(tabixReader!=null) tabixReader.close();
				rows.clear();
				}
			}
		}
	
	
	BamFrame(List<BamFileRef> BamFileRefs)
		{
		super((JFrame)null,"Bam View ("+BamFileRefs.size()+" files)",ModalityType.APPLICATION_MODAL);
		this.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		this.BamFileRefs=BamFileRefs;
		
		
		 
		
		addWindowListener(new WindowAdapter()
			{
			
			@Override
			public void windowOpened(WindowEvent e)
				{
				removeWindowListener(this);
				Dimension d=Toolkit.getDefaultToolkit().getScreenSize();
				d.width-=150;
				d.height-=150;
				for(BamFileRef vfr:BamFrame.this.BamFileRefs)
					{
					LOG.info("Reading "+vfr.bamFile);
					int w=(int)(d.width*0.8);
					int h=(int)(d.height*0.8);
					BamInternalFrame iFrame=new BamInternalFrame(vfr);
					iFrame.setBounds(
							Math.max((int)((d.width-w)*Math.random()),0),
							Math.max((int)((d.height-h)*Math.random()),0),
							w, h);
					desktopPane.add(iFrame);
					BamInternalFrames.add(iFrame);
					iFrame.setVisible(true);
					iFrame.jTable.addMouseListener(new MouseAdapter(){
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
		
		JMenu menu=new JMenu("File");
		bar.add(menu);
		menu.add(new AbstractAction("Quit") {
			@Override
			public void actionPerformed(ActionEvent e) {
				doMenuClose();
			}});
		
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
		
		JPanel contentPane=new JPanel(new BorderLayout(5,5));
		setContentPane(contentPane);
		this.desktopPane=new JDesktopPane();
		contentPane.add(this.desktopPane, BorderLayout.CENTER);
		
		JPanel top=new JPanel(new FlowLayout(FlowLayout.LEADING));
		contentPane.add(top, BorderLayout.NORTH);
		
		JLabel lbl=new JLabel("Max Rows:",JLabel.LEADING);
		JSpinner spinner=new JSpinner(this.numFetchModel=new SpinnerNumberModel(100, 1, 10000, 10));
		lbl.setLabelFor(spinner);
		top.add(lbl);
		top.add(spinner);
		
		lbl=new JLabel("Timeout (secs):",JLabel.LEADING);
		spinner=new JSpinner(this.numSecondsModel=new SpinnerNumberModel(2, 1, 10000, 1));
		lbl.setLabelFor(spinner);
		top.add(lbl);
		top.add(spinner);
		
		//lbl=new JLabel("JEXL:",JLabel.LEADING);
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
		
		
		
			
		}
	private synchronized void reloadFrameContent()
		{
		if(dataLoaderThread!=null)
			{
			try { dataLoaderThread.interrupt();}
			catch(Exception err) {}
			}
		dataLoaderThread=new Reload();
		
		for(SamFlag flag:SamFlag.values())
			{
			if(this.requiredFlags.get(flag.ordinal()).isSelected())
				{
				dataLoaderThread.g_flag_on |=flag.getFlag();
				}
			if(this.filteringFlags.get(flag.ordinal()).isSelected())
				{
				dataLoaderThread.g_flag_off.add(flag);
				}
			}
		
		dataLoaderThread.maxRows=numFetchModel.getNumber().intValue();
		dataLoaderThread.maxTimeSeconds=numSecondsModel.getNumber().intValue();
		dataLoaderThread.reg=parseOne(this.selectRgnField.getText());
		dataLoaderThread.start();
		}

	private void doMenuClose()
		{
		this.setVisible(false);
		this.dispose();
		}
	
	String igvIP=null;//"127.0.0.1"
	Integer igvPort=null;
	
	//http://plindenbaum.blogspot.fr/2011/07/controlling-igv-through-port-my.html
	private void showIgv(final Object chrom,final Object pos)
		{
		if(igvIP==null || igvPort==null) return;
		Thread thread=new Thread()
			{	
			@Override
			public void run()
				{
				PrintWriter out=null;
				BufferedReader in=null;
				Socket socket=null;

				try
					{
					socket = new Socket(igvIP, igvPort);
					out = new PrintWriter(socket.getOutputStream(), true);
					in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
					out.println("goto "+chrom+":"+pos);
					try{ Thread.sleep(5*1000);}
					catch(InterruptedException err2) {}
					 }
				catch(Exception err)
					{
					LOG.info(""+err.getMessage());
					}
				finally
					{
					if(in!=null) try {in.close();} catch(Exception err){}
					if(out!=null) try {out.close();} catch(Exception err){}
					if(socket!=null) try {socket.close();} catch(Exception err){}
					}
				}
			};
		thread.start();
		}
	
	private static Interval parseOne(String s)
		{
		s=s.trim();
		if(s.isEmpty()) return null;
		int colon=s.indexOf(':');
		if(colon==0) return null;
		if(colon==-1) return new Interval(s, 1, Integer.MAX_VALUE-10);
		String chrom=s.substring(0,colon);
		int hyphen=s.indexOf('-',colon+1);
		try {
			int start;
			int end;
			
			if(hyphen==-1)
				{
				start=Integer.parseInt(s.substring(colon+1));
				end= Integer.MAX_VALUE-10;
				}
			else
				{
				start=Integer.parseInt(s.substring(colon+1),hyphen);
				end=Integer.parseInt(s.substring(hyphen+1));
				}
			if(start<1 || end<start) return null;
			return new Interval(chrom,start,end);
		} catch (Exception e)
			{
			System.err.println(String.valueOf(e.getMessage()));
			return null;
			}
		
		}
	}

/**
 * 
 * BamViewGui
 *
 */
@Deprecated
public class BamViewGui
	extends AbstractCommandLineProgram
	{

    
    @Override
    public String getProgramDescription() {
    	return "Simple java-Swing-based BAM viewer.";
    	}
    @Override
    protected String getOnlineDocUrl() {
    	return "https://github.com/lindenb/jvarkit/wiki/BamViewGui";
    	}
    
    private BamFileRef create(File bamFile) throws IOException
    	{
    	BamFileRef vfr=new BamFileRef();
    	vfr.bamFile=bamFile;
    	SAMFileReader r=null;
    	try
    		{
    		r=new SAMFileReader(bamFile);
    		vfr.header=r.getFileHeader();
    		}
    	catch(Exception err)
    		{
    		vfr.header=new SAMFileHeader();
    		}
    	finally
    		{	
    		CloserUtil.close(r);
    		}
    	return vfr;
    	}
    
    private static boolean acceptBam(final File f)
		{
		if(f==null) return false;
		String name=f.getName();
		if(!name.endsWith(".bam")) return false;
		if(new File(f.getParentFile(),name+".bai").exists()) return true;
		return new File(
				f.getParentFile(),
				name.substring(0,name.length()-4)+".bai"
				).exists();
		}
    
    @Override
    public void printOptions(PrintStream out) {
    	out.println(" -H (host) IGV host example: '127.0.0.1' .Optional.");
    	out.println(" -P (port:integer) IGV port example: '60151' .Optional.");
    	super.printOptions(out);
    	}
    
    
	@Override
	public int doWork(String args[])
		{
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		JFrame.setDefaultLookAndFeelDecorated(true);
		JDialog.setDefaultLookAndFeelDecorated(true);
		List<BamFileRef> bams=new ArrayList<BamFileRef>();
	     String IGV_HOST=null;
	     Integer IGV_PORT=null;

		

		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"H:P:"))!=-1)
			{
			switch(c)
				{
				case 'H': IGV_HOST=opt.getOptArg();break;
				case 'P': IGV_PORT=Integer.parseInt(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		List<File> IN=new ArrayList<File>();
		
		if(opt.getOptInd()==args.length)
			{
			info("NO BAM provided; Opening dialog");
			JFileChooser chooser=new JFileChooser();
			chooser.setFileFilter(new FileFilter() {
				@Override
				public String getDescription() {
					return "Indexed BAM files.";
					}
				
				@Override
				public boolean accept(File f)
					{
					if(f.isDirectory()) return true;
					return acceptBam(f);
					};
				});
			chooser.setMultiSelectionEnabled(true);
			if(chooser.showOpenDialog(null)!=JFileChooser.APPROVE_OPTION)
				{
				info("user pressed cancel");
				return -1;
				}
			File fs[]=chooser.getSelectedFiles();
			if(fs!=null) IN.addAll(Arrays.asList(chooser.getSelectedFiles()));
			}
		else
			{
			for(int i=opt.getOptInd();i< args.length;++i)
				{
				File filename=new File(args[i]);
				if(!acceptBam(filename))
					{
					error("Cannot use "+filename+" as input Bam. bad extenstion ? index missing ?");
					return -1;
					}
				IN.add(filename);
				}
			}

		
		for(File in:IN)
			{
			try
				{
				bams.add(create(in));
				}
			catch(Exception err)
				{
				error(err);
				return -1;
				}
			}
		if(bams.isEmpty())
			{
			error("No Bam file");
			return -1;
			}
		info("showing BAM frame");
		final BamFrame frame=new BamFrame(bams);
		frame.igvIP=IGV_HOST;
		frame.igvPort=IGV_PORT;
		try
			{
			SwingUtilities.invokeAndWait(new Runnable()
				{
				@Override
				public void run()
					{
					Dimension screen=Toolkit.getDefaultToolkit().getScreenSize();
					frame.setBounds(50, 50, screen.width-100, screen.height-100);
					frame.setVisible(true);
					}
				});
			}
		catch(Exception err)
			{
			err.printStackTrace();
			System.exit(-1);
			}
		
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BamViewGui().instanceMainWithExit(args);
		}

	}
