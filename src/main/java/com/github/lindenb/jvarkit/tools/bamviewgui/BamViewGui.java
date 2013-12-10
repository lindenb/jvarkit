package com.github.lindenb.jvarkit.tools.bamviewgui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileFilter;

import com.github.lindenb.jvarkit.util.picard.IntervalUtils;
import com.github.lindenb.jvarkit.util.picard.SamFlag;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;

@SuppressWarnings("serial")
public class BamViewGui extends JFrame
	{
	private enum COLS{
		QNAME,FLAG,RNAME,POS,
		MAPQ,CIGAR,RNEXT,PNEXT,
		TLEN,SEQ,QUAL};
	private File bamFile;
	private JTextField selectField=null;
	private SamRecordModel tableModel;
	private List<JCheckBox> requiredFlags=new ArrayList<JCheckBox>();
	private List<JCheckBox> filteringFlags=new ArrayList<JCheckBox>();
	
	private static final int MAX_NUM_SAM_RECORDS=1000;
	
	private class SamRecordModel extends  AbstractGenericTable<SAMRecord>
		{
		
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
		
		@Override
		public Class<?> getColumnClass(int columnIndex)
			{
			switch(COLS.values()[columnIndex])
				{
				case FLAG:
				case POS: 
				case MAPQ: 
				case PNEXT: 
				case TLEN: return Integer.class;
				
				case CIGAR: 
				case RNEXT: 
				case QNAME:
				case RNAME: 
				case SEQ: //
				case QUAL: return String.class;
				}
			return Object.class;
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
				case TLEN: return o.getInferredInsertSize();
				case SEQ: return o.getReadString();
				case QUAL: return o.getBaseQualityString();
				}
			return null;
			}
		}
	
	
	private BamViewGui(File bamFile)
		{
		super(bamFile.getPath());
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		this.bamFile=bamFile;
		JPanel main=new JPanel(new BorderLayout(5,5));
		setContentPane(main);
		
		JMenuBar menuBar=new JMenuBar();
		setJMenuBar(menuBar);
		
		JPanel pane=new JPanel(new FlowLayout());
		main.add(pane,BorderLayout.NORTH);
		pane.add(new JLabel("select:"));
		this.selectField=new JTextField(20);
		pane.add(this.selectField);
		AbstractAction action=new AbstractAction("Go")
			{
			@Override
			public void actionPerformed(ActionEvent arg0)
				{
				doMenuSelect();
				}
			};
		this.selectField.addActionListener(action);
		pane.add(new JButton(action));
		this.tableModel=new SamRecordModel();
		JTable table=new JTable(this.tableModel);
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		table.setShowVerticalLines(false);
		JScrollPane scroll=new JScrollPane(table);
		main.add(scroll,BorderLayout.CENTER);
		this.addWindowListener(new WindowAdapter()
			{
		
			@Override
			public void windowClosing(WindowEvent e)
				{
				doMenuClose();
				}
			});
		
		this.addWindowListener(new WindowAdapter()
			{
			@Override
			public void windowOpened(WindowEvent e)
				{
				doMenuSelect();
				removeWindowListener(this);
				}
			});
		JMenu menu=new JMenu("File");
		menuBar.add(menu);
		menu.add(new AbstractAction("Quit")
			{
			@Override
			public void actionPerformed(ActionEvent arg0)
				{
				doMenuClose();
				}
			});
		menu=new JMenu("Flags");
		menuBar.add(menu);
		for(SamFlag flag:SamFlag.values())
			{
			JCheckBox cbox=new JCheckBox("Require "+flag);
			requiredFlags.add(cbox);
			menu.add(cbox);
			}
		menu.add(new JSeparator());
		for(SamFlag flag:SamFlag.values())
			{
			JCheckBox cbox=new JCheckBox("Filter "+flag);
			filteringFlags.add(cbox);
			menu.add(cbox);
			}
		}
	
	private void doMenuSelect()
		{
		String text=this.selectField.getText().trim();
		SAMFileReader sfr=null;
		SAMRecordIterator iter=null;
		try
			{
			sfr=new SAMFileReader(this.bamFile);
			if(text.isEmpty())
				{
				iter=sfr.iterator();
				}
			else
				{
				SAMSequenceDictionary dict=sfr.getFileHeader().getSequenceDictionary();
				Interval reg=IntervalUtils.parseOne(dict,text);
				if(reg==null) return;
				SAMSequenceRecord ssr=dict.getSequence(reg.getSequence());
				if(ssr==null) return;
				if(reg.getStart()>ssr.getSequenceLength()) return;
				if(reg.getEnd()>ssr.getSequenceLength()) 
					{
					reg=new Interval(
							reg.getSequence(),
							reg.getStart(),
							ssr.getSequenceLength()
							);
					}
				iter=sfr.queryOverlapping(reg.getSequence(), reg.getStart(), reg.getEnd());
				}
			
			int filteringFlag=0;
			int RequiredFlag=0;
			for(SamFlag flag:SamFlag.values())
				{
				if(this.requiredFlags.get(flag.ordinal()).isSelected())
					{
					filteringFlag |=flag.getFlag();
					}
				if(this.filteringFlags.get(flag.ordinal()).isSelected())
					{
					RequiredFlag |=flag.getFlag();
					}
				}
			
			sfr.setValidationStringency(ValidationStringency.SILENT);
			List<SAMRecord> buffer=new ArrayList<SAMRecord>(MAX_NUM_SAM_RECORDS);
			while(iter.hasNext() && buffer.size()< MAX_NUM_SAM_RECORDS)
				{
				SAMRecord rec=iter.next();
				if((rec.getFlags() & filteringFlag)!=0) continue;
				if((rec.getFlags() & RequiredFlag)!=RequiredFlag) continue;

				buffer.add(rec);
				}
			tableModel.setRows(buffer);
			}
		catch(Exception err)
			{
			err.printStackTrace();
			tableModel.clear();
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			}
		}
	
	private void doMenuClose()
		{
		this.setVisible(false);
		this.dispose();
		}
	
	private static boolean acceptBam(final File f)
		{
		if(f==null) return false;
		String name=f.getName();
		if(!name.endsWith(".bam")) return false;
		if(new File(f.getParentFile(),name+".bai").exists()) return true;
		return new File(f.getParentFile(),name.substring(name.length()-4)+".bai").exists();
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		
		File bam=null;
		if(args.length==1)
			{
			bam=new File(args[0]);
			}
		else
			{
			JFileChooser chooser=new JFileChooser();
			chooser.setFileFilter(new FileFilter()
				{
				@Override
				public String getDescription()
					{
					return "*.bam";
					}
				
				@Override
				public boolean accept(File f)
					{
					if(f.isDirectory()) return true;
					return acceptBam(f);
					}
				});
			if(chooser.showOpenDialog(null)!=JFileChooser.APPROVE_OPTION) return;
			bam=chooser.getSelectedFile();
			}
		
		
		if(!acceptBam(bam))
			{
			System.err.println("Bad bam file (missing index ?)");
			System.exit(-1);
			}
		final BamViewGui frame=new BamViewGui(bam);
		JFrame.setDefaultLookAndFeelDecorated(true);
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
		}

	}
