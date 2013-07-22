package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.net.Socket;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;
import java.util.Vector;
import java.util.regex.Pattern;

import javax.crypto.spec.IvParameterSpec;
import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JDesktopPane;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.AbstractTableModel;


import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.TabixReader;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.picard.IOUtils;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

class VCFFileRef
	{
	File vcfFile;
	VCFCodec codec=new VCFCodec();
	VCFHeader header;
	}

class VCFInternalFrame extends JInternalFrame
	{
	private static final long serialVersionUID = 1L;

	private JTable jTable;
	VCFTableModel tableModel;
	VCFFileRef ref;
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
		JScrollPane scroll=new JScrollPane(this.jTable);
		pane.add(scroll);
		
		pane=new JPanel(new BorderLayout(5,5));
		
		ByteArrayOutputStream baos=new ByteArrayOutputStream();
		VariantContextWriter w=VariantContextWriterFactory.create(baos, null, EnumSet.noneOf(Options.class));
		w.writeHeader(ref.header);
		
		tabbedPane.addTab("Header", pane);
		JTextArea area=new JTextArea(new String(baos.toByteArray()));
		area.setCaretPosition(0);
		area.setEditable(false);
		pane.add(new JScrollPane(area),BorderLayout.CENTER);
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
		return ref.header.getSampleNamesInOrder().get(column);
		}
	
	@Override
	public int getColumnCount()
		{
		int n=VCFHeader.HEADER_FIELDS.values().length;
		if(this.ref.header.getSampleNamesInOrder().isEmpty())
			{
			return n;
			}
		return n+1+ref.header.getSampleNamesInOrder().size();
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
				Dimension d=Toolkit.getDefaultToolkit().getScreenSize();
				d.width-=150;
				d.height-=150;
				LOG.info(d);
				for(VCFFileRef vfr:VCFFrame.this.vcfFileRefs)
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
					LOG.info(iFrame.getBounds());
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
	
	String igvIP=null;//"127.0.0.1"
	Integer igvPort=null;
	
	//http://plindenbaum.blogspot.fr/2011/07/controlling-igv-through-port-my.html
	private void showIgv(String chrom,int pos)
		{
		if(igvIP==null || igvPort==null) return;
		Socket socket=null;
		PrintWriter out=null;
		BufferedReader in=null;
		try
			{
			socket = new Socket(igvIP, igvPort);
			out = new PrintWriter(socket.getOutputStream(), true);
			in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

			out.println("goto "+chrom+":"+pos);
			//todo wait
			 }
		catch(Exception err)
			{
			LOG.error(err);
			}
		finally
			{
			if(in!=null) try {in.close();} catch(Exception err){}
			if(out!=null) try {out.close();} catch(Exception err){}
			if(socket!=null) try {socket.close();} catch(Exception err){}
			}
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
	extends CommandLineProgram
	{
	private static Log LOG=Log.getInstance(VcfViewGui.class);
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Simple java-Swing-based VCF viewer. ";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF files to process.",minElements=0)
	public List<File> IN=new ArrayList<File>();

    private VCFFileRef create(File vcfFile) throws IOException
    	{
    	VCFFileRef vfr=new VCFFileRef();
    	vfr.vcfFile=vcfFile;
    	AsciiLineReader r=null;
    	r=new AsciiLineReader(IOUtils.openFileForReading(vcfFile));
    	vfr.header=(VCFHeader)vfr.codec.readHeader(r);
    	r.close();
    	return vfr;
    	}
    
	@Override
	protected int doWork()
		{
		try
			{
			JFrame.setDefaultLookAndFeelDecorated(true);
			JDialog.setDefaultLookAndFeelDecorated(true);
			List<VCFFileRef> vfrs=new ArrayList<VCFFileRef>();
			if(IN.isEmpty())
				{
				LOG.info("NO VCF provided; Opening dialog");
				JFileChooser chooser=new JFileChooser();
				chooser.setFileFilter(new FileFilter() {
					@Override
					public String getDescription() {
						return "VCF indexed with tabix";
					}
					
					@Override
					public boolean accept(File f) {
						if(f.isDirectory()) return true;
						if(f.isFile() && f.canRead() && f.getName().endsWith(".vcf.gz"))
							{
							File tabix=new File(f.getParentFile(),f.getName()+".tbi");
							if(!tabix.exists()) return false;
							if(!tabix.canRead() || tabix.lastModified()< f.lastModified()) return false;
							return true;
							}
						return false;
						};
					});
				chooser.setMultiSelectionEnabled(true);
				if(chooser.showOpenDialog(null)!=JFileChooser.APPROVE_OPTION)
					{
					LOG.info("user pressed cancel");
					return -1;
					}
				File fs[]=chooser.getSelectedFiles();
				if(fs!=null) IN.addAll(Arrays.asList(chooser.getSelectedFiles()));
				}
			
			
			for(File in:IN)
				{
				vfrs.add(create(in));
				}
			if(vfrs.isEmpty())
				{
				return -1;
				}
			LOG.info("showing VCF frame");
			Dimension screen=Toolkit.getDefaultToolkit().getScreenSize();
			VCFFrame f=new VCFFrame(vfrs);
			f.setBounds(50, 50, screen.width-100, screen.height-100);
			f.setVisible(true);
			}
		catch(Exception err)
			{
			LOG.error(err);
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
