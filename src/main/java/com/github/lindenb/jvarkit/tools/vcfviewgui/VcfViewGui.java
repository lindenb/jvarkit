package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Vector;
import java.util.regex.Pattern;

import javax.swing.JDesktopPane;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.table.AbstractTableModel;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.TabixReader;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.picard.IOUtils;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;

class VCFFileRef
	{
	File vcfFile;
	VCFCodec codec=new VCFCodec();
	VCFHeader header;
	}

class VCFInternalFrame extends JInternalFrame
	{
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
		JScrollPane scroll=new JScrollPane();
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
	}

class VCFFrame extends JDialog
	{
	private JDesktopPane desktopPane;
	private List<VCFFileRef> vcfFileRefs;
	private List<VCFInternalFrame> vcfInternalFrames=new Vector<VCFInternalFrame>();
	class Reload extends Thread
		{
		int maxRows;
		String reg;
		@Override
		public void run()
			{
			TabixReader tabixReader=null;
			try
				{
				
				for(VCFInternalFrame cif:vcfInternalFrames)
					{
					String line;
					List<String> rows=new ArrayList<String>();
					if(reg!=null)
						{
						tabixReader=new TabixReader(cif.ref.vcfFile.getPath());
						TabixReader.Iterator iter=tabixReader.query(reg);
						while(iter!=null && (line=iter.next())!=null && rows.size()<maxRows)
							{
							rows.add(line);
							}
						tabixReader.close();
						tabixReader=null;
						
						}
					else
						{
						
						}
					cif.tableModel.setRows(rows);
					}
					
				}
			catch(Exception err)
				{
				
				}
			finally
				{
				if(tabixReader!=null) tabixReader.close();
				}
			}
		}
	
	
	VCFFrame(List<VCFFileRef> vcfFileRefs)
		{
		super((JFrame)null,"VCF View ("+vcfFileRefs.size()+" files)",ModalityType.APPLICATION_MODAL);
		this.vcfFileRefs=vcfFileRefs;
		addWindowListener(new WindowAdapter()
			{
			@Override
			public void windowOpened(WindowEvent e)
				{
				Dimension d=desktopPane.getMaximumSize();
				for(VCFFileRef vfr:VCFFrame.this.vcfFileRefs)
					{
					int w=(int)(d.width*0.8);
					int h=(int)(d.width*0.8);
					VCFInternalFrame iFrame=new VCFInternalFrame(vfr);
					iFrame.setBounds(
							(int)((d.width-w)*Math.random()),
							(int)((d.height-h)*Math.random()),
							w, h);
					desktopPane.add(iFrame);
					vcfInternalFrames.add(iFrame);
					iFrame.setVisible(true);
					}
				reloadFrameContent();
				}
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
		
		}
	private void reloadFrameContent()
		{
		
		}
	private void doMenuClose()
		{
		this.setVisible(false);
		this.dispose();
		}
	}

public class VcfViewGui extends CommandLineProgram
	{
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Creates the code to insert one or more VCF into a SQL database. ";
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
			List<VCFFileRef> vfrs=new ArrayList<VCFFileRef>();
			for(File in:IN)
				{
				vfrs.add(create(in));
				}
			if(vfrs.isEmpty())
				{
				return -1;
				}
			VCFFrame f=new VCFFrame(vfrs);
			f.setVisible(true);
			}
		catch(Exception err)
			{
			return -1;
			}
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		// TODO Auto-generated method stub

		}

	}
