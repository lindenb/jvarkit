package com.github.lindenb.jvarkit.tools.ws.rmi.client;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.image.BufferedImage;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDesktopPane;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.SwingUtilities;

import net.sf.picard.util.Interval;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMLineParser;
import net.sf.samtools.SAMRecord;

import com.github.lindenb.jvarkit.tools.ws.WSBam;
import com.github.lindenb.jvarkit.tools.ws.WSProject;
import com.github.lindenb.jvarkit.tools.ws.WSVcf;
import com.github.lindenb.jvarkit.tools.ws.rmi.NGSService;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

@SuppressWarnings("serial")
class ProjectsTableModel extends AbstractGenericTable<WSProject>
	{
	private static String COLS[]={"Label","Description"};
	@Override
	public int getColumnCount()
		{
		return COLS.length;
		}
	@Override
	public Object getValueOf(WSProject o, int columnIndex)
		{
		switch (columnIndex)
			{
			case 0: return o.getLabel();
			case 1: return o.getDescription();
			}
		return null;
		}
	}


@SuppressWarnings("serial")
class BamsTableModel extends AbstractGenericTable<WSBam>
	{
	private static String COLS[]={"Path","Reference"};
	BamsTableModel(List<WSBam> L)
		{
		super(L);
		}
	@Override
	public int getColumnCount()
		{
		return COLS.length;
		}
	@Override
	public Object getValueOf(WSBam o, int columnIndex)
		{
		switch (columnIndex)
			{
			case 0: return o.getPath();
			case 1: return o.getReferenceId();
			}
		return null;
		}
	}

@SuppressWarnings("serial")
class VcfTableModel extends AbstractGenericTable<WSVcf>
	{
	private static String COLS[]={"Path","Reference"};
	VcfTableModel(List<WSVcf> L)
		{
		super(L);
		}
	@Override
	public int getColumnCount()
		{
		return COLS.length;
		}
	@Override
	public Object getValueOf(WSVcf o, int columnIndex)
		{
		switch (columnIndex)
			{
			case 0: return o.getPath();
			case 1: return o.getReferenceId();
			}
		return null;
		}
	}

@SuppressWarnings("serial")
class BamFrame extends AbstractInternalFrame
	{
	private WSBam wsBam;
	private SAMFileHeader samFileHeader=null;
	private ImageIcon imageIcon;
	BamFrame(NGSServiceClientFrame owner,WSBam wsBam)
		{
		super(owner);
		this.wsBam=wsBam;
		setTitle(wsBam.getPath());
		
		JPanel pane=new JPanel(new BorderLayout(5,5));
		setContentPane(pane);
		
		this.imageIcon=new ImageIcon();
		pane.add(new JScrollPane(new JLabel(this.imageIcon)),BorderLayout.CENTER);
		}
	
	private void x() throws IOException
		{
		ByteArrayInputStream bais=new ByteArrayInputStream(null);//TODO
		SAMFileReader sfr=new SAMFileReader(bais);
		this.samFileHeader=sfr.getFileHeader();
		sfr.close();
		bais.close();
		}
	private List<SAMRecord> y()
		{
		List<SAMRecord> records=new ArrayList<SAMRecord>();
		SAMLineParser samLineParser=new SAMLineParser(this.samFileHeader);
		int lineNumber=0;
		for(String s:new String[0])
			{
			samLineParser.parseLine(s, lineNumber++);
			}
		samLineParser.getValidationStringency();
		return records;
		}
	
	public WSBam getWsBam()
		{
		return wsBam;
		}
	
	private class DataThread extends Thread
		{
		private Interval interval;
		private BufferedImage img=null;
		@Override
		public void run()
			{
			try {
				List<SAMRecord> records=new ArrayList<SAMRecord>();
				Collections.sort(records,new Comparator<SAMRecord>()
							{
							@Override
							public int compare(SAMRecord o1, SAMRecord o2) {
								return 0;//TODO
								}
							});
				//pileup
				
				} 
			catch (Exception e)
				{
				
				}
			}
		}
	
	}

@SuppressWarnings("serial")
class ProjectFrame extends AbstractInternalFrame
	{
	private WSProject project;
	private BamsTableModel bamTableModel;
	private VcfTableModel vcfTableModel;
	private JTable bamTable;
	private JTable vcfTable;
	
	ProjectFrame(NGSServiceClientFrame owner,WSProject project)
		{
		super(owner);
		this.project=project;
		setTitle("Project: "+project.getLabel());
		try
			{
			List<WSBam>  bams=new ArrayList<WSBam>();
			List<WSVcf>  vcfs=new ArrayList<WSVcf>();
			for(String s:project.getBamIds())
				{
				WSBam bam= owner.getService().getBamById(s);
				if(bam==null) continue;
				bams.add(bam);
				}
			for(String s:project.getVcfIds())
				{
				WSVcf v= owner.getService().getVcfById(s);
				if(v==null) continue;
				vcfs.add(v);
				}
			this.bamTableModel=new BamsTableModel(bams);
			this.vcfTableModel=new VcfTableModel(vcfs);
			
			JPanel pane=new JPanel(new BorderLayout(5,5));
			this.setContentPane(pane);
			
			JTabbedPane tabbedPane=new JTabbedPane();
			pane.add(tabbedPane, BorderLayout.CENTER);

			
			
			//BAM
			JPanel pane2=new JPanel(new BorderLayout(5,5));
			tabbedPane.addTab("BAMS", pane2);
			this.bamTable=new JTable(this.bamTableModel);
			JScrollPane scroll=new JScrollPane(this.bamTable);
			pane2.add(scroll,BorderLayout.CENTER);
			
			JPanel top=new JPanel(new FlowLayout(FlowLayout.LEADING));
			pane2.add(top,BorderLayout.SOUTH);
			AbstractAction action=new AbstractAction("View")
				{
				@Override
				public void actionPerformed(ActionEvent evt)
					{
					int row=bamTable.getSelectedRow();
					if(row==-1) return;
					WSBam w=bamTableModel.getElementAt(row);
					for(JInternalFrame jif:getOwnerWindow().getDesktop().getAllFrames())
						{
						if(!(jif instanceof BamFrame))
							{
							continue;
							}
						BamFrame bf=(BamFrame)jif;
						if(bf.getWsBam().getId().equals(w.getId()))
							{
							jif.moveToFront();
							return;
							}
						}
					
					}
				};
			this.getActionMap().put("action.bam.table", action);
			top.add(new JButton(action));
			
			//VCF
			pane2=new JPanel(new BorderLayout(5,5));
			tabbedPane.addTab("VCF", pane2);
			this.vcfTable=new JTable(this.vcfTableModel);
			scroll=new JScrollPane(this.vcfTable);
			pane2.add(scroll,BorderLayout.CENTER);

			
			}
		catch(Exception err)
			{
			
			}
		}
	
	public WSProject getWSProject()
		{
		return project;
		}
	}

@SuppressWarnings("serial")
class ProjectsFrame extends AbstractTableFrame<WSProject>
	{
	ProjectsFrame(NGSServiceClientFrame owner)
		{
		super(owner);
		}
	
	
	private void x()
		{
		new Thread()
			{
			public void run()
					{
					try
						{
						final List<? extends WSProject> L=getService().getProjects();
						SwingUtilities.invokeAndWait(new Runnable()
							{
							@Override
							public void run()
								{
								
								}
							});
						}
					catch (Exception e) {
						
						}
					}
			}.start();
		}
	}

abstract class AbstractTableFrame<T> extends AbstractInternalFrame
	{
	AbstractTableFrame(NGSServiceClientFrame owner)
		{
		super(owner);
		}
	}

@SuppressWarnings("serial")
abstract class AbstractInternalFrame
extends JInternalFrame
	{
	static final Log LOG=Log.getInstance(AbstractInternalFrame.class);
	private NGSServiceClientFrame owner;
	
	AbstractInternalFrame(NGSServiceClientFrame owner)
		{
		this.owner=owner;
		}
	
	NGSServiceClientFrame getOwnerWindow()
		{
		return this.owner;
		}
	
	NGSService getService()
		{
		return getOwnerWindow().getService();
		}
	}


@SuppressWarnings("serial")
class NGSServiceClientFrame extends JFrame
	{
	private NGSService service;
	private JDesktopPane desktop;
	static final Log LOG=Log.getInstance(NGSServiceClientFrame.class);
	
	NGSServiceClientFrame(NGSService service)
		{
		this.service=service;
		
		}
	
	public Interval getInterval()
		{
		return null;
		}
	
	public Log getLog()
		{	
		return LOG;
		}
	
	public JDesktopPane getDesktop() {
		return desktop;
	}
	
	NGSService getService()
		{
		return this.service;
		}
	}


@SuppressWarnings("serial")
public class NGSServiceClient extends AbstractCommandLineProgram
	{
	@Override
	protected int doWork() {
		return 0;
		}
	}
