package com.github.lindenb.jvarkit.tools.central;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.prefs.Preferences;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.text.Document;
import javax.swing.text.JTextComponent;

import org.slf4j.Logger;

import com.github.lindenb.jvarkit.tools.bam2svg.AbstractBamToSVG;
import com.github.lindenb.jvarkit.tools.bamstats04.AbstractBamStats05;
import com.github.lindenb.jvarkit.tools.bioalcidae.AbstractBioAlcidae;
import com.github.lindenb.jvarkit.tools.misc.AbstractAddLinearIndexToBed;
import com.github.lindenb.jvarkit.tools.misc.AbstractFindAVariation;
import com.github.lindenb.jvarkit.tools.misc.AbstractVcfHead;
import com.github.lindenb.jvarkit.tools.misc.AbstractVcfMultiToOneAllele;
import com.github.lindenb.jvarkit.tools.onesamplevcf.AbstractVcfMultiToOne;
import com.github.lindenb.jvarkit.tools.samjs.AbstractSamJavascript;
import com.github.lindenb.jvarkit.tools.vcf2sql.AbstractVcfToSql;
import com.github.lindenb.jvarkit.tools.vcfcmp.AbstractVcfCompareCallers;
import com.github.lindenb.jvarkit.tools.vcffilterjs.AbstractVCFFilterJS;
import com.github.lindenb.jvarkit.tools.vcffixindels.AbstractVCFFixIndels;
import com.github.lindenb.jvarkit.util.log.Logging;

@SuppressWarnings("serial")
public class Central extends JFrame {
	private static final Logger LOG = Logging.getLog(Central.class);
	private JMenu toolsMenu;
	protected JTabbedPane tabbedPane;
	protected Preferences preferences = Preferences.userNodeForPackage(Central.class);
	protected JTextComponent logArea;
	private Thread runningThread =null;
	protected Central() {
	super("JVARKIT");
	super.setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
	
	final JPanel contentPane = new JPanel(new BorderLayout(5,5));
	contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
	setContentPane(contentPane);
	
	
	
	final JMenuBar bar = new JMenuBar();
	this.setJMenuBar(bar);
	JMenu menu = new JMenu("File");
	bar.add(menu);
	menu.add(new JMenuItem(new AbstractAction("Quit")
		{
		@Override
		public void actionPerformed(ActionEvent e) {
			doMenuQuit();
			}
		}));
	
	this.toolsMenu = new JMenu("Tools");
	bar.add(this.toolsMenu);
	
	this.tabbedPane = new JTabbedPane(JTabbedPane.TOP);
	this.tabbedPane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
	contentPane.add(this.tabbedPane, BorderLayout.CENTER);
	
	this.addTab(AbstractAddLinearIndexToBed.AddLinearIndexToBedSwingUI.class);
	this.addTab(AbstractBamStats05.BamStats05SwingUI.class);
	this.addTab(AbstractBamToSVG.BamToSVGSwingUI.class);
	this.addTab(AbstractBioAlcidae.BioAlcidaeSwingUI.class);
	this.addTab(AbstractFindAVariation.FindAVariationSwingUI.class);
	this.addTab(AbstractSamJavascript.SamJavascriptSwingUI.class);
	this.addTab(AbstractVCFFixIndels.VCFFixIndelsSwingUI.class);
	this.addTab(AbstractVcfHead.VcfHeadSwingUI.class);
	this.addTab(AbstractVcfCompareCallers.VcfCompareCallersSwingUI.class);
	this.addTab(AbstractVCFFilterJS.VCFFilterJSSwingUI.class);
	this.addTab(AbstractVcfMultiToOneAllele.VcfMultiToOneAlleleSwingUI.class);
	this.addTab(AbstractVcfMultiToOne.VcfMultiToOneSwingUI.class);
	this.addTab(AbstractVcfToSql.VcfToSqlSwingUI.class);
	
	
	
	
	
	JPanel bottom = new JPanel(new BorderLayout(5,5));
	bottom.setBorder(new EmptyBorder(5, 5, 5, 5));
	contentPane.add(bottom, BorderLayout.SOUTH);
	this.logArea = new JTextArea(10, 40);
	this.logArea.setFont(new Font("Courier",Font.PLAIN,9));
	this.logArea.setEditable(false);
	JScrollPane scroll = new JScrollPane(this.logArea);
	bottom.add(scroll,BorderLayout.CENTER);
	
	
	this.addWindowListener(new WindowAdapter()
		{
		@Override
		public void windowClosing(WindowEvent e) {
			doMenuQuit();
			}
		});
	this.addWindowListener(new WindowAdapter()
		{
		@Override
		public void windowOpened(WindowEvent e)
			{
			System.setErr(new PrintStream(new ConsolePrintStream(128),true));
			System.setOut(new PrintStream(new ConsolePrintStream(2028),true));
			removeWindowListener(this);
			}
		});
	}

private void addTab(Class<? extends CentralPane> clazz)
	{
	final String label;
	final CentralPane pane;
	try {
		label = clazz.getSimpleName().replaceAll("SwingUI", "");
		pane = clazz.newInstance();
	} catch (Exception e) {
		e.printStackTrace();
		return;
		}
	
	final int tabCount =  tabbedPane.getTabCount();
	final ExcutePane executePane = new ExcutePane(pane);
	this.tabbedPane.addTab(label,executePane);
	final AbstractAction action  = new AbstractAction(label)
		{
		@Override
		public void actionPerformed(final java.awt.event.ActionEvent evt)
			{
			tabbedPane.setSelectedIndex(tabCount);
			}
		};
	action.putValue(AbstractAction.LONG_DESCRIPTION,pane.getDescription());		
	action.putValue(AbstractAction.SHORT_DESCRIPTION,pane.getLabel());		
	this.toolsMenu.add(new JMenuItem(action));
	}
	
private void doMenuQuit()
	{
	setVisible(false);
	dispose();
	try {
	preferences.sync();
	preferences.flush();
	} catch(Exception err) {}
	}

private class ExcutePane extends JPanel
	{
	private CentralPane delegate;
	private AbstractAction runAction;
	private AbstractAction cancelAction;
	
	
	

	
	ExcutePane(CentralPane pane)
		{
		super(new BorderLayout(5,5));
		this.delegate = pane;
		this.setBorder(new EmptyBorder(5, 5, 5, 5));
		JScrollPane scroll = new JScrollPane(pane);
		this.add(scroll,BorderLayout.CENTER);
		
		JPanel bottomPane = new JPanel(new FlowLayout(FlowLayout.TRAILING));
		this.add(bottomPane,BorderLayout.SOUTH);
		
		this.cancelAction = new AbstractAction("Cancel")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				if(Central.this.runningThread==null) return;
				try {
					Central.this.runningThread.interrupt();
					}
				catch (Exception e2) {
					}
				Central.this.runningThread=null;
				}
			};
		runAction = new AbstractAction(delegate.getLabel())
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				if(Central.this.runningThread!=null)
					{
					JOptionPane.showMessageDialog(Central.this, "Command already running");
					return;
					}
				final String msg = delegate.getValidationMessage();
				if(msg!=null)
					{
					JOptionPane.showMessageDialog(
							Central.this,
							msg);
					return;
					}
				final List<String> cmd = buildCommandLine();
				Central.this.runningThread = new CentralRunner(
						delegate.getMainClass(),
						cmd
						);
				Central.this.runningThread.start();
				}
			};
		Font font = new Font("Dialog", Font.BOLD, 18);
		JButton button = new JButton(cancelAction);
		button.setContentAreaFilled(true);
		button.setBackground(Color.ORANGE);
		button.setForeground(Color.WHITE);
		button.setFont(font);
		bottomPane.add(button);
		button = new JButton(runAction);
		button.setContentAreaFilled(true);
		button.setBackground(Color.GREEN);
		button.setForeground(Color.WHITE);
		button.setFont(font);
		bottomPane.add(button);
		}
	
	public List<String> buildCommandLine()
		{
		List<String> L = new ArrayList<>();
		delegate.fillCommandLine(L);
		return L;
		}

	}	
/**
 * GATKRunner
 *
 */
private  class CentralRunner extends Thread
	{
	private Class<? extends com.github.lindenb.jvarkit.util.command.Command > commandClass;
	private String args[];
	public CentralRunner(
				final  Class<? extends com.github.lindenb.jvarkit.util.command.Command > commandClass,
				final List<String> args)
		{
		this.commandClass = commandClass;
		this.args=args.toArray(new String[args.size()]);
		}
	@Override
	public void run()
		{
		LOG.info("starting "+Arrays.toString(args));
		

		try
			{
			com.github.lindenb.jvarkit.util.command.Command command  = commandClass.newInstance();
			//command.stderr( new PrintStream(new ConsolePrintStream(0),true));
			//command.stdout( new PrintStream( new ConsolePrintStream(1024)));
			
			List<Throwable> errors = new ArrayList<>();
			int ret=command.instanceMainWithExceptions(args, errors);
			if(ret == 0)
				{
				try {
					SwingUtilities.invokeAndWait(new Runnable() {
						@Override
						public void run() {
							if(CentralRunner.this != Central.this.runningThread) return;
							
							JOptionPane.showMessageDialog(Central.this,
							new JScrollPane(new JTextArea("Completed:"+Arrays.toString(args),5,20)),
							"Completed",
							JOptionPane.INFORMATION_MESSAGE,
							null);
							Central.this.runningThread=null;
						}
					});
				} catch (Exception e) {
					LOG.warn(e.getMessage());
					}
				}
			else
				{
				try {
					SwingUtilities.invokeAndWait(new Runnable() {
						@Override
						public void run() {
							if(CentralRunner.this != Central.this.runningThread) return;
							JOptionPane.showMessageDialog(
									Central.this,
							new JScrollPane(new JTextArea("Failure:"+Arrays.toString(args),5,20)),
							"Failure",
							JOptionPane.ERROR_MESSAGE,
							null);
							Central.this.runningThread=null;
						}
					});
				} catch (Exception e) {
					LOG.warn(e.getMessage());
					
					}
				}
			}
		catch(final Exception err)
			{
			try {
				SwingUtilities.invokeAndWait(new Runnable() {
					@Override
					public void run() {
						if(CentralRunner.this != Central.this.runningThread) return;
						JOptionPane.showMessageDialog(Central.this,
							new JScrollPane(new JTextArea("Failure:"+Arrays.toString(args),5,20)),
							"Failure",
							JOptionPane.ERROR_MESSAGE,
							null);
						Central.this.runningThread=null;
					}
				});
			} catch (Exception e) {
				LOG.warn(e.getMessage());
				
				}
			}
		}
	}

	private class ConsolePrintStream extends OutputStream
		{
		private int buffsize;
		private StringBuilder buffer= new StringBuilder();
		ConsolePrintStream(int buffsize)
			{
			this.buffsize = buffsize;
			}
		@Override
		public void write(int b) throws IOException {
			if(b==-1) flush();
			buffer.append((char)b);
			if(b=='\n' || buffer.length()>=this.buffsize) flush();
			}
		@Override
		public void flush() throws IOException {
			if(buffer.length()==0) return;
			writeToDevice(buffer.toString());
			buffer.setLength(0);
			}
		public void writeToDevice(final String logString)
			{
			try
			{
			SwingUtilities.invokeAndWait(new Runnable()
				{
				@Override
				public void run() {


					if(!logArea.isVisible())
						{
						//System.err.println(logString);
						return;
						}
					final Document doc = logArea.getDocument();	
					try
						{
						int L = doc.getLength();
						if(L>10000) doc.remove(0, L);
						L = doc.getLength();
						doc.insertString(doc.getLength(),logString, null );
						}
					catch(Exception err)
						{
						err.printStackTrace();
						}
					}
				});
			} catch(Exception err)
				{
				
				}
			}
		}
	public static void main(String[] args)
		{
		JFrame.setDefaultLookAndFeelDecorated(true);
		JDialog.setDefaultLookAndFeelDecorated(true);
		final Central app = new Central();
		try
			{
			SwingUtilities.invokeAndWait(new Runnable()
				{
				@Override
				public void run()
					{
					app.pack();
					final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
					final Dimension dim = app.getPreferredSize();
					Rectangle r = new Rectangle(										
							(screen.width-dim.width)/2,
							(screen.height-dim.height)/2,
							dim.width, dim.height
							);
					app.setBounds(
							Math.max(0,r.x),
							Math.max(0,r.y),
							Math.min(r.width,screen.width),
							Math.min(r.width,screen.height)
							);
					app.setVisible(true);
					}
				});
			}
		catch(Exception err)
			{
			LOG.error(err.getMessage());
			}
		}

}
