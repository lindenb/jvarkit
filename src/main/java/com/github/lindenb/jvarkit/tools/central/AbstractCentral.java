package com.github.lindenb.jvarkit.tools.central;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.prefs.Preferences;

import javax.swing.AbstractAction;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.border.EmptyBorder;
import javax.swing.text.JTextComponent;

import com.github.lindenb.jvarkit.tools.bam2svg.AbstractBamToSVG;
import com.github.lindenb.jvarkit.tools.misc.AbstractAddLinearIndexToBed;
import com.github.lindenb.jvarkit.tools.misc.AbstractVcfHead;
import com.github.lindenb.jvarkit.tools.misc.AddLinearIndexToBed;
import com.github.lindenb.jvarkit.tools.misc.VcfHead;
import com.github.lindenb.jvarkit.tools.vcfcmp.AbstractVcfCompareCallers;
import com.github.lindenb.jvarkit.util.swing.FormLayout;
import com.github.lindenb.jvarkit.util.swing.InputChooser;
import com.github.lindenb.jvarkit.util.swing.OutputChooser;

@SuppressWarnings("serial")
public abstract class AbstractCentral extends JFrame {
	protected JTabbedPane tabbedPane;
	protected Preferences preferences = Preferences.userNodeForPackage(AbstractCentral.class);
	protected JTextComponent logArea;
	
	protected AbstractCentral() {
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
	
	this.tabbedPane = new JTabbedPane(JTabbedPane.TOP);
	this.tabbedPane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
	contentPane.add(this.tabbedPane, BorderLayout.CENTER);
	
	JPanel configPane = new JPanel(new BorderLayout());
	this.tabbedPane.addTab("Config", configPane);
	JPanel pane2= new JPanel(new FormLayout());
	pane2.setBorder(new EmptyBorder(5, 5, 5, 5));
	pane2.add(new JLabel("REF"));
	pane2.add(new InputChooser());
	pane2.add(new JLabel("REF"));
	pane2.add(new OutputChooser());
	
	configPane.add(new JScrollPane(pane2),BorderLayout.CENTER);
	
	this.tabbedPane.addTab("VCFHead", new AbstractVcfHead.VcfHeadSwingUI());
	this.tabbedPane.addTab("VCFCompareCallers", new AbstractVcfCompareCallers.VcfCompareCallersSwingUI());
	this.tabbedPane.addTab("BamToSVG", new AbstractBamToSVG.BamToSVGSwingUI());
	this.tabbedPane.addTab("AbstractAddLinearIndexToBed", new AbstractAddLinearIndexToBed.AddLinearIndexToBedSwingUI());
	
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
			removeWindowListener(this);
			}
		});
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
	
}
