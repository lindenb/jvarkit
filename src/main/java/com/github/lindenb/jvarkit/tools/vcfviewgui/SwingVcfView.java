/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFFilterHeaderLineTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFFormatHeaderLineTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFGenotypesTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFInfoHeaderLineTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFInfoTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVariantsTableModel;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC

## Example



END_DOC
 */
@Program(name="swingvcfview",
description="VCFviewer using Java Swing UI",
keywords={"vcf","visualization","swing"},
creationDate="20210503",
modificationDate="20210503",
generate_doc=false
)
public class SwingVcfView extends Launcher
	{
	private static final Logger LOG = Logger.build(SwingVcfView.class).make();
	@Parameter(names={"-r","--regions","--interval"},description="default interval region on opening")
	private String defaultRegion="";

	
	@SuppressWarnings("serial")
	private static class XFrame extends JFrame {
		final SAMSequenceDictionary dict;
		final VCFReader vcfReader;
		final JTextField jtextFieldLocation;
		final SwingVariantsTableModel swingVariantsTableModel;
		
		
		XFrame(final Path vcfPath,String defaultLoc) {
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			setTitle(SwingVcfView.class.getSimpleName());
			final VCFReaderFactory vcfReaderFactory = VCFReaderFactory.makeDefault();
			vcfReader = vcfReaderFactory.open(vcfPath, true);
			final VCFHeader header = vcfReader.getHeader();
			
			this.dict = header.getSequenceDictionary();
			
			final JPanel mainPane = new JPanel(new BorderLayout(5, 5));
			mainPane.setBorder(new EmptyBorder(5, 5, 5, 5));
			this.setContentPane(mainPane);
			final Panel topPane = new Panel(new FlowLayout(FlowLayout.LEADING));
			mainPane.add(topPane,BorderLayout.NORTH);
			JLabel label = new JLabel("Location:", JLabel.RIGHT);
			topPane.add(label);
			this.jtextFieldLocation= new JTextField(defaultLoc,20);
			if(StringUtil.isBlank(defaultLoc) &&  !dict.isEmpty()) {
				final SAMSequenceRecord ssr= dict.getSequence(0);
				this.jtextFieldLocation.setText(ssr.getSequenceName()+":1-1000");
				}
			
			topPane.add(this.jtextFieldLocation);
			label.setLabelFor(this.jtextFieldLocation);
			final Action actionGo = new AbstractAction("Load")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					refreshInterval();
					}
				};
			this.jtextFieldLocation.addActionListener(actionGo);
			JButton button = new JButton(actionGo);
			topPane.add(button);
			topPane.add(new JSeparator());
			label = new JLabel("Cap:", JLabel.RIGHT);
			topPane.add(label);
			

			final Panel botPane = new Panel(new FlowLayout(FlowLayout.TRAILING));
			mainPane.add(botPane, BorderLayout.SOUTH);

						
			final JTabbedPane  tabPane = new JTabbedPane();
			mainPane.add(tabPane,BorderLayout.CENTER);
			
			final JPanel pane2 = new JPanel(new BorderLayout(5,5));
			this.swingVariantsTableModel = new SwingVariantsTableModel();
			final JTable variantTable =  new JTable(this.swingVariantsTableModel);
			pane2.add(new JScrollPane(variantTable),BorderLayout.WEST);
			final SwingVCFGenotypesTableModel swingVCFGenotypesTableModel = new SwingVCFGenotypesTableModel();
			final JTable genotypeTable =  new JTable(swingVCFGenotypesTableModel);
			pane2.add(new JScrollPane(genotypeTable),BorderLayout.CENTER);
			final SwingVCFInfoTableModel swingInfoTableModel = new SwingVCFInfoTableModel();
			final JTable infoTable =  new JTable(swingInfoTableModel);
			pane2.add(new JScrollPane(infoTable),BorderLayout.SOUTH);
			
			
			variantTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {				
				@Override
				public void valueChanged(ListSelectionEvent e) {
					int i = variantTable.getSelectedRow();
					VariantContext ctx =  i>=0? swingVariantsTableModel.getElementAt(i):null;
					swingVCFGenotypesTableModel.setVariant(ctx);
					swingInfoTableModel.setVariant(ctx);
				}
			});
			
			
			tabPane.addTab("Main", pane2);
			
			tabPane.addTab("INFO", new JScrollPane( new JTable(new SwingVCFInfoHeaderLineTableModel(header))));
			tabPane.addTab("FORMAT", new JScrollPane( new JTable(new SwingVCFFormatHeaderLineTableModel(header))));
			tabPane.addTab("FILTER", new JScrollPane( new JTable(new SwingVCFFilterHeaderLineTableModel(header))));

			if(this.dict!=null) {
				tabPane.addTab("REF", new JScrollPane( new JTable(new SwingSequenceDictionaryTableModel(this.dict))));
				}
			if(!StringUtil.isBlank(defaultLoc)) {
				this.addWindowListener(new WindowAdapter()
					{
					@Override
					public void windowOpened(final WindowEvent e)
						{
						refreshInterval();
						}
					public void windowClosed(WindowEvent e) {
						try {vcfReader.close(); }
						catch(Throwable err) {LOG.error(err);}
						};
					});
				}
			final JMenuBar menuBar= new JMenuBar();
			setJMenuBar(menuBar);
			final JMenu menu = new JMenu("File");
			menuBar.add(menu);
			menu.add(new JSeparator());
			menu.add(actionGo);
			menu.add(new JSeparator());
			menu.add(new JMenuItem(new AbstractAction("Quit")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					XFrame.this.setVisible(false);
					XFrame.this.dispose();
					}
				}));
			}
		
		final Optional<SimpleInterval> getUserInterval() {
			final String s = this.jtextFieldLocation.getText(); 
			if(StringUtil.isBlank(s)) return Optional.empty();
			return IntervalParserFactory.newInstance().
					dictionary(this.dict).
					make().
					apply(s);
			}
				
		private void refreshInterval() {
			final List<VariantContext> L;
			final Optional<SimpleInterval> location = getUserInterval();
			if(!location.isPresent()) {
				L = Collections.emptyList();
			}
			else
			{
			try(CloseableIterator<VariantContext> r=this.vcfReader.query(location.get())) {
				L = r.stream().collect(Collectors.toList());	
				}
			}
			this.swingVariantsTableModel.setRows(L);
		}
	}
	
	@Override
	public int doWork(final List<String> args)
		{
		try {
			final String input = this.oneAndOnlyOneFile(args);
			final Path vcfPath = Paths.get(input);
			IOUtil.assertFileIsReadable(vcfPath);
			
			JFrame.setDefaultLookAndFeelDecorated(true);
			final XFrame frame = new XFrame(vcfPath,defaultRegion);
			final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
			frame.setBounds(50, 50, screen.width-100, screen.height-100);

			
			SwingUtilities.invokeAndWait(()->{
				frame.setVisible(true);
				});
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args)
		{
		new SwingVcfView().instanceMain(args);//no exit
		}
	}
