/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.MouseAdapter;
import java.awt.event.WindowAdapter;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Vector;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.swing.SAMRecordPanel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.swing.PropertyChangeObserver;
import com.github.lindenb.jvarkit.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.iterator.LineIterators;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFGenotypesTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFInfoTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVariantsTableModel;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

## Example:

```
java -jar dist/jvarkit.jar swingbamview -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam
```

```
find dir -type f -name "*.bam" > out.list
java -jar dist/jvarkit.jar swingbamview -R src/test/resources/rotavirus_rf.fa out.list
```

END_DOC
 */
@Program(name="swingbamview",
description="Read viewer using Java Swing UI",
keywords={"bam","alignment","graphics","visualization","swing"},
creationDate = "20220503",
modificationDate="20230427",
jvarkit_amalgamion =  true,
menu="BAM Manipulation"
)
public class SwingBamView extends Launcher {
	private static final Logger LOG = Logger.build(SwingBamView.class).make();
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFile = null;
	@Parameter(names={"-V","--variant"},description="Indexed VCF File")
	private Path variantFile = null;
	@Parameter(names={"--gff","--gff3"},description="Tabix-Indexed GFF3 File")
	private String gffFile = null;

		
	
	
	
	private static class BamFile {
		final Path path;
		String sn;
		SAMFileHeader header= null;
		SamReader sr = null;
		BamFile(final Path path) {
			this.path = path;
			this.sn = IOUtils.getFilenameWithoutCommonSuffixes(path);
			}
		
		String getName() {
			return this.sn;
		}
		
		void open(final Path referenceFile) {
			if(sr!=null) return;
			final SamReaderFactory srf = SamReaderFactory.makeDefault();
			srf.referenceSequence(referenceFile);
			sr = srf.open(this.path);
			header = sr.getFileHeader();
			}
		
		void close() {
			if(sr!=null) try {
				sr.close();
				} catch(IOException er) {
					er.printStackTrace();
				}
			sr=null;
			header=null;
			}
		@Override
		public String toString() {
			return getName();
			}
		}

	
	@SuppressWarnings("serial")
	private static class XFrame extends JFrame {
		private final Path referenceFaidx;
		private final SAMSequenceDictionary dict;
		private final Path vcfFile;
		private final String gffFile;
		private List<VariantContext> all_variants_list = new Vector<>();
		private final JList<BamFile> jlistBamFilesJList;
		private final JTextField jtextFieldLocation;
		private final PropertyChangeObserver<Locatable> interval = new PropertyChangeObserver<>();
		private final PropertyChangeObserver<Integer> mappingQuality = new PropertyChangeObserver<>(1);
		private final JScrollBar vScrollBar;
		private final JPanel drawingArea;
		private final ReferenceSequenceFile referenceSequenceFile;
		private String referenceStr = null;
		private final List<List<Read>> rows = new Vector<>();
		private final List<Gff3Feature> gff3features = new Vector<>();
		private final JCheckBoxMenuItem jmenuShowClip;
		private final JCheckBoxMenuItem jmenuShowDuplicates;
		private final JCheckBoxMenuItem jmenuShowFailQuality;
		private final JCheckBoxMenuItem jmenuShowSecondary;
		private final JCheckBoxMenuItem jmenuShowSupplementary;
		private final JCheckBoxMenuItem jmenuShowCoverage;
		private final JCheckBoxMenuItem jmenuShowIntervalTitle;
		private final JCheckBoxMenuItem jmenuShowSampleTitle;
		private final JCheckBoxMenuItem jmenuShowReferenceSequence;
		private final JCheckBoxMenuItem jmenuShowBases;
		private final JCheckBoxMenuItem jmenuShowMismatches;
		private final JCheckBoxMenuItem jmenuShowInsertions;
		
		
		private final SwingVariantsTableModel swingVariantsTableModel;
		private final SwingVCFInfoTableModel swingVCFInfoTableModel;
		private final SwingVCFGenotypesTableModel swingVCFGenotypesTableModel;

		
		private final JMenuItem jmenuMappingQuality;
		private final JMenu jMenuHyperlinks;
		
		private class Read extends Rectangle2D.Float implements Locatable {
			private final SAMRecord rec;
			Read(final SAMRecord rec) {
				this.rec =  rec;
				}
			@Override
			public String getContig() {
				return rec.getContig();
				}
			@Override
			public int getStart() {
				return XFrame.this.isUsingClip()?rec.getUnclippedStart():rec.getStart();
				}
			@Override
			public int getEnd() {
				return XFrame.this.isUsingClip()?rec.getUnclippedEnd():rec.getEnd();
				}
			
			public String toHtml() {
				final StringBuilder sb = new StringBuilder("<HTML><BODY><DL>");
				sb.append("<DT>Name</DT><DD>"+ rec.getReadName()+ "</DD>");
				sb.append("<DT>Length</DT><DD>"+ rec.getReadLength()+ "</DD>");
				sb.append("<DT>Flags</DT><DD>"+ rec.getFlags() + " (" + Arrays.stream(SAMFlag.values()).filter(F->F.isSet(rec.getFlags())).map(F->F.name()).collect(Collectors.joining(",")) + ")</DD>");
				sb.append("<DT>Mapq</DT><DD>"+ rec.getMappingQuality() + "</DD>");
				if(!rec.getReadUnmappedFlag()) {
					sb.append("<DT>Contig</DT><DD>"+ rec.getContig()+ "</DD>");
					sb.append("<DT>Unclipped Start</DT><DD>"+ rec.getUnclippedStart()+ "</DD>");
					sb.append("<DT>Start</DT><DD>"+ rec.getUnclippedStart()+ "</DD>");
					sb.append("<DT>End</DT><DD>"+ rec.getUnclippedEnd()+ "</DD>");
					sb.append("<DT>Unclipped End</DT><DD>"+ rec.getUnclippedEnd()+ "</DD>");
					}
				else
					{
					sb.append("<DT>Contig</DT><DD><I>unmapped</I></DD>");
					}
				if(rec.getReadPairedFlag()) {
					if(!rec.getMateUnmappedFlag()) {
						sb.append("<DT>Mate Contig</DT><DD>"+ rec.getMateReferenceName()+ "</DD>");
						sb.append("<DT>Mate Start</DT><DD>"+ rec.getMateAlignmentStart()+ "</DD>");
						}
					else
						{
						sb.append("<DT>Mate Contig</DT><DD><I>unmapped</I></DD>");
						}
					}
				return sb.append("</DL></BODY></HTML>").toString();
				}
			
			@Override
			public String toString() {
				return rec.toString()+" ("+getX()+","+getY()+","+getWidth()+","+getHeight()+")";
				}
			}
		
		private abstract class ChangeViewAction extends AbstractAction {
			final double factor;
			ChangeViewAction(String title,double factor) {
				super(title);
				this.factor=factor;
				}
			@Override
			public Object getValue(String key)
				{
				if(key.equals(AbstractAction.SHORT_DESCRIPTION)) return getShortDesc();
				return super.getValue(key);
				}
			@Override
			public void actionPerformed(final ActionEvent e)
				{
				final Optional<SimpleInterval> optR = getUserInterval();
				if(!optR.isPresent()) return;
				Locatable r= optR.get();
				final SAMSequenceRecord ssr=dict.getSequence(r.getContig());
				if(ssr==null) return;
				r= change(ssr,r);
				if(r==null) return;
				XFrame.this.interval.setValue(r);
				}
			abstract Locatable change(SAMSequenceRecord ssr,Locatable loc);
			abstract String getShortDesc();
			@Override
			public String toString()
				{
				return String.valueOf(getValue(AbstractAction.NAME));
				}
			}
		
		
		private class ZoomAction extends ChangeViewAction {
			ZoomAction(double factor) {
				super("x"+factor,factor);
				}
			@Override
			Locatable change(SAMSequenceRecord ssr, Locatable r)
				{
				final double L= r.getLengthOnReference()/2.0*factor;
				final double mid = r.getStart()+r.getLengthOnReference()/2;
				int x1= (int)Math.max(1,mid-L);
				int x2= (int)Math.min(mid+L,ssr.getLengthOnReference());
				return new SimpleInterval(r.getContig(),x1,x2);
				}
			@Override
			String getShortDesc()
				{
				return "Scale by "+factor;
				}
			}
		
		private class ShiftAction extends ChangeViewAction {
			ShiftAction(String title,double factor) {
				super(title,factor);
				}
			@Override
			Locatable change(SAMSequenceRecord ssr, Locatable r)
				{
				final double L= r.getLengthOnReference()*Math.abs(factor);
				int x1,x2;
				if(factor<0) {
					x1 = Math.max(1,r.getStart() - (int)L);
					x2 = Math.min(ssr.getLengthOnReference(),x1 + r.getLengthOnReference());
					}
				else
					{
					x2 = Math.min(ssr.getLengthOnReference(),r.getEnd() + (int)L);
					x1 = Math.max(1,x2 -r.getLengthOnReference());
					}
				return new SimpleInterval(r.getContig(),x1,x2);
				}
			@Override
			String getShortDesc()
				{
				return "Shift by "+factor;
				}
			}

		
		
		XFrame(final Path referenceFaidx,
				final List<BamFile> bamFiles,
				final Path vcfFile,
				final String gffFile
				)
			{
			super(SwingBamView.class.getSimpleName());
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			this.vcfFile = vcfFile;
			this.gffFile = gffFile;
			final JMenuBar menuBar = new JMenuBar();
			setJMenuBar(menuBar);
			JMenu menu=new JMenu("File");
			menuBar.add(menu);
			menu.add(new JSeparator());
			menu.add(new JMenuItem(new AbstractAction("Quit") {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					XFrame.this.setVisible(false);
					XFrame.this.dispose();
					}
				}));
			
			this.referenceFaidx = referenceFaidx;
			this.referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referenceFaidx);
			this.dict = SequenceDictionaryUtils.extractRequired(this.referenceSequenceFile);
			
			final JPanel mainPanel=new JPanel(new BorderLayout());
			this.setContentPane(mainPanel);
			
			JTabbedPane jTabbedPane = new JTabbedPane();
			mainPanel.add(jTabbedPane,BorderLayout.CENTER);
			JPanel viewPane= new JPanel(new BorderLayout());
			jTabbedPane.addTab("View", viewPane);
			
			/** navigation */
			JPanel navPane = new JPanel(new FlowLayout(FlowLayout.LEADING));
			viewPane.add(navPane,BorderLayout.NORTH);
			
			JLabel label = new JLabel("Location:", JLabel.RIGHT);
			navPane.add(label);
			this.jtextFieldLocation= new JTextField("",20);
			navPane.add(this.jtextFieldLocation);
			label.setLabelFor(this.jtextFieldLocation);
			final Action actionGo = new AbstractAction("Load") {
				@Override
				public void actionPerformed(ActionEvent e)
					{
					final String s = jtextFieldLocation.getText().trim();
					if(StringUtil.isBlank(s)) {
						XFrame.this.interval.reset();
						return;
						}
					try {
						XFrame.this.interval.setValue(
							new IntervalParser(XFrame.this.dict).
							apply(s).
							orElse(null)
							);
						if(!XFrame.this.interval.isPresent()) {
							JOptionPane.showMessageDialog(XFrame.this, "Cannot parse "+s);
							}
						}
					catch(final Throwable err) {
						ThrowablePane.show(XFrame.this, err);
						}
					}
				};
			this.jtextFieldLocation.addActionListener(actionGo);
			navPane.add(new JButton(actionGo));
			navPane.add(new JSeparator());
			navPane.add(new JButton(new ZoomAction(0.1)));
			navPane.add(new JButton(new ZoomAction(0.333)));
			navPane.add(new JButton(new ZoomAction(0.5)));
			navPane.add(new JButton(new ZoomAction(1.5)));
			navPane.add(new JButton(new ZoomAction(3.0)));
			navPane.add(new JButton(new ZoomAction(10)));
			navPane.add(new JSeparator());
			navPane.add(new JButton(new ShiftAction("<<<",-0.9)));
			navPane.add(new JButton(new ShiftAction("<<",-0.5)));
			navPane.add(new JButton(new ShiftAction("<",-0.1)));
			navPane.add(new JButton(new ShiftAction(">",0.1)));
			navPane.add(new JButton(new ShiftAction(">>",0.5)));
			navPane.add(new JButton(new ShiftAction(">>>",0.9)));
			
			/** list of bam files */
			final DefaultListModel<BamFile> dm = new DefaultListModel<>();
			bamFiles.stream().
				sorted((A,B)->A.getName().compareTo(B.getName())).
				forEach(B->dm.addElement(B));
				
			this.jlistBamFilesJList = new JList<>(dm);
			this.jlistBamFilesJList.setPreferredSize(new Dimension(200,200));
			this.jlistBamFilesJList.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
			this.jlistBamFilesJList.setSelectedIndex(0);
			this.jlistBamFilesJList.getSelectionModel().addListSelectionListener(AE->{
				reloadReadsAndRef();
				});
			final JPanel bamPanel = new JPanel(new BorderLayout());
			
			bamPanel.setBorder(BorderFactory.createTitledBorder("BAMS N="+StringUtils.niceInt(bamFiles.size())));
			bamPanel.add(new JScrollPane(this.jlistBamFilesJList),BorderLayout.CENTER);
			viewPane.add(bamPanel,BorderLayout.WEST);
			
			
			final JPanel drawingPane = new JPanel(new BorderLayout());
			/** drawing area itself */
			this.drawingArea = new JPanel(null, true) {
				
				@Override
				public String getToolTipText(final java.awt.event.MouseEvent event) {
					final Optional<Read> r = findReadAt(event.getX(),event.getY());
					if(!r.isPresent()) return null;
					return r.get().toHtml();
					};
				@Override
				protected void paintComponent(final Graphics g) {
					paintDrawingArea(Graphics2D.class.cast(g));
					}
				};
			this.drawingArea.addMouseListener(new MouseAdapter() {
				public void mouseClicked(java.awt.event.MouseEvent e)
					{
					if(e.getClickCount()<2) return;
					final Optional<Read> r = findReadAt(e.getX(),e.getY());
					if(!r.isPresent()) return ;
					final JDialog dlg = new JDialog(XFrame.this);
					dlg.setTitle(r.get().rec.getReadName());
					dlg.setContentPane(new SAMRecordPanel(r.get().rec));
					dlg.pack();
					dlg.setVisible(true);
					}
				});
			this.drawingArea.setOpaque(true);
			this.drawingArea.setToolTipText("");
			drawingPane.add(this.drawingArea,BorderLayout.CENTER);
			
			/** scroll bar for vertical navigation */
			this.vScrollBar = new JScrollBar(JScrollBar.VERTICAL);
			this.drawingArea.addComponentListener(new ComponentAdapter() {
				public void componentResized(java.awt.event.ComponentEvent e) {
					adjustScrollBar();
					drawingArea.repaint();
				};
			});
			drawingPane.add(this.vScrollBar,BorderLayout.EAST);
			viewPane.add(drawingPane,BorderLayout.CENTER);
			
			/** ref tab */
			jTabbedPane.addTab("REF", new JScrollPane( new JTable(new SwingSequenceDictionaryTableModel(dict))));
			this.interval.addPropertyChangeListener(AE->{
				final Locatable loc=(Locatable)AE.getNewValue();
				this.jtextFieldLocation.setText(loc==null?"":loc.getContig()+":"+loc.getStart()+"-"+loc.getEnd());
				reloadReadsAndRef();
				this.drawingArea.repaint();
				});
			
			
			/** VCF tab */
			final JPanel vcfPane = new JPanel(new BorderLayout());
			final JTable variantTable = new JTable(this.swingVariantsTableModel= new SwingVariantsTableModel());
			JTable vcfInfoTable = new JTable(this.swingVCFInfoTableModel= new SwingVCFInfoTableModel());
			JTable vcfGenotypeTable = new JTable(this.swingVCFGenotypesTableModel=new SwingVCFGenotypesTableModel());
			JSplitPane split1 = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,new JScrollPane(variantTable),new JScrollPane(vcfGenotypeTable));
			JSplitPane split2 = new JSplitPane(JSplitPane.VERTICAL_SPLIT,split1,new JScrollPane(vcfInfoTable));
			vcfPane.add(split2,BorderLayout.CENTER);
			jTabbedPane.addTab("Variants", vcfPane);
			variantTable.getSelectionModel().addListSelectionListener(AE->{
				VariantContext ctx =null;
				int i = variantTable.getSelectedRow();
				if(i>=0) i = variantTable.convertRowIndexToModel(i);
				if(i>=0) {
					ctx = this.swingVariantsTableModel.getElementAt(i);
					}
				this.swingVCFInfoTableModel.setVariant(ctx);
				this.swingVCFGenotypesTableModel.setVariant(ctx);
				});
			
			
			final JMenu optionsMenu = new JMenu("Options");
			menuBar.add(optionsMenu);
		
			final ActionListener repaintAction = AE->reloadReadsAndRef();
			

			optionsMenu.add(this.jmenuShowBases=new JCheckBoxMenuItem("Show Bases", true));
			this.jmenuShowBases.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowMismatches=new JCheckBoxMenuItem("Show Mismatches", true));
			this.jmenuShowMismatches.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowInsertions=new JCheckBoxMenuItem("Show Insertions", true));
			this.jmenuShowInsertions.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowDuplicates=new JCheckBoxMenuItem("Show Duplicates", false));
			this.jmenuShowDuplicates.setToolTipText("Show/Hide Dup Reads");
			this.jmenuShowDuplicates.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowSecondary=new JCheckBoxMenuItem("Show Secondary", false));
			this.jmenuShowSecondary.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowSupplementary=new JCheckBoxMenuItem("Show Supplementary", false));
			this.jmenuShowSupplementary.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowFailQuality=new JCheckBoxMenuItem("Show Failing Quality", false));
			this.jmenuShowFailQuality.addActionListener(repaintAction);
			optionsMenu.add(new JSeparator());
			optionsMenu.add(this.jmenuShowClip=new JCheckBoxMenuItem("Show Clipping", false));
			this.jmenuShowClip.setToolTipText("Show/Hide Read Clipping");
			this.jmenuShowClip.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowCoverage=new JCheckBoxMenuItem("Show Coverage", true));
			this.jmenuShowCoverage.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowIntervalTitle=new JCheckBoxMenuItem("Show Interval Title", true));
			this.jmenuShowIntervalTitle.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowSampleTitle=new JCheckBoxMenuItem("Show Sample Title", true));
			this.jmenuShowSampleTitle.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuShowReferenceSequence=new JCheckBoxMenuItem("Show Reference Sequence", true));
			this.jmenuShowReferenceSequence.addActionListener(repaintAction);
			optionsMenu.add(this.jmenuMappingQuality=new JMenuItem("Mapping Quality ["+this.mappingQuality.orElse(0)+"]"));
			this.jmenuMappingQuality.addActionListener(AE->{
				String s= JOptionPane.showInputDialog(XFrame.this, "Mapping quality ?",mappingQuality.orElse(0));
				if(StringUtils.isBlank(s)) return;
				if(!StringUtils.isInteger(s)) {
					JOptionPane.showConfirmDialog(XFrame.this, "Not an integer: "+s);
					}
				this.mappingQuality.setValue(Integer.parseInt(s));
				});
			this.mappingQuality.addPropertyChangeListener(AE->{
				this.jmenuMappingQuality.setText("Mapping Quality ["+AE.getNewValue()+"]");
				reloadReadsAndRef();
			});
			
			
			
			


			this.jMenuHyperlinks = new JMenu("Hyperlinks");
			menuBar.add(this.jMenuHyperlinks);
			
			/* open close operations */
			this.addWindowListener(new WindowAdapter() {
				public void windowOpened(java.awt.event.WindowEvent e) {
					final SAMSequenceRecord rec = XFrame.this.dict.getSequence(0);
					interval.setValue(new SimpleInterval(rec.getContig(),1,Math.min(100, rec.getSequenceLength())));
					}
				public void windowClosed(java.awt.event.WindowEvent e)
					{
					runOnClose();
					}
				});
			}/* END CONSTRUCTOR ********************************/
		
		private Optional<Read> findReadAt(final double x,final double y) {
			for(List<Read> row: rows) {
				if(row.isEmpty()) continue;
				final Read first = row.get(0);
				if(first.getMinY() > y ) continue;
				if(first.getMaxY() < y ) continue;
				return row.stream().filter(R->R.contains(x, y)).findAny();
				}
			return Optional.empty();
			}
		
		
		private void runOnClose() {
			final DefaultListModel<BamFile> dm = (DefaultListModel<BamFile>)this.jlistBamFilesJList.getModel();
			for(int i=0;i< dm.size();i++) dm.getElementAt(i).close();
			try {
				this.referenceSequenceFile.close();
				}
			catch(final IOException err) {
				LOG.error(err);
				}
 			}
		
		private Shape createShape(double y,int position,int lengthOnReference,int arrow_type) {
			final Locatable rgn = this.interval.orElse(null);
			if(rgn==null) return null;
			if(position> rgn.getEnd()) return null;
			if(position+lengthOnReference< rgn.getStart()) return null;
			double x0 = pos2pixel(rgn,position);
			double x1 = pos2pixel(rgn,position+lengthOnReference);
			double featH = getReadRowHeight();
			double y0 = y + featH*0.05;
			double y1 = y + featH*0.95;
			double midy= y + featH*0.5;
			double da = getBasePixel()/2.0;
			switch(arrow_type) {
			case -1: {
				final GeneralPath g = new GeneralPath();
				g.moveTo(x0, midy);
				g.lineTo(x0+da, y0);
				g.lineTo(x1, y0);
				g.lineTo(x1, y1);
				g.lineTo(x0+da, y1);
				g.closePath();
				return g;
				}	
			case 1: {
				final GeneralPath g = new GeneralPath();
				g.moveTo(x1, midy);
				g.lineTo(x1-da, y0);
				g.lineTo(x0, y0);
				g.lineTo(x0, y1);
				g.lineTo(x1-da, y1);
				g.closePath();
				return g;
				}
			default:
				return new Rectangle2D.Double(x0, y0, (x1-x0), (y1-y0));
				}
			}
		
		private Color base2color(char c) {
			switch(c) {
			case 'a': case 'A': return Color.RED;
			case 't': case 'T': return Color.BLUE;
			case 'c': case 'C': return Color.GREEN;
			case 'g': case 'G': return Color.ORANGE;
			default: return Color.BLACK;
			}
		}

		private double getIntervalTitleHeight() { return jmenuShowIntervalTitle.isSelected()?12:0;}
		private double getSampleTitleHeight() { return jmenuShowSampleTitle.isSelected()?12:0;}
		private double getReferenceTitleHeight() {
			if(!jmenuShowReferenceSequence.isSelected()) return 0;
			final Locatable loc = this.interval.orElse(null);
			if(loc==null) return 0;
			final double base2pix = this.getWidth()/(double)loc.getLengthOnReference();
			if(base2pix< 5.0) return 0;
			return Math.min(12,base2pix);
			}
		private double getCoverageTrackHeight() {
			if(!jmenuShowCoverage.isSelected()) return 0;
			return 50;
			}
		
		private void paintDrawingArea(final Graphics2D g) {
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, drawingArea.getWidth(), drawingArea.getHeight());
			final Locatable rgn = this.interval.orElse(null);
			if(rgn==null) return ;
			final Hershey hershey=new Hershey();
			
			final BamFile bamFile = getSelectedBamFile();
			double y=0;
			
			
			double height = getIntervalTitleHeight();
			if(height>0.0 ) {
				g.setColor(Color.DARK_GRAY);
				final String title = new SimpleInterval(rgn).toNiceString()+" length:"+StringUtils.niceInt(rgn.getLengthOnReference())+" bp";
				hershey.paint(g,title,0,y,Math.min((double)this.drawingArea.getWidth(),title.length()*height),height*0.95);
				y+=height;
				}
			
			height = getSampleTitleHeight();
			if(height>0.0) {
				g.setColor(Color.DARK_GRAY);
				if(bamFile!=null) {
					final String title= bamFile.getName();
					hershey.paint(g,title,0,y,Math.min((double)this.drawingArea.getWidth(),title.length()*height),height*0.95);
					}
				y+=height;
				}
			
			final Function<Integer, Character> pos2ref = P->{
				int idx = P - rgn.getStart();
				if(idx<0) return 'N';
				if(idx>=referenceStr.length()) return 'N';
				return referenceStr.charAt(idx);
				};
			
			height = getReferenceTitleHeight();
			if(height>0.0) {
				//paint reference
				for(int x=rgn.getStart();x<=rgn.getEnd();++x) {
					final double x0 = pos2pixel(rgn, x);
					final double x1 = pos2pixel(rgn, x+1)-1.0;
					g.setColor(base2color(pos2ref.apply(x)));
					final char c = referenceStr.charAt(x-rgn.getStart());
					hershey.paint(g,String.valueOf(c),x0,y,(x1-x0),height*.95);
					}
				y+=height;
				}

			// plot GFF
			height = getGffTrackHeight();
			if(height>0.0) {
				final double midy = y+ height/2.0;
				g.setColor(Color.BLUE);
				for(final Gff3Feature gff3: this.gff3features ) {
					if(!gff3.getType().equals("gene")) continue;
					final double x0 = pos2pixel(rgn, Math.max(gff3.getStart(),rgn.getStart()));
					final double x1 = pos2pixel(rgn, Math.min(gff3.getEnd(),rgn.getEnd())+1);
					g.drawLine((int)x0, (int)midy,(int)x1, (int)midy);
					final String geneName = gff3.getAttributes().entrySet().stream().
							filter(KV->KV.getKey().equals("gene_name") || KV.getKey().equals("Name")).
							flatMap(KV->KV.getValue().stream()).
							findFirst().
							orElse(null);
					if(!StringUtils.isBlank(geneName)) {
						hershey.paint(g, geneName,new Rectangle2D.Double(
							x0,
							y +1,
							geneName.length()*10,
							10
							));
						}
					}
				g.setColor(Color.CYAN);
				for(final Gff3Feature gff3: this.gff3features ) {
					if(!(gff3.getType().equals("exon") || gff3.getType().equals("CDS"))) continue;
					final double x0 = pos2pixel(rgn, Math.max(gff3.getStart(),rgn.getStart()));
					final double x1 = pos2pixel(rgn, Math.min(gff3.getEnd(),rgn.getEnd())+1);
					final double h3 = gff3.equals("exon") ? 5 : 10;
					
					g.fillRect(
						(int)x0,
						(int)(midy - h3/2.0),
						(int)(x1-x0)+1,
						(int)h3
						);
					}
				y += height;
				}
			
			
			/* prepare coverage track */
			final double coverage_top = y;
			height= getCoverageTrackHeight();
			final Map<Character,int[]> base2coverage = new HashMap<>();
			if(height>1.0) {
				for(char base : new char[] {'a','c','g','t','n'}) {
					final int[] cov = new int[rgn.getLengthOnReference()];
					base2coverage.put(base, cov);
					}
				y+= height;
				}
			final double coverage_bottom = y;

			for(final VariantContext ctx:this.all_variants_list) {
				final Genotype gt = ctx.getGenotype(bamFile.sn);
				if(gt==null) continue;
				final  double x0 = pos2pixel(rgn,ctx.getStart());
				final double x1 = pos2pixel(rgn,ctx.getEnd()+1);
				g.setColor(Color.LIGHT_GRAY);
				g.fill(new Rectangle2D.Double(x0, y, (x1-x0), this.drawingArea.getHeight()-y));

				}
			
			
			final double readRowHeight = getReadRowHeight();
			
			final boolean show_bases = this.jmenuShowBases.isSelected();
			final boolean show_insertions = this.jmenuShowInsertions.isSelected();
			final boolean show_mismatches = this.jmenuShowMismatches.isSelected();

			for(int rowidx= 0;rowidx< this.rows.size();rowidx++) {
				final boolean row_is_visible=true;
				final List<Read> row = this.rows.get(rowidx);
				for(final Read read:row) {
					read.x =  (float)pos2pixel(rgn,read.getStart());
					read.y = (float)y;
					read.width = (float)pos2pixel(rgn,read.getEnd()+1) - read.x;
					read.height = (float)readRowHeight;					
					
					final Map<Integer,Integer> pos2insert_size = new HashMap<>();
					final byte[] readBases = read.rec.getReadBases();
					final Cigar cigar = read.rec.getCigar();
					// horizontal line for deletions
					if(row_is_visible && cigar.getCigarElements().stream().anyMatch(OP->OP.getOperator().equals(CigarOperator.D) || OP.equals(CigarOperator.N)) ) {
						final double x0 = pos2pixel(rgn,read.rec.getStart());
						final double x1 = pos2pixel(rgn,read.rec.getEnd());
						g.setColor(Color.DARK_GRAY);
						g.draw(new Line2D.Double(x0, y+readRowHeight/2.0, x1, y+readRowHeight/2.0));
						}
					
					
					int refpos = read.rec.getUnclippedStart();
					int readpos=0;
					int next_refpos=refpos;
					int next_readpos=readpos;
					for(int cx=0;cx<cigar.numCigarElements();++cx) {
						// should we display an arrow at start/end ?
						int arrow_type=0;
						if(cx==0 && read.rec.getReadNegativeStrandFlag()) arrow_type=-1;
						else if(cx+1==cigar.numCigarElements() && !read.rec.getReadNegativeStrandFlag()) arrow_type=1;
						
						final CigarElement ce = cigar.getCigarElement(cx);
						final CigarOperator op = ce.getOperator();
						final int len = ce.getLength();
						Shape shape=null;
						Color paper = Color.LIGHT_GRAY;
						Color pen = Color.DARK_GRAY;
						switch(op) {
							case P:break;
							case I: {
								if(row_is_visible) {
									pos2insert_size.put(refpos, len);
									}
								next_readpos+=len;
								break;
								}
							case D: case N: {
								if(row_is_visible) {
									double x0 = pos2pixel(rgn,Math.max(rgn.getStart(), refpos));								
									double x1 = pos2pixel(rgn,Math.min(rgn.getEnd(), refpos+len));
									g.setColor(Color.BLACK);
									String title = String.valueOf(len);
									double dx = Math.min(title.length()*readRowHeight,(x1-x0));
									hershey.paint(g, title, (x1-x0)/2.0-dx/2,y,dx,readRowHeight*0.95);
									}
								next_refpos+=len;
								break;
								}
							case S: case H:{
								if(row_is_visible && isUsingClip()) {
									paper = Color.YELLOW;
									shape = createShape(y,refpos,len,arrow_type);
									}
								if(op!=CigarOperator.H) next_readpos+= len;
								next_refpos+= len;
								break;
								}
							case M:case X: case EQ: {
								if(row_is_visible) {
									shape = createShape(y,refpos,len,arrow_type);
									}
								next_readpos+= len;
								next_refpos+= len;
								break;
								}
							default: throw new IllegalStateException();
							}
						
						if(row_is_visible && shape!=null && paper!=null) {
							g.setColor(paper);
							g.fill(shape);
							}
						if(row_is_visible && shape!=null && pen!=null) {
							g.setColor(pen);
							g.draw(shape);
							}
						if(
							((op.consumesReadBases() && !op.equals(CigarOperator.SOFT_CLIP)) ||
							 (op.equals(CigarOperator.SOFT_CLIP) && isUsingClip())) &&
							  readBases!=SAMRecord.NULL_SEQUENCE) {
							for(int x=0;x< len;++x) {
								if(refpos+x < rgn.getStart()) continue;
								if(refpos+x > rgn.getEnd()) break;
								char readBase = (char)readBases[readpos+x];
								char refBase = referenceStr.charAt(refpos+x-rgn.getStart());
								final boolean is_mimatch = Character.toUpperCase(readBase)!=Character.toUpperCase(refBase) && Character.toUpperCase(refBase)!='N';
								final int[] coverage = base2coverage.get(Character.toLowerCase(readBase));
								if(coverage!=null && refpos+x >= rgn.getStart() && refpos+x<=rgn.getEnd()) {
									coverage[refpos+x-rgn.getStart()]++;
									}
								if(row_is_visible) {
									if(is_mimatch && show_mismatches) {
										g.setColor(Color.YELLOW);
										double x0 = pos2pixel(rgn,refpos+x);								
										double x1 = pos2pixel(rgn,refpos+x+1);
										double y0 = y+ readRowHeight*0.05;
										double y1 = y+ readRowHeight*0.95;
										g.fill(new Rectangle2D.Double(x0, y0, (x1-x0), (y1-y0)));
										}
									if((is_mimatch && show_mismatches) || show_bases) {
										g.setColor(base2color(readBase));					
										hershey.paint(g,String.valueOf(readBase),
											pos2pixel(rgn,refpos+x),
											y,
											getBasePixel(),
											readRowHeight
											);
										}
									}
								}
							}
							
						/* plot insertions */
						if(row_is_visible && show_insertions) {
							for(Integer pos:pos2insert_size.keySet()) {
								final double x0 = pos2pixel(rgn,pos);
								g.setColor(Color.RED);
								g.fill(new Rectangle2D.Double(x0,y, 1, readRowHeight*0.95));
								}
							}
						readpos= next_readpos;
						refpos= next_refpos;
						}
					}
				y+= readRowHeight;
				}
			if(coverage_top < coverage_bottom) {
				double max_cov=0;
				double cov_height = coverage_bottom-coverage_top;
				for(int i=0;i< rgn.getLengthOnReference();i++) {
					int sum_cov=0;
					for(int[] cov:base2coverage.values()) {
						sum_cov+= cov[i];
						}
					max_cov=Math.max(sum_cov,max_cov);
					}
				for(int i=0;i< rgn.getLengthOnReference();i++) {
					final double x0 = pos2pixel(rgn, i+rgn.getStart()  );
					final double x1 = pos2pixel(rgn, i+rgn.getStart()+1);
					double y2= coverage_bottom;
					for(char base:new char[] {'a','c','g','t','n'}) {
						final double covi = base2coverage.get(base)[i];
						final double covh = covi/max_cov * cov_height;
						g.setColor(base2color(base));
						final Shape shape = new Rectangle2D.Double(x0,y2-covh,(x1-x0),covh);
						g.fill(shape);
						if(x1-x0>3) {
							g.setColor(Color.DARK_GRAY);
							g.draw(shape);
							}
						y2-=covh;
						}
					}
				}
			}

		
		private BamFile getSelectedBamFile() {
			return this.jlistBamFilesJList.getSelectedValue();
		}
		
		
		
		private void reloadReadsAndRef() {
			final Locatable rgn = this.interval.orElse(null);
			if(rgn!=null) {
				final ReferenceSequence refseq = this.referenceSequenceFile.getSubsequenceAt(
					rgn.getContig(),
					rgn.getStart(),
					rgn.getEnd()
					);
				this.referenceStr = refseq.getBaseString();
				
				final BamFile f = getSelectedBamFile();
				if(f==null) return;
				this.rows.clear();
				this.all_variants_list.clear();
				this.gff3features.clear();
				
				if(this.vcfFile!=null) {
					try(VCFReader vcfReader = VCFReaderFactory.makeDefault().open(this.vcfFile, true)) {
						try(CloseableIterator<VariantContext> iter = vcfReader.query(rgn)) {
							while(iter.hasNext()) {
								final VariantContext ctx = iter.next();
								this.all_variants_list.add(ctx);
								}
							}
						this.swingVariantsTableModel.setRows(this.all_variants_list);
						}
					catch(final IOException err) {
						LOG.error(err);
						}
					}

				if(!StringUtils.isBlank(this.gffFile)) {
					final Gff3Codec codec = new Gff3Codec(DecodeDepth.SHALLOW);
					try(TabixReader tabix= new TabixReader(this.gffFile)) {
						final TabixReader.Iterator iter = tabix.query(rgn.getContig(), rgn.getStart(), rgn.getEnd());
						
						for(;;) {
							final String line = iter.next();
							if(line==null) break;
							Gff3Feature gff3 = codec.decode(LineIterators.of(Collections.singletonList(line)));
							if(gff3==null) continue;
							if(!(gff3.getType().equals("gene") || gff3.getType().equals("exon") || gff3.getType().equals("CDS"))) continue;
							this.gff3features.add(gff3);
							}
						}
					catch(final IOException|TribbleException err) {
						LOG.error(err);
						}
					}
				
				
				final List<Read> buffer= new ArrayList<>(fetchReads());
				while(!buffer.isEmpty()) {
					final Read read = buffer.remove(0);
					int y = 0;
					for(y=0;y< this.rows.size();++y) {
						final List<Read> row = this.rows.get(y);
						final Read last = row.get(row.size()-1);
						if(last.getEnd() + 1  >= read.getStart()) continue;
						row.add(read);
						break;
						}
					if(y==this.rows.size()) {
						final List<Read> row = new Vector<>();
						row.add(read);
						this.rows.add(row);
						}
					}
				}
			else
				{
				this.referenceStr = "";
				this.rows.clear();
				this.gff3features.clear();
				}
			adjustScrollBar();
			this.drawingArea.repaint();
			}

		public void adjustScrollBar() {
						
			int n_visible = Math.max(1, (int)Math.floor((this.drawingArea.getHeight() - getMarginTop())/getReadRowHeight()));
			if(n_visible<=rows.size()) {
				this.vScrollBar.setValues(0, 0, 0, 0);
				}
			else
				{
				this.vScrollBar.setValues(
						0,n_visible,1, rows.size()-n_visible
						);
				}
			}
		
		final double getMarginTop() {
			return getReferenceTitleHeight() +
					getSampleTitleHeight() +
					getIntervalTitleHeight() +
					getCoverageTrackHeight() +
					getGffTrackHeight()
					;
			}
		
		final double getGffTrackHeight() {
			if(StringUtils.isBlank(this.gffFile)) return 0.0;
			return 40;
			}
		
		final double getReadRowHeight() {
			final double remain_height = (drawingArea.getHeight()-5) - getMarginTop();
			final double row_height = remain_height/ this.rows.size();
			return Math.min(row_height, 30);
			}
		
		final double pos2pixel(final Locatable loc,final int coord) {
			return (coord - loc.getStart())/((double)loc.getLengthOnReference())*this.drawingArea.getWidth();
			}
		final double getBasePixel() {
			final Locatable loc = this.interval.orElse(null);
			return (1.0/loc.getLengthOnReference())*this.drawingArea.getWidth();
			}
			
		final Optional<SimpleInterval> getUserInterval() {
			final String s = this.jtextFieldLocation.getText().trim();
			if(StringUtil.isBlank(s)) return Optional.empty();
			Optional<SimpleInterval> ret;

			try {
				ret = new IntervalParser(this.dict).
					apply(s);
				}
			catch(Throwable err) {
				return Optional.empty();
				}
		
			return ret;
			}
	
		/** fetch Read for current bam, and current interval **/
		private List<Read> fetchReads() {
			final Locatable loc = this.interval.orElse(null);
			if(loc==null) return Collections.emptyList();
			final BamFile bf = getSelectedBamFile();
			if(this.interval==null || bf==null) return Collections.emptyList();
			bf.open(this.referenceFaidx);
			final SamRecordFilter srf = createSamRecordFilter();
			final List<Read> L;
			try(CloseableIterator<SAMRecord> iter = bf.sr.query(
				loc.getContig(),
				loc.getStart(),
				loc.getEnd(),
				false)) {
				L =  iter.stream().
						filter(R->!srf.filterOut(R)).
						map(R->new Read(R)).
						sorted((A,B)->Integer.compare(A.getStart(),B.getStart())).
						collect(Collectors.toList());
				}
			bf.close();
			return L;
			}
		
		SamRecordFilter createSamRecordFilter() {
			final List<SamRecordFilter> filters = new ArrayList<>();
			filters.add(new AlignedFilter(true));
			if(!this.jmenuShowDuplicates.isSelected()) {
				filters.add(new DuplicateReadFilter());
				}
			if(!this.jmenuShowFailQuality.isSelected()) {
				filters.add(new FailsVendorReadQualityFilter());
				}
			if(!this.jmenuShowSecondary.isSelected()) {
				filters.add(new SecondaryAlignmentFilter());
				}
			if(!this.jmenuShowSupplementary.isSelected()) {
				filters.add(new SamRecordFilter() {
					@Override
					public boolean filterOut(SAMRecord first, SAMRecord second) {
						return filterOut(first) || filterOut(second);
					}@Override
					public boolean filterOut(SAMRecord record) {
						return record.getSupplementaryAlignmentFlag();
					}
				});
				}
			if(this.mappingQuality.isPresent()) {
				filters.add(new MappingQualityFilter(this.mappingQuality.orElse(0)));
				}
			return new AggregateFilter(filters);
			}
		
		private boolean isUsingClip() {
			return this.jmenuShowClip.isSelected();
			}
		}
	
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			if(IOUtils.unrollPaths(args).isEmpty()) {
				LOG.error("no BAM/CRAM was provided");
				return -1;
				}
			
			JFrame.setDefaultLookAndFeelDecorated(true);
		
			
			final List<Path> paths = IOUtils.unrollPaths(args);
			if(paths.isEmpty()) {
				LOG.warn("No File was provided");
				return -1 ;
				}
			final List<BamFile> bamFiles = new ArrayList<>(paths.size());
			for(final Path p:paths) {
				final BamFile bf = new BamFile(p);
				bf.open(this.referenceFile);
				bf.sn = bf.header.getReadGroups().
					stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtil.isBlank(S)).
					findFirst().
					orElse(IOUtils.getFilenameWithoutCommonSuffixes(p));
				bf.close();
				bamFiles.add(bf);
				}
			
			final XFrame frame = new XFrame(this.referenceFile,bamFiles,this.variantFile,this.gffFile);
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
	
	public static void main(String[] args) {
		new SwingBamView().instanceMain(args);
	}

}
