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
import java.awt.Component;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.net.URI;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.TableModel;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gff3.SwingGff3TableModel;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.UrlSupplier;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;
import com.github.lindenb.jvarkit.util.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.vcf.predictions.SmooveGenesParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffLofNmdParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffLofNmdParser.Prediction;
import com.github.lindenb.jvarkit.variant.swing.SwingAllelesTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingAnnPredictionTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingBcsqPredictionTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingPedigreeTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingTrioTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFFilterHeaderLineTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFFormatHeaderLineTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFGenotypesTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFInfoHeaderLineTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFInfoTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVariantsTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVepPredictionTableModel;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.tribble.readers.TabixIteratorLineReader;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC

## Example

```
java -jar dist/swingvcfview.jar input.vcf.gz
```

## Screenshot

 * https://twitter.com/yokofakun/status/1392173413079322625

END_DOC
 */
@Program(name="swingvcfview",
description="VCFviewer using Java Swing UI",
keywords={"vcf","visualization","swing"},
creationDate="20210503",
modificationDate="20210503",
generate_doc=true
)
public class SwingVcfView extends Launcher
	{
	private static final Logger LOG = Logger.build(SwingVcfView.class).make();
	@Parameter(names={"-r","--regions","--interval"},description="default interval region on opening")
	private String defaultRegion="";
	@Parameter(names={"--limit"},description="Limit number of variants. Ignore if < 0")
	private int limit_number_variant = -1;
	@Parameter(names={"--pedigree"},description=PedigreeParser.OPT_DESC)
	private Path pedigreePath = null;
	@Parameter(names={"--gff","--gff3"},description="GFF3 file Path used to find intervals by gene name")
	private Path gffPath = null;

	
	@SuppressWarnings("serial")
	private static class XFrame extends JFrame {
		final SAMSequenceDictionary dict;
		final VCFReader vcfReader;
		final JTextField jtextFieldLocation;
		final JTextField jtextFieldGEXL;
		final SwingVariantsTableModel swingVariantsTableModel;
		final int limit_number_variant;
		final JTable variantTable;
		final JTextField txtFieldSample;
		final Pedigree pedigree;
		final SwingVCFGenotypesTableModel swingVCFGenotypesTableModel;
		final SwingVCFInfoTableModel swingInfoTableModel;
		final SwingAnnPredictionTableModel swingAnnPredictionTableModel;
		final SwingVepPredictionTableModel swingVepPredictionTableModel;
		final SwingBcsqPredictionTableModel swingBcsqPredictionTableModel;
		final SnpEffNmdLOfTableModel lofSnpEffTableModel;
		final SnpEffNmdLOfTableModel nmdSnpEffTableModel;
		final SmooveGeneTableModel smooveGeneTableModel;
		final SwingAllelesTableModel swingAllelesTableModel;
		final SwingTrioTableModel swingTrioTableModel;
		final GenotypeTypeTableModel genotypeTypeTableModel;
		final Map<VariantContext.Type, JCheckBox> variantType2cbox = new HashMap<>();
		final JCheckBox filteredCtxcbox;
		final Map<GenotypeType, JCheckBox> genotypeType2cbox = new HashMap<>();
		final DefaultListModel<String> filterListModel;
		final GffTableModel gffTableModel;
		final Path gffPath;
		final JMenu menuOpenBrowser;
		
		XFrame(final Path vcfPath,String defaultLoc,
				final int limit_number_variant,
				final Pedigree pedigree,
				final Path gffPath) {
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			setTitle(SwingVcfView.class.getSimpleName()+": "+ vcfPath.getFileName().toString());
			this.limit_number_variant = limit_number_variant;
			this.pedigree = pedigree;
			this.gffPath = gffPath;
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
			if(this.limit_number_variant>=0) {
				topPane.add(new JLabel("(limit: "+StringUtils.niceInt(this.limit_number_variant)+")"));
				topPane.add(new JSeparator());
			}
			for(VariantContext.Type vtype: VariantContext.Type.values()) {
				final JCheckBox cb = new JCheckBox(vtype.name(),true);
				cb.setToolTipText("Filter Variants of type "+ vtype);
				topPane.add(cb);
				cb.addActionListener(AE->refreshInterval());
				this.variantType2cbox.put(vtype, cb);
				}
			this.filteredCtxcbox = new JCheckBox("FILTER",true);
			this.filteredCtxcbox.setToolTipText("FILTER-ed Variants");
			topPane.add(this.filteredCtxcbox);
			this.filteredCtxcbox.addActionListener(AE->refreshInterval());

			
			topPane.add(new JSeparator());
			label= new JLabel("JEXL:");
			this.jtextFieldGEXL = new JTextField(10);
			this.jtextFieldGEXL.setToolTipText("Filter variants using a JEXL expression.");
			label.setLabelFor(this.jtextFieldGEXL);
			topPane.add(label);
			topPane.add(this.jtextFieldGEXL);
			this.jtextFieldGEXL.addActionListener(AE->refreshInterval());
			
			
			
			final Panel botPane = new Panel(new FlowLayout(FlowLayout.TRAILING));
			mainPane.add(botPane, BorderLayout.SOUTH);

						
			final JTabbedPane  tabPane = new JTabbedPane();
			mainPane.add(tabPane,BorderLayout.CENTER);
			
			final JPanel pane2 = new JPanel(new BorderLayout(5,5));
			this.swingVariantsTableModel = new SwingVariantsTableModel();
			this.variantTable =  new JTable(this.swingVariantsTableModel);
			this.variantTable.addMouseListener(new MouseAdapter() {
				@Override
				public void mousePressed(java.awt.event.MouseEvent e) {
					if(!e.isPopupTrigger()) return;
			        int rowindex = variantTable.getSelectedRow();
			        if(rowindex<0) return;
			        rowindex = variantTable.convertRowIndexToModel(rowindex);
			        if(rowindex<0 || rowindex>= swingVariantsTableModel.getRowCount()) return;
			        final VariantContext vc  = swingVariantsTableModel.getElementAt(rowindex);
					final JPopupMenu pop=new JPopupMenu();
					fillMenuForVariant(pop,vc);
					pop.show(e.getComponent(), e.getX(), e.getY());
					}
				});
			
			this.swingVCFGenotypesTableModel = new SwingVCFGenotypesTableModel(this.pedigree);
			final JTable genotypeTable =  new JTable(this.swingVCFGenotypesTableModel);
			genotypeTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			
			final JPanel gtPane = new JPanel(new BorderLayout(1,1));
			final JPanel topPane2 = new JPanel(new FlowLayout(FlowLayout.LEADING));
			gtPane.add(topPane2,BorderLayout.NORTH);
			
			for(GenotypeType gttype: GenotypeType.values()) {
				final JCheckBox cb = new JCheckBox(gttype.name(),true);
				cb.setToolTipText("Filter Genotype of type "+gttype);
				topPane2.add(cb);
				cb.addActionListener(AE->refreshTables());
				this.genotypeType2cbox.put(gttype, cb);
				}
			topPane2.add(new JSeparator());
			label = new JLabel("Sample:");
			topPane2.add(label);
			this.txtFieldSample = new JTextField(5);
			label.setLabelFor(this.txtFieldSample);
			topPane2.add(this.txtFieldSample);
			this.txtFieldSample.addActionListener(AE->refreshTables());
			topPane2.add(new JSeparator());
			JButton save = new JButton("*");
			save.addActionListener(AE->saveTable("Genotypes",genotypeTable));
			save.setPreferredSize(new Dimension(5,5));
			topPane2.add(save);
			
			gtPane.add(new JScrollPane(genotypeTable),BorderLayout.CENTER);
			
			final JSplitPane split1 = new JSplitPane(
					JSplitPane.HORIZONTAL_SPLIT,
					wrapTable("Variants",variantTable),
					gtPane
					);
			split1.setResizeWeight(0.5);
			
			this.swingInfoTableModel = new SwingVCFInfoTableModel();
			
			this.swingAnnPredictionTableModel = new SwingAnnPredictionTableModel(header);
			this.swingVepPredictionTableModel = new SwingVepPredictionTableModel(header);
			this.swingBcsqPredictionTableModel = new SwingBcsqPredictionTableModel(header);
			this.genotypeTypeTableModel = new GenotypeTypeTableModel();
			this.swingAllelesTableModel = new SwingAllelesTableModel();
			
			final MouseAdapter showBrowser = new MouseAdapter() {
				@Override
				public void mousePressed(final MouseEvent e) {
					if(!e.isPopupTrigger()) return;
					if(!Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.BROWSE)) return;
					final Component component = e.getComponent();
					if(!(component instanceof JTable)) return;
					final JTable table = JTable.class.cast(component);
					final TableModel m = table.getModel();
			        int rowindex = table.getSelectedRow();
			        if(rowindex<0) return;
			        rowindex = table.convertRowIndexToModel(rowindex);
			        if(rowindex<0 || rowindex>= m.getRowCount()) return;
			        int colindex = table.getSelectedColumn();
			        colindex = table.convertColumnIndexToModel(colindex);
			        if(colindex<0 || colindex>= m.getColumnCount()) return;
	
			        Object o = m.getValueAt(rowindex, colindex);
			        if(o==null || !(o instanceof String)) return;
			        final String s= (String)o;
			        final String colName = m.getColumnName(colindex);
			        if(StringUtils.isBlank(s)) return;
					final JPopupMenu pop=new JPopupMenu();
					final UrlSupplier urlSupplier = new UrlSupplier(XFrame.this.dict);
					Stream.concat(
						urlSupplier.of(s).stream(),
						urlSupplier.of(colName,s).stream()).
						collect(Collectors.toSet()).
						stream().
						forEach(U->{
						final AbstractAction action = new AbstractAction(U.getLabel())
							{
							@Override
							public void actionPerformed(final ActionEvent e)
								{
								try {
									Desktop.getDesktop().browse(new URI(U.getUrl()));
									}
								catch(final Throwable err) {
									ThrowablePane.show(table, err);
									}
								}
							};
						action.putValue(AbstractAction.LONG_DESCRIPTION, U.getLabel());
						action.putValue(AbstractAction.SHORT_DESCRIPTION, U.getLabel());
						action.putValue(AbstractAction.NAME, U.getLabel());
						final JMenuItem mi = new JMenuItem(action);
						pop.add(mi);
						});
	
					pop.show(e.getComponent(), e.getX(), e.getY());
					}
				};
			
			final JTabbedPane tabbed2 = new JTabbedPane();
			tabbed2.addTab("INFO", wrapTable("INFO",new JTable(this.swingInfoTableModel)));
			tabbed2.addTab("FILTER", new JScrollPane(new JList<>(this.filterListModel = new DefaultListModel<>())));
			tabbed2.addTab("Types", wrapTable("Types",new JTable(this.genotypeTypeTableModel)));
			tabbed2.addTab("Alleles",wrapTable("Alleles",new JTable(this.swingAllelesTableModel)));
			final JTable snpEffTable  = new JTable(this.swingAnnPredictionTableModel);
			snpEffTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			snpEffTable.addMouseListener(showBrowser);
			tabbed2.addTab("SnpEff:ANN",wrapTable("ANN",snpEffTable));
			final JTable vepTable = new JTable(this.swingVepPredictionTableModel);
			vepTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			vepTable.addMouseListener(showBrowser);
			
			tabbed2.addTab("SnpEff:LOF",wrapTable("LOF",new JTable(this.lofSnpEffTableModel=new SnpEffNmdLOfTableModel(SnpEffLofNmdParser.LOF_TAG, header))));
			tabbed2.addTab("SnpEff:NMD",wrapTable("NMD",new JTable(this.nmdSnpEffTableModel=new SnpEffNmdLOfTableModel(SnpEffLofNmdParser.NMD_TAG, header))));
			tabbed2.addTab("VEP",wrapTable("VEP",vepTable));
			final JTable smooveTable = new JTable(this.smooveGeneTableModel=new SmooveGeneTableModel(header));
			vepTable.addMouseListener(showBrowser);
			tabbed2.addTab("Smoove",wrapTable("Smoove",smooveTable));
			if(this.gffPath==null) {
				this.gffTableModel = null;
				}	
			else {
				tabbed2.addTab("GTF",wrapTable("GTF",new JTable(this.gffTableModel=new GffTableModel(this.gffPath))));
			}
			
			final JTable bcfTable  = new JTable(this.swingBcsqPredictionTableModel);
			bcfTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			tabbed2.addTab("bcftools:BCSQ",wrapTable("BCSQ", bcfTable));
			if(this.pedigree!=null && !this.pedigree.getTrios().isEmpty()) {
				this.swingTrioTableModel = new SwingTrioTableModel(this.pedigree);
				tabbed2.addTab("TRIOS",new JScrollPane(new JTable(this.swingTrioTableModel)));

			} else {
				this.swingTrioTableModel = null;
			}
			
			
			final JSplitPane split2 = new JSplitPane(
					JSplitPane.VERTICAL_SPLIT,
					split1,
					tabbed2
					);
			split2.setResizeWeight(0.5);
			
			pane2.add(split2,BorderLayout.CENTER);
			this.variantTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
			this.variantTable.getSelectionModel().addListSelectionListener(E->refreshTables());
			
			tabPane.addTab("Main", pane2);
			
			tabPane.addTab("INFO", wrapTable("INFO", new JTable(new SwingVCFInfoHeaderLineTableModel(header))));
			tabPane.addTab("FORMAT",wrapTable("FILTER",new JTable(new SwingVCFFormatHeaderLineTableModel(header))));
			tabPane.addTab("FILTER",wrapTable("FILTER",new JTable(new SwingVCFFilterHeaderLineTableModel(header))));
			if(this.pedigree!=null) {
				tabPane.addTab("Pedigree", wrapTable("Pedigree", new JTable(new SwingPedigreeTableModel(this.pedigree))));
			}
			
			
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
			this.menuOpenBrowser = new JMenu("Browse");
			menuBar.add(this.menuOpenBrowser);
			}
		
		
		
		private void fillMenuForVariant(final JComponent menu,final VariantContext ctx) {
			if(ctx==null) return;
			if(!Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.BROWSE)) return;
			final UrlSupplier urlSupplier = new UrlSupplier(this.dict);
			urlSupplier.of(ctx).stream().forEach(U->{
				final AbstractAction action = new AbstractAction(U.getLabel())
					{
					@Override
					public void actionPerformed(final ActionEvent e)
						{
						try {
							Desktop.getDesktop().browse(new URI(U.getUrl()));
							}
						catch(final Throwable err) {
							ThrowablePane.show(menuOpenBrowser, err);
							}
						}
					};
				action.putValue(AbstractAction.LONG_DESCRIPTION, U.getLabel());
				action.putValue(AbstractAction.SHORT_DESCRIPTION, U.getLabel());
				action.putValue(AbstractAction.NAME, U.getLabel());
				final JMenuItem mi = new JMenuItem(action);
				menu.add(mi);
				});
			}
		
		
		/** get interval using GTF file */
		final Optional<SimpleInterval> getGtfInterval(final String s) {
			if(StringUtil.isBlank(s)) return Optional.empty();
			if(this.gffPath==null) return Optional.empty();
			
			try(BufferedReader br=IOUtils.openPathForBufferedReading(this.gffPath)) {
				final Gff3Codec codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
				final LineIterator lr = new LineIteratorImpl(new SynchronousLineReader(br));
				while(!codec.isDone(lr)) {
					final Gff3Feature G = codec.decode(lr);
					if(G==null) continue;
					boolean flag = false;
					if(G.getType().equals("gene")) {
						if(s.equals(G.getAttribute("gene_name"))) flag =  true;
						else if(s.equals(G.getAttribute("gene_id"))) flag = true;
						else if(s.equals(G.getAttribute("ID"))) flag = true;
						}
					else if(G.getType().equals("transcript")) {
						if(s.equals(G.getAttribute("transcript_name"))) flag = true;
						else if(s.equals(G.getAttribute("transcript_id"))) flag = true;
						else if(s.equals(G.getAttribute("ID"))) flag = true;
						}
					if(!flag) continue;
					return Optional.of(new SimpleInterval(G));
					}
				return Optional.empty();
				}
			catch(final Throwable err) {
				return Optional.empty();
				}
			}
		
		final Optional<SimpleInterval> getUserInterval() {
			final String s = this.jtextFieldLocation.getText();
			if(StringUtil.isBlank(s)) return Optional.empty();
			Optional<SimpleInterval> ret;

			try {
				ret = IntervalParserFactory.newInstance().
					dictionary(this.dict).
					make().
					apply(s);
				}
			catch(Throwable err) {
				return Optional.empty();
				}
			if(ret.isPresent()) return ret;
			ret  = getGtfInterval(s);
			if(ret.isPresent()) {
				this.jtextFieldLocation.setText(ret.get().toString());
				}
			return ret;
			}
				
		private void refreshInterval() {
			final List<VariantContext> L;
			final Optional<SimpleInterval> location = getUserInterval();
			if(!location.isPresent()) {
				L = Collections.emptyList();
				}
			else
				{
				final String jexlExpr = this.jtextFieldGEXL.getText();
				try(CloseableIterator<VariantContext> r=this.vcfReader.query(location.get())) {
					Stream<VariantContext> st = r.stream();
					st=st.filter(V->{
						final VariantContext.Type vtype = V.getType();
						if(!this.variantType2cbox.get(vtype).isSelected()) return false;
						if(V.isFiltered() && !this.filteredCtxcbox.isSelected()) return false;
						return true;
						});
					if(!StringUtils.isBlank(jexlExpr)) {
						try {
						final List<JexlVCMatchExp> exps= VariantContextUtils.initializeMatchExps(
							Collections.singletonList("CUSTOM_JEXL_FILTER"),
							Collections.singletonList(jexlExpr)
							);
							st = st.filter( V-> VariantContextUtils.match(V,exps.get(0)));
						} catch(Throwable err) {
							st = st.limit(0);
							LOG.error(err);
						}
					}
					
					if(this.limit_number_variant>=0) {
						st = st.limit(this.limit_number_variant);
						}
		
					L = st.collect(Collectors.toList());	
					}
				}
			this.swingVariantsTableModel.setRows(L);
			refreshTables();
		}
		
		private void refreshTables() {

		final int i = variantTable.getSelectedRow();
		final VariantContext ctx =  i>=0? swingVariantsTableModel.getElementAt(i):null;
		
		this.menuOpenBrowser.removeAll();
		fillMenuForVariant(this.menuOpenBrowser,ctx);
		this.filterListModel.clear();
		this.swingAnnPredictionTableModel.setVariant(ctx);
		this.swingVepPredictionTableModel.setVariant(ctx);
		this.swingBcsqPredictionTableModel.setVariant(ctx);
		this.swingAllelesTableModel.setVariant(ctx);
		this.genotypeTypeTableModel.setVariant(ctx);
		this.lofSnpEffTableModel.setVariant(ctx);
		this.nmdSnpEffTableModel.setVariant(ctx);
		this.smooveGeneTableModel.setVariant(ctx);
		if(this.gffTableModel!=null) {
			this.gffTableModel.setVariant(ctx);
		}
		if(this.swingTrioTableModel!=null) {
			this.swingTrioTableModel.setVariant(ctx);
		}
		
		
		final List<Genotype> genotypes;
		if(ctx==null) {
			genotypes = Collections.emptyList();
			}
		else {
			
			for(final String flt: ctx.getFilters()) {
				this.filterListModel.addElement(flt);
				}
			
			Pattern regexSample = null;
			if(!StringUtils.isBlank(txtFieldSample.getText())) {
				try {
					regexSample = Pattern.compile(txtFieldSample.getText());
					}
				catch(Throwable err) {
					LOG.error(err);
					regexSample = null;
					}
				}
			final Pattern finalRegexSample = regexSample;
			genotypes = ctx.getGenotypes().stream().
				filter(G->{
					if(finalRegexSample!=null &&
						!finalRegexSample.matcher(G.getSampleName()).matches()) {
						return false;
						}
					final GenotypeType gtype = G.getType();
					if(!this.genotypeType2cbox.get(gtype).isSelected()) return false;
					return true;
					}).
				collect(Collectors.toList());
			}
		
		
		this.swingVCFGenotypesTableModel.setRows(genotypes);
		final Map<String,Object> hash = ctx==null?
				new HashMap<>():
				new HashMap<>(ctx.getAttributes())
				;
		
		this.swingInfoTableModel.setAttributes(hash);
		}
		
	private JPanel wrapTable(final String title,final JTable t) {
		final JPanel pane = new JPanel(new BorderLayout(1,1));
		final JPanel top = new JPanel(new FlowLayout(FlowLayout.TRAILING,1,1));
		final JButton save = new JButton("[+]");
		save.setPreferredSize(new Dimension(7,7));
		save.setToolTipText("Save Table "+title+"...");
		save.addActionListener(A->saveTable(title,t));
		top.add(save);
		
		if(!StringUtils.isBlank(title)) {
			pane.setBorder(BorderFactory.createTitledBorder(title));
			}
		pane.add(top,BorderLayout.NORTH);
		pane.add(new JScrollPane(t),BorderLayout.CENTER);
		
		return pane;
	}
		
	private void saveTable(final String title,final JTable jtable) {
		if(jtable==null) return;
		try {
			final JFileChooser fc = new JFileChooser();
			if(!StringUtils.isBlank(title)) {
				fc.setSelectedFile(new File(title.replaceAll("[ \t]", "")+".tsv"));
				}
			if(fc.showSaveDialog(jtable)!=JFileChooser.APPROVE_OPTION) return;
			final File f = fc.getSelectedFile();
			if(f.exists() && JOptionPane.showConfirmDialog(jtable, "File "+f.getName()+" exists. Overwrite ?","Overwrite",JOptionPane.OK_CANCEL_OPTION)!=JOptionPane.OK_OPTION) return;
			try(PrintWriter pw = IOUtils.openFileForPrintWriter(f)) {
				final TableModel tm = jtable.getModel();
				for(int x=0;x< tm.getColumnCount();x++) {
					pw.print(x>0?"\t":"");
					pw.print(tm.getColumnName(x));
					}
				pw.println();
				for(int y=0;y< tm.getRowCount();y++) {
					for(int x=0;x< tm.getColumnCount();x++) {
						pw.print(x>0?"\t":"");
						pw.print(tm.getValueAt(y,x));
						}
					pw.println();
					}
				pw.flush();
				}
			}
		catch(Throwable err) {
			ThrowablePane.show(jtable, err);
			}
		}
	

	@SuppressWarnings("serial")
	private static class SnpEffNmdLOfTableModel extends AbstractGenericTable<SnpEffLofNmdParser.Prediction> {
		final SnpEffLofNmdParser parser;
		SnpEffNmdLOfTableModel(final String tag,VCFHeader header) {
			this.parser = new SnpEffLofNmdParser(tag,header);
 			}
		
		public void setVariant(final VariantContext ctx) {
			this.setRows(this.parser.parse(ctx));
		}
		
		@Override
		public int getColumnCount() { return 4;}
		@Override
		public String getColumnName(int column) {
			switch(column) {
				case 0: return "GeneName";
				case 1: return "GeneID";
				case 2: return "Number of Transcripts";
				case 3: return "% Transcripts Affected";
				default: throw new IllegalArgumentException();
				}
			}
		@Override
		public Class<?> getColumnClass(int columnIndex) {
			switch(columnIndex) {
				case 0: return String.class;
				case 1: return String.class;
				case 2: return Integer.class;
				case 3: return Float.class;
				default: throw new IllegalArgumentException();
				}
			}
		@Override
		public Object getValueOf(Prediction o, int columnIndex) {
			switch(columnIndex) {
				case 0: return o.getGeneName();
				case 1: return o.getGeneId();
				case 2: return o.getNumberOfTranscripts();
				case 3: return o.getPercentOfTranscriptsAffected();
				default: throw new IllegalArgumentException();
				}
			}
		}
	
	private static class GenotypeCount {
		GenotypeType gt;
		long count;
		int fraction;
		}
	@SuppressWarnings("serial")
	private static class GenotypeTypeTableModel extends AbstractGenericTable<GenotypeCount> {
		@Override
		public int getColumnCount() { return 3;}
		
		public void setVariant(final VariantContext ctx) {
			final List<GenotypeCount> counts = new ArrayList<>();
			if(ctx!=null && ctx.hasGenotypes()) {
				final Counter<GenotypeType> counter = new Counter<>();
				for(Genotype gt: ctx.getGenotypes()) {
					counter.incr(gt.getType());
				}
				for(GenotypeType gt:counter.keySet()) {
					final GenotypeCount c = new GenotypeCount();
					c.gt = gt;
					c.count = counter.count(gt);
					c.fraction = (int)(100.0*(counter.count(gt)/(double)counter.getTotal()));
					counts.add(c);
				}
			}
			setRows(counts);
		}
		
		@Override
		public Class<?> getColumnClass(int columnIndex)
			{
			switch(columnIndex) {
				case 0: return String.class;
				case 1: return Long.class;
				case 2: return Integer.class;
				default: throw new IllegalArgumentException();
				}
			}
		@Override
		public Object getValueOf(GenotypeCount F, int columnIndex)
			{
			switch(columnIndex) {
				case 0: return F.gt.name();
				case 1: return F.count;
				case 2: return F.fraction;
				default: throw new IllegalArgumentException();
				}
			}
		@Override
		public String getColumnName(int columnIndex)
			{
			switch(columnIndex) {
				case 0: return "Type";
				case 1: return "Count";
				case 2: return "Percent";
				default: throw new IllegalArgumentException();
				}
			}
		}
	
	
	@SuppressWarnings("serial")
	private static class SmooveGeneTableModel extends AbstractGenericTable<SmooveGenesParser.Prediction> {
		final SmooveGenesParser parser;
		SmooveGeneTableModel(final VCFHeader header) {
			this.parser = new SmooveGenesParser(header);
 			}
		
		public void setVariant(final VariantContext ctx) {
			this.setRows(this.parser.parse(ctx));
		}
		
		@Override
		public int getColumnCount() { return 4;}
		@Override
		public String getColumnName(int column) {
			switch(column) {
				case 0: return "GeneName";
				case 1: return "Feature";
				case 2: return "Features Count";
				case 3: return "Bases Count";
				default: throw new IllegalArgumentException();
				}
			}
		@Override
		public Class<?> getColumnClass(int columnIndex) {
			switch(columnIndex) {
				case 0: return String.class;
				case 1: return String.class;
				case 2: return Integer.class;
				case 3: return Integer.class;
				default: throw new IllegalArgumentException();
				}
			}
		@Override
		public Object getValueOf(SmooveGenesParser.Prediction o, int columnIndex) {
			switch(columnIndex) {
				case 0: return o.getGeneName();
				case 1: return o.getFeature();
				case 2: return o.getFeaturesCount();
				case 3: return o.getBasesCount();
				default: throw new IllegalArgumentException();
				}
			}
		}

	@SuppressWarnings("serial")
	private static class GffTableModel extends SwingGff3TableModel {
		final Path gffFile;
		GffTableModel(final Path gffFile) {
			this.gffFile = gffFile;
 			}
		
		public void setVariant(final VariantContext ctx) {
			if(ctx==null || this.gffFile==null) {
				setRows(Collections.emptyList());
				return;
				}
			final List<Gff3Feature> lines= new ArrayList<>();
			try(TabixReader tbr = new TabixReader(this.gffFile.toString())) {
				final Gff3Codec codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
				final ContigNameConverter ctgConverter = ContigNameConverter.fromContigSet(tbr.getChromosomes());
				final String ctg = ctgConverter.apply(ctx.getContig());
				if(!StringUtils.isBlank(ctg)) {
				final TabixReader.Iterator iter0 = tbr.query(ctg, ctx.getStart(), ctx.getEnd());
				final TabixIteratorLineReader iter1 = new TabixIteratorLineReader(iter0);
				final LineIterator iter = new LineIteratorImpl(iter1);
				while(!codec.isDone(iter)) {
						final Gff3Feature feature = codec.decode(iter);
						if(feature==null) continue;
						lines.add(feature);
						}
					}
				}
			catch(final Throwable err) {
				lines.clear();
				}
			setRows(lines);
			}
		}
	}
	
	@Override
	public int doWork(final List<String> args)
		{
		try {
			final List<Path> all_vcf_paths = new ArrayList<>(IOUtils.unrollPaths(args));
			if(all_vcf_paths.isEmpty()) {
				final JFileChooser jfc = new JFileChooser();
				jfc.setMultiSelectionEnabled(true);
				jfc.setFileFilter(new FileFilter()
					{
					@Override
					public String getDescription()
						{
						return "Indexed Variant File";
						}
						
					@Override
					public boolean accept(final File f)
						{
						if(f.isDirectory()) return true;
						if(!f.canRead()) return false;
						if(FileExtensions.VCF_LIST.stream().noneMatch(X->f.getName().endsWith(X))) return false;
						File idx =  Tribble.tabixIndexFile(f);
						if( idx.exists()) return true;
						idx = Tribble.indexFile(f);
						if( idx.exists()) return true;
						idx = new File(f.getParentFile(),f.getName()+ FileExtensions.CSI);
						if( idx.exists()) return true;
						return false;
						}
					});
				if(jfc.showOpenDialog(null)!=JFileChooser.APPROVE_OPTION) return -1;
				final File[] sel = jfc.getSelectedFiles();
				if(sel==null || sel.length==0) return -1;
				Arrays.asList(sel).stream().map(F->F.toPath()).forEach(P->all_vcf_paths.add(P));
				}
			
			
			for(final Path path: all_vcf_paths) {
				IOUtil.assertFileIsReadable(path);
				}
					
			final Pedigree ped;
			if(this.pedigreePath!=null) {
				ped = new PedigreeParser().parse(this.pedigreePath);
				}
			else
				{
				ped = null;
				}
			
			JFrame.setDefaultLookAndFeelDecorated(true);
			final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

			
			SwingUtilities.invokeAndWait(()->{
				int dx = 0;
				for(final Path vcfPath: all_vcf_paths) {
					final XFrame frame = new XFrame(vcfPath,defaultRegion,this.limit_number_variant,ped,this.gffPath);
					frame.setBounds(50 + dx , 50 + dx, screen.width-100, screen.height-100);
					frame.setVisible(true);
					dx += 2;
					}
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
