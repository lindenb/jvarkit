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
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
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
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;
import com.github.lindenb.jvarkit.util.swing.ThrowablePane;
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
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
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
	private Path pedigreePath;

	
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
		final SwingTrioTableModel swingTrioTableModel;
		final GenotypeTypeTableModel genotypeTypeTableModel;
		final Map<VariantContext.Type, JCheckBox> variantType2cbox = new HashMap<>();
		final Map<GenotypeType, JCheckBox> genotypeType2cbox = new HashMap<>();
		JTable lastSelectedJtable = null;
		
		XFrame(final Path vcfPath,String defaultLoc,
				final int limit_number_variant,
				final Pedigree pedigree) {
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			setTitle(SwingVcfView.class.getSimpleName());
			this.limit_number_variant = limit_number_variant;
			this.pedigree = pedigree;
			final VCFReaderFactory vcfReaderFactory = VCFReaderFactory.makeDefault();
			vcfReader = vcfReaderFactory.open(vcfPath, true);
			final VCFHeader header = vcfReader.getHeader();
			
			this.dict = header.getSequenceDictionary();
			
			final JPanel mainPane = new JPanel(new BorderLayout(5, 5));
			mainPane.setBorder(new EmptyBorder(5, 5, 5, 5));
			mainPane.addMouseListener(new MouseAdapter()
				{
				@Override
				public void mouseClicked(MouseEvent e)
					{
					Component c= e.getComponent();
					while(c!=null) {
						if(c instanceof JTable) {
							lastSelectedJtable = JTable.class.cast(c);
							break;
							}
						if(c.getParent()==null) break;
						if(!(c.getParent() instanceof Component)) break;
						c = Component.class.cast(c.getParent());
						}
					}
				});
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
				topPane.add(cb);
				cb.addActionListener(AE->refreshInterval());
				this.variantType2cbox.put(vtype, cb);
				}
			topPane.add(new JSeparator());
			label= new JLabel("JEXL:");
			this.jtextFieldGEXL = new JTextField(10);
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
			this.swingVCFGenotypesTableModel = new SwingVCFGenotypesTableModel(this.pedigree);
			final JTable genotypeTable =  new JTable(this.swingVCFGenotypesTableModel);
			
			final JPanel gtPane = new JPanel(new BorderLayout(1,1));
			gtPane.setBorder(BorderFactory.createTitledBorder("Genotypes"));
			final JPanel topPane2 = new JPanel(new FlowLayout(FlowLayout.LEADING));
			gtPane.add(topPane2,BorderLayout.NORTH);
			
			for(GenotypeType gttype: GenotypeType.values()) {
				final JCheckBox cb = new JCheckBox(gttype.name(),true);
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
			
			gtPane.add(new JScrollPane(genotypeTable),BorderLayout.CENTER);
			
			final JSplitPane split1 = new JSplitPane(
					JSplitPane.HORIZONTAL_SPLIT,
					new JScrollPane(variantTable),
					gtPane
					);
			
			
			this.swingInfoTableModel = new SwingVCFInfoTableModel();
			
			this.swingAnnPredictionTableModel = new SwingAnnPredictionTableModel(header);
			this.swingVepPredictionTableModel = new SwingVepPredictionTableModel(header);
			this.swingBcsqPredictionTableModel = new SwingBcsqPredictionTableModel(header);
			this.genotypeTypeTableModel = new GenotypeTypeTableModel();
			
			
			final JTabbedPane tabbed2 = new JTabbedPane();
			tabbed2.addTab("INFO",new JScrollPane(new JScrollPane( new JTable(this.swingInfoTableModel))));
			tabbed2.addTab("Types",new JScrollPane(new JScrollPane(new JTable(this.genotypeTypeTableModel))));
			tabbed2.addTab("SnpEff:ANN",new JScrollPane(new JTable(this.swingAnnPredictionTableModel)));
			tabbed2.addTab("VEP",new JScrollPane(new JTable(this.swingVepPredictionTableModel)));
			tabbed2.addTab("bcftools:BCSQ",new JScrollPane(new JTable(this.swingBcsqPredictionTableModel)));
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
			pane2.add(split2,BorderLayout.CENTER);
			this.variantTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
			this.variantTable.getSelectionModel().addListSelectionListener(E->refreshTables());
			
			tabPane.addTab("Main", pane2);
			
			tabPane.addTab("INFO", new JScrollPane( new JTable(new SwingVCFInfoHeaderLineTableModel(header))));
			tabPane.addTab("FORMAT", new JScrollPane( new JTable(new SwingVCFFormatHeaderLineTableModel(header))));
			tabPane.addTab("FILTER", new JScrollPane( new JTable(new SwingVCFFilterHeaderLineTableModel(header))));
			if(this.pedigree!=null) {
				tabPane.addTab("Pedigree", new JScrollPane( new JTable(new SwingPedigreeTableModel(this.pedigree))));
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
			menu.add(new JMenuItem(new AbstractAction("Save Table as...")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					saveTable(lastSelectedJtable);
					}
				}));
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
				final String jexlExpr = this.jtextFieldGEXL.getText();
				try(CloseableIterator<VariantContext> r=this.vcfReader.query(location.get())) {
					Stream<VariantContext> st = r.stream();
					st=st.filter(V->{
						final VariantContext.Type vtype = V.getType();
						if(!this.variantType2cbox.get(vtype).isSelected()) return false;
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
		
		this.swingAnnPredictionTableModel.setVariant(ctx);
		this.swingVepPredictionTableModel.setVariant(ctx);
		this.swingBcsqPredictionTableModel.setVariant(ctx);
		this.genotypeTypeTableModel.setVariant(ctx);
		if(this.swingTrioTableModel!=null) {
			this.swingTrioTableModel.setVariant(ctx);
		}
		
		final List<Genotype> genotypes;
		if(ctx==null) {
			genotypes = Collections.emptyList();
			}
		else {
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
	private void saveTable(JTable jtable) {
		if(jtable==null) return;
		try {
			JFileChooser fc = new JFileChooser();
			if(fc.showSaveDialog(jtable)!=JFileChooser.APPROVE_OPTION) return;
			final File f = fc.getSelectedFile();
			System.err.println("ok");
			}
		catch(Throwable err) {
			ThrowablePane.show(jtable, err);
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
	
	
	@Override
	public int doWork(final List<String> args)
		{
		try {
			final String input = this.oneAndOnlyOneFile(args);
			final Path vcfPath = Paths.get(input);
			IOUtil.assertFileIsReadable(vcfPath);
			
			final Pedigree ped;
			if(this.pedigreePath!=null) {
				ped = new PedigreeParser().parse(this.pedigreePath);
				}
			else
				{
				ped = null;
				}
			
			JFrame.setDefaultLookAndFeelDecorated(true);
			final XFrame frame = new XFrame(vcfPath,defaultRegion,this.limit_number_variant,ped);
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
