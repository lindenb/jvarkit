/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
import java.awt.Component;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.Vector;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileFilter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.swing.PreferredDirectory;
import com.github.lindenb.jvarkit.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import com.github.lindenb.jvarkit.variant.swing.ActionIndexVcf;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFFilterHeaderLineTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFFormatHeaderLineTableModel;
import com.github.lindenb.jvarkit.variant.swing.SwingVCFInfoHeaderLineTableModel;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.DynamicIndexCreator;
import htsjdk.tribble.index.IndexCreator;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.JavascriptVariantFilter;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC

## Example

```
java -jar dist/swingvcfjexl.jar src/test/resources/rotavirus_rf.vcf.gz
```

## Screenshot

https://twitter.com/yokofakun/status/1514287386410336259

![https://twitter.com/yokofakun/status/1514287386410336259](https://pbs.twimg.com/media/FQPUMpbWUAgJbCs?format=png&name=900x900)

END_DOC
 */
@Program(name="swingvcfjexl",
description="Filter VCF using Java Swing UI and JEXL/Javascript expression",
keywords={"vcf","visualization","swing","jexl","javascript"},
creationDate="20220413",
modificationDate="20220414",
generate_doc=true,
jvarkit_amalgamion =  true,
menu="VCF Manipulation"
)
public class SwingVcfJexlFilter extends Launcher {
	private static final Logger LOG = Logger.of(SwingVcfJexlFilter.class);
	private enum SavingState {IDLE,SAVING};
	private static final String[] languages = {"JEXL","JAVASCRIPT"};
	private static String ABOUT="";
	
	private static class Saver extends Thread {
		SAMSequenceDictionary dict;
		XFrame owner;
		File inputFile;
		File outputFile;
		String jexlExpr;
		Predicate<VariantContext> selVariant;
		Locatable interval;
		Long limit;
		boolean write_index=false;
		@Override
		public void run() {
			final File idx;
			
			if(this.inputFile.getName().toLowerCase().endsWith(".gz")){
				idx =  Tribble.tabixIndexFile(this.inputFile);
				}
			else if(this.inputFile.getName().toLowerCase().endsWith(".vcf")){
				idx =  Tribble.indexFile(this.inputFile);
				}
			else
				{
				idx = null;
				}
			final boolean require_tbi_index= this.interval!=null && idx!=null && idx.exists();
			IndexCreator indexCreator = null;
			// input is indexed
			if(this.write_index) {
				if(this.outputFile.getName().toLowerCase().endsWith(".gz")){
					indexCreator = this.dict == null ?
						new TabixIndexCreator(TabixFormat.VCF):
						new TabixIndexCreator(this.dict,TabixFormat.VCF);
					}
				else if(this.outputFile.getName().toLowerCase().endsWith(".vcf")){
					indexCreator = new DynamicIndexCreator(
							this.outputFile,
							 IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME
							);
					}
				}
			
			final VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
			vcwb.setCreateMD5(false);
			vcwb.setReferenceDictionary(this.dict);
			vcwb.setOutputFile(this.outputFile);
			if(indexCreator!=null) {
				vcwb.setIndexCreator(indexCreator);
				vcwb.setOption(Options.INDEX_ON_THE_FLY);
				}
			vcwb.setOutputFileType(outputFile.getName().toLowerCase().endsWith(".gz")?
					VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF:
					VariantContextWriterBuilder.OutputType.VCF
					);

			boolean show_success=true;
			final VCFReaderFactory vcfReaderFactory = VCFReaderFactory.makeDefault();			
			try (VCFReader vcfReader = vcfReaderFactory.open(inputFile, require_tbi_index);
				VariantContextWriter vcw = vcwb.build()) {
				final VCFHeader header = vcfReader.getHeader();
				header.addMetaDataLine(new VCFHeaderLine(SwingVcfJexlFilter.class.getSimpleName()+".jexl", StringUtils.normalizeSpaces(this.jexlExpr)));
				header.addMetaDataLine(new VCFHeaderLine(SwingVcfJexlFilter.class.getSimpleName()+".version", SwingVcfJexlFilter.ABOUT));
				vcw.writeHeader(header);
				try(CloseableIterator<VariantContext> iter = (require_tbi_index?
						vcfReader.query(this.interval):
						vcfReader.iterator()
						)) {
					long n_seen_variants = 0L;
					long n_saved_variants = 0L;
					long last_millisec = System.currentTimeMillis();
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						++n_seen_variants;
						if(this.owner.currentThead!=this) {
							LOG.error("Cancel (!thread)");
							show_success = false;
							break;
							}
						if(this.owner.state!=SavingState.SAVING) {
							LOG.error("Cancel (!state)");
							show_success = false;
							break;
							}
						long now = System.currentTimeMillis();
						if((now - last_millisec)/1_000 > 10) {
							final String title = ctx.getContig()+":"+ctx.getStart()+" N="+
										StringUtils.niceInt(n_saved_variants)+"/"+
										StringUtils.niceInt(n_seen_variants);
							SwingUtilities.invokeLater(()->{
								owner.jprogressbar.setToolTipText(title);
								});
							last_millisec = now;
							}
						if(this.interval!=null && !this.interval.overlaps(ctx)) {
							continue;
							}
						if(!this.selVariant.test(ctx)) continue;
						vcw.add(ctx);
						n_saved_variants++;
						if(this.limit!=null && this.limit.longValue()>=n_saved_variants) break;
						}
					if(show_success) {
						final String title = "Saved "+ this.outputFile.getName() +" "+
								StringUtils.niceInt(n_saved_variants)+" / "+
								StringUtils.niceInt(n_seen_variants)+" variants";
						SwingUtilities.invokeLater(()->{
							owner.stopSaving();
							JOptionPane.showMessageDialog(owner,title);
							});
						}
					}
				}
			catch(final Throwable err) {
				SwingUtilities.invokeLater(()->{
					owner.stopSaving();
					ThrowablePane.show(owner, err);
					});
				}
			finally {
				SwingUtilities.invokeLater(()->{
					owner.stopSaving();
				});
				}
			}
		}
	
	
	
	
	@SuppressWarnings("serial")
	private static class XFrame extends JFrame {
		final SAMSequenceDictionary dict;
		final VCFHeader vcfHeader;
		final File inputVcf;
		final JTextField jtextFieldLocation;
		final JTextField jtextFieldLimit;
		final JTextArea jtextAreaJEXL;
		final JList<String> jlistSamples;
		final JButton jbuttonRun;
		final JProgressBar jprogressbar;
		final JCheckBox jcboxWriteIndex;
		final JComboBox<String> jcomboxLang;
		Saver currentThead = null;
		transient SavingState state = SavingState.IDLE;
		
		/** insert example code in text area */
		private class ActionExample extends AbstractAction {
			final String code;
			ActionExample(final String title,final String code) {
			super(title);
			this.code = code;
			}
			@Override
			public void actionPerformed(ActionEvent e)
					{
					if(XFrame.this.state!=SavingState.IDLE) return;
					jtextAreaJEXL.setText(this.code);
					}
			}
		
		/** insert example code in text area */
		private class OpenWebPageAction extends AbstractAction {
			final URI url;
			OpenWebPageAction(final String title,final String url) {
			super(title);
			this.url = URI.create(url);
			}
			@Override
			public void actionPerformed(ActionEvent e)
					{
					final Desktop desktop = Desktop.isDesktopSupported() ? Desktop.getDesktop() : null;
				    if (desktop != null && desktop.isSupported(Desktop.Action.BROWSE)) {
				        try {
				            desktop.browse(this.url);
				        } catch (Exception err) {
				            ThrowablePane.show(XFrame.this, err);
				        }
				    }
					}
			}
		
		private class SampleCode extends AbstractAction {
			final Set<GenotypeType> types= new HashSet<>();
			final boolean andFlag;
			SampleCode(boolean andFlag,GenotypeType...types) {
				super((andFlag?"ALL":"Any")+" genotypes must be "+(types.length>1?"one of ":"")+Arrays.stream(types).map(GT->GT.name()).collect(Collectors.joining(" or ")));
				for(GenotypeType gt:types) this.types.add(gt);
				this.andFlag=andFlag;
				}
			String code(GenotypeType gt) {
				switch(gt) {
					case HET: return "isHet";
					case HOM_REF: return "isHomRef";
					case HOM_VAR: return "isHomVar";
					case NO_CALL: return "isNoCall";
					case MIXED: return "isMixed";
					case UNAVAILABLE: return "todo";
					}
				return "";
				}
			String one(final String sn) {
				String s= types.stream().
						map(GT-> "vc.getGenotype(\""+sn+"\")."+code(GT)+"()").
						collect(Collectors.joining(" || "));
				if(types.size()>1) s="("+s+")";
				return s;
				}
			@Override
			public void actionPerformed(ActionEvent e)
				{
				final List<String> samples;
				if(jlistSamples.isSelectionEmpty()) {
					samples = XFrame.this.vcfHeader.getSampleNamesInOrder();
					}
				else
					{
					samples = jlistSamples.getSelectedValuesList();
					}
				if(samples.isEmpty()) return;
				String txt = "("+
						samples.
						stream().
						map(S->one(S)).
						collect(Collectors.joining(andFlag?" &&  ":" || ")) +
						")"
						;
				jtextAreaJEXL.insert(txt, jtextAreaJEXL.getCaretPosition());
				}
			}
		
		
		XFrame(final File inputVcf) throws IOException {
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			this.inputVcf = inputVcf;
			setTitle(SwingVcfJexlFilter.class.getSimpleName()+": "+ inputVcf.getName());
			final JPanel mainPane = new JPanel(new BorderLayout(5, 5));
			mainPane.setBorder(new EmptyBorder(5, 5, 5, 5));
			this.setContentPane(mainPane);
			final JTabbedPane  tabPane = new JTabbedPane();
			mainPane.add(tabPane,BorderLayout.CENTER);

			final JPanel pane1 = new JPanel(new BorderLayout(5,5));
			pane1.add(new JScrollPane(this.jtextAreaJEXL=new JTextArea()),BorderLayout.CENTER);
			this.jtextAreaJEXL.setFont(this.jtextAreaJEXL.getFont().deriveFont(36f));
			tabPane.addTab("JEXL",pane1);
			JPanel pane2= new JPanel(new FlowLayout(FlowLayout.LEADING));
			JLabel label = new JLabel("Language:");
			pane2.add(label);
			this.jcomboxLang=new JComboBox<>(SwingVcfJexlFilter.languages);
			this.jcomboxLang.setSelectedIndex(0);
			label.setLabelFor(this.jcomboxLang);
			pane2.add(this.jcomboxLang);
			final String LABEL_JEXL = "JEXL filtering expressions. See https://gatk.broadinstitute.org/hc/en-us/articles/360035891011";
			final String LABEL_JS = "Javascript expression https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/filter/JavascriptVariantFilter.html .";
			final JLabel lblLang = new JLabel(LABEL_JEXL);
			pane2.add(lblLang);
			this.jcomboxLang.addActionListener(new ActionListener()
				{
				@Override
				public void actionPerformed(final ActionEvent e)
					{
					String msg="";
					int i = jcomboxLang.getSelectedIndex();
					if(i==0) {msg = LABEL_JEXL;}
					else if(i==1) {msg = LABEL_JS;}
					lblLang.setText(msg);
					lblLang.setToolTipText(msg);
					}
				});
			
			pane1.add(pane2,BorderLayout.NORTH);
			
			pane2= new JPanel(new FlowLayout(FlowLayout.TRAILING));
			pane1.add(pane2,BorderLayout.SOUTH);
			
			this.jprogressbar=new JProgressBar();
			this.jprogressbar.setPreferredSize(new Dimension(300, 20));
			pane2.add(this.jprogressbar);
			pane2.add(new JSeparator(SwingConstants.VERTICAL));
			//limit
			label = new JLabel("Limit:",JLabel.TRAILING);
			this.jtextFieldLimit = new JTextField(9);
			label.setLabelFor(this.jtextFieldLimit);
			label.setToolTipText("optional. limit number of variants.");
			pane2.add(label);
			pane2.add(this.jtextFieldLimit);
			pane2.add(new JSeparator(SwingConstants.VERTICAL));
			//region
			label = new JLabel("Region:",JLabel.TRAILING);
			this.jtextFieldLocation = new JTextField(15);
			label.setLabelFor(this.jtextFieldLocation);
			label.setToolTipText("optional. syntax: 'chrom:start-end'");
			pane2.add(label);
			pane2.add(this.jtextFieldLocation);
			pane2.add(new JSeparator(SwingConstants.VERTICAL));
			this.jcboxWriteIndex=new JCheckBox("Write Index", true);
			this.jcboxWriteIndex.setToolTipText("save tribble or tabix index");
			pane2.add(this.jcboxWriteIndex);
			pane2.add(new JSeparator(SwingConstants.VERTICAL));
			
			this.jbuttonRun = new JButton(new AbstractAction("Save...")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					if(XFrame.this.state==SavingState.IDLE) {
						startSaving(jbuttonRun);
						}
					else
						{
						stopSaving();
						}
					}
				});
			pane2.add(jbuttonRun);
			this.jbuttonRun.setOpaque(true);
			this.jbuttonRun.setFont(this.jbuttonRun.getFont().deriveFont(36f));
			this.jbuttonRun.setBackground(Color.GREEN);
			this.jbuttonRun.setForeground(Color.WHITE);
			final JMenu menuCode = new JMenu("Code");
			final VCFReaderFactory vcfReaderFactory = VCFReaderFactory.makeDefault();
			try (VCFReader vcfReader = vcfReaderFactory.open(inputVcf, false)) {
				this.vcfHeader = vcfReader.getHeader();
				this.dict = this.vcfHeader.getSequenceDictionary();
				


				tabPane.addTab("INFO", wrapTable("INFO", new JTable(new SwingVCFInfoHeaderLineTableModel(this.vcfHeader))));
				tabPane.addTab("FORMAT",wrapTable("FILTER",new JTable(new SwingVCFFormatHeaderLineTableModel(this.vcfHeader))));
				tabPane.addTab("FILTER",wrapTable("FILTER",new JTable(new SwingVCFFilterHeaderLineTableModel(this.vcfHeader))));
				
				if(this.vcfHeader.hasGenotypingData()) {
					jlistSamples =new JList<>(new Vector<String>(this.vcfHeader.getSampleNamesInOrder()));
					jlistSamples.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
					tabPane.addTab("Samples",new JScrollPane(jlistSamples));
					for(int orand=0;orand<2;++orand) {
						menuCode.add(new SampleCode(orand==0, GenotypeType.HET));
						menuCode.add(new SampleCode(orand==0, GenotypeType.HET,GenotypeType.HOM_VAR));
						menuCode.add(new SampleCode(orand==0, GenotypeType.HOM_REF));
						menuCode.add(new SampleCode(orand==0, GenotypeType.HOM_REF,GenotypeType.NO_CALL));
						menuCode.add(new SampleCode(orand==0, GenotypeType.HOM_REF,GenotypeType.NO_CALL,GenotypeType.HET));
						menuCode.add(new SampleCode(orand==0, GenotypeType.HOM_VAR));
						menuCode.add(new SampleCode(orand==0, GenotypeType.NO_CALL));
						if(orand==0) menuCode.add(new JSeparator());
						}
					
					
					menuCode.add(new JSeparator());
					}
				else
					{
					jlistSamples = null;
					}
				if(this.dict!=null) {
					tabPane.addTab("REF", new JScrollPane( new JTable(new SwingSequenceDictionaryTableModel(this.dict))));
					}

				}
			final JMenuBar menuBar= new JMenuBar();
			setJMenuBar(menuBar);
			JMenu menu = new JMenu("File");
			menuBar.add(menu);
			menu.add(new AbstractAction("Open...")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					doMenuOpen(XFrame.this);
					}
				});
			menu.add(new ActionIndexVcf());
			menu.add(new JSeparator());
			menu.add(new JMenuItem(new AbstractAction("Close")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					stopSaving();
					XFrame.this.setVisible(false);
					XFrame.this.dispose();
					}
				}));
			
			menuBar.add(menuCode);
			menuCode.add(new JMenuItem(new ActionExample("Example 1","!vc.isFiltered()")));
			menu = new JMenu("Documentation");
			menuBar.add(menu);
			menu.add(new JMenuItem(new OpenWebPageAction("VariantContext","https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html")));
			menu.add(new JMenuItem(new OpenWebPageAction("Genotype","https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Genotype.html")));
			menu.add(new JMenuItem(new OpenWebPageAction("Allele","https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Allele.html")));
			menu.add(new JMenuItem(new OpenWebPageAction("VCFHeader","https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html")));
			menu.add(new JMenuItem(new OpenWebPageAction("JEXL","https://gatk.broadinstitute.org/hc/en-us/articles/360035891011")));
			menu.add(new JMenuItem(new OpenWebPageAction("JS","https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/filter/JavascriptVariantFilter.html")));
			menu.add(new JMenuItem(new OpenWebPageAction("SwingVcfJexlFilter","http://lindenb.github.io/jvarkit/SwingVcfJexlFilter.html")));
			menu.add(new JSeparator());
			menu.add(new AbstractAction("About") {
				@Override
				public void actionPerformed(ActionEvent e) {
					JOptionPane.showMessageDialog(XFrame.this,SwingVcfJexlFilter.ABOUT);
				}
			});

			}
		

	
	private synchronized void stopSaving() {
		this.state = SavingState.IDLE;
		if(this.currentThead!=null) {
				try {
					this.currentThead.join(10_000);
				} catch(Throwable err) {
					LOG.error(err);
				}
			}
		this.currentThead=null;
		this.jbuttonRun.setBackground(Color.GREEN);
		this.jbuttonRun.setText("Save");
		this.jprogressbar.setToolTipText("");
		this.jprogressbar.setIndeterminate(false);
		}
	
	
	private synchronized void startSaving(Component parent) {
		stopSaving();
		final Saver saver = new Saver();
		saver.dict = dict;
		saver.owner = this;
		saver.inputFile = this.inputVcf;
		saver.write_index = this.jcboxWriteIndex.isSelected();
		saver.jexlExpr = this.jtextAreaJEXL.getText();
		boolean jexl_lang = this.jcomboxLang.getSelectedIndex()==0;
		if(StringUtil.isBlank(saver.jexlExpr)) {
			saver.selVariant = V->true;
			saver.jexlExpr = "ALL";
			}
		else if(jexl_lang){
			try {
				saver.selVariant = JexlVariantPredicate.create(saver.jexlExpr);
				}catch(Throwable err) {
				ThrowablePane.show(parent, err);
				return;
				}
			}
		else //javascript
			{
			try { 
				saver.selVariant = new JavascriptVariantFilter(saver.jexlExpr, this.vcfHeader) {
					@Override
					public String getRecordKey()
						{
						return "vc";
						}
					};
				} catch(Throwable err) {
				ThrowablePane.show(parent, err);
				return;
				}
			}
		
		final String s = this.jtextFieldLocation.getText().trim();
		if(StringUtil.isBlank(s)) {
			saver.interval = null;
		} else {
			Optional<SimpleInterval> ret;
			try {
				ret = new IntervalParser(this.dict).
					apply(s);
				}
			catch(final Throwable err) {
				ThrowablePane.show(parent, err);
				return;
				}
			saver.interval = ret.get();
			saver.jexlExpr+= " in interval "+saver.interval;
			}
		
		final String limitStr = this.jtextFieldLimit.getText().trim();
		if(StringUtil.isBlank(limitStr)) {
			saver.limit = null;
		} else
		{
			try {
				saver.limit = Long.parseLong(limitStr.replace(",", ""));
			} catch(final Throwable err) {
				ThrowablePane.show(parent, err);
				return;
				}
		}
		
		final JFileChooser fc = new JFileChooser();
		fc.setSelectedFile(new File("save.vcf.gz"));
		
		if(fc.showSaveDialog(parent)!=JFileChooser.APPROVE_OPTION) {
			return;
			}
		saver.outputFile =fc.getSelectedFile();
		if(saver.outputFile==null) return;
		if(saver.outputFile.equals(saver.inputFile)) {
			JOptionPane.showMessageDialog(parent, "input file == output file");
			return;
			}
		if(!(saver.outputFile.getName().toLowerCase().endsWith(".vcf") || saver.outputFile.getName().toLowerCase().endsWith(".vcf.gz"))) {
			JOptionPane.showMessageDialog(parent, "file suffix should be .vcf or vcf.gz ("+saver.outputFile.getName()+")");
			return;
			}
		
		if(saver.outputFile.exists() &&
				JOptionPane.showConfirmDialog(parent, saver.outputFile.toString()+" already exists. Overwrite ?", "Overwrite", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE, null)!=JOptionPane.OK_OPTION
				)
			{
			return;
			}

		
		this.state = SavingState.SAVING;
		this.jbuttonRun.setBackground(Color.RED);
		this.jbuttonRun.setText("Stop");
		this.jprogressbar.setToolTipText("");
		this.jprogressbar.setIndeterminate(true);
		this.currentThead = saver;
		saver.start();
		}

	private JPanel wrapTable(final String title,final JTable t) {
		final JPanel pane = new JPanel(new BorderLayout(1,1));
		final JPanel top = new JPanel(new FlowLayout(FlowLayout.TRAILING,1,1));
		pane.setBorder(BorderFactory.createTitledBorder(title));
		pane.add(top,BorderLayout.NORTH);
		pane.add(new JScrollPane(t),BorderLayout.CENTER);
		return pane;
		}
	}	
		
	private static List<File> selectVcfFiles(final Component parent) {
		final JFileChooser jfc = new JFileChooser(PreferredDirectory.get(SwingVcfJexlFilter.class));
		jfc.setMultiSelectionEnabled(true);
		jfc.setFileFilter(new FileFilter()
			{
			@Override
			public String getDescription()
				{
				return "Variant File";
				}
				
			@Override
			public boolean accept(final File f)
				{
				if(f.isDirectory()) return true;
				if(!f.canRead()) return false;
				if(FileExtensions.VCF_LIST.stream().noneMatch(X->f.getName().endsWith(X))) return false;
				return true;
				}
			});
		if(jfc.showOpenDialog(parent)!=JFileChooser.APPROVE_OPTION) return Collections.emptyList();
		final File[] sel = jfc.getSelectedFiles();
		if(sel==null || sel.length==0)  return Collections.emptyList();
		PreferredDirectory.update(SwingVcfJexlFilter.class, sel[0]);
		return Arrays.asList(sel);
		}
	
	private static void createNewFrame(Component parent,final List<File> all_vcf_paths ) {
		final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
		try {
			SwingUtilities.invokeAndWait(()->{
			int dx = 0;
			for(final File vcfPath: all_vcf_paths) {
				try {
					final XFrame frame = new XFrame(vcfPath);
					frame.setBounds(50 + dx , 50 + dx, screen.width-100, screen.height-100);
					frame.setVisible(true);
					dx += 2;
				} catch(final Throwable err) {
					ThrowablePane.show(null, err);
					System.exit(-1);
					}
				}
			});
		} catch(final Throwable err) {
			ThrowablePane.show(parent,err);
			}
		}
	
	private static void doMenuOpen(final Component parent) {
		final List<File> L = selectVcfFiles(parent);
		if(L.isEmpty()) return;
		createNewFrame(parent,L);
		}
	
	
	@Override
	public int doWork(final List<String> args)
		{
		try {
			SwingVcfJexlFilter.ABOUT =" Author: Pierre Lindenbaum "+
					" Version:"+getVersion()+
					" Date:"+getCompileDate()
					;

			
			final List<File> all_vcf_paths = 
					IOUtils.unrollPaths(args).
					stream().
					map(P->P.toFile()).
					collect(Collectors.toCollection(ArrayList::new))
					;
			if(all_vcf_paths.isEmpty()) {
				final List<File> L= selectVcfFiles(null);
				if(L.isEmpty()) return -1;
				all_vcf_paths.addAll(L);
				}
			
			
			for(final File path: all_vcf_paths) {
				IOUtil.assertFileIsReadable(path);
				}
					
			JFrame.setDefaultLookAndFeelDecorated(true);
			createNewFrame(null,all_vcf_paths);
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			ThrowablePane.show(null, err);
			return -1;
			}
		}
	
	public static void main(final String[] args)
		{
		new SwingVcfJexlFilter().instanceMain(args);//no exit
		}
	}
