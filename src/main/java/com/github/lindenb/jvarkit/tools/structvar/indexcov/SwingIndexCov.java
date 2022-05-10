/*
The MIT License (MIT)

Copyright (c) 2022 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.structvar.indexcov;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Shape;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentAdapter;
import java.awt.event.WindowAdapter;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
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
import java.util.stream.Collectors;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileFilter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.PreferredDirectory;
import com.github.lindenb.jvarkit.util.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;

public class SwingIndexCov extends Launcher {
	private static final Logger LOG = Logger.build(SwingIndexCov.class).make();
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFile = null;



	
	private static class XFrame extends JFrame {
		private final Path referenceFaidx;
		private final SAMSequenceDictionary dict;
		private final JList<String> jlistSamples;
		private final JList<String> jlistChromosomes;
		private final JTextField jtextFieldLocation;
		private Interval interval =null;
		private final JPanel drawingArea;
		private final ReferenceSequenceFile referenceSequenceFile;
		private TabixFileReader taxbixFileReader = null;
		
		
		
		XFrame(final Path referenceFaidx,final String tabixFileUri) throws IOException {
			super(SwingIndexCov.class.getSimpleName());
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			
			final JMenuBar menuBar = new JMenuBar();
			setJMenuBar(menuBar);
			JMenu menu=new JMenu("File");
			menuBar.add(menu);
			menu.add(new JMenuItem(new AbstractAction("Open...") {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					}
				}));
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
			
			
			this.taxbixFileReader = new TabixFileReader(tabixFileUri);
			//read samples
			
			String line =this.taxbixFileReader.readLine();
			if(line==null) {		
				throw new IOException( "Cannot read first line of input");
				}
			final String tokens[] = CharSplitter.TAB.split(line);
			if(tokens.length<4 ||
				!tokens[0].equals("#chrom") ||
				!tokens[1].equals("start") ||
				!tokens[2].equals("end")) {
				throw new IOException("bad first line "+line );
				}
			this.jlistSamples = new JList<String>(new Vector<>(
				Arrays.asList(tokens).subList(3,tokens.length))
				);
		
			// get chromosomes in the tabix
			final Vector<String> vecChroms = new Vector<>( this.taxbixFileReader.getChromosomes());
			//remove if not in dict
			if(vecChroms.removeIf(S->this.dict.getSequence(S)==null));
			if(vecChroms.isEmpty()) {
				throw new IllegalArgumentException("no common chromosomes between "+referenceFaidx+" and "+tabixFileUri);
				}
			//sort on dict
			vecChroms.sort((A,B)->Integer.compare(this.dict.getSequenceIndex(A), this.dict.getSequenceIndex(B)));
			this.jlistChromosomes = new JList<>(vecChroms);
			
			
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
					reloadReads();
					}
				};
			navPane.add(new JButton(actionGo));
			navPane.add(new JSeparator());
		
			
			final JPanel leftPane = new JPanel(new GridLayout(2,1));
			viewPane.add(leftPane,BorderLayout.WEST);
			JPanel pane2 = new JPanel(new BorderLayout(5, 5));
			pane2.setBorder(BorderFactory.createTitledBorder("Samples"));
			pane2.add(new JScrollPane(this.jlistSamples),BorderLayout.CENTER);
			leftPane.add(pane2);
			pane2 = new JPanel(new BorderLayout(5, 5));
			pane2.setBorder(BorderFactory.createTitledBorder("Chroms"));
			pane2.add(new JScrollPane(this.jlistChromosomes),BorderLayout.CENTER);
			leftPane.add(pane2);
			
			
			/** drawing area itself */
			this.drawingArea = new JPanel(null, true) {
				@Override
				protected void paintComponent(Graphics g) {
					paintDrawingArea(Graphics2D.class.cast(g));
					}
				};
			this.drawingArea.setOpaque(true);
			
			/** scroll bar for vertical navigation */
			this.drawingArea.addComponentListener(new ComponentAdapter() {
				public void componentResized(java.awt.event.ComponentEvent e) {
					drawingArea.repaint();
				};
			});
			
			/** ref tab */
			jTabbedPane.addTab("REF", new JScrollPane( new JTable(new SwingSequenceDictionaryTableModel(dict))));
			/* open close operations */
			this.addWindowListener(new WindowAdapter() {
				public void windowClosed(java.awt.event.WindowEvent e)
					{
					runOnClose();
					}
				});
			}
		
		private void runOnClose() {
			try {
				this.referenceSequenceFile.close();
				}
			catch(final IOException err) {
				LOG.error(err);
				}
			try {
				this.taxbixFileReader.close();
				}
			catch(final Throwable err) {
				LOG.error(err);
				}
 			}
		
		private Shape createShape(double y,int position,int lengthOnReference,int arrow_type) {
			if(position>this.interval.getEnd()) return null;
			if(position+lengthOnReference<this.interval.getStart()) return null;
			double x0 = pos2pixel(position);
			double x1 = pos2pixel(position+lengthOnReference);
			double featH = getFeatureHeight();
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
		
		private void paintDrawingArea(final Graphics2D g) {
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, drawingArea.getWidth(), drawingArea.getHeight());
			
			}

				
		private void reloadLocation() {
			final ReferenceSequence refseq = this.referenceSequenceFile.getSubsequenceAt(
					this.interval.getContig(),
					this.interval.getStart(),
					this.interval.getEnd()
					);
			
			}
		
		
		final double pos2pixel(int coord) {
			return (coord - this.interval.getStart())/((double)this.interval.getLengthOnReference())*this.drawingArea.getWidth();
			}
		final double getBasePixel() {
			return (1.0/this.interval.getLengthOnReference())*this.drawingArea.getWidth();
			}
			
		final Optional<SimpleInterval> getUserInterval() {
				final String s = this.jtextFieldLocation.getText().trim();
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
			
				return ret;
				}
	
		
		
		}
	
	
	private static void doMenuOpen(final Component parent,final Path reference) {
		final JFileChooser jfc = new JFileChooser(PreferredDirectory.get(SwingIndexCov.class));
		jfc.setFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "bed.gz";
				}
			
			@Override
			public boolean accept(final File f) {
				if(f.isFile() && f.canRead()) {
					final String fname= f.getName();
					if(fname.endsWith(".bed.gz")) return true;
					}
				return false;
			}
		});
		if(jfc.showOpenDialog(parent)!=JFileChooser.APPROVE_OPTION) return;
		final File file= jfc.getSelectedFile();
		if(file==null) return;
		openFrameForBamFiles(
			parent,
			reference,
			file
			);
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = oneAndOnlyOneFile(args);
			JFrame.setDefaultLookAndFeelDecorated(true);
			openFrameForBamFiles(null,this.referenceFile,input);
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new SwingIndexCov().instanceMain(args);
	}

}
