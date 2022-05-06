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
package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
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
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.PreferredDirectory;
import com.github.lindenb.jvarkit.util.swing.ThrowablePane;

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

public class SwingBamView extends Launcher {
	private static final Logger LOG = Logger.build(SwingBamView.class).make();
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFile = null;


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

	
	private static class XFrame extends JFrame {
		private final Path referenceFaidx;
		private final SAMSequenceDictionary dict;
		private List<BamFile> all_bam_files = new Vector<>();
		private final JList<BamFile> jlistBamFilesJList;
		private final JTextField jtextFieldLocation;
		private Interval interval =null;
		private BamFile currenBamFile = null;
		private boolean use_clip = false;
		private final JScrollBar vScrollBar;
		private final JPanel drawingArea;
		private final ReferenceSequenceFile referenceSequenceFile;
		private final List<List<Read>> rows = new Vector<>();
		private final JCheckBox cboxUseClip;
		private final JCheckBox cboxShowDup;
		private CharSequence referenceSequence;
		
		private class Read implements Locatable {
			private final SAMRecord rec;
			int row_index;
			Read(final SAMRecord rec) {
				this.rec =  rec;
				}
			@Override
			public String getContig() {
				return rec.getContig();
				}
			@Override
			public int getStart() {
				return isUsingClip()?rec.getUnclippedStart():rec.getStart();
				}
			@Override
			public int getEnd() {
				return isUsingClip()?rec.getUnclippedEnd():rec.getEnd();
				}
			}
		
		
		
		XFrame(final Path referenceFaidx,final List<BamFile> bamFiles) {
			super(SwingBamView.class.getSimpleName());
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			
			final JMenuBar menuBar = new JMenuBar();
			setJMenuBar(menuBar);
			JMenu menu=new JMenu("File");
			menuBar.add(menu);
			menu.add(new JMenuItem(new AbstractAction("Open...") {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					doMenuOpen(XFrame.this,XFrame.this.referenceFaidx);
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
			navPane.add(this.cboxUseClip=new JCheckBox("clip", false));
			this.cboxUseClip.setToolTipText("Show/Hide Read Clipping");
			this.cboxUseClip.addActionListener(AE->reloadReads());
			navPane.add(this.cboxShowDup=new JCheckBox("dup", false));
			this.cboxShowDup.setToolTipText("Show/Hide Dup Reads");
			this.cboxShowDup.addActionListener(AE->reloadReads());
			

			
			
			/** list of bam files */
			final DefaultListModel<BamFile> dm = new DefaultListModel<>();
			bamFiles.stream().
				sorted((A,B)->A.getName().compareTo(B.getName())).
				forEach(B->dm.addElement(B));
				
			this.jlistBamFilesJList = new JList<>(dm);
			this.jlistBamFilesJList.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
			this.jlistBamFilesJList.getSelectionModel().addListSelectionListener(AE->{
				
				});
			final JPanel bamPanel = new JPanel(new BorderLayout());
			bamPanel.setBorder(BorderFactory.createTitledBorder("BAMS"));
			bamPanel.add(new JScrollPane(this.jlistBamFilesJList),BorderLayout.CENTER);
			viewPane.add(bamPanel,BorderLayout.WEST);
			
			/** drawing area itself */
			drawingArea = new JPanel(null, true) {
				@Override
				protected void paintComponent(Graphics g) {
					paintDrawingArea(Graphics2D.class.cast(g));
					}
				};
			this.drawingArea.setOpaque(true);
			
			/** scroll bar for vertical navigation */
			this.vScrollBar = new JScrollBar(JScrollBar.VERTICAL);
			this.drawingArea.addComponentListener(new ComponentAdapter() {
				public void componentResized(java.awt.event.ComponentEvent e) {
					adjustScrollBar();
					drawingArea.repaint();
				};
			});
			viewPane.add(bamPanel,BorderLayout.EAST);
			
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
			final double dy = getFeatureHeight();
			double y=0;
			final Hershey hershey=new Hershey();
			y= dy;
			for(int rowidx= this.vScrollBar.getValue();rowidx< this.rows.size();rowidx++) {
				final List<Read> row = this.rows.get(rowidx);
				for(final Read read:row) {
					final Map<Integer,Integer> pos2insert_size = new HashMap<>();
					final byte[] readBases = read.rec.getReadBases();
					final Cigar cigar =read.rec.getCigar();
					// horizontal line for deletions
					if(cigar.getCigarElements().stream().anyMatch(OP->OP.equals(CigarOperator.D) || OP.equals(CigarOperator.N)) ) {
						double x0 = pos2pixel(read.rec.getStart());
						double x1 = pos2pixel(read.rec.getEnd());
						g.setColor(Color.DARK_GRAY);
						g.draw(new Line2D.Double(x0, y+dy/2.0, x1, y+dy/2.0));
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
						Color paper = Color.WHITE;
						Color pen = Color.DARK_GRAY;
						switch(op) {
							case P:break;
							case I: {
								pos2insert_size.put(refpos, len);
								next_readpos+=len;
								break;
								}
							case D: case N: {
								double x0 = pos2pixel(Math.max(interval.getStart(), refpos));								
								double x1 = pos2pixel(Math.min(interval.getEnd(), refpos+len));
								String title = String.valueOf(len); 
								next_refpos+=len;
								break;
								}
							case S: case H:{
								paper = Color.YELLOW;
								if(isUsingClip()) {
									shape = createShape(y,refpos,len,arrow_type);
									}
								if(op!=CigarOperator.H) next_readpos+= len;
								next_refpos+= len;
								break;
								}
							case M:case X: case EQ: {
								shape = createShape(y,refpos,len,arrow_type);
								next_readpos+= len;
								next_refpos+= len;
								break;
								}
							default: throw new IllegalStateException();
							}
						if(shape!=null) {
							if(paper!=null) {
								g.setColor(paper);
								g.fill(shape);
								}
							if(pen!=null) {
								g.setColor(pen);
								g.draw(shape);
								}
							if(op.consumesReadBases() && readBases!=SAMRecord.NULL_SEQUENCE) {
								for(int x=0;x< len;++x) {
									if(refpos+x < interval.getStart()) continue;
									if(refpos+x > interval.getEnd()) break;
									char readBase = (char)readBases[refpos+x];
									char refBase = referenceSequence.charAt(refpos+x-interval.getStart());
									hershey.paint(g,String.valueOf(readBase),
										pos2pixel(refpos+1),
										y,
										getBasePixel(),
										dy
										);
									}
								}
							}
						/* plot insertions */
						for(Integer pos:pos2insert_size.keySet()) {
							double x0 = pos2pixel(pos);
							g.setColor(Color.RED);
							
							}
						readpos= next_readpos;
						refpos= next_refpos;
						}
					}
				rowidx++;
				}
			}

		
		private BamFile getSelectedBamFile() {
			return this.jlistBamFilesJList.getSelectedValue();
		}
		
		private void reloadLocation() {
			final ReferenceSequence refseq = this.referenceSequenceFile.getSubsequenceAt(
					this.interval.getContig(),
					this.interval.getStart(),
					this.interval.getEnd()
					);
			this.referenceSequence = refseq.getBaseString();
			reloadReads();
		}
		
		private void reloadReads() {
			final BamFile f = getSelectedBamFile();
			if(f==null) return;
			this.rows.clear();
			final List<Read> buffer= new ArrayList<>(featchReads());
			while(!buffer.isEmpty()) {
				final Read read = buffer.remove(0);
				int y = 0;
				for(y=0;y< this.rows.size();++y) {
					final List<Read> row = this.rows.get(y);
					final Read last = row.get(row.size()-1);
					if(last.getEnd() >= read.getStart()) continue;
					read.row_index = y;
					row.add(read);
					break;
					}
				if(y==this.rows.size()) {
					final List<Read> row = new Vector<>();
					row.add(read);
					read.row_index = this.rows.size();
					this.rows.add(row);
					}
				}
			adjustScrollBar();
			this.drawingArea.repaint();
			}

		public void adjustScrollBar() {
			int dy=(int)getFeatureHeight();
			int  n_rows = 0;
			n_rows+= 1;//reference
			n_rows+= this.rows.size();
			
			int paperh= this.drawingArea.getHeight();
			int n_visible = Math.max(1,Math.min(paperh/dy,n_rows));
			if(paperh<=1) {
				this.vScrollBar.setValues(0, 0, 0, 0);
				}
			else
				{
				this.vScrollBar.setValues(
						0,1, 0, 0
						);
				}
			}
		
		final double getFeatureHeight() {
			return 20;
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
	
		/** fetch Read for current bam, and current interval **/
		private List<Read> featchReads() {
			final BamFile bf = getSelectedBamFile();
			if(this.interval==null || bf==null) return Collections.emptyList();
			bf.open(this.referenceFaidx);
			final SamRecordFilter srf = createSamRecordFilter();
			try(CloseableIterator<SAMRecord> iter = bf.sr.query(
				XFrame.this.interval.getContig(),
				XFrame.this.interval.getStart(),
				XFrame.this.interval.getEnd(),
				false)) {
				return iter.stream().
						filter(R->!srf.filterOut(R)).
						map(R->new Read(R)).
						sorted((A,B)->Integer.compare(A.getStart(),B.getStart())).
						collect(Collectors.toList());
				}
			}
		
		SamRecordFilter createSamRecordFilter() {
			final List<SamRecordFilter> filters = new ArrayList<>();
			filters.add(new AlignedFilter(true));
			if(!this.cboxShowDup.isSelected()) {
				filters.add(new DuplicateReadFilter());
				}
			return new AggregateFilter(filters);
			}
		
		private boolean isUsingClip() {
			return this.cboxUseClip.isSelected();
			}
		}
	
	private static void openFrameForBamFiles(final Component owner,final Path reference,final List<String> args) {
		try {
			final List<Path> paths = IOUtils.unrollPaths(args);
			if(paths.isEmpty()) {
				LOG.warn("No File was provided");
				return;
				}
			final List<BamFile> bamFiles = new ArrayList<>(paths.size());
			for(final Path p:paths) {
				final BamFile bf = new BamFile(p);
				bf.open(reference);
				bf.sn = bf.header.getReadGroups().
					stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtil.isBlank(S)).
					findFirst().
					orElse(IOUtils.getFilenameWithoutCommonSuffixes(p));
				bf.close();
				bamFiles.add(bf);
				}
			
			final XFrame frame = new XFrame(reference,bamFiles);
			final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
			frame.setBounds(50, 50, screen.width-100, screen.height-100);
	
			
			SwingUtilities.invokeAndWait(()->{
				frame.setVisible(true);
				});
			}
		catch(final Throwable err) {
			ThrowablePane.show(owner, err);
			}
		}
	
	private static void doMenuOpen(final Component parent,final Path reference) {
		final JFileChooser jfc = new JFileChooser(PreferredDirectory.get(SwingBamView.class));
		jfc.setMultiSelectionEnabled(true);
		jfc.setFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "Bam/Cram";
				}
			
			@Override
			public boolean accept(final File f) {
				if(f.isFile() && f.canRead()) {
					String fname= f.getName();
					if(fname.endsWith(".list")) return true;
					if(fname.endsWith(FileExtensions.SAM)) return true;
					if(fname.endsWith(FileExtensions.CRAM)) return true;
					}
				return false;
			}
		});
		if(jfc.showOpenDialog(parent)!=JFileChooser.APPROVE_OPTION) return;
		final File[] files= jfc.getSelectedFiles();
		if(files==null || files.length==0);
		openFrameForBamFiles(
			parent,
			reference,
			Arrays.stream(files).
			map(P->P.toString()).
			collect(Collectors.toList()));
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			if(IOUtils.unrollPaths(args).isEmpty()) {
				LOG.error("no BAM/CRAM was provided");
				return -1;
				}
			
			JFrame.setDefaultLookAndFeelDecorated(true);
			openFrameForBamFiles(null,this.referenceFile,args);
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
