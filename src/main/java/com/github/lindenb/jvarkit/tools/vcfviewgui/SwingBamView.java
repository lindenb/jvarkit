package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.Vector;
import java.util.stream.Collectors;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
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
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;

public class SwingBamView extends Launcher {
	private static final Logger LOG = Logger.build(SwingBamView.class).make();
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFile = null;

	
	
	private class XFrame extends JFrame {
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
		
		
		private class BamFile {
			final Path path;
			final String sn;
			SAMFileHeader header= null;
			SamReader sr = null;
			BamFile(final Path path) {
				this.path = path;
				open();
				this.sn = header.getReadGroups().
					stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtil.isBlank(S)).
					findFirst().
					orElse(IOUtils.getFilenameWithoutCommonSuffixes(path));
				close();
				}
			
			List<Read> query() {
				if(XFrame.this.interval==null) return Collections.emptyList();
				open();
				final SamRecordFilter srf = createSamRecordFilter();
				try(CloseableIterator<SAMRecord> iter = this.sr.query(
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
			
			void open() {
				if(sr!=null) return;
				final SamReaderFactory srf = SamReaderFactory.makeDefault();
				srf.referenceSequence(XFrame.this.referenceFaidx);
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
			}

		
		XFrame(final Path referenceFaidx,final List<Path> bamList) {
			
			
			
			this.referenceFaidx = referenceFaidx;
			this.referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referenceFaidx);
			this.dict = SequenceDictionaryUtils.extractRequired(this.referenceSequenceFile);
			
			JPanel mainPanel=new JPanel(new BorderLayout());
			
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
			final Action actionGo = new AbstractAction("Load")
				{
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
			for(final Path p:bamList) {
				dm.addElement(new BamFile(p));
				}
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
			viewPane.add(bamPanel,BorderLayout.EAST);
			
			
			jTabbedPane.addTab("REF", new JScrollPane( new JTable(new SwingSequenceDictionaryTableModel(dict))));
			}
		
		
		
		
		private void paintDrawingArea(Graphics2D g) {
			
			
			for(int y=0;y< this.rows.size();y++) {
				final List<Read> row = this.rows.get(y);
				for(final Read read:row) {
					final Cigar cigar =read.rec.getCigar();
					int refpos = read.rec.getUnclippedStart();
					int readpos;
					for(CigarElement ce:cigar) {
						
					}
				}
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
			BamFile f = getSelectedBamFile();
			if(f==null) return;
			this.rows.clear();
			final List<Read> buffer= f.query();
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
			this.drawingArea.repaint();
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
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final List<Path> paths = IOUtils.unrollPaths(args);
			if(paths.isEmpty()) {
				LOG.error("no BAM/CRAM was provided");
				return -1;
				}
			JFrame.setDefaultLookAndFeelDecorated(true);
			final XFrame frame = new XFrame(this.referenceFile,paths);
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
