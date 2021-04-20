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
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Panel;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.Vector;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.table.AbstractTableModel;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.CoverageFactory;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
/**
BEGIN_DOC

## Example



END_DOC
 */
@Program(name="swingbamcov",
description="Bam coverage viewer using Java Swing UI",
keywords={"bam","alignment","graphics","visualization","swing"},
creationDate="20210420",
modificationDate="20210420",
generate_doc=true
)
public class SwingBamCov extends Launcher
	{
	private static final Logger LOG = Logger.build(SwingBamCov.class).make();
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFile = null;
	@Parameter(names={"-r","--regions","--interval"},description="default interval region on opening")
	private String defaultRegion="";
	@Parameter(names={"-q","--mapq"},description="min mapq")
	private int minmapq=1;
	
	private static class BamInfo {
		final Path bamPath;
		String sample;
		BamInfo(final Path bamPath) {
			this.bamPath = bamPath;
			}
		}
	
	@SuppressWarnings("serial")
	private static class XFrame extends JFrame {
		final Path referenceFile;
		final SAMSequenceDictionary dict;
		final List<BamInfo> bamPaths;
		final JTable jtableBams;
		final JPanel drawingArea;
		final JTextField jtextFieldLocation;
		final int minMapq;
		final SamReaderFactory srf;
		
		private class ZoomAction extends AbstractAction {
			final double factor;
			ZoomAction(double factor) {
				super("x"+factor);
				this.factor=factor;
				}
			@Override
			public void actionPerformed(final ActionEvent e)
				{
				final Optional<SimpleInterval> optR = getUserInterval();
				if(!optR.isPresent()) return;
				final SimpleInterval r= optR.get();
				final SAMSequenceRecord ssr=dict.getSequence(r.getContig());
				if(ssr==null) return;
				final double L= r.getLengthOnReference()/2.0*factor;
				final double mid = r.getStart()+r.getLengthOnReference()/2;
				int x1= (int)Math.max(1,mid-L);
				int x2= (int)Math.min(mid+L,ssr.getLengthOnReference());
				jtextFieldLocation.setText(r.getContig()+":"+x1+"-"+x2);
				drawingArea.repaint();
				}
			}
		
		
		XFrame(final Path referenceFile,final List<Path> bamPaths,String defaultLoc,int minMapq) {
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			setTitle(SwingBamCov.class.getSimpleName());
			this.referenceFile = referenceFile;
			this.srf = SamReaderFactory.makeDefault().
					referenceSequence(this.referenceFile).
					validationStringency(ValidationStringency.LENIENT);
			this.dict = SequenceDictionaryUtils.extractRequired(this.referenceFile);
			this.bamPaths = bamPaths.stream().map(P->new BamInfo(P)).
					collect(Collectors.toCollection(Vector::new));
					new Vector<BamInfo>(bamPaths.size());
			for(final BamInfo bi:this.bamPaths) {
				try(SamReader sr=this.srf.open(bi.bamPath)) {
					final SAMFileHeader header = sr.getFileHeader();
					bi.sample = header.getReadGroups().stream().
							map(S->S.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(bi.bamPath));
					}
				catch(IOException err) {
					throw new RuntimeIOException(err);
					}
				}
			Collections.sort(this.bamPaths,(A,B)->A.sample.compareTo(B.sample));
					
			this.minMapq=minMapq;
			
			final JPanel mainPane = new JPanel(new BorderLayout(5, 5));
			mainPane.setBorder(new EmptyBorder(5, 5, 5, 5));
			this.setContentPane(mainPane);
			final Panel topPane = new Panel(new FlowLayout(FlowLayout.LEADING));
			mainPane.add(topPane,BorderLayout.NORTH);
			JLabel label = new JLabel("Location:", JLabel.RIGHT);
			topPane.add(label);
			jtextFieldLocation= new JTextField(defaultLoc,20);
			if(StringUtil.isBlank(defaultLoc) &&  !dict.isEmpty()) {
				final SAMSequenceRecord ssr= dict.getSequence(0);
				jtextFieldLocation.setText(ssr.getSequenceName()+":1-1000");
				}
			
			topPane.add(jtextFieldLocation);
			label.setLabelFor(jtextFieldLocation);
			final Action actionGo = new AbstractAction("Plot")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					drawingArea.repaint();
					}
				};
			this.jtextFieldLocation.addActionListener(actionGo);
			JButton button = new JButton(actionGo);
			topPane.add(button);
			topPane.add(new JSeparator());
			topPane.add(new JButton(new ZoomAction(0.5)));
			topPane.add(new JButton(new ZoomAction(1.5)));
			topPane.add(new JButton(new ZoomAction(2.0)));
			topPane.add(new JButton(new ZoomAction(10)));
			
			
			this.drawingArea = new JPanel(true) {
				@Override
				protected void paintComponent(Graphics g)
					{
					paintDrawingArea(Graphics2D.class.cast(g));
					}
				};
			this.drawingArea.setOpaque(true);
			this.drawingArea.setBackground(Color.WHITE);
			
			this.jtableBams =  new JTable(new AbstractTableModel()
				{
				public boolean isCellEditable(int rowIndex, int columnIndex) {return false;};
				@Override
				public Class<?> getColumnClass(int columnIndex)
					{
					return String.class;
					}
				@Override
				public String getColumnName(int column)
					{
					switch(column) {
						case 0: return "Sample";
						case 1: return "Path";
						default: throw new IllegalArgumentException();
						}
					}
				@Override
				public Object getValueAt(int rowIndex, int column)
					{
					switch(column) {
						case 0: return XFrame.this.bamPaths.get(rowIndex).sample;
						case 1: return XFrame.this.bamPaths.get(rowIndex).bamPath.toString();
						default: throw new IllegalArgumentException();
						}
					}
				@Override
				public int getRowCount() { return bamPaths.size(); }
				@Override
				public int getColumnCount() {return 2;}
				});//end new jtables
			this.jtableBams.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
			final JTabbedPane  tabPane = new JTabbedPane();
			mainPane.add(tabPane,BorderLayout.CENTER);
			tabPane.addTab("Paint", drawingArea);
			tabPane.addTab("BAMS", new JScrollPane(this.jtableBams));
			tabPane.addTab("REF", new JScrollPane( new JTable(new SwingSequenceDictionaryTableModel(dict))));
			if(!StringUtil.isBlank(defaultLoc)) {
				this.addWindowListener(new WindowAdapter()
					{
					@Override
					public void windowOpened(final WindowEvent e)
						{
						drawingArea.repaint();
						removeWindowListener(this);
						}
					});
				}
			final JMenuBar menuBar= new JMenuBar();
			setJMenuBar(menuBar);
			final JMenu menu = new JMenu("File");
			menuBar.add(menu);
			menu.add(new JSeparator());
			menu.add(new JMenuItem(new AbstractAction("Save as...")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					doMenuSaveAs();
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
		
		private void doMenuSaveAs() {
			JFileChooser chooser= new JFileChooser();
			if(chooser.showSaveDialog(drawingArea)!=JFileChooser.APPROVE_OPTION) return;
			final File f = chooser.getSelectedFile();
			if(f.exists() && JOptionPane.showConfirmDialog(this, "File "+f.getName()+" exists. Overwite ?", "Overwite ?", JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE, null)!=JOptionPane.OK_OPTION)
				{
				return;
				}
			try {
				BufferedImage img = new BufferedImage(drawingArea.getWidth(), drawingArea.getHeight(), BufferedImage.TYPE_INT_RGB);
				Graphics2D g = img.createGraphics();
				paintComponents(g);
				g.dispose();
				ImageIO.write(img,f.getName().toLowerCase().endsWith(".png")?"PNG":"JPG", f);
				}
			catch(Throwable err ) {
				LOG.error(err);
				return;
				}
			}
		
		private void paintDrawingArea(final Graphics2D g) {
			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, drawingArea.getWidth(), drawingArea.getHeight());
			int marginTop=12;
			
			final Optional<SimpleInterval> location = getUserInterval();
			if(!location.isPresent()) return;
			g.setColor(Color.DARK_GRAY);
			g.drawString(location.get().toNiceString()+" length:"+ StringUtils.niceInt(location.get().getLengthOnReference())+"bp " +
					SequenceDictionaryUtils.getBuildName(dict).orElse(referenceFile.getFileName().toString())+" MAPQ:"+this.minMapq, 3, marginTop-1);
			
			final List<BamInfo> paths;
			if(this.jtableBams.getSelectedRowCount()==0) {
				paths = this.bamPaths;
				} else {
				int[] selidx = this.jtableBams.getSelectedRows();
				paths = new Vector<>(selidx.length);
				for(int idx:selidx) {
					paths.add(bamPaths.get(idx));
					}
				}
			if(paths.isEmpty()) return;
			int ncols = (int)Math.max(1,Math.floor(Math.sqrt(paths.size())));
			int nrows = (int)Math.max(1, Math.ceil(paths.size()/(double)ncols));
			double wi  = drawingArea.getWidth()/ncols;
			double hi  = (drawingArea.getHeight()-marginTop)/nrows;
			int nr=0;
			int nc=0;
			
			for(final BamInfo bamPath: paths) {
				paintBam(g,srf,bamPath,location.get(),new Rectangle2D.Double(
						nc*wi,
						marginTop+nr*hi,
						wi,
						hi));
				nc++;
				if(nc==ncols) {
					nc=0;
					nr++;
					}
				}
			}
		
		private void paintBam(final Graphics2D g,
				final SamReaderFactory srf,
				final BamInfo bam,
				final Locatable loc
				,final Rectangle2D rect) {
			try {
				try(SamReader sr=this.srf.open(bam.bamPath)) {
					if(!sr.hasIndex()) return;					
					final CoverageFactory covFactory = new CoverageFactory().
							setMappingQuality(this.minMapq);
					final CoverageFactory.SimpleCoverage cov = covFactory.getSimpleCoverage(sr, loc, bam.sample);
					
					final double[] depths = cov.scaleMedian((int)rect.getWidth());
					final double maxDepth = Arrays.stream(depths).max().orElse(1.0);
					final GeneralPath gp = new GeneralPath();
					gp.moveTo(rect.getX(), rect.getMaxY());
					for(int i=0;i< depths.length;i++) {
						double y = rect.getMaxY() - (depths[i]/maxDepth)*rect.getHeight();
						gp.lineTo(rect.getX()+i, y);
						}
					gp.lineTo(rect.getMaxX(),rect.getMaxY());
					gp.closePath();
					g.setColor(Color.LIGHT_GRAY);
					g.fill(gp);
					
					OptionalDouble mean = cov.getMedian();
					if(mean.isPresent()) {
						g.setColor(Color.RED);
						double y = rect.getMaxY() -(mean.getAsDouble()/maxDepth)*rect.getHeight();
						g.draw(new Line2D.Double(rect.getX(), y, rect.getMaxX(), y));
						}
					
					mean = cov.getAverage();
					if(mean.isPresent()) {
						g.setColor(Color.GREEN);
						double y = rect.getMaxY() -(mean.getAsDouble()/maxDepth)*rect.getHeight();
						g.draw(new Line2D.Double(rect.getX(), y, rect.getMaxX(), y));
						}
					
					g.setColor(Color.BLUE);
					int fontSize=12;
					g.setFont(new Font("Courier",Font.PLAIN,fontSize));
					g.drawString(bam.sample+"(max DP: " + (int)maxDepth+")",(int)rect.getX()+5,(int)rect.getMaxY()-5);
					//frame
					g.setColor(Color.DARK_GRAY);
					g.draw(rect);
					}
				}
			catch(final Throwable err) {
				LOG.error(err);
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		try {
			final List<Path> paths = IOUtils.unrollPaths(args);
			if(paths.isEmpty()) {
				LOG.error("no BAM/CRAM was provided");
				return -1;
				}
			JFrame.setDefaultLookAndFeelDecorated(true);
			final XFrame frame = new XFrame(this.referenceFile,paths,defaultRegion,this.minmapq);
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
	public static void main(String[] args)
		{
		new SwingBamCov().instanceMain(args);//no exit
		}

	}
