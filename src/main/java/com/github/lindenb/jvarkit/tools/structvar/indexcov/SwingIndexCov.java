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
import java.awt.Cursor;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentAdapter;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Line2D.Double;
import java.awt.image.BufferedImage;
import java.awt.LinearGradientPaint;
import java.awt.Paint;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
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
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileFilter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.net.UrlSupplier;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.PreferredDirectory;
import com.github.lindenb.jvarkit.util.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.TabixIteratorLineReader;
import htsjdk.tribble.readers.TabixReader;

public class SwingIndexCov extends Launcher {
	private static final Logger LOG = Logger.build(SwingIndexCov.class).make();
	@Parameter(names={"-R","--reference"},description=DICTIONARY_SOURCE)
	private Path dictSource = null;
	@Parameter(names={"-t","--treshold"},description=IndexCovUtils.TRESHOLD_OPT_DESC+". " + FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double treshold = IndexCovUtils.DEFAULT_TRESHOLD;
	@Parameter(names={"--gtf","--gff"},description="GFF3 file indexed with tabix to plot the genes.")
	private String gff3Path = null;
	
	
	@SuppressWarnings("serial")
	private static class XFrame extends JFrame {
		private final List<Locatable> history = new ArrayList<>();
		private int index_in_history = -1;
		/** original srcDict */
		private final SAMSequenceDictionary srcDict;
		private final String tabixFileUri;
		private final JList<Sample> jlistSamples;
		private final JComboBox<String> jcomboxChromosomes;
		private final JPanel drawingArea;
		private BufferedImage offscreen = null;
		private Plotter currentThread = null;
		private final JTextField jtextFieldLocation;
		private final JCheckBox jcboxPerSample;
		private final JMenu jmenuHyperlinks;
		private final double treshold;
		private final String gff3Path;
		private static class Sample {
			final String name;
			final int column;
			Color color;
			Sample(final String name,int column) {
				this.name= name;
				this.column = column;
				this.color= Color.LIGHT_GRAY;
				}
			@Override
			public String toString() {
				return this.name;
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
				final SAMSequenceRecord ssr=srcDict.getSequence(r.getContig());
				if(ssr==null) return;
				r= change(ssr,r);
				if(r==null) return;
				jtextFieldLocation.setText(r.getContig()+":"+r.getStart()+"-"+r.getEnd());
				XFrame.this.offscreen=null;
				pushHistory();
				intervalChanged();
				drawingArea.repaint();
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
			Locatable change(final SAMSequenceRecord ssr, Locatable r)
				{
				final double L= Math.max(15_000/2.0,r.getLengthOnReference()/2.0*factor);
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
		
		private class Plotter extends Thread {
			int width;
			int height;
			String tabixUri;
			List<Sample> samples;
			Locatable interval;
			boolean per_sample=false;
			transient boolean abort_flag=false;
			final double min_cov = 0.0;
			final double max_cov = 3.0;
			final int top_margin=30;
			IndexCovUtils indexCovUtils=null;
			String gff3Path = null;
			
			private BufferedImage createImage() {
				BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
				return img;
			}
			
			
			private void updateDrawingArea(BufferedImage img) {
				if(this.abort_flag) return;
				SwingUtilities.invokeLater(()->{
					offscreen = img;
					drawingArea.repaint();
				});
			}
			
			
			private double pos2pixel(int pos) {
				return ((pos-(double)this.interval.getStart())/this.interval.getLengthOnReference())*this.width;
				}
			
			private void plotGenes(Graphics2D g,BufferedImage img) {
					if(this.gff3Path==null) return;
					if(this.interval==null || this.interval.getLengthOnReference()> 10_000_000) return;
					final double midy= 20;
					long now = System.currentTimeMillis();
					g.setColor(Color.MAGENTA);
					try (TabixReader tbr = new TabixReader(this.gff3Path)) {						
						final Gff3Codec codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
						final TabixReader.Iterator iter0=tbr.query(interval.getContig(),interval.getStart(),interval.getEnd());
						final LineIterator lr = new LineIteratorImpl( new TabixIteratorLineReader(iter0));
						while(!codec.isDone(lr) && !this.abort_flag) {
								final Gff3Feature gffline = codec.decode(lr);
								if(gffline==null) continue;
								final double h;
								if(gffline.getType().equals("transcript"))
									{
									h=1;
									}
								else if(gffline.getType().equals("exon")) {
									h=3;
									}
								else {
									continue;
									}
								double x0 = pos2pixel(gffline.getStart());
								double x1 = pos2pixel(gffline.getEnd());
								g.fill(new Rectangle2D.Double(x0,midy-h/2.0,(x1-x0),h));
								if(System.currentTimeMillis()-now>1000L) {
									updateDrawingArea(img);
									now = System.currentTimeMillis();
									}
								}
						codec.close(lr);
						}
					catch(Throwable err) {
						
						}
					}
			
			
			@Override
			public void run() {
				final BufferedImage img = createImage();
				final Graphics2D g = img.createGraphics();
				g.setColor(Color.WHITE);
				g.fillRect(0, 0, this.width, this.height);
				final double sample_height =
						this.per_sample?
						(this.height-top_margin)/(double)samples.size():
						(this.height-top_margin);
				final Hershey hershey = new Hershey();
				if(this.samples.isEmpty() || this.interval==null || this.width<=0 || sample_height<=1) {
					if(sample_height<=1) {
						g.setColor(Color.BLACK);
						final String title="Too many samples N="+this.samples.size()+" height:"+sample_height;
						hershey.paint(g,title,0,10,title.length()*12,12);
						}
					updateDrawingArea(img);
					return;
					}
				plotGenes(g,img);
				final String title = new SimpleInterval(this.interval).toNiceString()+" length: "+ StringUtils.niceInt(this.interval.getLengthOnReference())+" bp";
				g.setColor(Color.DARK_GRAY);
				hershey.paint(g, title, 0,0,title.length()*10,10);
				double prev_x=0;
				double[] prev_y=null;
				long now = System.currentTimeMillis();
				
				
				if(this.per_sample) {
					final double fontSize=Math.min(12,sample_height/10.0);
					for(int i=0;i< this.samples.size() && !abort_flag;i++) {
						final String sn = samples.get(i).name;
						g.setColor(samples.get(i).color);
						hershey.paint(g,sn,0,top_margin + sample_height*i+1.0, sn.length()*fontSize,fontSize);
						g.setColor(Color.GRAY);
						for(double z=0.0;z<=2.0;z+=0.5) {
							final double y = top_margin + (sample_height*(i+1)) - z * sample_height;
							
							g.draw(new Line2D.Double(0,y,this.width, y));
							}
						}
					}
				else
					{
					g.setColor(Color.GRAY);
					for(double z=0.0;z<=2.0;z+=0.5) {
						final String sn = String.valueOf(z);
						final double y = this.height - ((z - min_cov)/(max_cov-min_cov)) * sample_height;
						hershey.paint(g,sn,0,y, sn.length()*10,10);
						g.draw(new Line2D.Double(0,y,this.width, y));
						}
					}
				
				try(TabixFileReader tfr = new TabixFileReader(this.tabixUri)) {
					final Iterator<String> iter = tfr.iterator(this.interval.getContig(), this.interval.getStart(), this.interval.getEnd());
					while(!abort_flag && iter.hasNext()) {
						final String line = iter.next();
						if(line==null) break;
						if(line.startsWith("#")) continue;
						final String[] tokens = CharSplitter.TAB.split(line);
						final Locatable loc =  new SimpleInterval(tokens[0], Integer.parseInt(tokens[1]),  Integer.parseInt(tokens[2]));
						
						double x0 = pos2pixel(loc.getStart());
						double x1 = pos2pixel(loc.getEnd());
						double[] array_y = new double[samples.size()];;

						//special point ?
						
						int nHetDel=0;
						int nHomDel=0;
						int nHetDup=0;
						int nHomDup=0;
						int nOther=0;

						for(int i=0;i< samples.size() && !abort_flag;i++) {
							final Sample sn = samples.get(i);
							final float value = Float.parseFloat(tokens[sn.column]);
							if(indexCovUtils.isHomDel(value)) {
								nHomDel++;
								}
							else if(indexCovUtils.isHetDel(value)) {
								nHetDel++;
								}
							else if(indexCovUtils.isHetDup(value)) {
								nHetDup++;
								}
							else if(indexCovUtils.isHomDup(value)) {
								nHomDup++;
								}
							else {
								nOther++;
								}
							
							final double frac = ((Math.min(value,max_cov) - min_cov)/(max_cov-min_cov));
							final double y;
							
							if(this.per_sample) {
								final double bottom_y = top_margin + (sample_height*(i+1)); 
								y= bottom_y - frac*sample_height;
								
								final Point2D p_start = new Point2D.Double(0, top_margin + (sample_height*(i+1)));
     								final Point2D p_end = new Point2D.Double(0, top_margin + (sample_height*(i)));
							     	final float[] p_dist = {
							     		(float)((0.0-min_cov)/(max_cov-min_cov)),
							     		(float)((0.5-min_cov)/(max_cov-min_cov)),
							     		(float)((1.0-min_cov)/(max_cov-min_cov)),
							     		(float)((1.5-min_cov)/(max_cov-min_cov)),
							     		(float)((2.0-min_cov)/(max_cov-min_cov))
							     		};
							     	final Color[] p_colors = {Color.BLUE, Color.CYAN,Color.LIGHT_GRAY,Color.ORANGE, Color.RED};
							     	final LinearGradientPaint grad = new LinearGradientPaint(p_start,p_end, p_dist, p_colors);
							     	final Paint oldPaint = g.getPaint();
								g.setPaint(grad);
								final Rectangle2D rec = new Rectangle2D.Double(x0,y, (x1-x0),bottom_y-y);
								g.fill(rec);
								g.setPaint(oldPaint);
								if(x1-x0>3) {
									g.setColor(samples.get(i).color);
									g.draw(rec);
									}
								}
							else
								{
								g.setColor(samples.get(i).color);
								y = this.height - ((value - min_cov)/(max_cov-min_cov)) * sample_height;
								if(prev_y!=null) {
									g.draw(new Line2D.Double(prev_x, prev_y[i], prev_x, y));
									}
								g.draw(new Line2D.Double(x0, y, x1, y));
								array_y[i] = y;
								}
							}
						if(this.samples.size()>1) {
							Color c=null;
							if(nHetDel==1 && nHomDel==0 && nHetDup==0 && nHomDup==0) {
								c=Color.CYAN;
								}
							else if(nHetDel==0 && nHomDel==0 && nHetDup==1 && nHomDup==0) {
								c=Color.ORANGE;
								}
							if(c!=null) {
								final double mid = (x0+x1)/2.0;
								final double y2 = top_margin/2.0;
								g.setColor(c);
								g.fill(new Ellipse2D.Double(mid-2.5,y2-2.5,5,5));
								}
							}
						if(System.currentTimeMillis()-now>1000L) {
							updateDrawingArea(img);
							now = System.currentTimeMillis();
							}
						prev_x= x1;
						prev_y = array_y;
						}
					}
				catch(final IOException err) {
					err.printStackTrace();
					}
				updateDrawingArea(img);
			}
		}
		
		
		
		XFrame(final SAMSequenceDictionary srcDict,final String tabixFileUri,final double treshold,final String gff3path) throws IOException {
			super(SwingIndexCov.class.getSimpleName());
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			this.tabixFileUri = tabixFileUri;
			this.treshold = treshold;
			this.gff3Path = gff3path;
			final JMenuBar menuBar = new JMenuBar();
			setJMenuBar(menuBar);
			JMenu menu=new JMenu("File");
			menuBar.add(menu);
			menu.add(new JMenuItem(new AbstractAction("Open...") {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					}
				}));
			menu.add(new JMenuItem(new AbstractAction("Save As...") {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					doMenuSaveAs();
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
			menuBar.add(this.jmenuHyperlinks = new JMenu("Hyperlinks"));
			
			this.srcDict = srcDict;
			
			try(TabixFileReader tbr = new TabixFileReader(tabixFileUri)) {
				//read samples
				final String line = tbr.readLine();
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
				
				
				final Vector<Sample> vSamples = new Vector<>(tokens.length);
				for(int i=3;i< tokens.length;i++) {
					final Sample sn = new Sample(tokens[i], i);
					// generate colors : https://stackoverflow.com/questions/223971
					sn.color = Color.getHSBColor((float) (i-3) / (float)(tokens.length-3), 0.85f, 1.0f);
					vSamples.add(sn);
					}
				this.jlistSamples = new JList<Sample>(vSamples);
				this.jlistSamples.getSelectionModel().setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
				this.jlistSamples.addListSelectionListener(AE->samplesSelectionchanged());
				this.jlistSamples.setCellRenderer(new DefaultListCellRenderer() {
					public Component getListCellRendererComponent(javax.swing.JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
						final Component c = super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
						if(value!=null && value instanceof Sample) {
							final Sample sn =Sample.class.cast(value);
							if(c instanceof JComponent) {
								JComponent.class.cast(c).setToolTipText(sn.name);
								}
							if(!isSelected) {
								c.setForeground(sn.color);
								}
						}
						return c;
					};
				});
				
				// get chromosomes in the tabix
				final Vector<String> vecChroms = new Vector<>( tbr.getChromosomes());
				//remove if not in dict
				vecChroms.removeIf(S->this.srcDict.getSequence(S)==null);
				if(vecChroms.isEmpty()) {
					throw new IllegalArgumentException("no common chromosomes between REF and "+tabixFileUri);
					}
				//sort on dict
				vecChroms.sort((A,B)->Integer.compare(
						this.srcDict.getSequenceIndex(A),
						this.srcDict.getSequenceIndex(B))
						);
				this.jcomboxChromosomes = new JComboBox<>(vecChroms);
				this.jcomboxChromosomes.setSelectedIndex(0);
				this.jcomboxChromosomes.setEditable(false);
				this.jcomboxChromosomes.addActionListener((AE)->chromSelectionchanged());
				}
			
			final JPanel mainPanel=new JPanel(new BorderLayout());
			this.setContentPane(mainPanel);
			
			final JTabbedPane jTabbedPane = new JTabbedPane();
			mainPanel.add(jTabbedPane,BorderLayout.CENTER);
			final JPanel viewPane= new JPanel(new BorderLayout());
			jTabbedPane.addTab("View", viewPane);
			
			/** navigation */
			final JPanel navPane = new JPanel(new FlowLayout(FlowLayout.LEADING));
			viewPane.add(navPane,BorderLayout.NORTH);
			
			
			JLabel label = new JLabel("Location:", JLabel.RIGHT);
			navPane.add(label);
			this.history.add(this.srcDict.getSequence(0));
			this.jtextFieldLocation= new JTextField(this.history.get(0).getContig(),10);
			this.history.add(this.srcDict.getSequence(0));
			this.index_in_history=0;
			navPane.add(this.jtextFieldLocation);
			label.setLabelFor(this.jtextFieldLocation);
			final Action actionGo = new AbstractAction("Load") {
				@Override
				public void actionPerformed(ActionEvent e)
					{
					pushHistory();
					intervalChanged();
					offscreen=null;
					drawingArea.repaint();
					}
				};
			navPane.add(new JButton(actionGo));
			this.jtextFieldLocation.addActionListener(actionGo);
			navPane.add(new JSeparator());
			label = new JLabel("Chroms:", JLabel.RIGHT);
			navPane.add(label);
			navPane.add(this.jcomboxChromosomes);
			label.setLabelFor(this.jcomboxChromosomes);
			navPane.add(new JSeparator());
			this.jcboxPerSample = new JCheckBox("Per Sample", false);
			this.jcboxPerSample.addActionListener((AE)->samplesSelectionchanged());//will call repaint anyway...
			navPane.add(this.jcboxPerSample);
			navPane.add(new JSeparator());
			final JComboBox<ZoomAction> zoomComboBox = new JComboBox<>(
						new ZoomAction[] {
						new ZoomAction(0.1),
						new ZoomAction(0.25),
						new ZoomAction(0.5),
						new ZoomAction(0.75),
						new ZoomAction(2.0),
						new ZoomAction(3.0),
						new ZoomAction(5.0),
						new ZoomAction(10.0),
						new ZoomAction(20.0),
						new ZoomAction(50.0),
						new ZoomAction(100.0)
						}
					);
			zoomComboBox.setEditable(false);
			zoomComboBox.setSelectedIndex(0);
			navPane.add(zoomComboBox);
			final JButton zoomButton = new JButton(new AbstractAction("Zoom") {
				@Override
				public void actionPerformed(ActionEvent e)
					{
					final ZoomAction za = (ZoomAction)zoomComboBox.getSelectedItem();
					if(za==null) return;
					za.actionPerformed(e);
					}
				});
			navPane.add(zoomButton);
			navPane.add(new JSeparator());
			navPane.add(new JButton(new ShiftAction("<<<",-0.9)));
			navPane.add(new JButton(new ShiftAction("<<",-0.5)));
			navPane.add(new JButton(new ShiftAction("<",-0.1)));
			navPane.add(new JButton(new ShiftAction(">",0.1)));
			navPane.add(new JButton(new ShiftAction(">>",0.5)));
			navPane.add(new JButton(new ShiftAction(">>>",0.9)));
			navPane.add(new JSeparator());
			navPane.add(new JButton(new AbstractAction("Prev") {	
				@Override
				public void actionPerformed(ActionEvent arg0) {
					if(index_in_history>0) {
						index_in_history--;
						final Locatable loc = history.get(index_in_history);
						jtextFieldLocation.setText(loc.getContig()+":"+loc.getStart()+"-"+loc.getEnd());
						intervalChanged();
						offscreen=null;
						drawingArea.repaint();
						}
					}
				}));
			navPane.add(new JButton(new AbstractAction("Next") {	
				@Override
				public void actionPerformed(ActionEvent arg0) {
					if(index_in_history+1< history.size()) {
						index_in_history++;
						final Locatable loc = history.get(index_in_history);
						jtextFieldLocation.setText(loc.getContig()+":"+loc.getStart()+"-"+loc.getEnd());
						intervalChanged();
						offscreen=null;
						drawingArea.repaint();
						}
					}
				}));
			
			JPanel pane2 = new JPanel(new BorderLayout(5, 5));
			pane2.setBorder(BorderFactory.createTitledBorder("Samples"));
			pane2.add(new JScrollPane(this.jlistSamples),BorderLayout.CENTER);
			pane2.setPreferredSize(new Dimension(200,10));
			viewPane.add(pane2,BorderLayout.WEST);
			
			/** drawing area itself */
			this.drawingArea = new JPanel(null, true) {
				@Override
				public String getToolTipText(final MouseEvent event)	{
					final SimpleInterval loc =getUserInterval().orElse(null);
					if(loc==null) return "";
					int x1 = loc.getStart() + (int)((event.getX()/(double)drawingArea.getWidth())*loc.getLengthOnReference());
					return loc.getContig()+":"+StringUtils.niceInt(x1);
					}
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
			
			viewPane.add(this.drawingArea,BorderLayout.CENTER);
			
			/** ref tab */
			jTabbedPane.addTab("REF", new JScrollPane( new JTable(new SwingSequenceDictionaryTableModel(this.srcDict))));
			/* open close operations */
			this.addWindowListener(new WindowAdapter() {
				public void windowOpened(java.awt.event.WindowEvent e) {
					drawingArea.repaint();
					}
				public void windowClosed(java.awt.event.WindowEvent e)
					{
					runOnClose();
					}
				});
			/** mouse listener */
			final MouseAdapter mouse= new MouseAdapter() {
				Integer mouseStart = null;
				Integer mousePrev = null;
				public void mousePressed(final MouseEvent e) 
					{
					mouseStart = e.getX();
					mousePrev = null;
					e.getComponent().setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
					}
				public void mouseDragged(final MouseEvent e) {
					final Component c= e.getComponent();
					final Graphics2D g=(Graphics2D)c.getGraphics();
					g.setXORMode(g.getBackground());
					
					if(mousePrev!=null) {
							g.fillRect(
									Math.min(mouseStart, mousePrev),
									0,
									Math.abs(mouseStart-mousePrev),
									c.getHeight()
									);
							}
					mousePrev= e.getX();
					g.setPaintMode();
					}
				public void mouseReleased(final MouseEvent e){
					final Component c= e.getComponent();
					c.setCursor(Cursor.getDefaultCursor());
					if(mousePrev==null) return;
					Optional<SimpleInterval> optR = getUserInterval();
					if(!optR.isPresent())	 return ;
					int pixMin = Math.min(mouseStart,mousePrev);
					int pixMax = Math.max(mouseStart,mousePrev);
					final Locatable loc = optR.get();
					final int pos1 = (int)( loc.getStart()+ ((pixMin)/(double)c.getWidth())*loc.getLengthOnReference());
					int pos2 = (int)( loc.getStart()+ ((pixMax)/(double)c.getWidth())*loc.getLengthOnReference());
					if(pos1>=pos2) pos2=pos1+1;
					jtextFieldLocation.setText(loc.getContig()+":"+pos1+"-"+pos2);
					pushHistory();
					intervalChanged();
					XFrame.this.offscreen=null;
					mouseStart=null;
					mousePrev=null;
					XFrame.this.drawingArea.repaint();
					}
				};
			this.drawingArea.addMouseListener(mouse);
			this.drawingArea.addMouseMotionListener(mouse);
			this.drawingArea.setToolTipText("");
			}
		
		private void doMenuSaveAs() {
			final JFileChooser chooser= new JFileChooser(PreferredDirectory.get(SwingIndexCov.class));
			if(chooser.showSaveDialog(drawingArea)!=JFileChooser.APPROVE_OPTION) return;
			final File f = chooser.getSelectedFile();
			if(f.exists() && JOptionPane.showConfirmDialog(this, "File "+f.getName()+" exists. Overwite ?", "Overwite ?", JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE, null)!=JOptionPane.OK_OPTION)
				{
				return;
				}
			PreferredDirectory.update(f);
			try {
				final BufferedImage img = new BufferedImage(drawingArea.getWidth(), drawingArea.getHeight(), BufferedImage.TYPE_INT_RGB);
				final Graphics2D g = img.createGraphics();
				g.setColor(Color.WHITE);
				g.fillRect(0, 0, img.getWidth(), img.getHeight());
				drawingArea.paintComponents(g);
				g.dispose();
				ImageIO.write(img,f.getName().toLowerCase().endsWith(".png")?"PNG":"JPG", f);
				}
			catch(final Throwable err ) {
				ThrowablePane.show(XFrame.this, err);
				}
			}
		
		private void intervalChanged() {
			this.jmenuHyperlinks.removeAll();
			final Optional<SimpleInterval> optR = getUserInterval();
			if(optR.isPresent()) {
				if(!Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.BROWSE)) return;
				final UrlSupplier urlSupplier = new UrlSupplier(this.srcDict);
				urlSupplier.of(optR.get()).stream().forEach(U->{
					final AbstractAction action = new AbstractAction(U.getLabel())
						{
						@Override
						public void actionPerformed(final ActionEvent e)
							{
							try {
								Desktop.getDesktop().browse(new URI(U.getUrl()));
								}
							catch(final Throwable err) {
								ThrowablePane.show(XFrame.this, err);
								}
							}
						};
					action.putValue(AbstractAction.LONG_DESCRIPTION, U.getLabel());
					action.putValue(AbstractAction.SHORT_DESCRIPTION, U.getLabel());
					action.putValue(AbstractAction.NAME, U.getLabel());
					final JMenuItem mi = new JMenuItem(action);
					this.jmenuHyperlinks.add(mi);
					});
				}
			}
		
		private void samplesSelectionchanged() {
			this.offscreen=null;
			this.drawingArea.repaint();
			}
		private void chromSelectionchanged() {
			this.offscreen=null;
			int i=this.jcomboxChromosomes.getSelectedIndex();
			String loc="";
			if(i>=0) {
				final String ctg=this.jcomboxChromosomes.getModel().getElementAt(i);
				final SAMSequenceRecord ssr = this.srcDict.getSequence(ctg);
				loc=ssr.getContig()+":1-"+ssr.getSequenceLength();
				}
			this.jtextFieldLocation.setText(loc);
			pushHistory();
			intervalChanged();
			this.drawingArea.repaint();
			}
		
		private void runOnClose() {
			stopDrawingThread();
 			}
		
		private void paintDrawingArea(final Graphics2D g) {
			if(drawingArea.getWidth()==0) return;
			if(drawingArea.getHeight()==0) return;
			if(offscreen==null || offscreen.getWidth()!=drawingArea.getWidth() ||  offscreen.getWidth()!=drawingArea.getWidth()) {
				startDrawingThread();
				return;
				}
			g.drawImage(offscreen,0,0,null);
			}

		
		private synchronized void stopDrawingThread() {
			if(currentThread!=null) {
				currentThread.abort_flag=true;
				}
			currentThread=null;
			}
		
		private synchronized void startDrawingThread() {
			stopDrawingThread();
			final Plotter plotter = new Plotter();
			plotter.indexCovUtils = new IndexCovUtils(this.treshold);
			plotter.tabixUri = this.tabixFileUri;
			plotter.width = this.drawingArea.getWidth();
			plotter.height = this.drawingArea.getHeight();
			plotter.per_sample = this.jcboxPerSample.isSelected();
			plotter.interval = getUserInterval().orElse(null);
			plotter.gff3Path = this.gff3Path;
			final List<Sample> L2 = jlistSamples.getSelectedValuesList();
			plotter.samples = new ArrayList<>();
			if(L2==null || L2.isEmpty()) {
				for(int i=0;i<jlistSamples.getModel().getSize();i++) {
					plotter.samples.add(jlistSamples.getModel().getElementAt(i));
					}
				}
			else
				{
				final Set<String> set = L2.stream().map(S->S.name).collect(Collectors.toSet());
				for(int i=0;i<jlistSamples.getModel().getSize();i++) {
					final Sample sn = jlistSamples.getModel().getElementAt(i);
					if(!set.contains(sn.name)) continue;
					plotter.samples.add(sn);
					}
				}
			this.currentThread = plotter;
			this.currentThread.start();
			}
		
		
		private final Optional<SimpleInterval> getUserInterval() {
			final String s = this.jtextFieldLocation.getText().trim();
			if(StringUtil.isBlank(s)) return Optional.empty();			
			try {
				return IntervalParserFactory.newInstance().
					dictionary(this.srcDict).
					enableWholeContig().
					make().
					apply(s);
				}
			catch(final Throwable err) {
				return Optional.empty();
				}
			}
		
		
		private void pushHistory() {
			final Optional<SimpleInterval> optR=  getUserInterval();
			if(!optR.isPresent()) return;
			if(this.index_in_history+1<this.history.size()) {
				this.history.subList(this.index_in_history+1, this.history.size()).clear();
				}
			this.history.add(optR.get());
			this.index_in_history++;
			this.offscreen=null;
			this.drawingArea.repaint();
			}
		}
	
	
	
	private static SAMSequenceDictionary buildDefaultDict(final Component owner,final String uri) {
		LOG.warn("building default dict for "+uri);
		final Map<String, Integer> contig2len = new LinkedHashMap<>();
		try(TabixFileReader tbr= new TabixFileReader(uri)) {
			for(;;) {
				final String line = tbr.readLine();
				if(line==null) break;
				if(line.startsWith("#")) continue;
				final String[] tokens = CharSplitter.TAB.split(line, 4);
				final String ctg = tokens[0];
				final int len = Math.max(contig2len.getOrDefault(ctg, 0),Integer.parseInt(tokens[2]));
				contig2len.put(ctg, len);
				}
			return new SAMSequenceDictionary(
					contig2len.entrySet().
						stream().
						map(KV->new SAMSequenceRecord(KV.getKey(), KV.getValue())).
						collect(Collectors.toList())
					);
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = oneFileOrNull(args);
			JFrame.setDefaultLookAndFeelDecorated(true);
			final File inputFile;
			if(input==null)
				{
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
							if(!fname.endsWith(".bed.gz")) return false;
							return TabixFileReader.isValidTabixFile(f);
							}
						return false;
					}
				});
				if(jfc.showOpenDialog(null)!=JFileChooser.APPROVE_OPTION) return -1;
				inputFile = jfc.getSelectedFile();
				if(inputFile==null) return -1;
				}
			else
				{
				inputFile= new File(input);
				}
			if(!TabixFileReader.isValidTabixFile(inputFile)) {
				LOG.error("not a valid tabix file :"+inputFile);
				return -1;
				}
				
			final SAMSequenceDictionary dict = this.dictSource==null?
				buildDefaultDict(null, input):
				SequenceDictionaryUtils.extractRequired(this.dictSource)
				;
			try {
				final XFrame frame = new XFrame(dict,input.toString(),this.treshold,this.gff3Path);
				SwingUtilities.invokeLater(()->{
					final Dimension dim= Toolkit.getDefaultToolkit().getScreenSize();
					frame.setBounds(50, 50, dim.width-100, dim.height-100);
					frame.setVisible(true);
				});
			} catch(final Throwable err) {
				ThrowablePane.show(null, err);
				}
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
