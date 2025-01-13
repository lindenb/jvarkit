/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
import java.awt.LinearGradientPaint;
import java.awt.Paint;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ComponentAdapter;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.Vector;
import java.util.stream.Collectors;

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
import com.github.lindenb.jvarkit.gff3.SwingGff3TableModel;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.UrlSupplier;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.swing.PreferredDirectory;
import com.github.lindenb.jvarkit.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.TabixIteratorLineReader;
import htsjdk.tribble.readers.TabixReader;
/**
BEGIN_DOC

## Examples:

```
java -jar dist/swingindexcov.jar indexcov.bed.gz
```
or 
```
java -jar dist/swingindexcov.jar -R reference.fa indexcov.bed.gz --gff indexed.gff3.gz --helper "/TMP/myhelper.jar MyHelper"
```

## Screenshots


![https://pbs.twimg.com/media/FSkLDJdXwAEmZTO?format=png&name=small](https://pbs.twimg.com/media/FSkLDJdXwAEmZTO?format=png&name=small
)

![https://pbs.twimg.com/media/FSkLIGtXMAAizeK?format=jpg&name=small](https://pbs.twimg.com/media/FSkLIGtXMAAizeK?format=jpg&name=small)

## Helper


### Helper Source code

file `MyHelper.java`

```java
import com.github.lindenb.jvarkit.tools.structvar.indexcov.SwingIndexCov;
import java.util.*;
import java.awt.Color;

public class MyHelper extends SwingIndexCov.DefaultHelper {
	private final Set<String> cases = new HashSet<>(Arrays.asList("sample1,sample2,sample3,sample4".split("[,]")));
	private final Set<Integer> cases_indexes = new HashSet<>();
	@Override
	public void initialize(final String[] header) {
		super.initialize(header);
		for(int i=3;i< header.length;i++) {
			if(cases.contains(getNormName(header[i]))) cases_indexes.add(i);
			}
		}
	private String getNormName(String sn) {
		final String[] tokens = sn.split("[_]");
		for(int i=0;i< tokens.length;i++) {
			if(tokens[i].length()==7) {
				return tokens[i];
				}
			}
		return sn;
		}
	@Override
	public String getDisplayName(final String sn) {
		String s = getNormName(sn);
		if(cases.contains(s)) return "*"+s;
		return s;
		}
	@Override
	public Optional<Colors> getColor(int[] indexes,float[] values) {
		int n_cases = 0;
		int n_other = 0;
		float treshold =  0.45f;
		for(int i=0;i< indexes.length;i++) {
			final int idx = indexes[i];
			float v  = values[i];
			boolean is_cnv = v> (1.0f + treshold) || v< (1.0f - treshold);
			if(!is_cnv) continue;
			boolean is_case = cases_indexes.contains(idx);
			if(is_case) {n_cases++;}
			else {n_other++;}
			}

		if(n_cases>0 && n_other==0) return Optional.of(Colors.RED);
		return Optional.empty();
		}
	}
```

### Compiling the helper

```make
CP=/home/lindenb/src/jvarkit-git/dist/swingindexcov.jar

myhelper.jar: MyHelper.java $(CP)
	rm -rf tmp
	mkdir -p tmp
	javac -d tmp -cp "$(CP):." $<
	jar cvf $@ -C tmp .
	rm -rf tmp
```

END_DOC
*/

@Program(
		name="swingindexcov",
		description="indexcov visualization",
		keywords={"cnv","duplication","deletion","sv"},
		creationDate="2020511",
		modificationDate="2020512",
		jvarkit_amalgamion = true,
		menu="CNV/SV"
		)
public class SwingIndexCov extends Launcher {
	private static final Logger LOG = Logger.build(SwingIndexCov.class).make();
	private static final int MAX_GFF3_ROWS = 10_000;
	@Parameter(names={"-R","--reference"},description=DICTIONARY_SOURCE)
	private Path dictSource = null;
	@Parameter(names={"--gtf","--gff","--gff3"},description="GFF3 file indexed with tabix to plot the genes.")
	private String gff3Path = null;
	@Parameter(names={"--helper"},description="For expert users only. java archive implenting Helper. Syntax \"path/to/helper.jar package.helper.implementation.Name\"")
	private String helperPath = null;

	
	public static abstract class Helper {
		public abstract void initialize(final String[] header);
		public String getDisplayName(final String sn) {
			return sn;
			}
		public Color getColorForSample(final String srcName) {
			return Color.BLACK;
			}
		public Optional<Color> getColor(int[] indexes,float[] values) {
			return Optional.empty();
			}
		public void dispose() {
			}
		}
	
	public static class DefaultHelper extends Helper {
		private final Map<String, Integer> sample2index=new HashMap<>();
		protected final IndexCovUtils indexCovUtils = new IndexCovUtils(IndexCovUtils.DEFAULT_TRESHOLD);
		private List<String> samples = Collections.emptyList();
		@Override
		public void initialize(final String[] header) {
			this.samples = Arrays.asList(header);
			for(int i=3;i<header.length;i++) {
				this.sample2index.put(header[i], i);
				}
			}
		protected String getSampleByColumn(int index) {
			if(index<3) throw new IllegalArgumentException();
			if(index>=samples.size()) throw new IllegalArgumentException();
			return this.samples.get(index);
			}
		
		protected int getNSamples() {
			return sample2index.size();
			}
		@Override
		public String getDisplayName(final String sn) {
			return sn;
			}
		@Override
		public Color getColorForSample(final String srcName) {
			final Integer  i= sample2index.getOrDefault(srcName,null);
			if(i==null) throw new IllegalStateException("srcName:"+srcName);
			if(i<3) throw new IllegalArgumentException();
			return  Color.getHSBColor((float) (i-3) / (float)getNSamples(), 0.85f, 1.0f);
			}
		@Override
		public Optional<Color> getColor(int[] indexes,float[] values) {
			if(indexes.length<=1) return Optional.empty();
			if(indexes.length!=values.length) throw new IllegalArgumentException("indexes.size!=values.size");
			
			
			int nHetDel=0;
			int nHomDel=0;
			int nHetDup=0;
			int nHomDup=0;
			int nOther=0;

			
			for(float value: values) {
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
				}
			if(nHetDel==1 && nHomDel==0 && nHetDup==0 && nHomDup==0) {
				return Optional.of(Color.CYAN);
				}
			else if(nHetDel==0 && nHomDel==0 && nHetDup==1 && nHomDup==0) {
				return Optional.of(Color.ORANGE);
				}
			return Optional.empty();
			}
		@Override
		public void dispose() {
			}
		}
	
	
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
		private final JCheckBox jcboxShowName;
		private final JMenu jmenuHyperlinks;
		private final String gff3Path;
		private final Helper helper;
		private final SwingGff3TableModel gff3TableModel;
		
		
		private static class Sample {
			final String srcName;
			final String displayName;
			final int column;
			Color color;
			Sample(final String srcName,final String displayName,int column) {
				this.srcName= srcName;
				this.displayName= displayName;
				this.column = column;
				this.color= Color.LIGHT_GRAY;
				}
			public String getDisplayName() {
				return this.displayName;
				}
			@Override
			public String toString() {
				return this.getDisplayName();
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
				updateHyperlinks();
				updateGff3Table();
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
		
		private static class BestPoint {
			Sample sample;
			double x;
			double y;
			float v;
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
			final double max_cov = 2.5;
			final int top_margin=30;
			String gff3Path = null;
			final Hershey _hershey = new Hershey();
			Helper helper = null;
			boolean showNames = false;
			
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
			
			private void plotGenes(final Graphics2D g,final BufferedImage img) {
					if(this.gff3Path==null || this.interval==null) return;
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
								if(	gffline.getType().equals("gene") ||
									gffline.getType().equals("transcript") ||
									gffline.getType().equals("processed_transcript"))
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
			
			private void drawText(Graphics2D g,double x,double y,String title,double fontHeight,int side) {
				double textWidth = title.length()*fontHeight;
				double dx =0;
				if(side==0) dx = -textWidth/2.0;
				this._hershey.paint(g, title, x +dx ,y,textWidth,fontHeight);
				}
			
			@Override
			public void run() {
				final Map<String,BestPoint> sample2best = new HashMap<>();
				final BufferedImage img = new BufferedImage(this.width,this.height, BufferedImage.TYPE_INT_RGB);
				final Graphics2D g = img.createGraphics();
				g.setColor(Color.WHITE);
				g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
				g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
				g.fillRect(0, 0, this.width, this.height);
				final double sample_height =
						this.per_sample?
						(this.height-top_margin)/(double)samples.size():
						(this.height-top_margin);
				if(this.samples.isEmpty() || this.interval==null || this.width<=0 || sample_height<=1) {
					if(sample_height<=1) {
						g.setColor(Color.BLACK);
						drawText(g,0,10,"Too many samples N="+this.samples.size()+" height:"+sample_height,12,-1);
						}
					updateDrawingArea(img);
					return;
					}
				
				final String title = new SimpleInterval(this.interval).toNiceString()+" length: "+ StringUtils.niceInt(this.interval.getLengthOnReference())+" bp";
				g.setColor(Color.BLACK);
				drawText(g,0,0, title,14,-1);
				double prev_x=0;
				double[] prev_y=null;
				long now = System.currentTimeMillis();
				
				
				if(this.per_sample) {
					final double fontSize=Math.min(12,sample_height/10.0);
					for(int i=0;i< this.samples.size() && !abort_flag;i++) {
						g.setColor(samples.get(i).color);
						drawText(g,1,top_margin + sample_height*i+1.0, samples.get(i).getDisplayName(),fontSize,-1);
						g.setColor(new Color(240,240,240));
						for(double z=0.0;z<=2.0;z+=0.5) {
							final double y = top_margin + (sample_height*(i+1)) - ((z-min_cov)/(max_cov-min_cov) * sample_height);
							drawText(g,1,y,String.valueOf(z),9,-1);
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
						drawText(g,0,y,sn,10,-1);
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
						
						
						final float[] all_values = new float[samples.size()];
						final int[] all_indexes = new int[samples.size()];
						
						for(int i=0;i< samples.size() && !abort_flag;i++) {
							final Sample sn = samples.get(i);
							final float value = Float.parseFloat(tokens[sn.column]);
							all_values[i] = value;
							all_indexes[i] = sn.column;
							
							
							
							final double frac = (value - min_cov)/(max_cov-min_cov);
							final double y;
							
							if(this.per_sample) {
								final double bottom_y = top_margin + (sample_height*(i+1)); 
								final double mid_y = bottom_y - (1.0/(max_cov-min_cov))*sample_height;
								
								y= bottom_y - frac*sample_height;
								final Rectangle2D rec;
								final Paint oldPaint = g.getPaint();
								final LinearGradientPaint grad;
								if(value>=1.0)
									{
									final Color[] p_colors = {Color.LIGHT_GRAY,Color.CYAN,Color.BLUE};
							     	final float[] p_dist = {
								     		(float)((1.0-1.0)/(max_cov-1.0)),
								     		(float)((1.5-1.0)/(max_cov-1.0)),
								     		(float)((2.0-1.0)/(max_cov-1.0))
								     		};
									final double top_y = top_margin + (sample_height*(i)); 
									final Point2D p_start = new Point2D.Double(0,mid_y);
	 								final Point2D p_end = new Point2D.Double(0,top_y);
	 								grad = new LinearGradientPaint(p_start,p_end, p_dist, p_colors);
	 								rec = new Rectangle2D.Double(x0,y,(x1-x0),mid_y-y);
	 								
									}
								else
									{
									final Color[] p_colors = {Color.LIGHT_GRAY,Color.ORANGE, Color.RED};
									final float[] p_dist = { 0.0f,0.5f,1f};
								    
									final Point2D p_start = new Point2D.Double(0,mid_y);
	 								final Point2D p_end = new Point2D.Double(0,bottom_y);
							     	
						     		grad = new LinearGradientPaint(p_start,p_end, p_dist, p_colors);
									rec = new Rectangle2D.Double(x0,mid_y, (x1-x0),y-mid_y);
									}
								g.setPaint(grad);
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
								if(this.showNames && (value < 0.6 || value > 1.4)) {
									BestPoint best = sample2best.get(sn.srcName);
									if(best!=null && Math.abs(best.v) >  Math.abs(value)) continue;
									if(best==null) {
										best = new BestPoint();
										best.sample = sn;
										sample2best.put(sn.srcName,best);
										}
									best.x = (x0+x1)/2.0;
									best.y = y;
									best.v = value;
									}
								}
							}
						
						/* interesting points at top */
						final Optional<Color> interesting = this.helper.getColor(all_indexes, all_values);
						if(interesting.isPresent()) {
							final double mid = (x0+x1)/2.0;
							final double y2 = this.top_margin/2.0;
							g.setColor(interesting.get());
							g.fill(new Ellipse2D.Double(mid-2.5,y2-2.5,5,5));
							}
						
						
						if(System.currentTimeMillis()-now>1000L) {
							updateDrawingArea(img);
							now = System.currentTimeMillis();
							}
						prev_x= x1;
						prev_y = array_y;
						}
					//plot names
					for(String sn:sample2best.keySet()) {
						final BestPoint best = sample2best.get(sn);
						g.setColor(best.sample.color);
						drawText(g, best.x, best.y+2, best.sample.getDisplayName(), 13, 0);
						if(abort_flag) break;
						}
					plotGenes(g,img);
					}
				catch(final IOException err) {
					err.printStackTrace();
					}
				updateDrawingArea(img);
			}
		}
		
		
		
		XFrame(final SAMSequenceDictionary srcDict,final String tabixFileUri,final String gff3path,final Helper helper) throws IOException {
			super(SwingIndexCov.class.getSimpleName());
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			this.tabixFileUri = tabixFileUri;
			this.gff3Path = gff3path;
			this.helper = helper;
			final JMenuBar menuBar = new JMenuBar();
			setJMenuBar(menuBar);
			JMenu menu=new JMenu("File");
			menuBar.add(menu);
			menu.add(new JMenuItem(new AbstractAction("About...") {
				@Override
				public void actionPerformed(ActionEvent arg0) {
					JOptionPane.showMessageDialog(XFrame.this, SwingIndexCov.class.getSimpleName()+". Author : Pierre Lindenbaum PhD.");
					}
				}));
			menu.add(new JMenuItem(new AbstractAction("Save As...") {
				@Override
				public void actionPerformed(final ActionEvent arg0) {
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
				
				this.helper.initialize(tokens);
				
				final Vector<Sample> vSamples = new Vector<>(tokens.length);
				for(int i=3;i< tokens.length;i++) {
					final Sample sn = new Sample(tokens[i],helper.getDisplayName(tokens[i]), i);
					// generate colors : https://stackoverflow.com/questions/223971
					sn.color = helper.getColorForSample(sn.srcName);
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
								JComponent.class.cast(c).setToolTipText(sn.getDisplayName());
								}
							if(!isSelected) {
								c.setForeground(sn.color);
								}
						}
						return c;
					};
				});
				
				// get chromosomes in the tabix
				final Set<String> setChroms =  tbr.getChromosomes();
				final Vector<String> vecChroms = new Vector<>(setChroms);
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
					updateHyperlinks();
					updateGff3Table();
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
			
			this.jcboxShowName = new JCheckBox("Names", false);
			this.jcboxShowName.addActionListener((AE)->samplesSelectionchanged());//will call repaint anyway...
			navPane.add(this.jcboxShowName);
			
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
						updateHyperlinks();
						updateGff3Table();
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
						updateHyperlinks();
						updateGff3Table();
						offscreen=null;
						drawingArea.repaint();
						}
					}
				}));
			
			JPanel pane2 = new JPanel(new BorderLayout(5, 5));
			pane2.setBorder(BorderFactory.createTitledBorder("Samples N="+StringUtils.niceInt(this.jlistSamples.getModel().getSize())));
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
			
			this.drawingArea.addComponentListener(new ComponentAdapter() {
				public void componentResized(java.awt.event.ComponentEvent e) {
					drawingArea.repaint();
				};
			});
			
			viewPane.add(this.drawingArea,BorderLayout.CENTER);
			
			/** ref tab */
			jTabbedPane.addTab("REF", new JScrollPane( new JTable(new SwingSequenceDictionaryTableModel(this.srcDict))));
			/** gff3 tab */
			
			final JTable jtableGff3=  new JTable(this.gff3TableModel=new SwingGff3TableModel());
			jtableGff3.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			final JPanel jpanelgff3 = new JPanel(new BorderLayout());
			jpanelgff3.add(new JScrollPane(jtableGff3),BorderLayout.CENTER);
			jpanelgff3.setBorder(BorderFactory.createTitledBorder("MAX="+StringUtils.niceInt(MAX_GFF3_ROWS)));
			final String titleTabGff3="GFF";
			jTabbedPane.addTab(titleTabGff3,jpanelgff3 );
			

			
			
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
					updateHyperlinks();
					updateGff3Table();
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
			final File prefDir = PreferredDirectory.get(SwingIndexCov.class);
			final JFileChooser chooser= new JFileChooser(prefDir);
			final Optional<SimpleInterval> optR = getUserInterval();
			if(optR.isPresent()) {
				final Locatable loc = optR.get();
				final String fname = loc.getContig()+"_"+loc.getStart()+"_"+loc.getEnd()+".png";
				chooser.setSelectedFile(prefDir==null?new File(fname):new File(prefDir,fname));
				}
			
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
		
		private void updateHyperlinks() {
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
								JOptionPane.showInputDialog(XFrame.this,"URL",U.getUrl());
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
		
		private void updateGff3Table()  {
			final Locatable loc = getUserInterval().orElse(null);

			if(this.gff3Path==null || loc==null) {
				this.gff3TableModel.setRows(Collections.emptyList());
				return;
				}
			final List<Gff3Feature> lines= new ArrayList<>();
			try(TabixReader tbr = new TabixReader(this.gff3Path.toString())) {
				final Gff3Codec codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
				final TabixReader.Iterator iter0 = tbr.query(loc.getContig(), loc.getStart(), loc.getEnd());
				final TabixIteratorLineReader iter1 = new TabixIteratorLineReader(iter0);
				final LineIterator iter = new LineIteratorImpl(iter1);
				while(!codec.isDone(iter) && lines.size()< MAX_GFF3_ROWS) {
					final Gff3Feature feature = codec.decode(iter);
					if(feature==null) continue;
					lines.add(feature);
					}
				codec.close(iter);
				}
			catch(final Throwable err) {
				lines.clear();
				}
			this.gff3TableModel.setRows(lines);
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
			updateHyperlinks();
			updateGff3Table();
			this.drawingArea.repaint();
			}
		
		private void runOnClose() {
			this.helper.dispose();
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
			plotter.tabixUri = this.tabixFileUri;
			plotter.width = this.drawingArea.getWidth();
			plotter.height = this.drawingArea.getHeight();
			plotter.per_sample = this.jcboxPerSample.isSelected();
			plotter.interval = getUserInterval().orElse(null);
			plotter.gff3Path = this.gff3Path;
			plotter.helper = this.helper;
			plotter.showNames = this.jcboxShowName.isSelected();
			final List<Sample> L2 = jlistSamples.getSelectedValuesList();
			plotter.samples = new ArrayList<>();
			if(L2==null || L2.isEmpty()) {
				for(int i=0;i<jlistSamples.getModel().getSize();i++) {
					plotter.samples.add(jlistSamples.getModel().getElementAt(i));
					}
				}
			else
				{
				final Set<String> set = L2.stream().map(S->S.srcName).collect(Collectors.toSet());
				for(int i=0;i<jlistSamples.getModel().getSize();i++) {
					final Sample sn = jlistSamples.getModel().getElementAt(i);
					if(!set.contains(sn.srcName)) continue;
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
				return new IntervalParser(this.srcDict).
					enableWholeContig().
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
			final Helper helper;
			if(StringUtils.isBlank(this.helperPath)) {
				helper = new DefaultHelper();
			} else
				{
				try {
					final String[] tokens = this.helperPath.trim().split("\\s+");
					if(tokens.length<2) {
						LOG.error("Expected two words (file.jar java.class) in "+this.helperPath+" but got "+tokens.length);
						return -1;
						}
					final Path jarPath = Paths.get(tokens[0]);
					IOUtil.assertFileIsReadable(jarPath);
					final URLClassLoader child = new URLClassLoader(
					        new URL[] {jarPath.toUri().toURL()},
					        this.getClass().getClassLoader()
							);
					LOG.info("loading class \""+tokens[1]+"\" from "+jarPath);
					final Class<?> classToLoad = Class.forName(tokens[1], true, child);
					helper = (Helper)classToLoad.getConstructor().newInstance();
					}
				catch(final Throwable err) {
					LOG.error(err);
					return -1;
					}
				}
			
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
						if(f.isDirectory()) return true;
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
				final XFrame frame = new XFrame(dict,input.toString(),this.gff3Path,helper);
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
