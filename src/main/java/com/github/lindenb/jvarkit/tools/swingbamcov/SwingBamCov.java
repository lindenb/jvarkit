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
package com.github.lindenb.jvarkit.tools.swingbamcov;

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.LinearGradientPaint;
import java.awt.Paint;
import java.awt.Panel;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Stroke;
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
import java.io.PrintWriter;
import java.net.URI;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.OptionalInt;
import java.util.Set;
import java.util.Vector;
import java.util.function.IntFunction;
import java.util.function.IntToDoubleFunction;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
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
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.table.AbstractTableModel;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.gff3.SwingGff3TableModel;
import com.github.lindenb.jvarkit.hershey.Hershey;
import com.github.lindenb.jvarkit.illumina.IlluminaConstants;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.net.UrlSupplier;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.swing.PreferredDirectory;
import com.github.lindenb.jvarkit.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.TabixIteratorLineReader;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
/**
BEGIN_DOC

## Example

```
java -jar dist/jvarkit.jar swingbamcov -R ref.fa *.bam
```


```
find dir -type f -name "*.bam" > all.list
java -jar dist/jvarkit.jar swingbamcov -R ref.fa all.list
```

## TruePath

TruePath reads processed with Dragen are splitted in 3 parts: phase 1/ phase undetermined / phase 2.


END_DOC
 */
@Program(name="swingbamcov",
description="Bam coverage viewer using Java Swing UI",
keywords={"bam","alignment","graphics","visualization","swing"},
creationDate="20210420",
modificationDate="20261007",
biostars= {9603151},
jvarkit_amalgamion =  true,
menu="BAM Manipulation"
)
public class SwingBamCov extends Launcher
	{
	private static final int MAX_GFF3_ROWS = 10_000;
	private static final Color[] HAP_COLORS = new Color[] {
				new Color(245,249,252),
				new Color(251,243,231),
				new Color(237,246,241)
				};
	
	private static final Logger LOG = Logger.of(SwingBamCov.class);
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFile = null;
	@Parameter(names={"-r","--regions","--interval"},description="default interval region on opening")
	private String defaultRegion="";
	@Parameter(names={"-q","--mapq"},description="min mapq")
	private int minmapq=1;
	@Parameter(names={"--gtf","--gff"},description="GFF3 file indexed with tabix to plot the genes.")
	private String gffPath = null;
	@Parameter(names={"--small"},description="Display the reads when the region is small than 'x' bp. " + DistanceParser.OPT_DESCRIPTION,splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int smallRegionLength = 200;
	@Parameter(names={"--bed"},description="Load this bed file and use the intervals as a set of menus to jump to a specific location.")
	private Path intervalBed = null;
	@Parameter(names={"--disable-truepath"},description="Disable detection of phasing with illumina truepath.")
	private boolean disable_illumina_truepath_detection = false;

	private static class BamInfo {
		final Path bamPath;
		String sample;
		/** activate splitting for phasing */
		boolean is_illumina_truepath = false;
		BamInfo(final Path bamPath) {
			this.bamPath = bamPath;
			}
		}
	
	private static class HaplotypeData {
		Rectangle2D rect;
		final int[] coverage;
		double[] norm_cov;
		Color color;
		boolean got_read=false;
		HaplotypeData(int genomic_len) {
			this.coverage=new int[genomic_len];
			Arrays.fill(this.coverage, 0);
		}
	}
	
	private static class HaplotypePileup {
		Rectangle2D rect;
		Pileup<SAMRecord> pileup;
		boolean got_read= false;
		Color color;
	}
		
	@SuppressWarnings("serial")
	private static class XFrame extends JFrame {
		final int smallRegionLength;
		final Path referenceFile;
		final SAMSequenceDictionary dict;
		final List<BamInfo> bamPaths;
		final JTable jtableBams;
		final JPanel drawingArea;
		final JTextField jtextFieldLocation;
		final JTextField jtextFieldCap;
		final JCheckBoxMenuItem jcboxDisablePhasing;
		final int minMapq;
		Thread drawingThread = null;
		BufferedImage offScreenImage = null;
		JProgressBar progressBar = null;
		final String gffPath;
		private final JMenu jmenuHyperlinks;
		private final SwingGff3TableModel gff3TableModel;

		/** abstract action to change the genomic view (zooming/shifting) */
		private abstract class ChangeViewAction extends AbstractAction {
			final double factor;
			ChangeViewAction(String title,double factor) {
				super(title);
				this.factor=factor;
				}
			@Override
			public Object getValue(final String key)
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
				jtextFieldLocation.setText(r.getContig()+":"+r.getStart()+"-"+r.getEnd());
				XFrame.this.offScreenImage=null;
				XFrame.this.drawingThread=null;
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

		/** action to zoom */
		private class ZoomAction extends ChangeViewAction {
			ZoomAction(final double factor) {
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
		/** action to shift genomic view */
		private class ShiftAction extends ChangeViewAction {
			ShiftAction(String title,double factor) {
				super(title,factor);
				}
			@Override
			Locatable change(final SAMSequenceRecord ssr, Locatable r)
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
		
		
	private class DrawingThread extends Thread {
			final BufferedImage img;
			final SimpleInterval location;
			final OptionalInt optCap;
			final List<BamInfo> bamInfos;
			final SamReaderFactory srf;
			final boolean disable_phasing;
			DrawingThread(
					BufferedImage img,
					final SimpleInterval location,
					final OptionalInt optCap,
					final List<BamInfo> bamInfos,
					boolean disable_phasing
					) {
				this.img = img;
				this.location=location;
				this.optCap = optCap;
				this.bamInfos = bamInfos;
				this.disable_phasing = disable_phasing;
				this.srf = SamReaderFactory.makeDefault()
						.referenceSequence(XFrame.this.referenceFile)
						.validationStringency(ValidationStringency.LENIENT);
				}
			@Override
			public void run()
				{
				final Graphics2D g = this.img.createGraphics();
				paintDrawingArea(g);
				g.dispose();
				}
			
			private boolean isCurrentThread() {
				return XFrame.this.drawingThread==this &&
					XFrame.this.offScreenImage==this.img;
				}
			private void repaintDrawingArea() {
				if(!isCurrentThread()) return;
				try {
					SwingUtilities.invokeAndWait(()->{
						XFrame.this.repaint();
						});
					} catch(final Throwable err) {LOG.error(err);}
				}
			
			private void paintDrawingArea(final Graphics2D g) {
				g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
				g.setColor(Color.WHITE);
				g.fillRect(0, 0, this.img.getWidth(), this.img.getHeight());
				if(location==null) return;
				final  int marginTop=12;
				g.setColor(Color.DARK_GRAY);
				g.drawString(location.toNiceString()+" length:"+ StringUtils.niceInt(location.getLengthOnReference())+"bp " +
						SequenceDictionaryUtils.getBuildName(XFrame.this.dict).orElse(XFrame.this.referenceFile.getFileName().toString())+" MAPQ:"+XFrame.this.minMapq, 3, marginTop-1);
				
				final int ncols = (int)Math.max(1,Math.floor(Math.sqrt(this.bamInfos.size())));
				final int nrows = (int)Math.max(1, Math.ceil(this.bamInfos.size()/(double)ncols));
				final double wi  = this.img.getWidth()/ncols;
				final double hi  = (this.img.getHeight()-marginTop)/nrows;
				int nr=0;
				int nc=0;
				
				for(final BamInfo bamPath: this.bamInfos) {
					if(!isCurrentThread()) break;
					paintBam(g,this.srf,bamPath,this.location,new Rectangle2D.Double(
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
				if(!isCurrentThread()) return;
				try {
					SwingUtilities.invokeAndWait(()->{
						XFrame.this.progressBar.setIndeterminate(false);
						});
					} catch(final Throwable err) {LOG.error(err);}
				}
		
		private void paintBam(final Graphics2D g,
				final SamReaderFactory srf,
				final BamInfo bam,
				final Locatable loc,
				final Rectangle2D rect
				) {
			if(loc.getLengthOnReference()<= XFrame.this.smallRegionLength ) {
				paintShortBamSection(g,srf,bam,loc,rect,this.disable_phasing);
				}
			else
				{
				paintLargeBamSection(g,srf,bam,loc,rect,this.disable_phasing);
				}
			}
		
		private boolean acceptRead(final SAMRecord rec) {
			 if(
				rec.getReadUnmappedFlag() ||
				rec.getDuplicateReadFlag() ||
		 		rec.getReadFailsVendorQualityCheckFlag() ||
		 		rec.isSecondaryOrSupplementary() ||
		 		rec.getMappingQuality()< XFrame.this.minMapq) return false;
			 return true;
			}
		
		private Color baseToColor(final char c) {
			return  AcidNucleics.baseToColor(c);
			}
		
		
		/** print BAM for small interval, displaying reads */
		private void paintShortBamSection(
				final Graphics2D g,
				final SamReaderFactory srf,
				final BamInfo bam,
				final Locatable region,
				final Rectangle2D rect_frame,
				final boolean disable_phasing
				){
			final boolean is_phased= (!disable_phasing && bam.is_illumina_truepath);
			final List<HaplotypePileup> hapData= new ArrayList<>();
			final Shape oldClip = g.getClip();
			final Stroke oldStroke = g.getStroke();	
			g.setClip(rect_frame); 
		    final int margin_top=12;
			try {
				final Hershey hershey = new Hershey();
				final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*rect_frame.getWidth() + rect_frame.getX();
				for(int i=0;i< (is_phased?3:1);i++) {
					final HaplotypePileup hap = new HaplotypePileup();
					hap.pileup = new Pileup<>((L,R)->position2pixel.applyAsDouble(L.getUnclippedEnd()+1) +1  < position2pixel.applyAsDouble(R.getUnclippedStart()));
					hap.color=(is_phased?HAP_COLORS[i]:Color.WHITE);
					hapData.add(hap);
					}
				
				try(SamReader sr=srf.open(bam.bamPath)) {
					 try(CloseableIterator<SAMRecord> iter=sr.query(
							 region.getContig(),
							 Math.max(0,region.getStart()-XFrame.this.smallRegionLength), //extend to get clipR
							 region.getEnd()+XFrame.this.smallRegionLength,false)) {
						 while(iter.hasNext()) {
							 final SAMRecord rec=iter.next();
							 if(!acceptRead(rec)) continue;
							 if(rec.getUnclippedEnd() < region.getStart()) continue;
							 if(rec.getUnclippedStart() > region.getEnd()) continue;
							 final Cigar cigar = rec.getCigar();
							 if(cigar==null || cigar.isEmpty()) continue;
							 final int phase_idx;
							 if(is_phased) {
								 switch(IlluminaConstants.getTruePathPhase(rec)) {
								 	case Phase_1: phase_idx = 0; break;//top
								 	case Phase_2: phase_idx = 2; break;//bottom
								 	default: phase_idx = 1; break;//middle
								 	}
							 	}
							 else
								{
								phase_idx = 0;
								}
							 final HaplotypePileup hap=hapData.get(phase_idx);
							 hap.pileup.add(rec);
							 hap.got_read=true;
						 	 }
						 }//end iterator
					}//end samreader
				
				// phase 1 and 2 are empty, keep middle
				if(hapData.size()==3 && !hapData.get(0).got_read  && !hapData.get(2).got_read) {
					hapData.get(1).color = Color.WHITE;
					hapData.remove(2);
					hapData.remove(0);
					}
				/** resize rectangles */

				final double total_weight = hapData
						.stream()
						.mapToDouble(H->Math.max(1.0, (double)H.pileup.getRowCount()))
						.max().orElse(1.0);								
				double y = rect_frame.getY() + margin_top;
				for(HaplotypePileup hap: hapData) {
					// weight data
					final double w = Math.max(1.0, (double)hap.pileup.getRowCount());
					final double adjusted_height = (w/total_weight)*(rect_frame.getHeight()-margin_top);
					hap.rect= new Rectangle2D.Double(
							rect_frame.getX(),
							y,
							rect_frame.getWidth(),
							adjusted_height*0.9
							);
					y+=adjusted_height;
					}
				
				ReferenceSequence refInInterval=null;
				 try (ReferenceSequenceFile refseq=ReferenceSequenceFileFactory.getReferenceSequenceFile(XFrame.this.referenceFile)) {
					 final SAMSequenceRecord ssr = XFrame.this.dict.getSequence(region.getContig());
					 if(region.getStart()<=ssr.getSequenceLength()) {
						 refInInterval = refseq.getSubsequenceAt(
								 region.getContig(),
								 region.getStart(),
								 Math.min(region.getEnd(),ssr.getSequenceLength())
								 );
					 	}
				 	}
				 /* clear rect */
						     
			     
			     
			     double featureHeight= 20;
			     for(HaplotypePileup hap: hapData) {
			    	final double h= hap.rect.getHeight()/Math.max(1.0, hap.pileup.getRowCount());
			    	featureHeight = Math.min(featureHeight, h);
			     	}
			     
			     for(int hap_idx=0;hap_idx< hapData.size();++hap_idx) {
			    	 final HaplotypePileup hap = hapData.get(hap_idx);
			    	 g.setColor(hap.color);
					 g.fill(hap.rect);
				     
				     g.setStroke(new BasicStroke(0.5f));
				     g.setColor(Color.WHITE);
				     
				     y= hap.rect.getMaxY()-featureHeight;
				     
				     for(final List<SAMRecord> row: hap.pileup) {
				    	final double h2= Math.min(featureHeight*0.9,featureHeight-2);
		
				    	
				    	for(final SAMRecord rec: row) {
				    		final Cigar cigar=rec.getCigar();
				    		if(cigar==null || cigar.isEmpty()) continue;
				    		
				    		/* draw rec itself */
				    		final double midy=y+h2/2.0;
				    		g.setColor(Color.DARK_GRAY);
				    		g.draw(new Line2D.Double(
				    				position2pixel.applyAsDouble(rec.getUnclippedStart()),
				    				midy,
				    				position2pixel.applyAsDouble(rec.getUnclippedEnd()),
				    				midy));
				    		final int unclipped_start = rec.getUnclippedStart();
				    		int ref1 = unclipped_start;
				    		final List<Double> insertions = new ArrayList<>();
				    		for(final CigarElement ce: cigar.getCigarElements()) {
				    			if(ref1> region.getEnd()) break;
				    			final CigarOperator op=ce.getOperator();
				    			Shape shape = null;
				    			Color fill=null;
				    			Color stroke=Color.DARK_GRAY;
				    			switch(op) {
				    				case P: break;
				    				case M://through
				    				case X://through
				    				case EQ://through
				    				case S: //through
				    				case H: 
				    						final double x1=position2pixel.applyAsDouble(ref1);
				    						
				    						shape = new Rectangle2D.Double(
				    						x1, y,
				    						position2pixel.applyAsDouble(ref1+ce.getLength())-x1,h2
				    						);
				    						
				    						ref1+=ce.getLength();
				    						switch(op) {
				    							case H: case S: fill=Color.YELLOW;break;
				    							case X: fill=Color.RED;break;
				    							case EQ: case M: fill=Color.LIGHT_GRAY;break;
				    							default:break;
				    							}
				    						break;
				    				case N://through
				    				case D: shape=null;fill=null;stroke=null;ref1+=ce.getLength();break;
				    				case I: shape=null;fill=null;stroke=null;insertions.add(position2pixel.applyAsDouble(ref1));break;
				    				default: throw new IllegalStateException(""+op);
				    				}
				    			if(ref1 < region.getStart()) continue;
				    			
				    			if(shape!=null) {
				    				if(fill!=null) {g.setColor(fill);g.fill(shape);}
				    				if(stroke!=null && h2>4)  {g.setColor(stroke);g.draw(shape);}
				    				}
				    			} // end loop cigar
				    		
				    		
				    		
				    		 /* draw mismatched bases */
				   	     	if(refInInterval!=null && rec.getReadBases()!=null && rec.getReadBases()!=SAMRecord.NULL_SEQUENCE) {
				   	     		final byte bases[]=rec.getReadBases();
				   	     		final IntFunction<Character> baseRead= IDX-> IDX<0 || IDX>=bases.length || bases==SAMRecord.NULL_SEQUENCE?'N':(char)Character.toUpperCase(bases[IDX]);
				   	     		int read0=0;
				   	     		ref1 = rec.getAlignmentStart();
					   	     	for(CigarElement ce: cigar.getCigarElements()) {
					    			if(ref1> region.getEnd()) break;
					    			final CigarOperator op=ce.getOperator();
					    			switch(op) {
						    			case P:break;
						    			case H:break;
						    			case D: case N: ref1+=ce.getLength(); break;
						    			case S: case I: read0+=ce.getLength(); break;
						    			case EQ:case M: case X:
						    				{
						    				for(int j=0;j< ce.getLength();j++) {
						    					if(ref1+j< region.getStart()) continue;
						    					if(ref1+j>=region.getStart()+refInInterval.length()) break;
						    					final int ref_base_idx = ref1-region.getStart()+j;
						    					char ctgBase =(char)(ref_base_idx<0 || ref_base_idx>=refInInterval.length()?'N':Character.toUpperCase(refInInterval.getBases()[ref_base_idx]));
						    					if(ctgBase=='N') continue;
						    					char readBase = baseRead.apply(read0+j);
						    					if(readBase=='N') continue;
						    					if(readBase==ctgBase) continue;
						    					g.setColor(Color.ORANGE);
						    					final double x1 = position2pixel.applyAsDouble(ref1+j);
						    					final double x2 = position2pixel.applyAsDouble(ref1+j+1);
						    					g.fill( new Rectangle2D.Double( x1, y,x2-x1,h2));
						    					}
						    				read0+=ce.getLength();
						    				ref1+=ce.getLength();
						    				break;
						    				}
						    			default:break;
					    				}
					   	     		}
				   	     		
				   	     		/** draw bases */
					    		ref1 = unclipped_start;
					    		read0=0;
						    	for(final CigarElement ce: cigar.getCigarElements()) {
					    			if(ref1> region.getEnd()) break;
					    			final CigarOperator op=ce.getOperator();
					    			switch(op) {
					    				case P: break;
					    				case H: ref1+= ce.getLength();break;
					    				case N:  case D: ref1+= ce.getLength(); break;
					    				case I: read0+=ce.getLength();break;
					    				case M: case X: case EQ: case S:
							    			{
							    			for(int j=0;j< ce.getLength();j++) {
							    				final char readBase = baseRead.apply(read0+j);
							    				double h3 = h2*0.9;
						    					final double x1 = position2pixel.applyAsDouble(ref1+j);
						    					final double x2 = position2pixel.applyAsDouble(ref1+j+1);
						    					g.setColor(baseToColor(readBase));
						    					hershey.paint(g, String.valueOf(readBase), x1,y+(h2-h3)/2.0,(x2-x1)*0.95,h3);
								    			}
						    				read0+=ce.getLength();
						    				ref1+=ce.getLength();
							    			break;
							    			}
					    				default: throw new IllegalStateException(""+op);
					    				}
					    			if(ref1 < region.getStart()) continue;
					    			
					    			} // end loop cigar
				   	     		}
				    		
				    		for(double px:insertions) {
				    			g.setColor(Color.RED);
				    			g.draw(new Line2D.Double(px,y-0.5,px,y+h2+0.5));
				    			}
				    		}
				    	y-=featureHeight;
				     	}
				     if(is_phased) {
				    	 g.setColor(Color.MAGENTA);
						 g.draw(hap.rect);
				     	}
				     }
			     
			     g.setColor(Color.DARK_GRAY);
			     g.drawString("Sample:"+ bam.sample +" "+new SimpleInterval(region).toNiceString() , (int)rect_frame.getX()+10,(int)rect_frame.getY()+ 10);
			     g.setColor(Color.LIGHT_GRAY);
			     g.draw(rect_frame);
				}
			catch(final IOException err) {
			     g.drawString("Sample:"+ bam.sample +" Error "+err.getMessage() , 10, 10);
				LOG.error(err);
				}
			g.setStroke(oldStroke);
			g.setClip(oldClip);
			repaintDrawingArea();
			}
		
		private double[] scaleMedian(final int[] coverage,final int width_pixel) {
			final double[] array = new double[width_pixel];
			final DiscreteMedian<Integer> med = new DiscreteMedian<>();
			for(int i=0;i< array.length;i++) {
				int idx0 = (int)(((i+0)/(double)width_pixel)*coverage.length);
				int idx1 = (int)(((i+1)/(double)width_pixel)*coverage.length);
				if(idx0==idx1) idx1++;
				while(idx0<idx1 && idx0 < coverage.length) {
					med.add(coverage[idx0]);
					idx0++;
					}
				array[i]=med.getMedian().orElse(0.0);
				med.clear();
				}
			return array;
			}
		
		private void paintLargeBamSection(final Graphics2D g,
				final SamReaderFactory srf,
				final BamInfo bam,
				final Locatable loc,
				final Rectangle2D rect_frame,
				final boolean disable_phasing
				) {
			final boolean is_phased  = bam.is_illumina_truepath && !disable_phasing;
			final Shape oldClip = g.getClip();
			g.setClip(rect_frame);
			try {
				
				try(SamReader sr=this.srf.open(bam.bamPath)) {
					if(!sr.hasIndex()) return;
					final List<HaplotypeData> hapData = new ArrayList<>();
					final double hap_height= rect_frame.getHeight()/(is_phased?3.0:1.0);
					// top , middle, bottom
					for(int i=0;i< 3;i++) {
						final HaplotypeData hap=new HaplotypeData(loc.getLengthOnReference());
						
						hap.rect = new Rectangle2D.Double(
								rect_frame.getX(),
								rect_frame.getY()+(hap_height*i),
								rect_frame.getWidth(),
								hap_height
								);
						hap.color = (is_phased? HAP_COLORS[i]: Color.WHITE);
						
						hapData.add(hap);
						
						
						
						if(!is_phased) break;
						}
					
					// create two arrays of depth for each phase (if any ) */
					final Set<Integer> haplotype_blocks_positions= new HashSet<>();
					Object prev_block_tag=null;

					
					try(CloseableIterator<SAMRecord> iter=sr.query(loc.getContig(),loc.getStart(),loc.getEnd(),false)) {
						while(iter.hasNext()) {
							final SAMRecord rec= iter.next();

							if(!acceptRead(rec)) continue;
							final Cigar cigar = rec.getCigar();
							if(cigar==null || cigar.isEmpty()) continue;
							final int phase_idx;
							if(is_phased) {
								// detect haplotype block boundaries
								final Object block_tag =  rec.getAttribute(IlluminaConstants.TRUEPATH_PS);
								if(block_tag!=null ) {
									if(prev_block_tag!=null && !prev_block_tag.equals(block_tag)) {
										final int p = rec.getStart();
										if(p>=loc.getStart() && p<=loc.getEnd()) {
											haplotype_blocks_positions.add(p);
											}
										}
									prev_block_tag = block_tag;
									}
								
								
								final Object hp_tag = rec.getAttribute(IlluminaConstants.TRUEPATH_HP);
								if(hp_tag==null || !(hp_tag instanceof Number)) {
									phase_idx = 1;
									}
								else
									{
									final int hp_phase= Number.class.cast(hp_tag).intValue();
									switch(hp_phase) {
										case 1: phase_idx = 0; break;//top
										case 2: phase_idx = 2; break;//bottom
										default: phase_idx= 1 ;break;//middle
										}
									}
								}
							else
								{
								phase_idx = 0;
								}
							// loop over cigar string
							final int[] coverage= hapData.get(phase_idx).coverage;
							boolean got_read=false;
							for(AlignmentBlock ab:rec.getAlignmentBlocks()) {
								for(int x=0;x< ab.getLength();++x) {
									final int ref_pos1 = ab.getReferenceStart()+x;
									final int array_pos0= ref_pos1 - loc.getStart();
									if(array_pos0<0) continue;
									if(array_pos0>=coverage.length) break;
									coverage[array_pos0]++;
									got_read = true;
									}
								}
							if(got_read) {
								hapData.get(phase_idx).got_read=true;
								}
							} // end loop iterator		
						} //end open iterator
					
					// no read phased in 0 or 1
					if(hapData.size()==3 && !hapData.get(0).got_read && !hapData.get(1).got_read) {
						hapData.remove(2);
						hapData.remove(0);
						hapData.get(0).rect=rect_frame;
						hapData.get(0).color=Color.WHITE;
						}
					
					
					/* get median and average depth */
					final DiscreteMedian<Integer> discrete_median = new DiscreteMedian<Integer>();
					hapData.stream()
						.flatMapToInt(H->Arrays.stream(H.coverage))
						.forEach(DP->discrete_median.add(DP));
					final int max_depth_any_phase= hapData.stream()
							.flatMapToInt(H->Arrays.stream(H.coverage))
							.max()
							.orElse(1);

					
					final double median_cov = discrete_median.getMedian().orElse(0);
					final double average_cov = discrete_median.getAverage().orElse(0);
					//System.err.println(median_cov);
					//System.err.println(average_cov);
					//System.err.println(max_depth_any_phase);
					
					for(HaplotypeData hap: hapData) {
						hap.norm_cov =  scaleMedian(hap.coverage,(int)rect_frame.getWidth());
					}
					
					final double maxDepth0 =  hapData.stream()
							.flatMapToDouble(H->Arrays.stream(H.norm_cov))
							.max()
							.orElse(1.0);
					
					final double maxDepth = this.optCap.isPresent()?
							Math.min(optCap.getAsInt(), maxDepth0):
							maxDepth0;
					
					/* loop over each phase */
					for(int phasing_side=0;phasing_side<hapData.size() ; phasing_side++) {
						final HaplotypeData hap = hapData.get(phasing_side);
						/* rect to display that phase */
						final Rectangle2D rect = hap.rect;
						final double[] depths = hap.norm_cov;
						g.setClip(hap.rect);
						g.setColor(hap.color);
						g.fill(hap.rect);
						g.setColor(new Color(240,240,240));
						g.draw(hap.rect);
						
						
						final ToDoubleFunction<Double> toYPixel = DP ->{
							double y = rect.getMaxY() - (DP/maxDepth)*rect.getHeight();
							return y;
							};
						@SuppressWarnings("unused")
						final ToDoubleFunction<Integer> toXPixel = POS ->{
							return rect.getX()+((POS-loc.getStart())/(double)loc.getLengthOnReference())*rect.getWidth();
							};
						
							/*
						for(int hap_block_pos: haplotype_blocks_positions) {
							g.setColor(Color.GREEN);
							final double x = toXPixel.applyAsDouble(hap_block_pos);
							g.fillRect((int)x,(int)rect.getY(),1,(int)rect.getHeight());
							}*/
						
						
						final GeneralPath gp = new GeneralPath();
						gp.moveTo(rect.getX(), rect.getMaxY());
						for(int i=0;i< depths.length;i++) {
							double y = toYPixel.applyAsDouble(depths[i]);
							if(y<rect.getY()) y = rect.getY();//capping
							gp.lineTo(rect.getX()+i, y );
							}
						gp.lineTo(rect.getMaxX(),rect.getMaxY());
						gp.closePath();
						g.setColor(Color.LIGHT_GRAY);
						Color colors[]=new Color[]{Color.LIGHT_GRAY,Color.GRAY};
						final Paint oldPaint = g.getPaint();
						g.setPaint(new LinearGradientPaint(0,(float)rect.getY(),0,(float)rect.getMaxY(),new float[]{0.0f,1.0f},colors));
						g.fill(gp);
						g.setPaint(oldPaint);
						
						
						if(average_cov>0) {
							g.setColor(Color.RED);
							final double y = toYPixel.applyAsDouble(average_cov);
							if(y>=rect.getY()) g.draw(new Line2D.Double(rect.getX(), y, rect.getMaxX(), y));
							}
						
						if(median_cov>0) {
							g.setColor(Color.GREEN);
							final double y = toYPixel.applyAsDouble(median_cov);
							if(y>=rect.getY()) g.draw(new Line2D.Double(rect.getX(), y, rect.getMaxX(), y));
							}
						
						// just draw genes and label for phase 1
						if(phasing_side==0) {
							writeGenes(g,loc,rect);
							}
						// draw label at the botom
						if(phasing_side+1==hapData.size()) {
							g.setColor(Color.BLUE);
							final int fontSize=Math.min(12,(int)(rect_frame.getHeight()/10.0));
							g.setFont(new Font("Courier",Font.PLAIN,fontSize));
							g.drawString(bam.sample+" (max DP: " + max_depth_any_phase+") "+ getGenes(loc).
									filter(G->G.getType().equals("gene")).
									map(G->G.getName()).
									filter(S->!StringUtils.isBlank(S)).
									collect(Collectors.toSet()).
									stream().
									collect(Collectors.joining(",")),
									(int)rect.getX()+5,(int)rect.getMaxY()-1);
							}
						if(!isCurrentThread()) break;
						} // end phasing
					//frame
					g.setClip(rect_frame);
					g.setColor(Color.DARK_GRAY);
					g.draw(rect_frame);
					}
				repaintDrawingArea();
				}
			catch(final Throwable err) {
				LOG.error(err);
				}
			finally {
				g.setClip(oldClip);
				}
			}
		
		
		private Stream<Gff3Feature> getGenes(final Locatable region) {
			if(XFrame.this.gffPath==null) return Stream.empty();
			TabixReader tbr = null;//DO NOT USE AUTOCLOSE RESOURCE
			try {
				tbr = new TabixReader(XFrame.this.gffPath);
				final TabixReader tbfinal = tbr;
				final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
				final String ctg = cvt.apply(region.getContig());
				if(StringUtils.isBlank(ctg)) {
					return Stream.empty();
					}
				
				final Gff3Codec codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
				final TabixReader.Iterator iter0=tbr.query(ctg,region.getStart(),region.getEnd());
				final LineIterator lr = new LineIteratorImpl( new TabixIteratorLineReader(iter0));
				final AbstractIterator<Gff3Feature> iter2= new AbstractIterator<Gff3Feature>() {
					@Override
					protected Gff3Feature advance() {
						try {
							while(!codec.isDone(lr)) {
								final Gff3Feature gffline = codec.decode(lr);
								if(gffline==null) continue;
								return gffline;
								}
							return null;
							} catch (final IOException e) {
							LOG.error(e);
							return null;
							}
						}
					};
				return StreamSupport.stream(new IterableAdapter<Gff3Feature>(iter2).spliterator(),false).
						onClose(()->{ tbfinal.close(); });
				}
			catch(Throwable err) {
				return Stream.empty();
				}
			}

		
		private void writeGenes(final Graphics2D g,final Locatable region,final Rectangle2D rect) {
			if(XFrame.this.gffPath==null) return;
			final IntToDoubleFunction position2pixel = X->rect.getX() + ((X-region.getStart())/(double)region.getLengthOnReference())*rect.getWidth();
			final Composite oldComposite = g.getComposite();
			final Stroke oldStroke = g.getStroke();
			g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.4f));
			g.setColor(Color.ORANGE);
			final double y= rect.getMaxY()-4.0;
			try {
				getGenes(region).
					filter(G->G.getType().equals("exon") || G.getType().equals("transcript")).
					forEach(feature->{
						final double x1 = Math.max(rect.getX(),position2pixel.applyAsDouble(feature.getStart()));
						final double x2 = Math.min(rect.getMaxX(),position2pixel.applyAsDouble(feature.getEnd()));
						if(feature.getType().equals("exon") ) {
							g.draw(new Rectangle2D.Double(x1, y-1, (x2-x1), 3));
							}
						else if(feature.getType().equals("transcript") ) {
							g.draw(new Line2D.Double(x1, y, x2, y));
							}
						});

				}
			catch(final Throwable err) {
				
				}
			finally {
				g.setComposite(oldComposite);
				g.setStroke(oldStroke);
				}
			}
		}
		
		
		XFrame(final Path referenceFile,
				final List<BamInfo> bamPaths,
				String defaultLoc,
				int minMapq,
				final String gffPath,
				final int smallRegionLength,
				final Path intervalBed
				) {
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			setTitle(SwingBamCov.class.getSimpleName());
			this.referenceFile = referenceFile;
			this.gffPath = gffPath;
			this.smallRegionLength = smallRegionLength;
			this.dict = SequenceDictionaryUtils.extractRequired(this.referenceFile);
			this.bamPaths = bamPaths;
			this.minMapq=minMapq;
			
			
			
			
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
			final Action actionGo = new AbstractAction("Plot")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					offScreenImage=null;
					drawingThread=null;
					updateHyperlinks();
					updateGff3Table();
					drawingArea.repaint();
					}
				};
			this.jtextFieldLocation.addActionListener(actionGo);
			JButton button = new JButton(actionGo);
			topPane.add(button);
			topPane.add(new JSeparator());
			label = new JLabel("Cap:", JLabel.RIGHT);
			topPane.add(label);
			this.jtextFieldCap= new JTextField("",3);
			topPane.add(this.jtextFieldCap);
			label.setLabelFor(this.jtextFieldCap);
			this.jtextFieldCap.addActionListener(actionGo);

			
			topPane.add(new JSeparator());
			final JComboBox<ZoomAction> zoomComboBox = new JComboBox<>(
						new ZoomAction[] {
						new ZoomAction(0.1),
						new ZoomAction(0.25),
						new ZoomAction(0.5),
						new ZoomAction(0.75),
						new ZoomAction(2.0),
						new ZoomAction(3.0),
						new ZoomAction(5.0),
						new ZoomAction(10.0)
						}
					);
			zoomComboBox.setEditable(false);
			topPane.add(zoomComboBox);
			final JButton zoomButton = new JButton(new AbstractAction("Zoom")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					final ZoomAction za = (ZoomAction)zoomComboBox.getSelectedItem();
					if(za==null) return;
					za.actionPerformed(e);
					}
				});
			topPane.add(zoomButton);
			
	
			topPane.add(new JSeparator());
			topPane.add(new JButton(new ShiftAction("<<<",-0.9)));
			topPane.add(new JButton(new ShiftAction("<<",-0.5)));
			topPane.add(new JButton(new ShiftAction("<",-0.1)));
			topPane.add(new JButton(new ShiftAction(">",0.1)));
			topPane.add(new JButton(new ShiftAction(">>",0.5)));
			topPane.add(new JButton(new ShiftAction(">>>",0.9)));

			final Panel botPane = new Panel(new FlowLayout(FlowLayout.TRAILING));
			mainPane.add(botPane, BorderLayout.SOUTH);
			this.progressBar=new JProgressBar();
			botPane.add(this.progressBar);

			
			
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
						updateHyperlinks();
						updateGff3Table();
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
			menu.add(new JMenuItem(new AbstractAction("Save IGV Session...")
				{
				@Override
				public void actionPerformed(ActionEvent e)
					{
					doMenuSaveIgvSession();
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
			
			this.jmenuHyperlinks=new JMenu("Hyperlinks");
			menuBar.add(this.jmenuHyperlinks);
			
			
			if(intervalBed!=null) {
				JMenu rMenu = null;
				try(BedLineReader blr = new BedLineReader(intervalBed)) {
					blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict)).
						setValidationStringency(ValidationStringency.LENIENT);
					while(blr.hasNext()) {
						final BedLine rec=blr.next();
						final Interval r = new Interval(rec.getContig(),rec.getStart(),rec.getEnd(),false,rec.getOrDefault(3, "."));
						
						final JMenuItem rItem = new JMenuItem(new AbstractAction(new SimpleInterval(r).toNiceString()+" ("+ StringUtils.niceInt(r.getLengthOnReference())+" bp.)") {
								@Override
								public void actionPerformed(ActionEvent e) {
									jtextFieldLocation.setText(r.getContig()+":"+r.getStart()+"-"+r.getEnd());
									XFrame.this.offScreenImage=null;
									XFrame.this.drawingThread=null;
									updateHyperlinks();
									updateGff3Table();
									drawingArea.repaint();
								}
							});
						if(rMenu==null) {
							rMenu = new JMenu("Intervals");
							menuBar.add(rMenu);
							}
						rMenu.add(rItem);
						}
					}
				}
			
			final JMenu toolMenu=new JMenu("Tools");
			menuBar.add(toolMenu);
			this.jcboxDisablePhasing = new JCheckBoxMenuItem("Disable Phasing", false);
			this.jcboxDisablePhasing.addActionListener(actionGo);
			toolMenu.add(this.jcboxDisablePhasing);
			
			
			/** gff3 tab */
			
			final JTable jtableGff3=  new JTable(this.gff3TableModel=new SwingGff3TableModel());
			jtableGff3.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			final JPanel jpanelgff3 = new JPanel(new BorderLayout());
			jpanelgff3.add(new JScrollPane(jtableGff3),BorderLayout.CENTER);
			jpanelgff3.setBorder(BorderFactory.createTitledBorder("MAX="+StringUtils.niceInt(MAX_GFF3_ROWS)));
			final String titleTabGff3="GFF";
			tabPane.addTab(titleTabGff3,jpanelgff3 );

		}/* end constructor */
		
		
		private void updateHyperlinks() {
			this.jmenuHyperlinks.removeAll();
			final Optional<SimpleInterval> optR = getUserInterval();
			if(optR.isPresent()) {
				if(!Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.BROWSE)) return;
				final UrlSupplier urlSupplier = new UrlSupplier(this.dict);
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
								JOptionPane.showInputDialog(
										XFrame.this,
										U.getDomain(),
										U.getUrl()
										);
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
			if(this.gffPath==null || loc==null) {
				this.gff3TableModel.setRows(Collections.emptyList());
				return;
				}
			final List<Gff3Feature> lines= new ArrayList<>();
			try(TabixReader tbr = new TabixReader(this.gffPath.toString())) {
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
		
		final Optional<SimpleInterval> getUserInterval() {
			final String s = this.jtextFieldLocation.getText().trim().replaceAll(",", "");
			if(StringUtil.isBlank(s)) return Optional.empty();
			try {
				final String tokens[]=s.split("\\s+");
				if(tokens.length==3 && StringUtils.isInteger(tokens[1]) && StringUtils.isInteger(tokens[2])) {
					return  new IntervalParser(this.dict).
							apply(tokens[0]+":"+tokens[1]+"-"+tokens[2]);
					}
				
				return new IntervalParser(this.dict).apply(s);
				}
			catch(Throwable err) {
				err.printStackTrace();
				return Optional.empty();
				}
			}
		
		private void doMenuSaveIgvSession() {
			final Optional<SimpleInterval> optinterval = getUserInterval();
			if(!optinterval.isPresent()) return;
			final File prefDir = PreferredDirectory.get(SwingBamCov.class);
			final SimpleInterval interval = optinterval.get();
			final JFileChooser chooser= new JFileChooser(prefDir);
			final String fname = interval.getContig()+"_"+interval.getStart()+"_"+interval.getEnd()+".igv_session.xml";
			chooser.setSelectedFile(prefDir==null?new File(fname):new File(prefDir,fname));
			if(chooser.showSaveDialog(drawingArea)!=JFileChooser.APPROVE_OPTION) return;
			final File f = chooser.getSelectedFile();
			if(f.exists() && JOptionPane.showConfirmDialog(this, "File "+f.getName()+" exists. Overwite ?", "Overwite ?", JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE, null)!=JOptionPane.OK_OPTION)
				{
				return;
				}
			PreferredDirectory.update(f);
			try {
				
				final XMLOutputFactory xof = XMLOutputFactory.newFactory();
				try(final PrintWriter pw = new PrintWriter(f,"UTF-8")) {
					final XMLStreamWriter w=xof.createXMLStreamWriter(pw);
					w.writeStartDocument("UTF-8","1.0");
					w.writeStartElement("Session");
					w.writeAttribute("genome", this.referenceFile.toAbsolutePath().toString());
					w.writeAttribute("hasSequenceTrack","true");
					w.writeAttribute("hasGeneTrack","false");
					w.writeAttribute("version","8");
					w.writeAttribute("locus",interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd());
					w.writeStartElement("Resources");
					for(final BamInfo bamInfo : this.bamPaths) {
						w.writeEmptyElement("Resource");
						w.writeAttribute("name",bamInfo.sample);
						w.writeAttribute("path",bamInfo.bamPath.toAbsolutePath().toString());
						}
					w.writeEndElement();//Resources
					w.writeEndElement();//Session
					w.writeEndDocument();
					w.flush();
					w.close();
					pw.flush();
					}
				}
			catch(final Throwable err) {
				LOG.error(err);
				return;
				}
			}
		private void doMenuSaveAs() {
			final File prefDir = PreferredDirectory.get(SwingBamCov.class);
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
			PreferredDirectory.update(SwingBamCov.class,f);
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
				LOG.error(err);
				return;
				}
			}
		
		private void paintDrawingArea(final Graphics2D g) {
			try {
				if(offScreenImage==null || offScreenImage.getWidth()!=this.drawingArea.getWidth() ||
					offScreenImage.getHeight()!=this.drawingArea.getHeight())
					{
					final Optional<SimpleInterval> location = getUserInterval();
					String capStr = this.jtextFieldCap.getText();
					OptionalInt optCap = OptionalInt.empty();
					if(StringUtils.isInteger(capStr)) {
						optCap = OptionalInt.of(Integer.parseInt(capStr));
						}
					
					
					
					offScreenImage = new BufferedImage(this.drawingArea.getWidth(), this.drawingArea.getHeight(), BufferedImage.TYPE_INT_RGB);
					
					final List<BamInfo> paths;
					if(this.jtableBams.getSelectedRowCount()==0) {
						paths = this.bamPaths;
						} 
					else {
						final int[] selidx = this.jtableBams.getSelectedRows();
						paths = new Vector<>(selidx.length);
						for(int idx:selidx) {
							paths.add(bamPaths.get(idx));
							}
						}
					drawingThread = new DrawingThread(this.offScreenImage,
							location.orElse(null),
							optCap,
							paths,
							this.jcboxDisablePhasing.isSelected()
							);
					this.progressBar.setIndeterminate(true);
					drawingThread.start();
					}
				else
					{
					g.drawImage(this.offScreenImage,0,0,null);
					}
				}
			catch(final Throwable err) {
				ThrowablePane.show(this.drawingArea, err);
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		try {
			final SamReaderFactory srf = super.createSamReaderFactory()
					.validationStringency(ValidationStringency.LENIENT)
					.referenceSequence(this.referenceFile);
			final List<BamInfo> bams = IOUtils.unrollPaths(args)
				.stream()
				.map(P->new BamInfo(P))
				.collect(Collectors.toCollection(Vector::new));
			
			for(final BamInfo bi: bams) {
				IOUtil.assertFileIsReadable(bi.bamPath);
				try(SamReader sr= srf.open(bi.bamPath)) {
					if(!sr.hasIndex()) {
						LOG.error("file is not indexed: "+bi.bamPath);
						return -1;
						}
					final SAMFileHeader header = sr.getFileHeader();
					bi.sample = header.getReadGroups().stream().
							map(S->S.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(bi.bamPath));
					
					
					
					// detect illumina truepath
					if(!this.disable_illumina_truepath_detection &&  IlluminaConstants.isDragenProximityEnabled(header)) {
						bi.is_illumina_truepath = true;
						LOG.info("TruePath detected for "+bi.bamPath);	
						}
					}
				catch(final IOException err) {
					LOG.error(err);
					return -1;
					}
				}
			Collections.sort(bams,(A,B)->A.sample.compareTo(B.sample));
				
			
			if(bams.isEmpty()) {
				LOG.error("no BAM/CRAM was provided");
				return -1;
				}
			JFrame.setDefaultLookAndFeelDecorated(true);
			final XFrame frame = new XFrame(
					this.referenceFile,
					bams,
					this.defaultRegion,
					this.minmapq,
					this.gffPath,
					this.smallRegionLength,
					this.intervalBed
					);
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
		new SwingBamCov().instanceMain(args);//no exit
		}

	}
