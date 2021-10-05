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

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Composite;
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
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.OptionalInt;
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
import javax.swing.JButton;
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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.CoverageFactory;
import com.github.lindenb.jvarkit.samtools.reference.SwingSequenceDictionaryTableModel;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

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
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
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
java -jar dist/swingbamcov.jar -R ref.fa *.bam
```

## Screenshot

 * https://twitter.com/yokofakun/status/1392173415684100105
 * https://twitter.com/yokofakun/status/1443187891480502279
 * https://twitter.com/yokofakun/status/1443208754229563392

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
	@Parameter(names={"--gtf","--gff"},description="GFF3 file indexed with tabix to plot the genes.")
	private String gffPath = null;
	@Parameter(names={"--small"},description="Display the reads when the region is small than 'x' bp. " + DistanceParser.OPT_DESCRIPTION,splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int smallRegionLength = 200;
	
	private static class BamInfo {
		final Path bamPath;
		String sample;
		BamInfo(final Path bamPath) {
			this.bamPath = bamPath;
			}
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
		final int minMapq;
		Thread drawingThread = null;
		BufferedImage offScreenImage = null;
		JProgressBar progressBar = null;
		final String gffPath;
		
		
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
				final SAMSequenceRecord ssr=dict.getSequence(r.getContig());
				if(ssr==null) return;
				r= change(ssr,r);
				if(r==null) return;
				jtextFieldLocation.setText(r.getContig()+":"+r.getStart()+"-"+r.getEnd());
				XFrame.this.offScreenImage=null;
				XFrame.this.drawingThread=null;
				drawingArea.repaint();
				}
			abstract Locatable change(SAMSequenceRecord ssr,Locatable loc);
			abstract String getShortDesc();
			}
		
		
		private class ZoomAction extends ChangeViewAction {
			ZoomAction(double factor) {
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
		
		
	private class DrawingThread extends Thread {
			final BufferedImage img;
			final SimpleInterval location;
			final OptionalInt optCap;
			final List<BamInfo> bamInfos;
			final SamReaderFactory srf;
			DrawingThread(
					BufferedImage img,
					final SimpleInterval location,
					final OptionalInt optCap,
					final List<BamInfo> bamInfos
					) {
				this.img = img;
				this.location=location;
				this.optCap = optCap;
				this.bamInfos = bamInfos;
				this.srf = SamReaderFactory.makeDefault().
						referenceSequence(XFrame.this.referenceFile).
						validationStringency(ValidationStringency.LENIENT);
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
				paintShortBamSection(g,srf,bam,loc,rect);
				}
			else
				{
				paintLargeBamSection(g,srf,bam,loc,rect);
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
			switch(c)
				{
				case 'A':case 'a': return Color.BLUE;
				case 'T': case 't' : return Color.GREEN;
				case 'C': case 'c': return Color.ORANGE;
				case 'G': case 'g': return Color.RED;
				default: break;
				}
			return Color.BLACK;
			}
		
		
		/** print BAM for small interval, displaying reads */
		private void paintShortBamSection(
				final Graphics2D g,
				final SamReaderFactory srf,
				final BamInfo bam,
				final Locatable region,
				final Rectangle2D rect
				){
			final Shape oldClip = g.getClip();
			final Stroke oldStroke = g.getStroke();	
			g.setClip(rect);
			try {
				final Hershey hershey = new Hershey();
				final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*rect.getWidth() + rect.getX();
				final Pileup<SAMRecord> pileup = new Pileup<>((L,R)->position2pixel.applyAsDouble(L.getUnclippedEnd()+1) +1  < position2pixel.applyAsDouble(R.getUnclippedStart()));
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
							 pileup.add(rec);
						 	 }
						 }//end iterator
					}//end samreader
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
				 g.setColor(new Color(240,240,240));
				 g.fill(rect);
			     
			     g.setStroke(new BasicStroke(0.5f));
			     g.setColor(Color.WHITE);		     
			     
			     final int margin_top=12;
			     
			     final double featureHeight= Math.min(20,(rect.getHeight()-margin_top)/Math.max(1.0,(double)pileup.getRowCount()));
			     
			     double y= rect.getMaxY()-featureHeight;
			     
			     for(final List<SAMRecord> row:pileup) {
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
			    
			     
			     g.setColor(Color.DARK_GRAY);
			     g.drawString("Sample:"+ bam.sample +" "+new SimpleInterval(region).toNiceString() , (int)rect.getX()+10,(int)rect.getY()+ 10);
			     g.setColor(Color.LIGHT_GRAY);
			     g.draw(rect);
				}
			catch(final IOException err) {
			     g.drawString("Sample:"+ bam.sample +" Error "+err.getMessage() , 10, 10);
				LOG.error(err);
				}
			g.setStroke(oldStroke);
			g.setClip(oldClip);
			repaintDrawingArea();
			}
		
		private void paintLargeBamSection(final Graphics2D g,
				final SamReaderFactory srf,
				final BamInfo bam,
				final Locatable loc,
				final Rectangle2D rect
				) {
			
			final Shape oldClip = g.getClip();
			g.setClip(rect);
			try {
				try(SamReader sr=this.srf.open(bam.bamPath)) {
					if(!sr.hasIndex()) return;			
					final CoverageFactory covFactory = new CoverageFactory().
							setMappingQuality(XFrame.this.minMapq);
					final CoverageFactory.SimpleCoverage cov = covFactory.getSimpleCoverage(sr, loc, bam.sample);
					
					final double[] depths = cov.scaleMedian((int)rect.getWidth());
					final double maxDepth0 = Arrays.stream(depths).max().orElse(1.0);
					final double maxDepth = this.optCap.isPresent()?
							Math.min(optCap.getAsInt(), maxDepth0):
							maxDepth0;
					final ToDoubleFunction<Double> toYPixel = DP ->{
						double y = rect.getMaxY() - (DP/maxDepth)*rect.getHeight();
						return y;
					};
					
					
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
					
					
					OptionalDouble mean = cov.getMedian();
					if(mean.isPresent()) {
						g.setColor(Color.RED);
						final double y = toYPixel.applyAsDouble(mean.getAsDouble());
						if(y>=rect.getY()) g.draw(new Line2D.Double(rect.getX(), y, rect.getMaxX(), y));
						}
					
					mean = cov.getAverage();
					if(mean.isPresent()) {
						g.setColor(Color.GREEN);
						final double y = toYPixel.applyAsDouble(mean.getAsDouble());
						if(y>=rect.getY()) g.draw(new Line2D.Double(rect.getX(), y, rect.getMaxX(), y));
						}
					
					writeGenes(g,loc,rect);
					
					g.setColor(Color.BLUE);
					final int fontSize=Math.min(12,(int)(rect.getHeight()/10.0));
					g.setFont(new Font("Courier",Font.PLAIN,fontSize));
					g.drawString(bam.sample+" (max DP: " + (int)maxDepth0+") "+ getGenes(loc).
							filter(G->G.getType().equals("gene")).
							map(G->G.getName()).
							filter(S->!StringUtils.isBlank(S)).
							collect(Collectors.toSet()).
							stream().
							collect(Collectors.joining(",")),
							(int)rect.getX()+5,(int)rect.getMaxY()-1);
					//frame
					g.setColor(Color.DARK_GRAY);
					g.draw(rect);
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
				final List<Path> bamPaths,
				String defaultLoc,
				int minMapq,
				final String gffPath,
				final int smallRegionLength
				) {
			super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			setTitle(SwingBamCov.class.getSimpleName());
			this.referenceFile = referenceFile;
			this.gffPath = gffPath;
			this.smallRegionLength = smallRegionLength;
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					referenceSequence(this.referenceFile).
					validationStringency(ValidationStringency.LENIENT);
			this.dict = SequenceDictionaryUtils.extractRequired(this.referenceFile);
			this.bamPaths = bamPaths.stream().map(P->new BamInfo(P)).
					collect(Collectors.toCollection(Vector::new));
					new Vector<BamInfo>(bamPaths.size());
			for(final BamInfo bi:this.bamPaths) {
				try(SamReader sr= srf.open(bi.bamPath)) {
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
			topPane.add(new JButton(new ZoomAction(0.5)));
			topPane.add(new JButton(new ZoomAction(1.5)));
			topPane.add(new JButton(new ZoomAction(2.0)));
			topPane.add(new JButton(new ZoomAction(10)));
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
					apply(s.trim());
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
						paths);
				this.progressBar.setIndeterminate(true);
				drawingThread.start();
				}
			else
				{
				g.drawImage(this.offScreenImage,0,0,null);
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
			final XFrame frame = new XFrame(this.referenceFile,paths,defaultRegion,this.minmapq,this.gffPath,this.smallRegionLength);
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
