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
package com.github.lindenb.jvarkit.tools.hic;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.hic.HicReader;
import com.github.lindenb.jvarkit.hic.HicReaderFactory;
import com.github.lindenb.jvarkit.hic.Normalization;
import com.github.lindenb.jvarkit.hic.Unit;
import com.github.lindenb.jvarkit.io.CustomSeekableStreamFactory;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;

/**
BEGIN_DOC

### Example


```

END_DOC
 */

@Program(name="hicstraw",
	description="Query a Hi-C file",
	keywords={"hic"},
	creationDate="20190613",
	modificationDate="20190614",
	generate_doc=false
	)
public class HicStraw  extends Launcher {
	private static final Logger LOG = Logger.build(HicStraw.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT+" If filename ends with '.png' , '.jpg' or '.svg' the output will be an image")
	private Path outputFile = null;
	@Parameter(names={"-i","--interval1"},description="Interval 1",required=true)
	private String interval1Str = null;
	@Parameter(names={"-j","--interval2"},description="Interval 2. Use '*' to map all the chromosomes.",required=true)
	private String interval2Str = null;
	@Parameter(names={"-u","--unit"},description="Unit")
	private Unit unit = Unit.BP;
	@Parameter(names={"-n","--normalization"},description="normalization")
	private Normalization norm = Normalization.VC;
	@Parameter(names={"-b","--bin"},description="bin size")
	private int binSize = 2_500_000;
	@Parameter(names={"-min-distance"},description="min distance between two intervals on the same chromosome. Don' print the value if they're closer than this value")
	private Integer minCisDistance = null;
	@Parameter(names={"-min-value"},description="Don't print the value if it's lower than 'v'")
	private Float minValue = null;
	@Parameter(names={"-max-value"},description="Don't print the value if it's greater than 'v'")
	private Float maxValue = null;

	private abstract class AbstractCallBack implements HicReader.QueryCallBack {
		PrintStream pw = null;
		boolean first = true;
		String source;
		void finish() {
			if(this.pw!=null) {
			this.pw.flush();
			this.pw.close();
			}
		}
	}
	
	private static class XYV {
		final int x;
		final int y;
		final float v;
		XYV(int x,int y,float v) {
		this.x=x;
		this.y=y;
		this.v=v;
		}
	}
	
	private abstract class ImageCallBack extends AbstractCallBack {
		protected int binsize = 1;
		protected String contig1;
		protected String contig2;
		protected  final List<XYV> contacts = new ArrayList<>(100_000);
		@Override
		public void reportContact(
				String contig1,int start1,int end1,
				String contig2,int start2,int end2,
				final Normalization norm,
				final Unit unit,
				final int binsize, 
				final float value
				)
			{
			if(this.first) {
				first=false;
				this.contig1 = contig1;
				this.contig2 = contig2;
				this.binsize = binsize;
				}
			this.contacts.add(new XYV(start1,start2,value));
			}
		
	}
	
	private class RasterCallBack extends ImageCallBack {
	
		@Override
		void finish() {
			
			final int margin = 20;
			final int imageSize=1000;
			final int drawingSize = imageSize + margin;
			
			final int minX = this.contacts.stream().mapToInt(P->P.x).min().orElse(0);
			final int maxX = this.contacts.stream().mapToInt(P->P.x).max().orElse(0);
			final int distanceX = maxX-minX;
			
			
			final int minY = this.contacts.stream().mapToInt(P->P.y).min().orElse(0);
			final int maxY = this.contacts.stream().mapToInt(P->P.y).max().orElse(0);
			final int distanceY = maxY-minY;
			
			final int distance = Math.max(distanceX, distanceY);

			final float maxV =(float)this.contacts.stream().mapToDouble(P->P.v).max().orElse(10.0);
			final double logMaxV =Math.log(maxV); 
			
	
			final BufferedImage img = new BufferedImage(drawingSize, drawingSize, BufferedImage.TYPE_INT_RGB);
			final Graphics2D g = img.createGraphics();
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, drawingSize, drawingSize);
			g.setColor(Color.GRAY);
			g.drawPolygon(
					new int[] {margin,margin+imageSize,margin},
					new int[] {0,imageSize,imageSize},
					3);
			
			for(final XYV contact:this.contacts) {
				if(contact.v< 1 ) continue;
				final int gray = 255-(int)(255*((Math.log(contact.v))/logMaxV));
				double x  =  margin + ((contact.x-minX)/(double)(distance))*imageSize;
				double y  = ((contact.y-minY)/(double)(distance))*imageSize;
				double w  = (this.binsize/(double)distance)*imageSize;
				double h  = w;
				final Rectangle2D rect = new Rectangle2D.Double(x, y, w, h);
				g.setColor(new Color(gray,gray,gray));
				g.fill(rect);
				if(w > 5) {
					g.setColor(Color.DARK_GRAY);
					g.draw(rect);
					}
				}
			g.setColor(Color.DARK_GRAY);
			g.drawString(this.contig2+":"+StringUtils.niceInt(minY)+"-"+StringUtils.niceInt(maxY),
					margin,
					drawingSize -2
					);
			final AffineTransform oldtr = g.getTransform();
			g.setTransform(AffineTransform.getRotateInstance(Math.PI/2.001));
			g.drawString(this.contig1+":"+StringUtils.niceInt(minX)+"-"+StringUtils.niceInt(maxX),
					2,
					-5
					);
			g.setTransform(oldtr);
			g.dispose();
			try {
			ImageIO.write(img,outputFile==null || outputFile.getFileName().toString().endsWith(".png")?"PNG":"JPG", this.pw);
			} catch(final IOException err)
			{
				throw new RuntimeIOException(err);
			}
		
			super.finish();
			}
	}
	
	private class SVGCallBack extends ImageCallBack {
	
		private String format(double v) {
			return String.valueOf(v);
		}
		@Override
		void finish() {
			final int margin = 50;
			final int imageSize=1000;
			final int drawingSize = imageSize + margin;
			
			final int minX = this.contacts.stream().mapToInt(P->P.x).min().orElse(0);
			final int maxX = this.contacts.stream().mapToInt(P->P.x).max().orElse(0);
			final int distanceX = maxX-minX;
			
			
			final int minY = this.contacts.stream().mapToInt(P->P.y).min().orElse(0);
			final int maxY = this.contacts.stream().mapToInt(P->P.y).max().orElse(0);
			final int distanceY = maxY-minY;
			
			final int distance = Math.max(distanceX, distanceY);

			final float maxV =(float)this.contacts.stream().mapToDouble(P->P.v).max().orElse(10.0);
			final double logMaxV =Math.log(maxV); 
	

			final double dw  = (this.binsize/(double)distance)*imageSize;
			final double dh  = dw;
			try 
			{
				final XMLOutputFactory xof = XMLOutputFactory.newFactory();
				final XMLStreamWriter w  = xof.createXMLStreamWriter(this.pw);
				w.writeStartDocument("1.0");
				w.writeStartElement("svg");
				w.writeDefaultNamespace(SVG.NS);
				w.writeAttribute("width", String.valueOf(drawingSize+1));
				w.writeAttribute("height",  String.valueOf(drawingSize+1));
				w.writeAttribute("style", "");
				
				w.writeStartElement("title");
				w.writeCharacters(this.source);
				w.writeEndElement();//defs
				
				w.writeStartElement("defs");
				
				w.writeEmptyElement("rect");
				w.writeAttribute("id","Q");
				w.writeAttribute("width", format(dw));
				w.writeAttribute("height", format(dh));
				w.writeAttribute("style","stroke:lightgray");
				
				w.writeEndElement();//defs

				
				w.writeStartElement("g");
				w.writeAttribute("style", "stroke:lightgray;");

				for(final XYV contact:this.contacts) {
					if(contact.v< 1 ) continue;
					final int gray = 255-(int)(255*((Math.log(contact.v))/logMaxV));
					double x  =  margin + ((contact.x-minX)/(double)(distance))*imageSize;
					double y  = ((contact.y-minY)/(double)(distance))*imageSize;
					String grayC="rgb("+gray+","+gray+","+gray+")";
					w.writeEmptyElement("use");
					w.writeAttribute("x", format(x));
					w.writeAttribute("y", format(y));
					w.writeAttribute("href","#Q");
					w.writeAttribute("style","fill:"+grayC);
				}
				
				w.writeStartElement("g");
				w.writeAttribute("style", "stroke:darkgray;");
				
				w.writeStartElement("text");
				w.writeAttribute("x", format(margin));
				w.writeAttribute("y", format(drawingSize));
				w.writeCharacters(this.contig2+":"+StringUtils.niceInt(minY)+"-"+StringUtils.niceInt(maxY));
				w.writeEndElement();//text
				
				w.writeStartElement("text");
				w.writeAttribute("x", format(2));
				w.writeAttribute("y", format(-5));
				w.writeAttribute("transform","rotate(90)");
				w.writeCharacters(this.contig1+":"+StringUtils.niceInt(minX)+"-"+StringUtils.niceInt(maxX));
				w.writeEndElement();//text

				
				w.writeEndElement();//g
				
				w.writeEndElement();//g
				w.writeEndElement();//svg
				w.writeEndDocument();
				w.flush();
				w.close();
			}
			catch(final XMLStreamException err)
			{
				throw new RuntimeIOException(err);
			}
			
		
			super.finish();
			}
	}
	
	private class DefaultCallBack extends AbstractCallBack {
		@Override
		public void reportContact(
				String contig1,int start1,int end1,
				String contig2,int start2,int end2,
				final Normalization norm,
				final Unit unit,
				final int binsize, 
				final float value
				)
			{
			if(this.first) {
				pw.println("##source="+this.source);
				pw.println("##unit="+unit);
				pw.println("##normalisation="+norm);
				pw.println("##bin-size="+binsize);
				pw.println("#CHROM1\tSTART1\tEND1\tCHROM2\tSTART2\tEND2\tVALUE");
				this.first = false;
				}
			if(minValue!=null && value < minValue.floatValue()) return;
			if(maxValue!=null && value > maxValue.floatValue()) return;
			
			if(minCisDistance!=null && contig1.equals(contig2)) {
				final int distance;
				if(CoordMath.overlaps(start1, end1, start2, end2)) {
					distance = 0;
					}
				else if(end1 < start2) {
					distance = start2 - end1;
					}
				else
					{
					distance = start1 - end2;
					}
				if(distance < minCisDistance) return;
				}
			pw.print(contig1);
			pw.print("\t");
			pw.print(start1);
			pw.print("\t");
			pw.print(end1);
			pw.print("\t");
			pw.print(contig2);
			pw.print("\t");
			pw.print(start2);
			pw.print("\t");
			pw.print(end2);
			pw.print(contig1);
			pw.print("\t");
			pw.print(start1);
			pw.print("\t");
			pw.print(value);
			pw.println();
			}
		};
	
		
		
	@Override
	public int doWork(final List<String> args) {
		try
			{
			final ISeekableStreamFactory seekableStreamFactory = new CustomSeekableStreamFactory().
					setUsingHttpHead(false).
					setNormalizeURI(false);
			final AbstractCallBack callback;
			
			
			if(this.outputFile==null) {
				callback = new DefaultCallBack();
			} else if(this.outputFile.getFileName().toString().endsWith(".svg") ||
					this.outputFile.getFileName().toString().endsWith(".svg.gz")) {
				callback = new SVGCallBack();
			}
			else if(this.outputFile.getFileName().toString().endsWith(".png") ||
					this.outputFile.getFileName().toString().endsWith(".jpeg") ||
					this.outputFile.getFileName().toString().endsWith(".jpg")) {
				callback = new RasterCallBack();
				}
			else  {
				callback = new DefaultCallBack();
				}

			callback.pw= super.openPathOrStdoutAsPrintStream(outputFile);
			
			for(final String input :args) {
				callback.source = input;
				
				try(final HicReader hicReader = new HicReaderFactory().
							setSeekableStreamFactory(seekableStreamFactory).
							open(input)) { 
				
					final Function<String,Locatable > parseInterval = (S)->{
						final Optional<Locatable> loc = hicReader.parseInterval(S);
						if(!loc.isPresent()) {
							LOG.error("bad interval : \""+S+"\" available are "+ hicReader.
									getDictionary().getSequences().stream().
									map(SR->SR.getSequenceName()).collect(Collectors.joining(" ; ")));
							return null;
							}
						return loc.get();
						};
					
					if(!hicReader.getBasePairResolutions().contains(this.binSize)) {
						LOG.error("bad binSize : \""+this.binSize+"\" available are "+ hicReader.getBasePairResolutions().stream().map(S->String.valueOf(S)).collect(Collectors.joining(" ; ")));
						return -1;
						}
						
					final Locatable loc1 = parseInterval.apply(this.interval1Str);
					if(loc1==null) return -1;
					
					final List<Locatable> loc2list;
					if(!("*".equals(this.interval2Str))) {
						final Locatable loc2 = parseInterval.apply(this.interval2Str);
						if(loc2==null) return -1 ;
						loc2list = java.util.Collections.singletonList(loc2);
						}
					else
						{
						loc2list = hicReader.getDictionary().
								getSequences().stream().
								map(SR->new Interval(SR)).
								collect(Collectors.toList());
						}
					
					for(final Locatable loc2:loc2list) {
						callback.first = true;
						hicReader.query(loc1, loc2,norm, this.binSize, this.unit,callback);
						}
					}
				}
			callback.finish();
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args) {
		new HicStraw().instanceMainWithExit(args);
		}
}
