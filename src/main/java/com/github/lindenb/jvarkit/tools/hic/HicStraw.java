/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.hic.HicReader;
import com.github.lindenb.jvarkit.hic.HicReaderFactory;
import com.github.lindenb.jvarkit.hic.Normalization;
import com.github.lindenb.jvarkit.hic.Unit;
import com.github.lindenb.jvarkit.io.CustomSeekableStreamFactory;
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

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
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
	@Parameter(names={"-min-value"},description="Don' print the value if it's lower than 'v'")
	private Float minValue = null;
	@Parameter(names={"-max-value"},description="Don' print the value if it's greater than 'v'")
	private Float maxValue = null;

	private abstract class AbstractCallBack implements HicReader.QueryCallBack {
		PrintWriter pw = null;
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
	
	private class SVGCallBack extends AbstractCallBack {
		private int binsize = 1;
		private  final List<XYV> contacts = new ArrayList<>(100_000);
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
				this.binsize = binsize;
				}
			this.contacts.add(new XYV(start1,start2,value));
			}
		
		private String format(double v) {
			return String.valueOf(v);
		}
		@Override
		void finish() {
			final int minX = this.contacts.stream().mapToInt(P->P.x).min().orElse(0);
			final int maxX = this.contacts.stream().mapToInt(P->P.x).max().orElse(0);
			final int distanceX = maxX-minX;
			final int count_items_x = distanceX/this.binsize;
			
			final int minY = this.contacts.stream().mapToInt(P->P.y).min().orElse(0);
			final int maxY = this.contacts.stream().mapToInt(P->P.y).max().orElse(0);
			final int distanceY = maxY-minY;
			final int count_items_y = distanceY/this.binsize;
			
			final int pageSize=800;
			final int distance = Math.max(distanceX, distanceY);
			final int dPix=pageSize/Math.max(count_items_x, count_items_y);


			final float maxV =(float)this.contacts.stream().mapToDouble(P->P.v).max().orElse(1.0);
			
			final Function<Integer, Double> convertXToPixel = (X)->{
				return ((X-minX)/(double)(maxX-minX))*pageSize;
			};
			
			final Function<Integer, Double> convertYToPixel = (Y)->{
				return ((Y-minY)/(double)(maxY-minY))*pageSize;
			};
	

			
			try 
			{
				final XMLOutputFactory xof = XMLOutputFactory.newFactory();
				final XMLStreamWriter w  = xof.createXMLStreamWriter(this.pw);
				w.writeStartDocument("1.0");
				w.writeStartElement("svg");
				w.writeDefaultNamespace(SVG.NS);
				w.writeAttribute("width", String.valueOf(pageSize));
				w.writeAttribute("height",  String.valueOf(pageSize));
				w.writeAttribute("style", "");
				
				w.writeStartElement("title");
				w.writeCharacters(this.source);
				w.writeEndElement();//defs
				
				w.writeStartElement("defs");
				
				w.writeEmptyElement("rect");
				w.writeAttribute("id","Q");
				w.writeAttribute("width", format(dPix));
				w.writeAttribute("height", format(dPix));
				w.writeAttribute("style","stroke:lightgray");
				
				w.writeEndElement();//defs

				
				w.writeStartElement("g");
				w.writeAttribute("style", "stroke:lightgray;");

				for(final XYV contact:this.contacts) {
					final int g = (int)(255*(contact.v/(double)maxV));
					String gray="rgb("+g+","+g+","+g+")";
					w.writeEmptyElement("use");
					w.writeAttribute("x", format(convertXToPixel.apply(contact.x)));
					w.writeAttribute("y", format(convertYToPixel.apply(contact.y)));
					w.writeAttribute("href","#Q");
					w.writeAttribute("style","fill:"+gray);
				}
				
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
	
	private class MyQueryCallBack extends AbstractCallBack {
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
			final ISeekableStreamFactory seekableStreamFactory = new CustomSeekableStreamFactory();
			final AbstractCallBack callback = new SVGCallBack();
			
			callback.pw= super.openPathOrStdoutAsPrintWriter(outputFile);
			
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
								map(SR->new Interval(SR.getSequenceName(),1,SR.getSequenceLength())).
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
