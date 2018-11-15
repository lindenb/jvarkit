/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bam2svg;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Result;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.DocumentFragment;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.tools.misc.ConcatSam;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Example

```bash
```

## Gallery

https://twitter.com/yokofakun/status/523031098541232128

![bam2svg](https://pbs.twimg.com/media/B0IuAw2IgAAYfNM.jpg)

https://twitter.com/yokofakun/status/522415314425090048

![bam2svg-2](https://pbs.twimg.com/media/Bz_99ayIMAAK57s.jpg)



END_DOC
 */
@Program(name="bam2svg",
description="BAM to raster graphics",
keywords={"bam","alignment","graphics","visualization","svg"},
generate_doc=false
)
public class SvToSVG extends Launcher
	{
	private static final int HEIGHT_MAIN_TITLE=100;
	private static final int HEIGHT_SAMPLE_NAME=50;
	
	
	private static final Logger LOG = Logger.build(SvToSVG.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-r","-i","--interval","--region"},description="interval CHROM:START-END",required=true)
	private String intervalStr = null;
	@Parameter(names={"-w","--width"},description="Page width")
	private int drawinAreaWidth = 1000 ;
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition samRecordPartition = SAMRecordPartition.sample;

	private Interval interval=null;
	private Map<String, Sample> sampleHash=new TreeMap<String, Sample>();
	private double featureHeight = 10;
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private final DecimalFormat niceIntFormat = new DecimalFormat("###,###");
	private Document document = null;
	
	private class Sample
		{
		String name;
		List<List<SAMRecord>> lines=new ArrayList<List<SAMRecord>>();
		
		public double getHeight()
			{
			double dim_height= HEIGHT_SAMPLE_NAME;
			dim_height+= SvToSVG.this.featureHeight;//ref seq
			dim_height+= SvToSVG.this.featureHeight;//consensus
			dim_height+= this.lines.size()*SvToSVG.this.featureHeight;//consensus
			return dim_height;
			}
		}
	
	
	/** convert double to string */
	private String format(double v)
		{
		return this.decimalFormater.format(v);
		}
		
	private int trim(int pos0)
		{
		return Math.min(Math.max(pos0,this.interval.getStart()),this.interval.getEnd());
		}


	private int left(final SAMRecord rec)
		{
		return rec.getUnclippedStart();
		}

	private int right(final SAMRecord rec)
		{
		return rec.getUnclippedEnd();
		}
	private double baseToPixel(int pos)
		{
		return  ((pos - this.interval.getStart())/(double)(1+this.interval.getEnd()-this.interval.getStart()))*(this.drawinAreaWidth)
				;
		}
	private void readBamStream(final ConcatSam.ConcatSamIterator iter) throws IOException
		{
		while(iter.hasNext())
			{
			final SAMRecord rec = iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.getReadFailsVendorQualityCheckFlag()) continue;
			if(rec.isSecondaryOrSupplementary()) continue;
			if(rec.getDuplicateReadFlag()) continue;
			final Cigar cigar = rec.getCigar();
			if(cigar==null || cigar.isEmpty()) continue;
			
			if( !rec.getContig().equals(this.interval.getContig())) continue;
			if( right(rec)  < this.interval.getStart()) continue;
			if( left(rec)  > this.interval.getEnd())continue;
			
			String sampleName=this.samRecordPartition.getPartion(rec);
			if(sampleName==null) sampleName = "undefined sample name";
			
			Sample sample= this.sampleHash.get(sampleName);
			if(sample==null)
				{
				sample = new Sample();
				sample.name=sampleName;
				this.sampleHash.put(sampleName,sample);
				}
			
			int y=0;
			for(y=0;y< sample.lines.size();++y)
				{
				final List<SAMRecord> line= sample.lines.get(y);
				final SAMRecord last = line.get(line.size()-1);
				if( right(last)+1 < left(rec) )
					{
					line.add(rec);
					break;
					}
				}
			
			if( y == sample.lines.size())
				{
				final List<SAMRecord> line= new ArrayList<SAMRecord>();
				line.add(rec);
				sample.lines.add(line);
				}
			}
		}
		
		
		
		private Element buildSample(final Sample sample)
			{
			double y=0;
			final Element sampleRoot = element("g");
			
			
			
			final Element sampleLabel= element("text",sample.name);
			sampleLabel.setAttribute("x", "5");
			sampleLabel.setAttribute("y", "12");
			sampleLabel.setAttribute("class", "samplename");
			sampleRoot.appendChild(sampleLabel);
			y+= 20;
				
			
			/* print all lines */
			for(int nLine=0;nLine< sample.lines.size();++nLine)
				{
				final Element g = element("g");
				sampleRoot.appendChild(g);
				g.setAttribute("transform", "translate(0,"+format(y)+")");
				final DocumentFragment fragBack = this.document.createDocumentFragment();
				final DocumentFragment fragFor = this.document.createDocumentFragment();
				final double midy =  this.featureHeight/2.0;
				final double maxy = this.featureHeight;
				final List<SAMRecord> line= sample.lines.get(nLine);
				//loop over records on that line
				for(final SAMRecord record: line)
					{
					final List<SAMRecord> others = SAMUtils.getOtherCanonicalAlignments(record);
					int readpos=0;
					int ref=record.getUnclippedStart();
					final Cigar cigar = record.getCigar();
					for(int i=0;i< cigar.numCigarElements();i++) {
						final CigarElement ce = cigar.getCigarElement(i);
						final CigarOperator op = ce.getOperator();
						int next_ref = ref;
						int next_read = readpos;
						switch(op)
							{
							case I: readpos+= ce.getLength();continue;
							case P: continue;
							case S:case H:
							case M:case X: case EQ:
								next_ref+= ce.getLength();
								next_read+= ce.getLength();
								final double distance_pix = baseToPixel(next_ref)-baseToPixel(ref);
								
								
								final StringBuilder sb=new StringBuilder();
								final double arrow_w = 5;
								final Element path = element("path");
								path.setAttribute("class", "op"+op.name());
								if(i==0 && record.getReadNegativeStrandFlag()) {
									sb.append( "M ").append(format(baseToPixel(ref)+arrow_w)).append(',').append(0);
									sb.append(" h ").append(format(distance_pix-arrow_w));
									sb.append(" v ").append(format(maxy));
									sb.append(" h ").append(format(-(distance_pix-arrow_w)));
									sb.append(" L ").append(format(baseToPixel(ref))).append(',').append(midy);
									sb.append(" Z");
									}
								else if(i+1== cigar.numCigarElements() && !record.getReadNegativeStrandFlag()) {
									sb.append( "M ").append(format(baseToPixel(next_ref)-arrow_w)).append(',').append(0);
									sb.append(" h ").append(format(-(distance_pix-arrow_w)));
									sb.append(" v ").append(format(maxy));
									sb.append(" h ").append(format(distance_pix-arrow_w));
									sb.append(" L ").append(format(baseToPixel(next_ref))).append(',').append(midy);
									sb.append(" Z");
									
									}
								else
									{
									sb.append( "M ").append(format(baseToPixel(ref))).append(',').append(0);
									sb.append(" h ").append(format(distance_pix));
									sb.append(" v ").append(format(maxy));
									sb.append(" h ").append(format(-(distance_pix)));
									sb.append(" Z");
									}
								path.setAttribute("d", sb.toString());
								fragBack.appendChild(path);
								
								if(op.isClipping()) {
									for(final SAMRecord rec2:others) {
										for(AlignmentBlock ab:rec2.getAlignmentBlocks()) {
											if(!CoordMath.overlaps(ab.getReadStart()-1,ab.getReadStart()-1+ab.getLength(),
												readpos,next_read)) continue;
											final Element x= element("rect");
											x.setAttribute("style", "fill:orange;stroke:red;");
											x.setAttribute("x", format(baseToPixel(ref)));
											x.setAttribute("y", format(midy));
											x.setAttribute("width", format(distance_pix));
											x.setAttribute("height", format(this.featureHeight));
											
											final String repeatCount="1";
											final float stepVisibleStart=0f;
											final float visibleDuration=30f;
									
											if(rec2.getContig().equals(record.getContig()))
												{
												Element a = element("animate");
												x.appendChild(a);
												a.setAttribute("attributeType", "XML");
												a.setAttribute("attribute", "y");
												a.setAttribute("from",format(midy));
												a.setAttribute("to",format(-1000));
												}
											else
												{	
												/* move rectangle.x  to destination.x */
												Element w=element("animate");
												x.appendChild(w);
												
												w.setAttribute("attributeType","XML");
												w.setAttribute("attributeName","x");
												w.setAttribute("begin",String.valueOf(stepVisibleStart));
												w.setAttribute("dur",String.valueOf(visibleDuration));
												w.setAttribute("from",format(baseToPixel(ref)));
												w.setAttribute("to",format(baseToPixel(ab.getReferenceStart())));
												w.setAttribute("repeatCount",repeatCount);
												w.setAttribute("fill","freeze");
												
												/* move rectangle.y  to destination.y */
												w=element("animate");
												x.appendChild(w);
												w.setAttribute("attributeType","XML");
												w.setAttribute("attributeName","y");
												w.setAttribute("begin",String.valueOf(stepVisibleStart));
												w.setAttribute("dur",String.valueOf(visibleDuration));
												w.setAttribute("from",format(y));
												w.setAttribute("to",format(y));
												w.setAttribute("repeatCount",repeatCount);
												w.setAttribute("fill","freeze");
												
												if(record.getReadNegativeStrandFlag()!=rec2.getReadNegativeStrandFlag()) {
													w=element("animateTransform");
													x.appendChild(w);
													w.setAttribute("attributeType","XML");
													w.setAttribute("attributeName","transform");
													w.setAttribute("type","rotate");
													w.setAttribute("begin",String.valueOf(stepVisibleStart));
													w.setAttribute("dur",String.valueOf(visibleDuration));
													w.setAttribute("repeatCount",repeatCount);
													//TODO
													//w.setAttribute("from","0 "+rect.getCenterX()+" "+rect.getCenterY());
													//w.setAttribute("to","180 "+next.getCenterX()+" "+next.getCenterY());
													w.setAttribute("fill","freeze");
													}
												}
											}
										}
									}
								
								break;
							case D: case N:
								next_ref+= ce.getLength();
								final Element lineE = element("line");
								lineE.setAttribute("class", "opD");
								lineE.setAttribute("x1", format(baseToPixel(ref)));
								lineE.setAttribute("y1", format(midy));
								lineE.setAttribute("x2", format(baseToPixel(next_ref)));
								lineE.setAttribute("y2", format(midy));
								fragBack.insertBefore(lineE, fragBack.getFirstChild());
								break;
							}
						ref = next_ref;
						readpos = next_read; 
						}
					for(final SAMRecord rec2:SAMUtils.getOtherCanonicalAlignments(record))
						{
						int ref2 = rec2.getAlignmentStart();
						final Cigar cigar2= rec2.getCigar();
						if(cigar2==null || cigar2.isEmpty()) continue;
						
						}
					}
				g.appendChild(fragBack);
				g.appendChild(fragFor);
				
				y+= (this.featureHeight+2);
				}
			
			
			final Element frame= element("rect");
			frame.setAttribute("class","frame");
			frame.setAttribute("x",format(0));
			frame.setAttribute("y",format(0));
			frame.setAttribute("width",format(drawinAreaWidth));
			frame.setAttribute("height",format(y));
			sampleRoot.appendChild(frame);
			
			sampleRoot.setAttribute("height",format(y));
			return sampleRoot;
			}
				
		private Element element(final String tag) {
			return this.document.createElementNS(SVG.NS, tag);
			}
		private Text text(final Object o) {
			return this.document.createTextNode(o==null?"":String.valueOf(o));
			}
		private Element element(final String tag,final Object content) {
			final Element E = element(tag);
			E.appendChild(text(content));
			return E;
			}
		private void buildDocument() 
			{
			final Element svgRoot = element("svg");
			this.document.appendChild(svgRoot);
			svgRoot.setAttribute("width",format(this.drawinAreaWidth));
			
			int doc_height =0;
			
			final Element title = element("title");
			svgRoot.appendChild(title);
			title.appendChild(text(this.intervalStr));
			
			final Element descr = element("description");
			svgRoot.appendChild(descr);
			descr.appendChild(text("Author: Pierre Lindenbaum\n"));
			
			final Element style = element("style");
			svgRoot.appendChild(style);
			style.appendChild(text(
					"g.maing {stroke:black;stroke-width:0.5px;fill:none;}\n"+
					".maintitle {stroke:blue;fill:none;}\n"+
					"rect.frame {stroke:darkgray;fill:none;}\n" + 
					"path.opEQ {stroke:black;fill:gainsboro;}\n" + 
					"path.opX {stroke:black;fill:tomato;}\n" +
					"path.opM {stroke:black;fill:gainsboro;}\n" +
					"path.opS {stroke:black;fill:yellow;}\n" + 
					"path.opH {stroke:black;fill:yellow;}\n" + 
					"line.opN {stroke:black;}\n"+
					"line.opD {stroke:black;}\n"+
					"text.samplename {stroke:none;fill:black;stroke-width:1px;}\n"+
					""
					));

			final Element mainG = element("g");
			mainG.setAttribute("class","maing");
			svgRoot.appendChild(mainG);
			
			//loop over each sample
			for(final Sample sample:this.sampleHash.values())
				{
				final Element div = buildSample(sample);
				final Attr att = div.getAttributeNode("height");
				div.removeAttributeNode(att);
				div.setAttribute("transform", "translate(0,"+format(doc_height)+")");
				doc_height += Double.parseDouble(att.getValue());
				mainG.appendChild(div);
				}
			
			svgRoot.setAttribute("height",format(doc_height));
			}
		
		
		
		@Override
		public int doWork(final List<String> args) {
			/* parse interval */
			if(StringUtil.isBlank(this.intervalStr))
				{
				LOG.error("interval undefined");
				return -1;
				}	
			this.drawinAreaWidth = Math.max(100,this.drawinAreaWidth );
			ConcatSam.ConcatSamIterator iter = null;
			FileOutputStream fout=null;
			try
				{
				final ConcatSam.Factory concat = new ConcatSam.Factory();
				concat.addInterval(this.intervalStr);
				concat.setEnableUnrollList(true);
				iter = concat.open(args);
				final IntervalParser parser=new IntervalParser(iter.getFileHeader().getSequenceDictionary());
				this.interval = parser.parse(this.intervalStr);
				
				readBamStream(iter);
				iter.close(); iter = null;
				
				
				
				
				final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
				dbf.setNamespaceAware(true);
				final DocumentBuilder db = dbf.newDocumentBuilder();
				this.document=db.newDocument();
				
				buildDocument();
				
				final Transformer tr = TransformerFactory.newInstance().newTransformer();
				final Result result;
				
				if(this.outputFile!=null)
					{
					result = new StreamResult(this.outputFile);
					}
				else
					{
					result = new StreamResult(stdout());
					}
				tr.transform(new DOMSource(this.document),result);
				
				return RETURN_OK;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(iter);
				CloserUtil.close(fout);
				this.interval=null;
				this.document = null;
				}
			}
		
	
	public static void main(final String[] args)
		{
		new SvToSVG().instanceMainWithExit(args);
		}
	}
