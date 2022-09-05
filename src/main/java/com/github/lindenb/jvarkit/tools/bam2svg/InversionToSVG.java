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
package com.github.lindenb.jvarkit.tools.bam2svg;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Result;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StaticCodeExtractor;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Example

```bash
$ samtools view -b "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA20845/high_coverage_alignment/NA20845.wgs.ILLUMINA.bwa.GIH.high_cov_pcr_free.20140203.bam" "7:8352157-8356597" > jeter.bam && samtools index jeter.bam
$ java -jar dist/svg2svg.jar jeter.bam > jeter.svg
```

### A translocation


Translocation described in [PMC5932280](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5932280/)

> Identification of Balanced Chromosomal Rearrangements Previously Unknown Among Participants in the 1000 Genomes Project: Implications for Interpretation of Structural Variation in Genomes and the Future of Clinical Cytogenetics

```
$ samtools view -b "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02260/alignment/HG02260.mapped.ILLUMINA.bwa.PEL.low_coverage.20130415.bam" "9:137229907-137231907" > jeter1.bam
$ samtools view -b "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02260/alignment/HG02260.mapped.ILLUMINA.bwa.PEL.low_coverage.20130415.bam" "14:79838174-79840174" > jeter2.bam
$ samtools merge jeter3.bam jeter1.bam jeter2.bam
$ samtools index jeter3.bam
$ java -jar dist/sv2svg.jar -r "9:137229907-137231907" -r "14:79838174-79840174"  jeter3.bam > jeter.svg
```


## Gallery

[https://gist.github.com/lindenb/bf48989b8da31eeafdc2caa0694361eb](https://gist.github.com/lindenb/bf48989b8da31eeafdc2caa0694361eb)

[https://twitter.com/yokofakun/status/1063369955406725120](https://twitter.com/yokofakun/status/1063369955406725120)

![https://imgur.com/EVOrXuc](https://i.imgur.com/EVOrXuc.gif)

[https://twitter.com/yokofakun/status/1063511215832539136](https://twitter.com/yokofakun/status/1063511215832539136)

![https://pbs.twimg.com/media/DsJZKRrWsAE_QqA.jpg](https://pbs.twimg.com/media/DsJZKRrWsAE_QqA.jpg)

[https://gist.github.com/lindenb/877d1d00d9f19c618f2d8505a2fe5614](https://gist.github.com/lindenb/877d1d00d9f19c618f2d8505a2fe5614)

[https://twitter.com/yokofakun/status/1064484537059684355](https://twitter.com/yokofakun/status/1064484537059684355)

![https://twitter.com/yokofakun/status/1064484537059684355](https://pbs.twimg.com/media/DsXObwuXoAAepSw.jpg)

[https://twitter.com/yokofakun/status/1064503996285681666](https://twitter.com/yokofakun/status/1064503996285681666)

![https://pbs.twimg.com/media/DsXgnKrWwAAYiBS.jpg](https://pbs.twimg.com/media/DsXgnKrWwAAYiBS.jpg)


[https://gist.github.com/lindenb/88bb702478cafe732d00e2694e77bc09](https://gist.github.com/lindenb/88bb702478cafe732d00e2694e77bc09)

![https://twitter.com/yokofakun/status/1448993550746759169](https://pbs.twimg.com/media/FBva8sTWUAMNApn?format=jpg&name=large)

END_DOC
 */
@Program(name="sv2svg",
description="BAM to SVG. Used to display inversions.",
keywords={"bam","inversions","alignment","graphics","visualization","svg","structural","variant"},
modificationDate="20220905",
creationDate="20220905",
generate_doc=false
)
public class InversionToSVG extends Launcher
	{
	private static final Logger LOG = Logger.build(InversionToSVG.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-r","-i","--interval","--region"},description="interval CHROM:START-END",required=true)
	private String intervalStr = "";
	@Parameter(names={"-w","--width"},description="Page width")
	private int drawinAreaWidth = 1000 ;
	@Parameter(names={"-R","--reference"},description= "For CRAM. " + INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidxFile = null ;
	@Parameter(names={"-d","--duration"},description="Animation duration, in secs. <=0 disable animation.")
	private int svgDuration=10;
	@Parameter(names= {"--repeat-count"},description="SVG animation repeat count")
	private String svgRepeatCount="indefinite";
	@Parameter(names= {"--svg"},description="Write SVG only document. Default is to write a XHTML+SVG document.")
	private boolean write_svg_only = false;

	
	private double featureHeight = 10;
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private Document document = null;
	private final double arrow_w = 5;
	
	
	private class ReadName implements Locatable {
		final List<SAMRecord> records = new ArrayList<>();

		@Override
		public String getContig() {
			return records.get(0).getContig();
		}
		@Override
		public int getStart() {
			return records.stream().mapToInt(R->R.getUnclippedStart()).min().getAsInt();
			}
		@Override
		public int getEnd() {
			return records.stream().mapToInt(R->R.getUnclippedEnd()).max().getAsInt();
			}
	}
	
	private class Window implements Locatable {
		Locatable loc;
		Element clipRect;
		Element g;
		final double x;
		final double width;
		double height = 0;
		Window(Locatable loc,double x,double width) {
			this.loc = loc;
			this.x = x;
			this.width = width;
			}
		double baseToPixel(double pos) {
			return ((pos - getStart())/(double)getLengthOnReference())*width;
			}
		@Override
		public String getContig() {
			return loc.getContig();
		}
		@Override
		public int getStart() {
			return loc.getStart();
			}
		@Override
		public int getEnd() {
			return loc.getEnd();
			}
		public String getClipId(){
			return StringUtils.md5(loc.toString());
		}
	}

		
	/** convert double to string */
		private String format(final double v)
			{
			return this.decimalFormater.format(v);
			}
		
		private Element html(final String tag) {
			return this.document.createElement(tag);
			}

		private Element html(final String tag,final String content) {
			final Element E = html(tag);
			E.appendChild(text(content));
			return E;
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
		
		
		
		@Override
		public int doWork(final List<String> args) {
			
			this.drawinAreaWidth = Math.max(100,this.drawinAreaWidth );
			try
				{
				
				
				final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
				dbf.setNamespaceAware(true);
				final DocumentBuilder db = dbf.newDocumentBuilder();
				this.document=db.newDocument();
				int extend = 50;
				
				
				final Path bamFile = Paths.get(this.oneAndOnlyOneFile(args));
				IOUtil.assertFileIsReadable(bamFile);
				
				final SamReaderFactory srf = super.createSamReaderFactory().referenceSequence(this.faidxFile);
				
				final Map<String,ReadName> readNameMap = new HashMap<>();
				/* loop over each bam file */
					try(final SamReader sr = srf.open(bamFile)) {
						final SAMFileHeader header = sr.getFileHeader();
						final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
						final Function<String,Optional<SimpleInterval>> intervalParser = IntervalParserFactory.newInstance().dictionary(dict).make();
						final SimpleInterval interval = intervalParser.apply(intervalStr).orElseThrow(IntervalParserFactory.exception(intervalStr));
	
						final SimpleInterval leftInterval = new SimpleInterval(
								interval.getContig(),
								interval.getStart()-extend,
								interval.getStart()+extend
								);
						
						final SimpleInterval rightInterval = new SimpleInterval(
								interval.getContig(),
								interval.getEnd()-extend,
								interval.getEnd()+extend
								);
						final Window[] windows;
						if(leftInterval.overlaps(rightInterval)) {
							windows= new Window[]{
									new Window(new SimpleInterval(leftInterval.getContig(),leftInterval.getStart(),rightInterval.getEnd()),0,0)
									};
						} else {
							double w2 = this.drawinAreaWidth/2.0;
							windows= new Window[]{
								new Window(leftInterval,0,w2),
								new Window(rightInterval,w2,w2)
								};
						}
	
						
						final String sampleName = header.getReadGroups().stream().
								map(G->G.getSample()).
								filter(S->!StringUtil.isBlank(S)).
								findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamFile));
						try(CloseableIterator<SAMRecord> iter = sr.query(interval.getContig(), interval.getStart(), interval.getEnd(), false) ) {
							while(iter.hasNext()) {
								final SAMRecord record = iter.next();
								if(!SAMRecordDefaultFilter.accept(record)) continue;
								if(
									!CoordMath.overlaps(record.getUnclippedStart(), record.getUnclippedEnd(), leftInterval.getStart(), leftInterval.getEnd()) &&
									!CoordMath.overlaps(record.getUnclippedStart(), record.getUnclippedEnd(), rightInterval.getStart(), rightInterval.getEnd())
									)
									{
									continue;
									}
								final Cigar cigar = record.getCigar();
								if(cigar==null || cigar.isEmpty()) continue;
								ReadName rn = readNameMap.get(record.getReadName());
								if(rn==null) {
									rn = new ReadName();
									readNameMap.put(record.getReadName(), rn);
									}
								rn.records.add(record);
								
								for(SAMRecord rec2 : SAMUtils.getOtherCanonicalAlignments(record)) {
									if(!rec2.contigsMatch(record)) continue;
									if(
										!CoordMath.overlaps(rec2.getUnclippedStart(), rec2.getUnclippedEnd(), leftInterval.getStart(), leftInterval.getEnd()) &&
										!CoordMath.overlaps(rec2.getUnclippedStart(), rec2.getUnclippedEnd(), rightInterval.getStart(), rightInterval.getEnd())
										)
										{
										continue;
										}
									rn.records.add(rec2);
									}
								
								}/* end iterator */
							} /* end loop interval */
						
					for(String rn : readNameMap.keySet()) {
						Collections.sort(readNameMap.get(rn).records,(A,B)->Integer.compare(A.getStart(),B.getStart()));
						}
					
	
					final Element svgRoot = element("svg");
					final Element metadata = element("metadata");
					svgRoot.appendChild(metadata);
					metadata.setAttribute("id", "metadata");
					
					if(write_svg_only) {
						this.document.appendChild(svgRoot);
						}
					else
						{
						final Element html = html("html");
						this.document.appendChild(html);
						final Element head = html("head");
						html.appendChild(head);
						final Element title = html("title");
						title.appendChild(text(interval.toNiceString()));
	
						head.appendChild(title);
	
						Element meta = html("meta");
						meta.setAttribute("charset", "UTF-8");
						head.appendChild(meta);
	
						Element style=html("style");
						head.appendChild(style);
						
						style.appendChild(text(
								"body { font-family: sans-serif;}\n" 
								+ "table {  border-collapse: collapse; margin: 25px 0;  font-size: 0.9em; font-family: sans-serif;min-width: 400px;box-shadow: 0 0 20px rgba(0, 0, 0, 0.15);}\n"
								+"thead{ background-color: #009879; color: #ffffff; text-align: left;}\n"
								+"th,td { padding: 12px 15px}\n"
								+"tbody tr {border-bottom: 1px solid #dddddd;}\n"
								+"tbody tr:nth-of-type(even) {background-color: #f3f3f3;}\n"
								+"tbody tr:last-of-type { border-bottom: 2px solid #009879;}\n"
								+"tbody tr.active-row { font-weight: bold;color: #009879;}\n"
								));
						final Element body = html("body");
						html.appendChild(body);
	
						final Element h1 = html("h1");
						body.appendChild(h1);
						h1.appendChild(text(interval.toNiceString()));				
						Element htmlDiv1 = html("div");
						body.appendChild(htmlDiv1);
						
						htmlDiv1.setAttribute("style","display:none;");
						htmlDiv1.setAttribute("id","__PLACEHOLDER__");
						htmlDiv1.appendChild(document.createComment("Use this div to insert things later. "));
						htmlDiv1.appendChild(text("__PLACEHOLDER__"));
	
						htmlDiv1 = html("div");
						body.appendChild(htmlDiv1);
	
						
						htmlDiv1.appendChild(svgRoot);
						
						htmlDiv1.appendChild(html("hr"));
						Element p = html("p");
						htmlDiv1.appendChild(p);
						p.appendChild(text("Pierre Lindenbaum PhD 2022. Made with "+getProgramName()+" version:"+getVersion()+
							". Command was: " + getProgramCommandLine()
							));
						}
					
					svgRoot.setAttribute("width",format(this.drawinAreaWidth+1));
					
					int doc_height =0;
					
					final Element title = element("title");
					svgRoot.appendChild(title);
					title.appendChild(text(interval.toString()));
					
					final Element defs = element("defs");
					svgRoot.appendChild(defs);
	
					
					final Element descr = element("desc");
					svgRoot.appendChild(descr);
					descr.appendChild(text("Author: Pierre Lindenbaum"));
					
					final Element style = element("style");
					svgRoot.appendChild(style);
					style.appendChild(text(StaticCodeExtractor.forClass(InversionToSVG.class).extract("CSS").get()));
					
					final Element mainG = element("g");
					mainG.setAttribute("class","maing");
					svgRoot.appendChild(mainG);
					
					
					for(Window win: windows) {
						final Element clipPath = element("clipPath");
						clipPath.setAttribute("id",win.getClipId());
						defs.appendChild(clipPath);
						win.clipRect = element("rect");
						clipPath.appendChild(win.clipRect);
						win.clipRect.setAttribute("x",format(0));
						win.clipRect.setAttribute("y",format(0));
						win.clipRect.setAttribute("width",format(win.width));
						win.clipRect.setAttribute("height",format(0));
						
						win.g = element("g");
						win.g.setAttribute("transform", "translate("+win.x+",0)");
						win.g.setAttribute("clip-path","url(#"+win.getClipId()+")");
						mainG.appendChild(win.g);
					}
					
					//loop over each sample
					double y1=0;
					int record_index = 0;
					for(ReadName readName: readNameMap.values()) {
						double y_sample = y1;
						List<Element> backgrounds = new ArrayList<>(readName.records.size());
						
						for(SAMRecord rec:readName.records) {
							double y_rec = y_sample;
							final String readTitle = rec.getReadName()+" "+(rec.getReadNegativeStrandFlag()?"-":"+")+" "+
									(rec.getSupplementaryAlignmentFlag()?"supplementary":"");
							
							for(Window win: windows) {
									final Element g = element("g");
									g.setAttribute("transform", "translate(0,"+y_rec+")");
									win.g.appendChild(g);
								
									final Element background = element("rect");
									backgrounds.add(background);
									
									win.g.appendChild(background);
									background.setAttribute("x", "0");
									background.setAttribute("y", "0");
									background.setAttribute("height", "0");
									background.setAttribute("width", format(win.width));
									
									if(rec.overlaps(win)) {
									
									final double midy =  this.featureHeight/2.0;
									final double maxy = this.featureHeight;
		
									
									
									Element line = element("line");
									line.setAttribute("x1",format(win.baseToPixel(rec.getUnclippedStart())));
									line.setAttribute("x2",format(win.baseToPixel(rec.getUnclippedEnd()+1)));
									line.setAttribute("y1",format(midy));
									line.setAttribute("y2",format(midy));
									line.appendChild(element("title",readTitle));
									g.appendChild(line);
									
									int ref1 = rec.getUnclippedStart();
									final Cigar cigar = rec.getCigar();
									final List<Element> rectInserts = new ArrayList<>();
									for(int cigarIdx=0;cigarIdx<cigar.numCigarElements();++cigarIdx) {
										final CigarElement ce = cigar.getCigarElement(cigarIdx);
										final CigarOperator op = ce.getOperator();
										final double leftX = win.baseToPixel(ref1);
										int next_ref = ref1;
										switch(op)
											{
											case I: 
												{
												if(CoordMath.overlaps(ref1, ref1+ce.getLength(),win.getStart(),win.getEnd())) {
													final Element rectInsert = element("rect");
													rectInsert.setAttribute("class", "insert");
													rectInsert.setAttribute("x", format(leftX));
													rectInsert.setAttribute("y", format(0));
													rectInsert.setAttribute("width", format(1));
													rectInsert.setAttribute("height", format(this.featureHeight));
													rectInsert.appendChild(element("title",readTitle+" op:"+op.name()+" "+ format(ce.getLength())));
													rectInserts.add(rectInsert);
													}
												break;
												}
											case P: continue;
											case S:case H:
											case M:case X: case EQ:
												{
												next_ref += ce.getLength();
													
												if(!(ref1 > win.getEnd() || next_ref < win.getStart())) {
													final double distance_pix = win.baseToPixel(next_ref+1)-leftX;
													
													final StringBuilder sb=new StringBuilder();
													
													final Element path = element("path");
													path.setAttribute("class", "op"+op.name()+(rec.getSupplementaryAlignmentFlag()?"x":""));
													path.appendChild(element("title",readTitle+" op:"+op.name()));
													
													// arrow <--
													if(cigarIdx==0 && rec.getReadNegativeStrandFlag()) {
														sb.append( "M ").append(format(leftX)).append(',').append(0);
														sb.append(" h ").append(format(distance_pix));
														sb.append(" v ").append(format(maxy));
														sb.append(" h ").append(format(-(distance_pix)));
														sb.append(" l ").append(format(-arrow_w)).append(',').append(-featureHeight/2.0);
														sb.append(" Z");
														}
													// arrow -->
													else if(cigarIdx+1==cigar.numCigarElements() && !rec.getReadNegativeStrandFlag()) {
														sb.append( "M ").append(format(leftX+distance_pix)).append(',').append(0);
														sb.append(" h ").append(format(-(distance_pix)));
														sb.append(" v ").append(format(maxy));
														sb.append(" h ").append(format(distance_pix));
														sb.append(" l ").append(format(arrow_w)).append(',').append(-featureHeight/2.0);
														sb.append(" Z");
														
														}
													else
														{
														sb.append( "M ").append(format(leftX)).append(',').append(0);
														sb.append(" h ").append(format(distance_pix));
														sb.append(" v ").append(format(maxy));
														sb.append(" h ").append(format(-(distance_pix)));
														sb.append(" Z");
														}
													path.setAttribute("d", sb.toString());
													g.appendChild(path);
													}
												break;
												}
											case D: case N:
												{
												next_ref+= ce.getLength();
												break;
												}
											default: throw new IllegalStateException(op.name());
											}
										ref1 = next_ref;
										} // end cigar
									for(Element insertE:rectInserts) {
										g.appendChild(insertE);
										}
									} //end loop over record
								y_rec += this.featureHeight;
								y_rec += 3;
								}
							y_sample = Math.max(y_rec,y_sample);
							}
						y1 = y_sample;
						}
					svgRoot.setAttribute("height", format(y1+10));
					for(Window win:windows) {
						win.height = y1;
						// frame
						final Element frame = element("rect");
						frame.setAttribute("class", "frame");
						frame.setAttribute("x", format(win.x));
						frame.setAttribute("y", "0");
						frame.setAttribute("width",format(win.width));
						frame.setAttribute("height",format(win.height));
						win.g.appendChild(frame);
						
						win.clipRect.setAttribute("height", format(win.height));
						
						mainG.appendChild(frame);
						}
					}/* end open SAM */
				
				final Transformer tr = TransformerFactory.newInstance().newTransformer();
				
				if(!this.write_svg_only) tr.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
				tr.setOutputProperty(OutputKeys.METHOD, "xml");

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
				
				return 0;
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				this.document = null;
				}
			}
		
	
	public static void main(final String[] args)
		{
		new InversionToSVG().instanceMainWithExit(args);
		}
	}

/**
BEGIN_CSS

g.maing {stroke:black;stroke-width:0.5px;fill:none;}
.maintitle {stroke:blue;fill:none;}
rect.frame {stroke:darkgray;fill:none;} 
path.opEQ {stroke:black;fill:gainsboro;} 
path.opX {stroke:black;fill:tomato;}
path.opM {stroke:black;fill:gainsboro;}
path.opS {stroke:black;fill:yellow;} 
path.opH {stroke:black;fill:yellow;} 
line.opN {stroke:black;}
line.opD {stroke:black;}
path.opEQx {stroke:yellow;fill:salmon;opacity:0.8;} 
path.opXx {stroke:yellow;fill:tomato;opacity:0.8;}
path.opMx {stroke:yellow;fill:salmon;opacity:0.8;}
path.opSx {stroke:yellow;fill:yellow;opacity:0.8;} 
path.opHx {stroke:yellow;fill:yellow;opacity:0.8;} 
line.opNx {stroke:yellow;opacity:0.8;}
line.opDx {stroke:yellow;opacity:0.8;}
rect.mismatch {fill:red;opacity:0.7;}
rect.insert {fill:none;stroke:red;stroke-width:1.5px}
text.samplename {stroke:none;fill:black;stroke-width:1px;}
.discordant {stroke:darkred; fill:none;stroke-dasharray:2;}
line.ruler {stroke:darkgray;stroke-dasharray:4;fill:none;stroke-width:1;}
rect.variant  {stroke:none;fill:orange;stroke-width:1px;opacity:0.8;}
rect.variantHOM_REF  {stroke:none;fill:green;stroke-width:1px;opacity:0.8;}
rect.variantHOM_VAR  {stroke:none;fill:red;stroke-width:1px;opacity:0.8;}
rect.variantHET  {stroke:none;fill:orange;stroke-width:1px;opacity:0.8;}
rect.variantNO_CALL  {stroke:none;fill:blue;stroke-width:1px;opacity:0.8;}
path.coverage {stroke:darkslateblue;fill:darkseaGreen}


END_CSS
*/
