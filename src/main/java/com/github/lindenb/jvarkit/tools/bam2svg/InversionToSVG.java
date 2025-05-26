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
import java.util.function.IntToDoubleFunction;
import java.util.function.ToIntFunction;

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
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.svg.SVG;

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
@Program(name="inversion2svg",
description="BAM to SVG. Used to display inversions.",
keywords={"bam","inversions","alignment","graphics","visualization","svg","structural","variant"},
modificationDate="20220905",
creationDate="20220905",
generate_doc=false
)
public class InversionToSVG extends Launcher
	{
	private static final Logger LOG = Logger.of(InversionToSVG.class);


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
	private static final double arrow_w = 5;
	final Function<SAMRecord, Integer> mateEnd = REC->SAMUtils.getMateCigar(REC)!=null?
			SAMUtils.getMateAlignmentEnd(REC):
			REC.getMateAlignmentStart()
			;

	
	private class Segment implements Locatable {
		final String readName;
		final String contig;
		int startA;
		int endA;
		boolean revA;
		int startB;
		int endB;
		boolean revB;

		Segment(SAMRecord rec) {
			this.readName = rec.getReadName();
			this.contig = rec.getContig();
			this.startA = rec.getStart();
			this.endA = rec.getEnd();
			this.startB = rec.getMateAlignmentStart();
			this.endB = mateEnd.apply(rec);
			}
		@Override
		public String getContig()
			{
			return contig;
			}
		@Override
		public int getStart()
			{
			return Math.min(startA, startB);
			}
		@Override
		public int getEnd()
			{
			return Math.max(endA, endB);
			}
		}
	
	
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
		public String getReadName() {
			return records.get(0).getReadName();
			}
		boolean isInversion(final Locatable r) {

			for(int i=0;i< records.size();i++) {
				final SAMRecord R1 = records.get(i);
				for(int side=0;side<2;side++) {
					int pos1 =side==0?r.getStart():r.getEnd();
					int pos2=side==0?r.getEnd():r.getStart();
					if(CoordMath.overlaps(pos1,pos1, R1.getUnclippedStart(), R1.getUnclippedEnd())) {
						for(SAMRecord RS : SAMUtils.getOtherCanonicalAlignments(R1)) {
							if(!R1.contigsMatch(RS)) continue;
							if(R1.getReadNegativeStrandFlag()==RS.getReadNegativeStrandFlag()) continue;
							if(CoordMath.overlaps(pos2,pos2, RS.getUnclippedStart(), RS.getUnclippedEnd())) {
								return true;
								}
							}
						}
					}
				if(R1.getReadPairedFlag() && 
					!R1.getMateUnmappedFlag()  &&
					R1.getReferenceName().equals(R1.getMateReferenceName()) &&
					R1.getReadNegativeStrandFlag()==R1.getMateNegativeStrandFlag()
					)
					{
					int mpos = R1.getMateAlignmentStart();
					int mend = mateEnd.apply(R1);
					//R1 before/after breakpoint, forward and R2
					if( (R1.getStart()/*oui, pas end */ < r.getStart() || R1.getEnd()/*oui*/ > r.getEnd() ) &&
						mpos > r.getStart() && mend < r.getEnd()
						) {
						return true;
						}
					
					}
				
				}
			return false;
			}
		}
	
	private class Window implements Locatable {
		final Locatable loc;
		final double x;
		final double width;
		double y = 0;
		double height = 0;
		Window(Locatable loc,double x,double width) {
			this.loc = loc;
			this.x = x;
			this.width = width;
			}
		public double getWidth() {
			return this.width;
			}
		public double getMinX() {
			return this.x;
			}
		public double getMaxX() {
			return getMinX()+ this.getWidth();
			}
		
		double baseToPixel(double pos) {
			return getMinX() + ((pos - getStart())/(double)getLengthOnReference())*getWidth();
			}
		
		double trim(double x) {
			return Math.min(Math.max(x, getMinX()),getMaxX());
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
					final Function<String,Optional<SimpleInterval>> intervalParser = new IntervalParser(dict);
					final SimpleInterval interval = intervalParser.apply(intervalStr).orElseThrow(IntervalParser.exception(intervalStr));

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
					
					final List<Segment> segments = new ArrayList<>();
					
					final String sampleName = header.getReadGroups().stream().
							map(G->G.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamFile));
					try(CloseableIterator<SAMRecord> iter = sr.query(interval.getContig(), interval.getStart(), interval.getEnd(), false) ) {
						while(iter.hasNext()) {
							final SAMRecord record = iter.next();
							if(!SAMRecordDefaultFilter.accept(record)) continue;
							
							if(record.getReadPairedFlag() &&
								!record.getMateUnmappedFlag() && 
								record.getReadNegativeStrandFlag()==record.getMateNegativeStrandFlag() &&
								record.getReferenceName().equals(record.getMateReferenceName())
								) {
								if(segments.stream().noneMatch(S->S.readName.equals(record.getReadName()))) {
									segments.add(new Segment(record));
									}
								}
							
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
					title.appendChild(text(sampleName +" "+interval.toNiceString()));

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
					h1.appendChild(text(sampleName+" "+interval.toNiceString()));				
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
				
				
				final Element title = element("title");
				svgRoot.appendChild(title);
				title.appendChild(text(interval.toString()));
				
				final Element defs = element("defs");
				svgRoot.appendChild(defs);
				for(int suppl=0;suppl<2;++suppl) {
					for(CigarOperator op:CigarOperator.values()) {
						if(!(op.isAlignment() || op.isClipping())) continue;
						// pattern
						for(int side=0;side<2;++side) {
							final double w = this.featureHeight;
							final double w2 = (w*2.0)*0.95;
							final double h2 = w/2.0;
							final Element patt = element("pattern");
							defs.appendChild(patt);
							patt.setAttribute("id",op.name()+(side==0?"F":"R")+(suppl==1?"x":""));
							patt.setAttribute("x","0");
							patt.setAttribute("y","0");
							patt.setAttribute("width",format(w*2));
							patt.setAttribute("height",format(w));
							patt.setAttribute("patternUnits","userSpaceOnUse");
							
							final Element rect = element("rect");
							patt.appendChild(rect);
							rect.setAttribute("x","0");
							rect.setAttribute("y","0");
							rect.setAttribute("width",format(w*2));
							rect.setAttribute("height",format(w));
							String css= "stroke:none;fill:";
							switch(op) {
								case EQ: css+="gainsboro";break;
								case M: css+="gainsboro";break;
								case X: css+="tomato";break;
								case S: css+="yellow";break;
								case H: css+="yellow";break;
								default: css+="green";break;
								}
							css+=";";
							rect.setAttribute("style",css);
							
							final Element arrow = element("path");
							arrow.setAttribute("style", "fill:darkgray;stroke:"+(suppl==1?"pink":"blue")+";");
							patt.appendChild(arrow);
							final StringBuilder sb = new StringBuilder();
							if(side==0) {
								sb.append("M ").append(0).append(",").append(h2);
								sb.append(" l ").append(w2).append(",0");
								sb.append(" l ").append(-h2).append(",").append(-h2);
								sb.append(" l ").append(0).append(",").append(2*h2);
								sb.append(" l ").append(h2).append(",").append(-h2);
								}
							else {
								sb.append("M ").append(w2).append(",").append(h2);
								sb.append(" l ").append(-w2).append(",0");
								sb.append(" l ").append(h2).append(",").append(-h2);
								sb.append(" l ").append(0).append(",").append(2*h2);
								sb.append(" l ").append(-h2).append(",").append(-h2);
								}
							sb.append(" Z");
							arrow.setAttribute("d", sb.toString());
							}
						} // end operator
					} // end suppl
				
				
				final Element descr = element("desc");
				svgRoot.appendChild(descr);
				descr.appendChild(text("Author: Pierre Lindenbaum"));
				
				final Element style = element("style");
				svgRoot.appendChild(style);
				style.appendChild(text(""
						 + "g.maing {stroke:black;stroke-width:0.5px;fill:none;} "
						 + ".maintitle {stroke:blue;fill:none;} "
						 + "rect.frame {stroke:darkgray;fill:none;} "
						 + "path.opEQ {stroke:black;fill:gainsboro;} "
						 + "path.opX {stroke:black;fill:tomato;} "
						 + "path.opM {stroke:black;fill:gainsboro;} "
						 + "path.opS {stroke:black;fill:yellow;} "
						 + "path.opH {stroke:black;fill:yellow;} "
						 + "line.opN {stroke:black;} "
						 + "line.opD {stroke:black;} "
						 + "path.opEQx {stroke:yellow;fill:salmon;opacity:0.8;} "
						 + "path.opXx {stroke:yellow;fill:tomato;opacity:0.8;} "
						 + "path.opMx {stroke:yellow;fill:salmon;opacity:0.8;} "
						 + "path.opSx {stroke:yellow;fill:yellow;opacity:0.8;} "
						 + "path.opHx {stroke:yellow;fill:yellow;opacity:0.8;} "
						 + "line.opNx {stroke:yellow;opacity:0.8;} "
						 + "line.opDx {stroke:yellow;opacity:0.8;} "
						 + "rect.insert {fill:none;stroke:red;stroke-width:1.5px} "
						 + "text.samplename {stroke:none;fill:black;stroke-width:1px;} "
						 + "text.intervalLbl {stroke:none;fill:black;stroke-width:1px;} "
						 + ".discordant {stroke:darkred; fill:none;stroke-dasharray:2;} "
						 + "line.ruler {stroke:darkgray;stroke-dasharray:4;fill:none;stroke-width:1;} "
						 + "rect.variant  {stroke:none;fill:orange;stroke-width:1px;opacity:0.8;} "
						 + "rect.variantHOM_REF  {stroke:none;fill:green;stroke-width:1px;opacity:0.8;} "
						 + "rect.variantHOM_VAR  {stroke:none;fill:red;stroke-width:1px;opacity:0.8;} "
						 + "rect.variantHET  {stroke:none;fill:orange;stroke-width:1px;opacity:0.8;} "
						 + "rect.variantNO_CALL  {stroke:none;fill:blue;stroke-width:1px;opacity:0.8;} "
						 + "path.coverage {stroke:darkslateblue;fill:darkseaGreen} "
						 + "rect.readrect0 {stroke:none;fill:lavender;} "
						 + "rect.readrect1 {stroke:none;fill:lavenderblush;} "
						 + "path.rR {fill:url(#rR)} "
						 + "path.rF {fill:url(#rF)} "
						 + "line.breakpoint {stroke:black;stroke-dasharray:5;stroke-width:2px;opacity:0.8;} "
						 + "rect.frame {stroke:black;fill:none;} "
						));
				
				final Element mainG = element("g");
				mainG.setAttribute("class","maing");
				svgRoot.appendChild(mainG);
				
				
				
				//loop over each sample
				double y1=20;
				
				for(Window win: windows) {
					win.y = y1;
					final Element txt = element("text",win.loc.toString());
					txt.setAttribute("class", "intervalLbl");
					txt.setAttribute("x", format(win.getMinX()+3));
					txt.setAttribute("y", format(y1-3));
					mainG.appendChild(txt);
					}
				
				int read_index=0;
				List<ReadName> sortedReadNames = new ArrayList<>(readNameMap.values());
				// remove non-informative reads.
				sortedReadNames.removeIf(R->R.records.stream().allMatch(RR->RR.getUnclippedStart() > interval.getStart() && RR.getUnclippedEnd() < interval.getEnd()));
				
				
				final ToIntFunction<SAMRecord> recToScore=A->{
					int n=0;
					if(A.getCigar().isClipped()) n+=10;
					for(Window win:windows) {
						if(win.overlaps(A)) n+=10;
						}
					if(A.getReadPairedFlag() && !A.getProperPairFlag()) {
						n+=10;
						}
					return n;
					};
				Collections.sort(sortedReadNames,(A,B)->{
					int scoreA = A.records.stream().mapToInt(recToScore).sum();
					int scoreB = B.records.stream().mapToInt(recToScore).sum();
					return Integer.compare(scoreB,scoreA);
					});
				for(ReadName readName: sortedReadNames) {
					read_index++;
					final Element readG = element("g");
					readG.setAttribute("transform", "translate(0,"+y1+")");
					mainG.appendChild(readG);
					final Element readRect = element("rect");
					readG.appendChild(readRect);
					readRect.setAttribute("x", "0");
					readRect.setAttribute("y", format(0));
					readRect.setAttribute("width", format(this.drawinAreaWidth));
					readRect.setAttribute("height", format(0));
					readRect.setAttribute("class", "readrect"+(read_index%2));
					readRect.appendChild(element("title", readName.getReadName()));
					
					y1+=3;
					double y_read_name = y1;
					
					double y_rec=0;
					for(SAMRecord rec:readName.records) {
						final String readTitle = rec.getReadName()+" "+
								(rec.getReadNegativeStrandFlag()?"-":"+")+" "+
								(rec.getSupplementaryAlignmentFlag()?"supplementary":"")+
								(rec.getReadPairedFlag()?(rec.getFirstOfPairFlag()?"R1":"R2"):"");
						
						for(Window win: windows) {
								if(!rec.overlaps(win)) continue;

								final Element g = element("g");
								g.setAttribute("transform", "translate(0,"+format(y_rec)+")");
								readG.appendChild(g);
															
								
								final double midy =  this.featureHeight/2.0;
								final double maxy = this.featureHeight;
	
								
								
								final Element line = element("line");
								line.setAttribute("class", "del");
								line.setAttribute("x1",format(win.trim(win.baseToPixel(rec.getUnclippedStart()))));
								line.setAttribute("x2",format(win.trim(win.baseToPixel(rec.getUnclippedEnd()+1))));
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
									int next_ref = ref1;
									switch(op)
										{
										case I: 
											{
											final double x1 = win.baseToPixel(ref1);
											final double x2 = win.baseToPixel(ref1+ce.getLength());
											
											if(!(x2 < win.getMinX() || x1 > win.getMaxX())) {
												final Element rectInsert = element("rect");
												rectInsert.setAttribute("class", "insert");
												rectInsert.setAttribute("x", format(win.trim(x1)));
												rectInsert.setAttribute("y", format(0));
												rectInsert.setAttribute("width", format(win.trim(x2)-win.trim(x1)));
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
											
											double x1 = win.baseToPixel(ref1);
											double x2 = win.baseToPixel(ref1+ce.getLength());

											
											if(!(x2 < win.getMinX() || x1 > win.getMaxX())) {
												x1 = win.trim(x1);
												x2 = win.trim(x2);
												final double distance_pix = x2-x1;
												
												final StringBuilder sb=new StringBuilder();
												
												final Element path = element("path");
												path.setAttribute("style",
														(readName.isInversion(interval)?"stroke:green;stroke-width:2;":"stroke:gray;")+
														"fill:url(#"+op.name()+(rec.getReadNegativeStrandFlag()?"R":"F")+(rec.isSecondaryOrSupplementary()?"x":"")+")" );
												path.appendChild(element("title",readTitle+" op:"+op.name()+" "+StringUtils.niceInt(ref1)+":"+StringUtils.niceInt(ref1+ce.getLength())));
												
												// arrow <--
												if(cigarIdx==0 && rec.getReadNegativeStrandFlag()) {
													sb.append( "M ").append(format(x1)).append(',').append(0);
													sb.append(" h ").append(format(distance_pix));
													sb.append(" v ").append(format(maxy));
													sb.append(" h ").append(format(-(distance_pix)));
													sb.append(" l ").append(format(-arrow_w)).append(',').append(-featureHeight/2.0);
													sb.append(" Z");
													}
												// arrow -->
												else if(cigarIdx+1==cigar.numCigarElements() && !rec.getReadNegativeStrandFlag()) {
													sb.append( "M ").append(format(x2)).append(',').append(0);
													sb.append(" h ").append(format(-(distance_pix)));
													sb.append(" v ").append(format(maxy));
													sb.append(" h ").append(format(distance_pix));
													sb.append(" l ").append(format(arrow_w)).append(',').append(-featureHeight/2.0);
													sb.append(" Z");
													
													}
												else
													{
													sb.append( "M ").append(format(x1)).append(',').append(0);
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
								} //end window
							y_rec += this.featureHeight;
							y_rec += 1;
							y_read_name += this.featureHeight;
							y_read_name += 1;

							}
					readRect.setAttribute("height", format(y_read_name - y1));
					y1 = y_read_name;
					}// end loop over each readName
					
					
				
				for(Window win:windows) {
					win.height = y1 - win.y;
					// frame
					final Element frame = element("rect");
					frame.setAttribute("class", "frame");
					frame.setAttribute("x", format(win.x));
					frame.setAttribute("y",format(win.y));
					frame.setAttribute("width",format(win.width));
					frame.setAttribute("height",format(win.height));
					mainG.appendChild(frame);
					
					for(int side=0;side<2;++side) {
						int pos = side==0?interval.getStart():interval.getEnd();
						double x = win.baseToPixel(pos);
						if(x< win.getMinX() || x> win.getMaxX()) continue;
						final Element line = element("line");
						line.setAttribute("class", "breakpoint");
						line.setAttribute("x1", format(x));
						line.setAttribute("x2", format(x));
						line.setAttribute("y1",format(win.y));
						line.setAttribute("y2",format(win.y+win.height));
						mainG.appendChild(line);
						}
					
					} //end window
				y1+=10;
				if(!segments.isEmpty()) {
					final int x0 = segments.stream().mapToInt(R->R.getStart()-10).min().orElse(0);
					final int x1 = segments.stream().mapToInt(R->R.getEnd()+10).max().orElse(0);
					final IntToDoubleFunction base2pos = POS-> ((POS-x0)/(double)(x1-x0))*this.drawinAreaWidth;
					final Pileup<Segment> pileup = new Pileup<>();
					pileup.addAll(segments);
					final double y0=y1;
					final Element g = element("g");
					mainG.appendChild(g);
					for(List<Segment> row:pileup.getRows()) {
						for(Segment seg:row) {
							final Element g2 = element("g");
							g.appendChild(g2);
							final Element path = element("path");
							path.appendChild(element("title",seg.readName));
							StringBuilder sb = new StringBuilder();
							sb.append("M" ).append(base2pos.applyAsDouble(seg.getStart())).append(",").append(y1+this.featureHeight/2.0);
							sb.append(" L " ).append(base2pos.applyAsDouble(seg.getEnd())).append(",").append(y1+this.featureHeight/2.0);
							sb.append(" Z ");
							
							sb.append(" M" ).append(base2pos.applyAsDouble(seg.startA)).append(",").append(y1);
							sb.append(" L " ).append(base2pos.applyAsDouble(seg.endA)).append(",").append(y1);
							sb.append(" L " ).append(base2pos.applyAsDouble(seg.endA)).append(",").append(y1+this.featureHeight);
							sb.append(" L " ).append(base2pos.applyAsDouble(seg.startA)).append(",").append(y1+this.featureHeight);
							sb.append(" Z ");
							
							sb.append(" M" ).append(base2pos.applyAsDouble(seg.startB)).append(",").append(y1);
							sb.append(" L " ).append(base2pos.applyAsDouble(seg.endB)).append(",").append(y1);
							sb.append(" L " ).append(base2pos.applyAsDouble(seg.endB)).append(",").append(y1+this.featureHeight);
							sb.append(" L " ).append(base2pos.applyAsDouble(seg.startB)).append(",").append(y1+this.featureHeight);
							sb.append(" Z ");

							
							path.setAttribute("d", sb.toString());
							g2.appendChild(path);
							}
						y1+= this.featureHeight;
						}
					
					for(int side=0;side<2;++side) {
						double x= base2pos.applyAsDouble(side==0?interval.getStart():interval.getEnd());
						Element line=element("line");
						line.setAttribute("class", "breakpoint");
						line.setAttribute("x1",format(x));
						line.setAttribute("x2",format(x));
						line.setAttribute("y1",format(y0));
						line.setAttribute("y2",format(y1));
						g.appendChild(line);
						}
					final Element frame = element("rect");
					frame.setAttribute("class", "frame");
					frame.setAttribute("x", "0");
					frame.setAttribute("y0",format(y0));
					frame.setAttribute("width", format(this.drawinAreaWidth-1));
					frame.setAttribute("height", format(y1-y0-1));
					g.appendChild(frame);
					y1+=10;
					}
				
				
				svgRoot.setAttribute("height", format(y1+10));
				
				
				
				} // end samReader
				
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

