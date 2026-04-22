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
package com.github.lindenb.jvarkit.tools.bed2hilbert;

import java.awt.Color;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.RectangularShape;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.hilbert.HilbertCurve;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.DoubleRounder;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
BEGIN_DOC

## Motivation

Creates a Hilbert SVG Graph  for a BED file.


## Example

```
 java -jar dist/jvarkit.jar bed2hilbert  -R src/test/resources/rotavirus_rf.dict input.bed > out.svg
```

END_DOC
*/
@Program(
		name="bed2hilbert",
		description="BED to hilbert curve as SVG.",
		keywords={"bed","xml","svg","hilbert"},
		creationDate = "20260417",
		modificationDate = "20260417",
		jvarkit_amalgamion = true
		)
public class BedToHilbert extends Launcher {
	private static final Logger LOG=Logger.of(BedToHilbert.class);
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;

	@Parameter(names={"--regex"},description="keep chromosomes matching that regular expression")
	private String contig_regex=null;
	@Parameter(names={"--columns"},description="Optional Comma separated list names of column starting from the 3rd column (after 'end'). Use '.' to ignore.")
	private String column_names_str="";
	@Parameter(names={"--min-contig-length"},description="keep chromosomes which length is greater than 'x'")
	private int min_contig_length=0;
	@Parameter(names={"-R","--reference"},description=DICTIONARY_SOURCE,required = true)
	private Path faidPath=null;
    /** with/height of the final picture */
    @Parameter(names={"-w","--width"},description="Image width")
    private int imageWidth=1000;
    /** level of recursion */
    @Parameter(names={"-r","--recursion"},description="Hilbert Curve level of recursion")
    private int recursionLevel=6;
	@Parameter(names={"--omit-xml-desclaration"},description="Don't print XM declaration.")
	private boolean omit_xml_decl = false;
	@Parameter(names={"--first-line-header"},description="The first line of the bed file is a header")
	private boolean first_line_is_header = false;
	@Parameter(names={"--css"},description="load custom CSS style file")
	private Path cssPath=null;

    private  final DoubleRounder dblRounder= new DoubleRounder(2);
	
	
    private static long toIndex(final SAMSequenceDictionary dict,final String contig,final int pos) {
    	long n = 0L;
    	for(SAMSequenceRecord ssr: dict.getSequences()) {
    		if(ssr.getContig().equals(contig)) {
    			return n+pos;
    			}
    		n+= ssr.getLengthOnReference();
    		}
    	throw new IllegalArgumentException(contig+":"+pos);
    	}
    
    private String format(double v) {
    	return dblRounder.format(v);
    	}
    
    
    private String toString(final List<Point2D.Double> points) {
    	final StringBuilder sb = new StringBuilder();
		for(int i=0;i< points.size();i++) {
			if(i>0) sb.append(" ");
			sb.append(format(points.get(i).getX()));
			sb.append(",");
			sb.append(format(points.get(i).getY()));
			}
		return sb.toString();
    	}
    private RectangularShape toShape(final List<Point2D.Double> points) {
    	final double x = points.stream().mapToDouble(P->P.getX()).min().orElse(0.0);
    	final double y = points.stream().mapToDouble(P->P.getY()).min().orElse(0.0);
    	final double w = points.stream().mapToDouble(P->P.getX()-x).max().orElse(1.0);
    	final double h = points.stream().mapToDouble(P->P.getY()-y).max().orElse(1.0);
    	return new Rectangle2D.Double(x, y, w, h);
    	}
    private Point2D.Double toCenter(final List<Point2D.Double> points) {
    	final RectangularShape r = toShape(points);
    	return new Point2D.Double(r.getCenterX(),r.getCenterY());
    	}
    
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = oneFileOrNull(args);
			final Pattern contigRegex = (this.contig_regex==null?null:Pattern.compile(this.contig_regex));
			final SAMSequenceDictionary dict0 = new SequenceDictionaryExtractor().extractRequiredDictionary(faidPath);
			
			final SAMSequenceDictionary dict = new SAMSequenceDictionary(
					dict0.getSequences().stream().
					filter(SSR->min_contig_length<=0 || SSR.getSequenceLength()>=min_contig_length).
					filter(SSR->contigRegex==null?true:contigRegex.matcher(SSR.getSequenceName()).matches()).
					collect(Collectors.toList())
					);
			if(dict.isEmpty()) {
				LOG.warn("empty dictionary");
				}
			final Hyperlink hyperlink = Hyperlink.compile(dict0);
			final long genomeLength= dict.getReferenceLength();
			if(genomeLength<=0L) {
				LOG.error("no empty dict");
				return -1;
				}
			
			
			
			final HilbertCurve hilbertCurve = new HilbertCurve(this.imageWidth,genomeLength, this.recursionLevel);
			try(BufferedReader br = (input==null?IOUtils.openStreamForBufferedReader(stdin()): IOUtils.openURIForBufferedReading(input))) {
				final ContigNameConverter ctgConvert = ContigNameConverter.fromOneDictionary(dict);
				String line;
				final List<String> header_tokens;
				if(first_line_is_header) {
					line = br.readLine();
					if(line==null) throw new IOException("Cannot read first line of input.");
					header_tokens  = new ArrayList<>(CharSplitter.TAB.splitAsStringList(line));
					}
				else
					{
					header_tokens  = new ArrayList<>(Arrays.asList("chromosome","start","end"));
					}
				if(header_tokens.size()<3) {
					throw new JvarkitException.TokenErrors(3, header_tokens.toArray(new String[header_tokens.size()]));
					}
				if(!StringUtils.isBlank(this.column_names_str)) {
					String[] tokens = CharSplitter.COMMA.split(this.column_names_str);
					for(int i=0;i< tokens.length;i++) {
						if(i+3< header_tokens.size()) {
							header_tokens.set(i+3, tokens[i]);
							}
						else
							{
							header_tokens.add(tokens[i]);
							}
						}
					}
				FileHeader fileHeader= null;
				

				
				try(OutputStream w0 =( this.outputFile==null?stdout():IOUtils.openPathForWriting(this.outputFile))) {
					final XMLOutputFactory xof = XMLOutputFactory.newFactory();
					final XMLStreamWriter w = xof.createXMLStreamWriter(w0, "UTF-8");
					if(!omit_xml_decl) w.writeStartDocument("UTF-8", "1.0");
					w.writeStartElement("svg");
					w.writeDefaultNamespace(SVG.NS);
					w.writeAttribute("width", String.valueOf(this.imageWidth+1));
					w.writeAttribute("height", String.valueOf(this.imageWidth+1));
					
					w.writeStartElement("title");
					w.writeCharacters(input==null?"bed2hilbert":input);
					w.writeEndElement();

					
					w.writeStartElement("style");
					w.writeCharacters(
						 ".bckg {stroke:darkgray;fill:whitesmoke;}\n"
						+".ka {fill:none;;stroke-dasharray:1,1;stroke-width:2px;}\n"
						+".k1 {stroke:darkslateblue;fill:none;}\n"
						+".kX {stroke:blue;fill:none;}\n"
						+".kY {stroke:pink;fill:none;}\n"
						+".rec {stroke:yellow;fill:none;stroke-width:5;opacity:0.5;}\n"
						+".rec:hover,.rec:focus {stroke:darkgreen;stroke-width:4px;}\n"
						+".ctgLabel {stroke:none;opacity:0.9;text-anchor:middle;font-size:10px;}\n"
						+"circle.edge {fill:black;stroke:black;opacity:0.8;}\n"
						);
					w.writeEndElement();
					if(this.cssPath!=null) {
						w.writeStartElement("style");
						w.writeCharacters(IOUtils.slurpPath(cssPath));
						w.writeEndElement();
						}
					
					w.writeStartElement("defs");
					w.writeStartElement("g");
					w.writeAttribute("id", "genome");
					
					
					
					for(SAMSequenceRecord ssr: dict.getSequences()) {
						final Color c = Color.getHSBColor((float) ssr.getSequenceIndex() / (float)dict.size(), 0.7f, 0.5f);
						final long index = toIndex(dict,ssr.getContig(),0); 
						final List<Point2D.Double> points = hilbertCurve.getPoints( index, index+ssr.getLengthOnReference());
						if(points.isEmpty()) continue;
						w.writeStartElement("g");
						w.writeStartElement("polyline"); //path
						String stroke=ColorUtils.toRGB(c)+";";
						String className =  "k"+(ssr.getSequenceIndex()%2);
						if(ssr.getSequenceName().matches("(chr)?X")) {
							stroke = "blue;";
							className = "kX";
							}
						else if(ssr.getSequenceName().matches("(chr)?Y")) {
							stroke = "pink;";
							className = "kY";
							}
						
						w.writeAttribute("class","ka "+className);
						if(!StringUtils.isBlank(stroke)) {
							w.writeAttribute("style", "stroke:"+stroke);
							}
						w.writeAttribute("points", toString(points));
						
						w.writeStartElement("title");
						w.writeCharacters(ssr.getContig()+":"+StringUtils.niceInt(ssr.getSequenceLength())+" bp");
						w.writeEndElement();
						
						w.writeEndElement(); // path
						
						
						Point2D.Double center = toCenter(points);
						w.writeStartElement("text");
						w.writeAttribute("class", "ctgLabel");
						w.writeAttribute("style", "fill:"+stroke);
						w.writeAttribute("x", format(center.getX()));
						w.writeAttribute("y", format(center.getY()));
						w.writeCharacters(ssr.getContig());
						w.writeEndElement();
						
						if(ssr.getSequenceIndex()+1 < dict.size()) {
							w.writeEmptyElement("circle");
							w.writeAttribute("class", "edge");
							w.writeAttribute("r", "3");
							w.writeAttribute("cx", format(points.get(points.size()-1).getX()));
							w.writeAttribute("cy", format(points.get(points.size()-1).getY()));
							}
						
						
						w.writeEndElement(); // g
						}
					
					w.writeEndElement();//g
					w.writeEndElement();//defs
					
					
					
					w.writeStartElement("g");
					w.writeAttribute("id", "bed_records");
					
					
					w.writeEmptyElement("rect");
					w.writeAttribute("class", "bckg");
					w.writeAttribute("x", String.valueOf(0));
					w.writeAttribute("y", String.valueOf(0));
					w.writeAttribute("width", String.valueOf(this.imageWidth));
					w.writeAttribute("height", String.valueOf(this.imageWidth));
					
					w.writeEmptyElement("use");
					w.writeAttribute("x","0");
					w.writeAttribute("y","0");
					w.writeAttribute("href","#genome");
					
					
					while((line=br.readLine())!=null) {
						if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
						final String[] tokens = CharSplitter.TAB.split(line);
						if(fileHeader==null) {
							while(header_tokens.size()< tokens.length) {
								header_tokens.add("$"+(header_tokens.size()+1));
								}
							fileHeader = new FileHeader(header_tokens, CharSplitter.TAB);
							}
						
						
						final FileHeader.RowMap row= fileHeader.toMap(line);
						final String  ctg  = ctgConvert.apply(row.at(0));
						if(StringUtils.isBlank(ctg)) continue;
						final SAMSequenceRecord ssr = dict.getSequence(ctg);
						final int chromStart0  = Integer.parseInt(row.at(1));
						if(chromStart0<0) throw new IllegalArgumentException("negative start in "+row);
						if(chromStart0> ssr.getSequenceLength()) continue;
						final int chromEnd  = Math.min(ssr.getLengthOnReference(), Integer.parseInt(row.at(2)));
						if(chromEnd<chromStart0) throw new IllegalArgumentException("end < start in "+row);
						if(chromStart0== chromEnd) continue;

						long index = toIndex(dict, ctg, chromStart0); 
						final List<Point2D.Double> points = hilbertCurve.getPoints( index, index+(chromEnd-chromStart0));
						
						
						
						if(points.isEmpty()) continue;
						
						
						
						w.writeStartElement("g");
						
						String href = row.getOrDefault("href", hyperlink.apply(new SimpleInterval(ctg,chromStart0+1,chromEnd)).orElse(""));
						String title = row.getOrDefault("title", ctg+":"+StringUtils.niceInt(chromStart0+1)+"-"+StringUtils.niceInt(chromEnd)+" len:"+StringUtils.niceInt(chromEnd-chromStart0));
						if(!StringUtils.isBlank(href)) {
							w.writeStartElement("a");
							w.writeAttribute("href", href);
							w.writeAttribute("target", "_blank");
							}
						
						w.writeStartElement("polyline"); //path
						
						w.writeAttribute("class", row.getOrDefault("class","rec" ));
							
						if(!StringUtils.isBlank( row.getOrDefault("style","" ))) {
							w.writeAttribute("style", row.getOrDefault("style","" ));
							}
						
						w.writeAttribute("points", toString(points));
						
						if(!StringUtils.isBlank(title)) {
							w.writeStartElement("title");
							w.writeCharacters(title);
							w.writeEndElement();
							}

						
						w.writeEndElement(); // path
						
						final String label = row.getOrDefault("label","");
						if(!StringUtils.isBlank(label)) {
							double wd = Math.min(toShape(points).getWidth()*0.8,label.length()*7)/label.length();
							Point2D.Double center = toCenter(points);
							w.writeStartElement("text");
							w.writeAttribute("style","stroke:none;text-anchor;middle;font-size:"+ wd);
							w.writeAttribute("x", format(center.getX()));
							w.writeAttribute("y", format(center.getY()));
							w.writeCharacters(label);
							w.writeEndElement();
							}
						
						if(!StringUtils.isBlank(href)) {
							w.writeEndElement();// anchor
							}
						w.writeEndElement();// g
						}
					/*
					w.writeEmptyElement("use");
					w.writeAttribute("x","0");
					w.writeAttribute("y","0");
					w.writeAttribute("href","#genome");
					*/
					w.writeEndElement();//g
					w.writeEndElement();
					if(!omit_xml_decl) w.writeEndDocument();
					w.flush();
					w.close();
					w0.flush();
					}
				}
			return 0;
			}
		catch(Throwable err ) {
			LOG.error(err);
			return -1;
			}
		finally {
			
			}
		
		}
	
	public static void main(String[] args) {
		new BedToHilbert().instanceMainWithExit(args);
	}
}
