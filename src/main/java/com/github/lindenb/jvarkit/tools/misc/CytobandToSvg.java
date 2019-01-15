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
package com.github.lindenb.jvarkit.tools.misc;

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.XLINK;
import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
/**
BEGIN_DOC

## Input

input is a tab delimited file similar to http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz

no header

columns are

* chrom
* chromStart
* chromEnd
* name
* gieStain

Chromosomes will be ordered occording to the first time they're seen in the file, so you'd better sort them.

```
$ wget -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" |\
	gunzip -c |\
	sort -t $'\t' -k1,1V -k2,2n > cytoBand.txt
```

END_DOC
 */
@Program(name="cytoband2svg",
	description="Creates a svg karyotype .",
	keywords={"karyotype","svg","ideogram"}
	)
public class CytobandToSvg extends Launcher {
private static final Logger LOG = Logger.build(CytobandToSvg.class).make();


@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private File outputFile = null;
@Parameter(names={"-C","--cytobands"},description="Cytoband URI. tab delimited file. no header. chrom/chromStart/chromEnd/name/gieStain. '-'=stdin")
private String cytobandsUri = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz";
@Parameter(names={"-t","--title"},description="title")
private String title = "";
@Parameter(names={"-W","-width","--width"},description="Image width")
private int width = 1000;
@Parameter(names={"-H","--height"},description="Image height")
private int height = 700;

private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
private final int contig_title_size = 12;
private final int round_rect = 5;

private class Cytoband
	implements Locatable
	{
	private final String contig;
	private final int chromStart;
	private final int chromEnd;
	private final String name;
	private final String gieStain;
	public Cytoband(
			final String contig,
			final int chromStart,
			final int chromEnd,
			final String name,
			final String gieStain
			) {
		this.contig = contig;
		this.chromStart = chromStart;
		this.chromEnd = chromEnd;
		this.name = name;
		this.gieStain = gieStain;
		}
	@Override
	public String getContig() { return this.contig;}
	@Override
	public int getStart() { return this.chromStart;}
	@Override
	public int getEnd() { return this.chromEnd;}
	public String getName() {
		return name;
		}
	public String getStain() {
		return gieStain;
		}
	/* https://github.com/ENCODE-DCC/kentUtils/blob/master/src/hg/lib/hCytoBand.c#L43 */
	public String getFillColor() {
		final String stain = this.getStain();
		if (stain.startsWith("gneg"))
		    {
		    return "lightblue";
		    }
		else if (stain.startsWith("gpos") && stain.length()>4)
		    {
		    int percentage;
		    try {
		    	percentage = Math.max(0, Math.min(100,Integer.parseInt(stain.substring(4))));
		    	}
		    catch(final NumberFormatException err) {
		    	percentage = 100;
		    	}
		    final int g = 40 + (int)(215.0*(percentage/100.0));
		    return "rgb("+g+","+g+","+g+")";
		    }
		else if (stain.startsWith("gvar"))
		    {
		    return "slategray";
		    }
		else 
		    {
		    return "honeydew";
		    }
		}
	}

private class Contig
	{
	final Rectangle2D.Double bounds = new Rectangle2D.Double();
	private Rectangle2D.Double _contig_bounds = null;
	final String name;
	final List<Cytoband> cytobands;
	final int length ;
	int telomerePos =0;
	Contig(final List<Cytoband> cytobands) {
		this.name = cytobands.get(0).getContig();
		this.cytobands = cytobands;
		this.length = cytobands.stream().mapToInt(C->C.getEnd()).max().getAsInt();
		for(int i=0;i+1<cytobands.size();i++)
			{
			if(cytobands.get(i).getName().contains("p") &&
			   cytobands.get(i+1).getName().contains("q"))
				{
				this.telomerePos = cytobands.get(i).getEnd();
				}
			}
		}
	public String getContig() {
		return this.name;
		}
	int getSequenceLength() {
		return this.length;
		}
	Rectangle2D.Double getContigBounds() {
		if(this._contig_bounds==null) {
			double top = (contig_title_size)+(StringUtil.isBlank(title)?0:12);
			final double w = this.bounds.getWidth()/3.0;
			this._contig_bounds =  new Rectangle2D.Double(
				this.bounds.getX()+w,
				this.bounds.getY()+top,
				w,
				(this.bounds.getHeight()-top)*(this.getSequenceLength()/(double)CytobandToSvg.this.max_contig_length)
				);
			}
		return this._contig_bounds;
		}
	
	double pos2pixel(int pos0, Rectangle2D.Double r) {
		final double y = r.getY() + (pos0/(double)this.getSequenceLength())*r.getHeight();
		return y;
		}
	
	
	void def(final XMLStreamWriter w) throws XMLStreamException {
		w.writeStartElement("clipPath");
		w.writeAttribute("id", "clip_"+getContig());
		if(this.telomerePos>0) {
			w.writeStartElement("rect");
			writecontigPAttrs(w);
			w.writeEndElement();
			}
		w.writeStartElement("rect");
			writecontigQAttrs(w);
		w.writeEndElement();//rect
		
		w.writeEndElement();//clipPath
		}
	void paint(final XMLStreamWriter w) throws XMLStreamException {
		final Rectangle2D.Double r= getContigBounds();
		w.writeStartElement("g");
		w.writeAttribute("id", getContig());
		
		w.writeStartElement("text");
		w.writeAttribute("class","ctglabel");
		w.writeAttribute("x",format(this.bounds.getCenterX()));
		w.writeAttribute("y",format(this.bounds.getY()+contig_title_size - 2 + (StringUtil.isBlank(title)?0:12)));
		w.writeCharacters(this.getContig());//title
		w.writeEndElement();
		
		w.writeStartElement("g");
		
		
		/* background shape */
		w.writeStartElement("g");
		w.writeAttribute("transform","translate(5,5) ");
		if(this.telomerePos>0) {
			w.writeStartElement("rect");
			w.writeAttribute("class","ctgback");
			w.writeAttribute("filter","url(#filter01)");
			writecontigPAttrs(w);
			w.writeEndElement();
			}
		
		w.writeStartElement("rect");
		w.writeAttribute("class","ctgback");
		w.writeAttribute("filter","url(#filter01)");
			writecontigQAttrs(w);
		w.writeEndElement();//rect
		w.writeEndElement();//g
		
		
		for(final  Cytoband c:this.cytobands) {
			w.writeStartElement("rect");
			w.writeAttribute("class","cytoband");
			w.writeAttribute("style","fill:"+c.getFillColor()+";");
			w.writeAttribute("clip-path","url(#clip_"+getContig()+")");
			w.writeAttribute("x",format(r.getX()));
			w.writeAttribute("y",format(pos2pixel(c.getStart(),r)));
			w.writeAttribute("width",format(r.getWidth()));
			w.writeAttribute("height",format((pos2pixel(c.getEnd(),r)-pos2pixel(c.getStart(),r))));
			
			if(!StringUtil.isBlank(c.getName())) {
				w.writeStartElement("title");
				w.writeCharacters(c.getName());
				w.writeEndElement();
				}
			
			w.writeEndElement();
			}
		
		if(this.telomerePos>0) {
			w.writeStartElement("rect");
			w.writeAttribute("class","ctgborderp");
			writecontigPAttrs(w);
			w.writeEndElement();
			}
		
		
		w.writeStartElement("rect");
			w.writeAttribute("class","ctgborderq");
			writecontigQAttrs(w);
		w.writeEndElement();//rect
		
		w.writeEndElement();//g
		
		w.writeEndElement();//g
		}
	
	private void writecontigPAttrs(final XMLStreamWriter w) throws XMLStreamException {
		final Rectangle2D.Double r= getContigBounds();
		w.writeAttribute("x",format(r.getX()));
		w.writeAttribute("y",format(pos2pixel(0,r)));
		w.writeAttribute("width",format(r.getWidth()));
		w.writeAttribute("height",format(pos2pixel(this.telomerePos,r)-pos2pixel(0,r)));
		w.writeAttribute("rx",format(round_rect));
		w.writeAttribute("ry",format(round_rect));
		}
	private void writecontigQAttrs(final XMLStreamWriter w)  throws XMLStreamException {
		final Rectangle2D.Double r= getContigBounds();
		w.writeAttribute("x",format(r.getX()));
		w.writeAttribute("y",format(pos2pixel(this.telomerePos,r)));
		w.writeAttribute("width",format(r.getWidth()));
		w.writeAttribute("height",format(pos2pixel(this.getSequenceLength(),r) - pos2pixel(this.telomerePos,r)));
		w.writeAttribute("rx",format(round_rect));
		w.writeAttribute("ry",format(round_rect));
		}

	}

private final Map<String,Contig> name2contig = new LinkedHashMap<>();
private long max_contig_length = 0L;



private List<Cytoband>  parseCytobands(final BufferedReader r) throws IOException
	{
	final Pattern tab=Pattern.compile("[\t]");
	return r.lines().
			map(L->tab.split(L)).
			map(T->new Cytoband(
					T[0],
					Integer.parseInt(T[1]),
					Integer.parseInt(T[2]),
					T[3],
					T[4]
					)).
			collect(Collectors.toList())
			;
	}

/** convert double to string */
private String format(final double v)
	{
	return this.decimalFormater.format(v);
	}


@Override
public int doWork(final List<String> args) {
	XMLStreamWriter w=null;
	FileOutputStream fout=null;

	try {
		final List<Cytoband> cytobands;
		if(StringUtil.isBlank(this.cytobandsUri) || this.cytobandsUri.equals("-")) {
			cytobands = parseCytobands(IOUtils.openStdinForBufferedReader());
			}
		else
			{
			cytobands = parseCytobands(IOUtils.openURIForBufferedReading(this.cytobandsUri));
			}
		
		cytobands.stream().
			map(C->C.getContig()).
			collect(Collectors.toCollection(LinkedHashSet::new)).
			stream().
			map(C->new Contig(cytobands.stream().filter(T->T.getContig().equals(C)).sorted((A,B)->Integer.compare(A.getStart(),B.getStart())).collect(Collectors.toList()))).
			forEach(C->this.name2contig.put(C.getContig(), C));
			;
		
		this.max_contig_length = this.name2contig.values().stream().mapToLong(C->C.getSequenceLength()).max().orElse(0L);	
			
		if(this.name2contig.isEmpty() || this.max_contig_length<=0)
			{
			LOG.error("no contig found");
			return -1;
			}
		
		final Rectangle drawingArea = new Rectangle(0, 0, this.width, this.height);
		
		
		final XMLOutputFactory xof=XMLOutputFactory.newFactory();
		if(this.outputFile==null)
			{
			w=xof.createXMLStreamWriter(stdout(), "UTF-8");
			}
		else
			{
			fout = new FileOutputStream(this.outputFile);
			w=xof.createXMLStreamWriter(fout, "UTF-8");
			}
		w.writeStartDocument("UTF-8", "1.0");
		
		w.writeStartElement("svg");
		w.writeAttribute("width", format(drawingArea.width));
		w.writeAttribute("height", format(drawingArea.height));
		w.writeAttribute("version","1.1");
		w.writeDefaultNamespace(SVG.NS);
		w.writeNamespace("xlink", XLINK.NS);
		
		
		w.writeStartElement("title");
		w.writeCharacters(StringUtil.isBlank(this.title)?CytobandToSvg.class.getName():this.title);
		w.writeEndElement();
		
		w.writeStartElement("description");
		w.writeCharacters("Cmd:"+getProgramCommandLine()+"\n");
		w.writeCharacters("Version:"+getVersion()+"\n");
		w.writeCharacters("Author: Pierre Lindenbaum\n");
		w.writeEndElement();

		
		w.writeStartElement("defs");
		
		w.writeStartElement("linearGradient");
		w.writeAttribute("id","grad01");
		w.writeAttribute("x1","0%");
		w.writeAttribute("x2","100%");
		w.writeAttribute("y1","50%");
		w.writeAttribute("y2","50%");
		w.writeEmptyElement("stop");
			w.writeAttribute("offset","0%");
			w.writeAttribute("style","stop-color:white;stop-opacity:1;");
		w.writeEmptyElement("stop");
			w.writeAttribute("offset","50%");
			w.writeAttribute("style","stop-color:white;stop-opacity:0;");
		w.writeEmptyElement("stop");
			w.writeAttribute("offset","75%");
			w.writeAttribute("style","stop-color:black;stop-opacity:0;");
		w.writeEmptyElement("stop");
			w.writeAttribute("offset","100%");
			w.writeAttribute("style","stop-color:black;stop-opacity:0.9;");
		
		w.writeEndElement();
		
		w.writeStartElement("filter");
			w.writeAttribute("id","filter01");
			w.writeEmptyElement("feGaussianBlur");
				w.writeAttribute("in","SourceGraphic");
				w.writeAttribute("stdDeviation","10");
			
		w.writeEndElement();//filter
		
		
		
		double x = drawingArea.getX();
		double contig_width = drawingArea.getWidth()/(double)this.name2contig.size();
		for(final Contig contig:this.name2contig.values()) {
			contig.bounds.x = x;
			contig.bounds.y = drawingArea.getY();
			contig.bounds.width = contig_width;
			contig.bounds.height = drawingArea.getHeight();
			contig.def(w);
			x+= contig_width;
		}
		w.writeEndElement();
		
		w.writeStartElement("style");
		w.writeCharacters(
				".ctgborderp {fill:url(#grad01);stroke:green;}" +
				".ctgborderq {fill:url(#grad01);stroke:green;}" +
				".ctglabel {text-anchor:middle;stroke:none;fill:darkgrey;font: bold 10px Verdana, Helvetica, Arial, sans-serif;}" +
				".cytoband {fill:silver;stroke:none;}" +
				".bedlabel {stroke:red;fill:none;text-anchor:start;font: normal 7px Verdana, Helvetica, Arial, sans-serif;}" +
				".maintitle {stroke:none;fill:darkgrey;text-anchor:middle;font: normal 12px Verdana, Helvetica, Arial, sans-serif;}" +
				".ctgback {fill:gainsboro;stroke:none;filter:url(#filter01);}"
				);
		
		w.writeEndElement();//style

		w.writeStartElement("g");
		w.writeAttribute("id", "body");
		
		
		/* title */
		if(!StringUtil.isBlank(this.title)) {
			w.writeStartElement("text");
			w.writeAttribute("class", "maintitle");
			w.writeAttribute("x", format(this.width/2.0));
			w.writeAttribute("y", format(12));
			w.writeCharacters(this.title);
			w.writeEndElement();
			}
		
		
		w.writeStartElement("g");
		w.writeAttribute("id", "karyotype");
		
		for(final Contig contig:this.name2contig.values()) {
			contig.paint(w);
			}
		
		w.writeEndElement();//g@id=karyotype
		
		
		final ContigNameConverter ctgNameConverter = ContigNameConverter.fromContigSet(this.name2contig.keySet());
		for(final String filename: args) {
			final Pattern tab = Pattern.compile("[\t]");
			w.writeStartElement("g");
			w.writeAttribute("id", "input"+filename);
			BufferedReader br = IOUtils.openURIForBufferedReading(filename);
			String line;
			while((line=br.readLine())!=null)
				{
				if(line==null || StringUtil.isBlank(line) || BedLine.isBedHeader(line)) continue;
				final String tokens[] = tab.split(line,4);
				if(tokens.length<3) continue;
				final String ctgName = ctgNameConverter.apply(tokens[0]);
				if(StringUtil.isBlank(ctgName)) continue;
				final Contig ctg = this.name2contig.get(ctgName);
				if(ctg==null) continue;
				final Rectangle2D.Double r = ctg.getContigBounds();
				final int chromStart = Integer.parseInt(tokens[1]);
				final int chromEnd = Integer.parseInt(tokens[2]);
				
				final Map<String,String> attributes = new HashMap<>();
				if(tokens.length>3)
					{
					final StreamTokenizer st=new StreamTokenizer(new StringReader(tokens[3]));
					st.wordChars('_', '_');
					String key=null;
					while(st.nextToken() != StreamTokenizer.TT_EOF)
						{
						String s=null;
						switch(st.ttype)
							{
							case StreamTokenizer.TT_NUMBER: s=String.valueOf(st.nval);break;
							case '"': case '\'' : case StreamTokenizer.TT_WORD: s=st.sval;break;
							case ';':break;
							default:break;
							}
						if(s==null) continue;
						if(key==null)
							{
							key=s;
							}
						else
							{
							attributes.put(key,s);
							key=null;
							}
						}
					}
				w.writeStartElement("g");
				if(attributes.containsKey("href"))
					{
					w.writeStartElement("a");
					w.writeAttribute("href", attributes.get("href"));
					}
				w.writeStartElement("rect");
				w.writeAttribute("style","stroke:red;");
				w.writeAttribute("x",format(r.getMaxX()+2));
				w.writeAttribute("width",format(2));
				w.writeAttribute("y",format(ctg.pos2pixel(chromStart,r)));
				w.writeAttribute("height",format(ctg.pos2pixel(chromEnd,r) - ctg.pos2pixel(chromStart,r)));
				if(attributes.containsKey("title"))
					{
					w.writeStartElement("title");
					w.writeCharacters(attributes.get("title"));
					w.writeEndElement();
					}
				w.writeEndElement();
				if(attributes.containsKey("label"))
					{
					w.writeStartElement("text");
					w.writeAttribute("class","bedlabel");
					w.writeAttribute("x",format(r.getMaxX()+4));
					w.writeAttribute("y",format(3 + ctg.pos2pixel(chromStart,r)));
					w.writeCharacters(attributes.get("label"));
					w.writeEndElement();
					}
				
				
				if(attributes.containsKey("href"))
					{
					w.writeEndElement();
					}
				w.writeEndElement();//g
				}
			
			
			br.close();
			
			w.writeEndElement();
		}
		
		w.writeEndElement();//g@id=body
		
		w.writeEndElement();//svg
		w.writeEndDocument();
		w.flush();
		w.close();
		if(fout!=null) {
			fout.flush();
			fout.close();
			fout=null;
			}
		return 0;
		}
	catch (final Exception e) {
		LOG.error(e);
		return -1;
		} 
	finally {
		CloserUtil.close(fout);
		}
	}
public static void main(final String[] args) {
	new CytobandToSvg().instanceMain(args);
	}
}
