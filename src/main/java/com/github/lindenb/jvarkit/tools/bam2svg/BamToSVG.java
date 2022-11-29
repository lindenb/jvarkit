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

import java.awt.Dimension;
import java.awt.Insets;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
//import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.RDF;
import com.github.lindenb.jvarkit.util.ns.XLINK;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;


/**
BEGIN_DOC

## Input

input is a set of indexed BAM/CRAM files or a file with the path to the bam with the '.list' suffix.

## Example

```bash

$ find dir dir2 -type f -name "*.bam" > file.list

$ java -jar dist/bam2svg.jar \
    -R human_g1k_v37.fasta \
    -i "19:252-260" \
    -S variants.vcf.gz \
    -o output.zip
    file.list
```

## Gallery

https://twitter.com/yokofakun/status/523031098541232128

![bam2svg](https://pbs.twimg.com/media/B0IuAw2IgAAYfNM.jpg)

https://twitter.com/yokofakun/status/522415314425090048

![bam2svg-2](https://pbs.twimg.com/media/Bz_99ayIMAAK57s.jpg)



END_DOC
 */
@Program(name="bam2svg",
description="BAM to Scalar Vector Graphics (SVG)",
keywords={"bam","alignment","graphics","visualization","svg"},
creationDate="20141013",
modificationDate="20210728"
)
@SuppressWarnings("fallthrough")
public class BamToSVG extends Launcher {
	private static final int HEIGHT_MAIN_TITLE=30;
	private static final int HEIGHT_SAMPLE_NAME=20;
	
	
	private static final Logger LOG = Logger.build(BamToSVG.class).make();

	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private File outputFile = null;
	@Parameter(names={"-i","--interval","--region"},description=IntervalParserFactory.OPT_DESC,required=true)
	private String intervalStr = null;
	@Parameter(names={"-w","--width"},description="Page width")
	private int drawinAreaWidth = 1000 ;
	@Parameter(names={"-c","--showclipping","--clip"},description="Show clipping")
	private boolean showClipping = false;
	@Parameter(names={"-S","--vcf"},description="Indexed VCF. the Samples's name must be the same than in the BAM")
	private Path vcf = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFile = null;
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition samRecordPartition = SAMRecordPartition.sample;
	//@Parameter(names={"--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	//private SamRecordFilter samRecordFilter = SamRecordFilterFactory.getDefault();
	@Parameter(names={"--mapq"},description="min mapping quality")
	private int mappingQuality =  1;
	@Parameter(names={"--prefix"},description="file prefix")
	private String prefix="";
	@DynamicParameter(names = "-D", description = "custom parameters. '-Dkey=value'. Undocumented.")
	private Map<String, String> dynaParams = new HashMap<>();


	private final Hershey hershey=new Hershey();
	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private SAMSequenceDictionary referenceDict=null;
	private GenomicSequence genomicSequence=null;
	private SimpleInterval interval=null;
	
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
		
	private class Context
		{
		String sampleName;
		Path path;
		List<List<SAMRecord>> lines= null;//filled by pileup
		final List<VariantContext> pos2variant=new ArrayList<VariantContext>();
		int HEIGHT_RULER=200;
		double featureHeight =1;
		double featureWidth =1;
		
		public double getHeight()
			{
			double dim_height= HEIGHT_SAMPLE_NAME;
			dim_height+= this.featureHeight;//ref seq
			dim_height+= this.featureHeight;//consensus
			if(!pos2variant.isEmpty())dim_height+=  this.featureHeight;//variants
			dim_height+= this.lines.size()* this.featureHeight;//consensus
			return dim_height;
			}
		
		private void readVariantFile(final Path vcf)  throws IOException
			{
			try(VCFReader r= VCFReaderFactory.makeDefault().open(vcf, true)) {
				r.query(BamToSVG.this.interval).stream().
					filter(V->V.getGenotype(this.sampleName)!=null).
					filter(V->!(V.getGenotype(this.sampleName).isNoCall() || V.getGenotype(this.sampleName).isHomRef())).
					forEach(ctx->this.pos2variant.add(ctx));
				}
			}
		
		}
		
		
		
		/** scale a rect by ratio */
		private static Rectangle2D.Double scaleRect(final Rectangle2D.Double r,double ratio)
			{
			double w=r.getWidth()*ratio;
			double h=r.getHeight()*ratio;
			double x=r.getX()+(r.getWidth()-w)/2.0;
			double y=r.getY()+(r.getHeight()-h)/2.0;
			return new Rectangle2D.Double(x,y,w,h);
			}
		/** scale a rect by default ratio */
		private static Rectangle2D.Double scaleRect(final Rectangle2D.Double r)
			{
			return scaleRect(r,0.95);
			}
	
		
		/** convert double to string */
		private String format(final double v)
			{
			return this.decimalFormater.format(v);
			}
		
	
		/** trim position in this.interval */
		private int trim(final int pos0)
			{
			return Math.min(Math.max(pos0, interval.getStart()),interval.getEnd()+1);
			}
		
		private int left(final SAMRecord rec)
			{
			return (this.showClipping?rec.getUnclippedStart():rec.getAlignmentStart());
			}
	
		private int right(final SAMRecord rec)
			{
			return (this.showClipping?rec.getUnclippedEnd():rec.getAlignmentEnd());
			}
		
		private double baseToPixel(int pos)
			{
			return  ((pos - this.interval.getStart())/(double)(this.interval.length()))*(this.drawinAreaWidth);
			}
		
		private void writeTitle(final XMLStreamWriter w,String title) throws XMLStreamException {
			if(StringUtils.isBlank(title)) return;
			w.writeStartElement("title");
			w.writeCharacters(title);
			w.writeEndElement();
			}
		
		/** get read name that is displayed in the popup window */
		private String getReadNameForPopup(final SAMRecord rec) {
			final StringBuilder sb = new StringBuilder(rec.getReadName());
			
			sb.append(rec.getReadNegativeStrandFlag()?" -":" +"); 
			
			if(rec.getReadPairedFlag()) {
				if(rec.getFirstOfPairFlag()) sb.append(" 1");
			if(rec.getSecondOfPairFlag()) sb.append(" 2");
			sb.append("/2");
			sb.append(" f"+rec.getFlags());
			if(rec.getMateUnmappedFlag())  {
				sb.append(" mate:unmapped");
				}
			else
				{
				sb.append(" mate:"+rec.getMateReferenceName()+":"+rec.getMateAlignmentStart());
					}
				}
			if(rec.getReadFailsVendorQualityCheckFlag())  sb.append(" fails-quality");
			if(rec.getDuplicateReadFlag()) sb.append(" duplicate");
			if(rec.isSecondaryAlignment())  sb.append(" secondary");
			if(rec.getSupplementaryAlignmentFlag())  sb.append("supplementary");
			sb.append(" MAPQ:"+rec.getMappingQuality());
			return sb.toString();
			}
		
		private void printGradientDef(
				final XMLStreamWriter w,
				final String gradId,
				final String styleTop,//e.g: "stop-color:black;stop-opacity:1;"
				final String styleMid //e.g: "stop-color:white;stop-opacity:1;"
				) throws XMLStreamException
			{
			w.writeStartElement("linearGradient");
			w.writeAttribute("id",gradId);
			w.writeAttribute("x1","50%");
			w.writeAttribute("x2","50%");
			w.writeAttribute("y1","0%");
			w.writeAttribute("y2","100%");
			w.writeEmptyElement("stop");
				w.writeAttribute("offset","0%");
				w.writeAttribute("style",styleTop);
			w.writeEmptyElement("stop");
				w.writeAttribute("offset","50%");
				w.writeAttribute("style",styleMid);
			w.writeEmptyElement("stop");
				w.writeAttribute("offset","100%");
				w.writeAttribute("style",styleTop);
			w.writeEndElement();
			}
		
		private void printSample(final XMLStreamWriter w,double top_y,final Context context)
				throws XMLStreamException
			{
			w.writeComment("START SAMPLE: "+context.sampleName);
			w.writeStartElement("g");
			w.writeAttribute("transform", "translate(0,"+top_y+")");
			double y=0;
			
			/* write title */
			w.writeStartElement("path");
			w.writeAttribute("class","samplename");
			
			w.writeAttribute("d", this.hershey.svgPath(
					context.sampleName,
					scaleRect(new Rectangle2D.Double(
							0,
							y,
							this.drawinAreaWidth,
							HEIGHT_SAMPLE_NAME
							)))
					);
			writeTitle(w,context.sampleName);
			w.writeEndElement();//path
			y+=HEIGHT_SAMPLE_NAME;
	
			/* write REFERENCE */
			w.writeComment("REFERENCE");
			for(int pos=this.interval.getStart();
					context.featureWidth>5 && //ignore if too small
					pos<=this.interval.getEnd();++pos)
				{
				char c=(this.genomicSequence==null?'N':this.genomicSequence.charAt(pos-1));
				double x0  = baseToPixel(pos);
				w.writeEmptyElement("use");
				w.writeAttribute("x",format(x0));
				w.writeAttribute("y",format(y));
				w.writeAttribute("xlink", XLINK.NS, "href", "#b"+c);
				}
			y+= context.featureHeight;
			
			/* skip line for later consensus */
			double consensus_y=y;
			y+=context.featureHeight;
			final Map<Integer,Counter<Character>> pos2consensus=new HashMap<Integer,Counter<Character>>();
			
			/* print variants */
			if(!context.pos2variant.isEmpty())
				{
				for(VariantContext ctx:context.pos2variant)
					{
					Genotype g=ctx.getGenotype(context.sampleName);
					if(g==null || !g.isCalled() || g.isHomRef()) continue;
					if(ctx.getEnd()< this.interval.getStart()) continue;
					if(ctx.getStart()> this.interval.getEnd()) continue;
					double x0  = baseToPixel(ctx.getStart());
					w.writeStartElement("use");
					w.writeAttribute("x",format(x0));
					w.writeAttribute("y",format(y));
					
					if(g.isHomVar())
						{
						w.writeAttribute("xlink", XLINK.NS, "href", "#homvar");				
						}
					else if (g.isHet())
						{
						w.writeAttribute("xlink", XLINK.NS, "href", "#het");
						}
					writeTitle(w,ctx.getAltAlleleWithHighestAlleleCount().getDisplayString());
					w.writeEndElement();//use
					}
				y+=context.featureHeight;
				}
			
			/* print all lines */
			w.writeComment("Alignments");
			for(int nLine=0;nLine< context.lines.size();++nLine)
				{
				w.writeStartElement("g");
				w.writeAttribute("transform", "translate(0,"+format(y+nLine*context.featureHeight)+")");
				List<SAMRecord> line= context.lines.get(nLine);
				//loop over records on that line
				for(SAMRecord record: line)
					{
					printSamRecord(w, context,record,pos2consensus);
					}
				w.writeEndElement();//g
				}
			
			/* write consensus */
			w.writeComment("Consensus");
			for(int pos=this.interval.getStart();
					context.featureWidth>5 && //ignore if too small
					pos<=this.interval.getEnd();++pos)
				{
				Counter<Character> cons = pos2consensus.get(pos);
				if(cons==null) continue;
				char c=cons.getMostFrequent();
				
				double x0  = baseToPixel(pos);
				w.writeEmptyElement("use");
				w.writeAttribute("x",format(x0));
				w.writeAttribute("y",format(consensus_y));
				w.writeAttribute("xlink", XLINK.NS, "href", "#b"+c);
					
				}
			
			
			/* surrounding frame for that sample */
			w.writeEmptyElement("rect");
			w.writeAttribute("class","frame");
			w.writeAttribute("x",format(0));
			w.writeAttribute("y",format(0));
			w.writeAttribute("width",format(drawinAreaWidth));
			w.writeAttribute("height",format(context.getHeight()));
	
			w.writeEndElement();//g for sample
			w.writeComment("END SAMPLE: "+context.sampleName);
	
			}
		
		private void printDocument(
			final XMLStreamWriter w,
			final SAMSequenceDictionary dict,	
			final Context context
			) 
			throws XMLStreamException
			{
			final String IDT_XNS="https://umr1087.univ-nantes.fr/";
			final Insets insets=new Insets(20, 20, 20, 20);
			final Dimension dim=new Dimension(
					insets.left+insets.right+this.drawinAreaWidth
					,insets.top+insets.bottom);
			dim.height+=HEIGHT_MAIN_TITLE;
			dim.height+=context.HEIGHT_RULER;
			dim.height+=context.getHeight();
			//dim.height+=100;
			
			w.writeStartElement("svg");
			w.writeAttribute("width", format(dim.width));
			w.writeAttribute("height", format(dim.height));
			w.writeDefaultNamespace(SVG.NS);
			w.writeNamespace("u",IDT_XNS);
			w.writeNamespace("rdf",RDF.NS);
			w.writeNamespace("xlink", XLINK.NS);
			
			/** write meta data */
			w.writeStartElement("meta");
			w.writeStartElement("rdf", "RDF", RDF.NS);
			w.writeStartElement("rdf", "Description", RDF.NS);
			w.writeAttribute("rdf", RDF.NS,"ID","#");
			w.writeStartElement("u", "contig", IDT_XNS);w.writeCharacters(this.interval.getContig());w.writeEndElement();
			w.writeStartElement("u", "start", IDT_XNS);w.writeCharacters(String.valueOf(this.interval.getStart()));w.writeEndElement();
			w.writeStartElement("u", "end", IDT_XNS);w.writeCharacters(String.valueOf(this.interval.getEnd()));w.writeEndElement();
			w.writeStartElement("u", "sample", IDT_XNS);w.writeCharacters(context.sampleName);w.writeEndElement();
			w.writeStartElement("u", "bam", IDT_XNS);w.writeCharacters(context.path.toString());w.writeEndElement();
			w.writeStartElement("u", "reference", IDT_XNS);w.writeCharacters(this.referenceFile.toString());w.writeEndElement();
			w.writeEndElement();//Description
			w.writeEndElement();//rdf
			w.writeEndElement();

			
			w.writeStartElement("title");
			w.writeCharacters(context.sampleName+" : "+this.interval.toNiceString());
			w.writeEndElement();
			
			w.writeStartElement("description");
			w.writeCharacters("Cmd:"+getProgramCommandLine()+"\n");
			w.writeCharacters("Version:"+getVersion()+"\n");
			w.writeCharacters("Author: Pierre Lindenbaum\n");
			w.writeEndElement();
			
			
			w.writeStartElement("style");
			w.writeCharacters(
					"g.maing {stroke:black;stroke-width:0.5px;fill:none;}\n"+
					".maintitle {stroke:blue;fill:none;}\n"+
					".bA {stroke:green;}\n" + 
					".bC {stroke:blue;}\n" +
					".bG {stroke:black;}\n" +
					".bT {stroke:red;}\n" +
					".bN {stroke:gray;}\n" +
					".homvar {stroke:blue;fill:blue;}\n"+
					".het {stroke:blue;fill:blue;}\n"+
					"line.insert {stroke:red;stroke-width:4px;}\n"+
					"line.rulerline {stroke:lightgray;stroke-width:0.5px;}\n"+
					"line.rulerlabel {stroke:gray;stroke-width:2px;}\n"+
					"text.maintitle {fill:darkgray;stroke:black;text-anchor:middle;font-size:30px;}\n"+
					"path.samplename {stroke:black;stroke-width:3px;}\n"+
					"circle.mutation {stroke:red;fill:none;stroke-width:"+(context.featureWidth<5?0.5:3.0)+"px;}"+
					""
					);
			w.writeEndElement();//style
			
			w.writeStartElement("defs");
			
			//mutations
			if(!context.pos2variant.isEmpty())
				{
				//hom var
				w.writeEmptyElement("rect");
				w.writeAttribute("id","homvar");
				w.writeAttribute("class","homvar");
				w.writeAttribute("x", "0");
				w.writeAttribute("y", "0");
				w.writeAttribute("width", format(context.featureWidth));
				w.writeAttribute("height",  format(context.featureHeight));
				
				//het
				w.writeStartElement("g");
				w.writeAttribute("id","het");
				w.writeAttribute("class","het");
					
				
					w.writeEmptyElement("rect");
					w.writeAttribute("x", "0");
					w.writeAttribute("y", "0");
					w.writeAttribute("width", format(context.featureWidth/2.0));
					w.writeAttribute("height",  format(context.featureHeight/2.0));
					w.writeEmptyElement("rect");
					w.writeAttribute("x", format(context.featureWidth/2.0));
					w.writeAttribute("y",  format(context.featureHeight/2.0));
					w.writeAttribute("width", format(context.featureWidth/2.0));
					w.writeAttribute("height",  format(context.featureHeight/2.0));
					w.writeEmptyElement("rect");
					w.writeAttribute("style", "fill:none;");
					w.writeAttribute("x", "0");
					w.writeAttribute("y", "0");
					w.writeAttribute("width", format(context.featureWidth));
					w.writeAttribute("height",  format(context.featureHeight));
	
				w.writeEndElement();
				}
			
			
	
			//base
			for(String base:new String[]{"a","A","t","T","g","G","c","C","n","N"})
				{
				w.writeStartElement("path");
				w.writeAttribute("id","b"+base);
				
				w.writeAttribute("class","b"+base.toUpperCase());
				w.writeAttribute("d",this.hershey.svgPath(
						base,
						0,
						0,
						context.featureWidth*0.95,
						context.featureHeight*0.95
						));
				writeTitle(w,base);
				w.writeEndElement();
				}
			//mutation
			w.writeStartElement("circle");
			w.writeAttribute("id","mut");
			w.writeAttribute("class","mutation");
			
			w.writeAttribute("cx",format(context.featureWidth/2.0));
			w.writeAttribute("cy",format(context.featureHeight/2.0));
			w.writeAttribute("r",format(context.featureHeight/2.0));
			writeTitle(w,"mutation");
			w.writeEndElement();
						
			printGradientDef(w,
					"fail",//fail qc
					"stop-color:#DDA0DD;stop-opacity:1;",
					"stop-color:#DA70D6;stop-opacity:1;"
					);
			/* properly paired */
			printGradientDef(w,
					"fprop",
					"stop-color:#90EE90;stop-opacity:1;",
					"stop-color:#98FB98;stop-opacity:1;"
					);
			/* discordant reads */
			printGradientDef(w,
					"fdisc",
					"stop-color:#D8BFD8;stop-opacity:1;",
					"stop-color:#F5F5DC;stop-opacity:1;"
					);
			/* paired mate unmapped */
			printGradientDef(w,
					"fum",
					"stop-color:#add8e6;stop-opacity:1;", /* light blue */
					"stop-color:#e0ffff;stop-opacity:1;" /* light cyan */
					);
			/* MAPQ=0 */
			printGradientDef(w,
					"fq0",
					"stop-color:#ffe4e1;stop-opacity:1;",/* misty rose */
					"stop-color:#fffff;stop-opacity:1;"
					);
			
			printGradientDef(w,
					"f0",
					"stop-color:#dcdcdc;stop-opacity:1;",/* gainsboro */
					"stop-color:#f5f5f5;stop-opacity:1;" /* white smoke */
					);
			
			w.writeEndElement();//defs
			
			w.writeStartElement("g");
			w.writeAttribute("transform", "translate("+insets.left+","+insets.top+")");
			w.writeAttribute("class","maing");
			
			int y=insets.top;
			/* write title */
			
			final Hyperlink hyperlink = Hyperlink.compile(dict);
			if(hyperlink.apply(this.interval).isPresent()) {
				w.writeStartElement("a");
				w.writeAttribute("href",hyperlink.apply(this.interval).orElse(""));
				}
			
			w.writeStartElement("text");
			w.writeAttribute("class","maintitle");
			w.writeAttribute("x",String.valueOf(this.drawinAreaWidth/2.0));
			w.writeAttribute("y",String.valueOf(y));
			w.writeCharacters(this.interval.toNiceString());
			w.writeEndElement();
			if(hyperlink.apply(this.interval).isPresent()) {
				w.writeEndElement();//a
				}
			y+=HEIGHT_MAIN_TITLE;
			
			/* write ruler */
			int prev_ruler_printed=-1;
			double prev_ruler_printed_x=-1;
			w.writeStartElement("g");
			for(int pos=this.interval.getStart(); pos<=this.interval.getEnd();++pos)
				{
				double x= this.baseToPixel(pos);
				if(prev_ruler_printed_x!=-1 &&  prev_ruler_printed_x+ Math.max(20.0,context.featureHeight) >= x)continue;
				prev_ruler_printed_x = x;
				
				w.writeStartElement("line");
				w.writeAttribute("class","rulerline");
				
				w.writeAttribute("x1",format(x));
				w.writeAttribute("x2",format(x));
				w.writeAttribute("y1",format(y+context.HEIGHT_RULER-5));
				w.writeAttribute("y2",format(dim.height));
				if(prev_ruler_printed==-1 ||  this.baseToPixel(prev_ruler_printed)+context.featureHeight < x)
					{
					x= (this.baseToPixel(pos)+this.baseToPixel(pos+1))/2.0;
					w.writeStartElement("path");
					w.writeAttribute("class","rulerlabel");
					w.writeAttribute("d",this.hershey.svgPath(StringUtils.niceInt(pos),0,0,Math.min(StringUtils.niceInt(pos).length()*context.featureHeight,context.HEIGHT_RULER)-20,context.featureHeight));
					w.writeAttribute("transform","translate("+(x-context.featureHeight/2.0)+","+ (y+(context.HEIGHT_RULER-10)) +") rotate(-90) ");
					writeTitle(w,StringUtils.niceInt(pos));
					w.writeEndElement();
					prev_ruler_printed=pos;
					}
				writeTitle(w,StringUtils.niceInt(pos));
				w.writeEndElement();//line
				}
			w.writeEndElement();//g for ruler
			y+=context.HEIGHT_RULER;
			
			printSample(w,y,context);
			y+=context.getHeight();			
			w.writeEndElement();//g
			w.writeEndElement();//svg
			}
		
		private void printSamRecord(
				final XMLStreamWriter w,
				final Context context,
				final SAMRecord record,
				final Map<Integer,Counter<Character>> consensus
				)
			throws XMLStreamException
			{
			double mid_y= context.featureHeight/2.0;
			final double y_h95= context.featureHeight*0.90;
			final double y_top5=(context.featureHeight-y_h95)/2.0;
			final double y_bot5= y_top5+y_h95;
			final double arrow_w= context.featureHeight/3.0;
	
			
			/* print that sam record */
			final int unclipped_start= record.getUnclippedStart();
			Cigar cigar = record.getCigar();
			if(cigar==null) return;
			byte bases[]=record.getReadBases();
			if(bases==null || bases.equals(SAMRecord.NULL_SEQUENCE)) return;
			byte qualities[]=record.getBaseQualities();
			if(qualities==null||qualities.equals(SAMRecord.NULL_QUALS)) return;
			
			w.writeStartElement("g");
			writeTitle(w,getReadNameForPopup(record));

			
			int readPos=0;
			Map<Integer,String> pos2insertions=new HashMap<Integer,String>();
			List<CigarElement> cigarElements= cigar.getCigarElements();
			
			/* find position of arrow */
			int arrow_cigar_index=-1;
			for(int cidx=0; cidx< cigarElements.size(); cidx++ )
				{
				final CigarElement ce = cigarElements.get(cidx);
				final CigarOperator op=ce.getOperator();
				
				switch(op)
					{
					case H:case S: if(!this.showClipping) break;//threw
					case M:case EQ: case X:
						{
						arrow_cigar_index=cidx;
						}
					default:break;
					}
				if(record.getReadNegativeStrandFlag() && arrow_cigar_index!=-1)
					{
					break;
					}
				}
	
			
			int refPos=unclipped_start;
			
			/* loop over cigar string */
			for(int cidx=0; cidx< cigarElements.size(); cidx++ )
				{
				final CigarElement ce = cigarElements.get(cidx);
				final CigarOperator op=ce.getOperator();
				boolean in_clip=false;
				switch(ce.getOperator())
					{
					case D://cont
					case N:
						{
						int c_start = trim(refPos);
						int c_end   = trim(refPos + ce.getLength());
						if(c_start<c_end)
							{
							w.writeStartElement("line");
							w.writeAttribute("class","indel");
							w.writeAttribute("x1", format(baseToPixel(c_start)));
							w.writeAttribute("x2", format(baseToPixel(c_end)));
							w.writeAttribute("y1", format(mid_y));
							w.writeAttribute("y2", format(mid_y));
							writeTitle(w,op.name()+ce.getLength());
							w.writeEndElement();
							}
						
						refPos += ce.getLength();
						break;
						}
					case I: 
						{
						final StringBuilder sb=new StringBuilder();
						for(int i=0;i< ce.getLength();++i)
							{
							sb.append((char)bases[readPos++]);
							}
						pos2insertions.put(refPos, sb.toString());
						break;
						}
					case H:
						{
						if(!this.showClipping)
							{
							refPos+=ce.getLength();
							break;
							}
						in_clip=true;
						//NO break;
						}
					case S:
						{
						if(!this.showClipping)
							{
							readPos+=ce.getLength();
							refPos+=ce.getLength();
							break;
							}
						in_clip=true;
						//NO break;
						}
					case X:
					case EQ:
					case M:
						{
						int match_start = refPos;
						int match_end = refPos + ce.getLength();
						
						//print sam background
						StringBuilder sb=new StringBuilder();
						if(record.getReadNegativeStrandFlag() &&
							match_start >= this.interval.getStart() && 
							match_start <= this.interval.getEnd() &&
							cidx==arrow_cigar_index)
							{
							sb.append(" M ").append(format(baseToPixel(match_start)+arrow_w)).append(',').append(y_top5);
							sb.append(" h ").append(format(baseToPixel(trim(match_end))-baseToPixel(match_start)-arrow_w));
							sb.append(" v ").append(format(y_h95));
							sb.append(" h ").append(format(-(baseToPixel(trim(match_end))-baseToPixel(match_start)-arrow_w)));
							sb.append(" L ").append(format(baseToPixel(match_start))).append(',').append(mid_y);
							sb.append(" Z");
							}
						else if(!record.getReadNegativeStrandFlag() &&
								match_end >= this.interval.getStart() && 
								match_end <= this.interval.getEnd() &&
								cidx==arrow_cigar_index
								)
							{
							sb.append(" M ").append(format(baseToPixel(match_end)-arrow_w)).append(',').append(y_top5);
							sb.append(" h ").append(format(-(baseToPixel(match_end)-baseToPixel(trim(match_start))-arrow_w)));
							sb.append(" v ").append(format(y_h95));
							sb.append(" h ").append(format(baseToPixel(match_end)-baseToPixel(trim(match_start))-arrow_w));
							sb.append(" L ").append(format(baseToPixel(match_end))).append(',').append(mid_y);
							sb.append(" Z");
							}
						else
							{
							sb.append(" M ").append(format(baseToPixel(trim(match_start)))).append(',').append(y_top5);
							sb.append(" h ").append(format(baseToPixel(trim(match_end))-baseToPixel(trim(match_start))));
							sb.append(" v ").append(format(y_h95));
							sb.append(" h ").append(format(-(baseToPixel(trim(match_end))-baseToPixel(trim(match_start)))));
							sb.append(" Z");
							}
						w.writeEmptyElement("path");
						w.writeAttribute("d",sb.toString().trim());
						if( ce.getOperator()==CigarOperator.H || ce.getOperator()==CigarOperator.S)
							{
							w.writeAttribute("style","fill:yellow;");
							}
						else
							{
							final String rec_style;
							if(record.getReadFailsVendorQualityCheckFlag()) {
								rec_style = "fill:url(#fail);";
								}
							else if(record.getReadPairedFlag() && !record.getMateUnmappedFlag() && !record.getMateReferenceName().equals(this.interval.getContig())) {
								rec_style = "fill:url(#fdisc);";
								}
							else if(record.getReadPairedFlag() && record.getMateUnmappedFlag()) {
								rec_style = "fill:url(#fum);";
								}
							else if(record.getReadPairedFlag() && !record.getProperPairFlag()) {
								rec_style = "fill:url(#fprop);";
								}
							else if(record.getMappingQuality()==0) {
								rec_style = "fill:url(#fq0);";
								}
							else
								{
								rec_style = "fill:url(#f0);";
								}
							w.writeAttribute("style",rec_style);
							}
						
						if(op.consumesReadBases())
							{
							for(int i=0;i< ce.getLength();++i)
								{
								char ca=(char)bases[readPos];
								if(!in_clip)
									{
									Counter<Character> count= consensus.get(refPos+i) ;
									if(count==null)
										{
										count=new Counter<>();
										consensus.put(refPos+i,count) ;
										}
									count.incr(ca);
									}
								char cb='N';
								if(genomicSequence!=null)
									{
									cb=Character.toUpperCase(genomicSequence.charAt(refPos+i-1));
									}
								
								if(cb!='N' && ca!='N' && cb!=ca && !in_clip && this.interval.contains(refPos+i))
									{
									w.writeEmptyElement("use");
									w.writeAttribute("x",format( baseToPixel(refPos+i)));
									//w.writeAttribute("y",format(y_top));never needed
									w.writeAttribute("xlink",XLINK.NS,"href", "#mut");
									}
								
								if(this.interval.contains(refPos+i) && 
									context.featureWidth >  5)
									{
									//int qual = qualities[readPos];
									w.writeEmptyElement("use");
									w.writeAttribute("x",format( baseToPixel(refPos+i)));
									//w.writeAttribute("y",format(y_top));never needed
									w.writeAttribute("xlink",XLINK.NS,"href", "#b"+ca);
									}
								readPos++;
								}
							}
						
	
						
						refPos+=ce.getLength();
						break;
						}
					default:
						{
						throw new RuntimeException("Unknown SAM op: "+ce.getOperator());
						}
					}
				}
			
			for(Integer pos:pos2insertions.keySet())
				{
				if(pos < this.interval.getStart())  continue;
				if(pos > this.interval.getEnd())  continue;
				final String insertion=pos2insertions.get(pos);
				w.writeStartElement("line");
				final double x= baseToPixel(pos);
				
				w.writeAttribute("class","insert");
				w.writeAttribute("x1",format(x));
				w.writeAttribute("x2",format(x));
				w.writeAttribute("y1",format(y_top5));
				w.writeAttribute("y2",format(y_bot5));
				writeTitle(w,"Insertion "+insertion);
				w.writeEndElement();//line
				}
			w.writeEndElement();//g
			}
		
		private void plotSVG(
				final ArchiveFactory archive,
				final PrintWriter manifest,
				final Path path
				) throws IOException,XMLStreamException {
			XMLStreamWriter w=null;
			
			final SamReaderFactory sfrf= SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT).
					referenceSequence(this.referenceFile);

			final Context context = new Context();
			context.path = path;
			
			try(SamReader samReader = sfrf.open(path)) {
				final SAMFileHeader header= samReader.getFileHeader();
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
				SequenceUtil.assertSequenceDictionariesEqual(dict, this.referenceDict);
				context.sampleName = header.getReadGroups().stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtils.isBlank(S)).
						findFirst().
						orElse(IOUtils.getFilenameWithoutCommonSuffixes(path));
				
				final int extend_for_clip = 100;//get a chance to get clipping close to bounds
				// read SamRecord
				try(CloseableIterator<SAMRecord> iter= samReader.query(this.interval.getContig(),
						Math.max(1, this.interval.getStart() - extend_for_clip),this.interval.getEnd() + extend_for_clip,false)) {
					final Pileup<SAMRecord> pileup = new Pileup<>((A,B)->!CoordMath.overlaps(left(A)-1,right(A)+1,left(B),right(B)));
					while(iter.hasNext())
						{
						final SAMRecord rec = iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						
						//if( this.samRecordFilter.filterOut(rec)) continue;
						if(!SAMRecordDefaultFilter.accept(rec, this.mappingQuality)) continue;
						if( !rec.getReferenceName().equals(this.interval.getContig())) continue;
						if( right(rec)  < this.interval.getStart()) continue;
						if( left(rec)  > this.interval.getEnd())continue;
						pileup.add(rec);
						}
					context.lines = pileup.getRows();
					}
				
				
				if(this.vcf!=null)
					{
					context.readVariantFile(vcf);
					}
				
				
				
				context.featureWidth=  this.drawinAreaWidth/(double)((this.interval.getEnd() - this.interval.getStart())+1); 
				context.featureHeight= Math.min(Math.max(5.0,context.featureWidth),30); 
				context.HEIGHT_RULER=(int)(StringUtils.niceInt(this.interval.getEnd()).length()*context.featureHeight+5);
				LOG.info("Feature height:"+context.featureHeight);
				
				final String buildName= SequenceDictionaryUtils.getBuildName(this.referenceDict).orElse("");
				final String filename= this.prefix+buildName+(buildName.isEmpty()?"":".")+this.interval.getContig()+"_"+this.interval.getStart()+"_"+this.interval.getEnd()+"."+context.sampleName+".svg";
				

				LOG.info("writing "+ filename);
				manifest.print(context.sampleName);
				manifest.print("\t");
				manifest.print(path);
				manifest.print("\t");
				manifest.print(this.interval);
				manifest.print("\t");
				manifest.print(this.referenceFile);
				manifest.print("\t");
				manifest.print(filename);
				manifest.println();

				
				final XMLOutputFactory xof=XMLOutputFactory.newFactory();
				try(OutputStream fout= archive.openOuputStream(filename)) {
					w=xof.createXMLStreamWriter(fout, "UTF-8");
					w.writeStartDocument("UTF-8", "1.0");
					printDocument(w,dict,context);
					w.writeEndDocument();
					w.flush();
					w.close();
					fout.flush();
					}
				}
			}
		
		
		@Override
		public int doWork(final List<String> args) {
			/* parse interval */
			if(StringUtils.isBlank(this.intervalStr))
				{
				LOG.error("Interval undefined");
				return -1;
				}
			
			
			
			try
				{
				this.indexedFastaSequenceFile =  Objects.requireNonNull(ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referenceFile));
				this.referenceDict = SequenceDictionaryUtils.extractRequired(indexedFastaSequenceFile);
				this.interval  =  IntervalParserFactory.
					newInstance(referenceDict).
					make().
					apply(this.intervalStr).
					orElseThrow(IntervalParserFactory.exception(this.intervalStr))
					;
				
				this.genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, this.interval.getContig());
				this.drawinAreaWidth = Math.max(100,this.drawinAreaWidth );
				
				
				
				final List<Path> inputFiles = IOUtils.unrollPaths(args);
				/* read SAM data */
				if(inputFiles.isEmpty())
					{
					LOG.error("No BAM");
					return -1;
					}
				
				try(ArchiveFactory archive = ArchiveFactory.open(this.outputFile)) {
					try(PrintWriter mf = archive.openWriter(this.prefix+"manifest.tsv")) {
						mf.print("sample");
						mf.print("\t");
						mf.print("bam");
						mf.print("\t");
						mf.print("interval");
						mf.print("\t");
						mf.print("reference");
						mf.print("\t");
						mf.print("svg");
						mf.println();
						
						for(final Path filename: inputFiles)
							{
							LOG.info("Reading from "+filename);
							plotSVG(archive,mf,filename);
							}
						mf.flush();
						}
					}				
				return 0;
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(this.indexedFastaSequenceFile);
				this.indexedFastaSequenceFile=null;
				this.interval=null;
				}
			}
		
	
	
	public static void main(final String[] args)
		{
		new BamToSVG().instanceMainWithExit(args);
		}
	}
