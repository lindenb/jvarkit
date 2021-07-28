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
package com.github.lindenb.jvarkit.tools.bam2svg;

import java.awt.Dimension;
import java.awt.Insets;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
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
public class BamToSVG extends Launcher {
	private static final int HEIGHT_MAIN_TITLE=100;
	private static final int HEIGHT_SAMPLE_NAME=50;
	
	
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
		private String format(double v)
			{
			return this.decimalFormater.format(v);
			}
		
	
		
		private int trim(int pos0)
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
		
		private String getReadNameForPopup(final SAMRecord rec) {
			 final StringBuilder sb = new StringBuilder(rec.getReadName());
			
			 sb.append(rec.getReadNegativeStrandFlag()?" -":" +"); 
			
			if(rec.getReadPairedFlag()) {
				if(rec.getFirstOfPairFlag()) sb.append(" 1");
				if(rec.getSecondOfPairFlag()) sb.append(" 2");
				sb.append("/2");
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
				XMLStreamWriter w,
				String gradId,
				String styleTop,//e.g: "stop-color:black;stop-opacity:1;"
				String styleMid //e.g: "stop-color:white;stop-opacity:1;"
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
		
		private void printSample(XMLStreamWriter w,double top_y,final Context context)
				throws XMLStreamException
			{
			w.writeComment("START SAMPLE: "+context.sampleName);
			w.writeStartElement(SVG.NS,"g");
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
		
		private void printDocument(final XMLStreamWriter w,final Context context) 
			throws XMLStreamException
			{
			Insets insets=new Insets(20, 20, 20, 20);
			Dimension dim=new Dimension(
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
			w.writeNamespace("xlink", XLINK.NS);
			
	
			
			w.writeStartElement(SVG.NS,"title");
			w.writeCharacters(intervalStr);
			w.writeEndElement();
			
			w.writeStartElement(SVG.NS,"description");
			w.writeCharacters("Cmd:"+getProgramCommandLine()+"\n");
			w.writeCharacters("Version:"+getVersion()+"\n");
			w.writeCharacters("Author: Pierre Lindenbaum\n");
			w.writeEndElement();
			
			
			w.writeStartElement(SVG.NS,"style");
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
					"path.maintitle {stroke:blue;stroke-width:5px;}\n"+
					"path.samplename {stroke:black;stroke-width:3px;}\n"+
					"circle.mutation {stroke:red;fill:none;stroke-width:"+(context.featureWidth<5?0.5:3.0)+"px;}"+
					""
					);
			w.writeEndElement();//style
			
			w.writeStartElement(SVG.NS,"defs");
			
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
			
			
			String pastels[]={"fff8dc", "cd5c5c", "708090", "ff4500", "4682b4", "ffa500", "d8bfd8", "fa8072", "00ffff", "3cb371", "000000", "1e90ff", "cd853f", "4169e1", "f0ffff", "5f9ea0", "eeeee0", "a020f0", "bc8f8f", "ff69b4", "00ff7f", "000080", "e9967a", "daa520", "32cd32", "ee82ee", "7b68ee", "ffffe0", "8b8378", "db7093", "a0522d", "ffe4b5", "9370db", "afeeee", "eee5de", "8fbc8f", "8b8682", "8470ff", "8b8b83", "ff1493", "eedfcc", "f0f8ff", "ff6347", "faebd7", "add8e6", "fffafa", "c1cdc1", "cdc0b0", "c71585", "7fffd4", "7fff00", "cdaf95", "e0eee0", "556b2f", "dcdcdc", "ffebcd", "cdc8b1", "fff0f5", "ffe4e1", "9acd32", "ffc0cb", "8b8989", "ffff00", "cdb79e", "f5f5f5", "2f4f4f", "eee9e9", "cdcdc1", "eee8aa", "bebebe", "ffffff", "2e8b57", "fffff0", "b0e0e6", "fff5ee", "778899", "fffacd", "191970", "d2691e", "eecbad", "40e0d0", "ffd700", "eed5b7", "fffaf0", "ffefd5", "98fb98", "b22222", "87ceeb", "483d8b", "adff2f", "b8860b", "66cdaa", "f5fffa", "ff8c00", "00fa9a", "f4a460", "dda0dd", "fafad2", "f5f5dc", "9400d3", "006400", "cdc5bf", "a52a2a", "8b7d6b", "0000ff", "d02090", "ffb6c1", "48d1cc", "e0ffff", "f8f8ff", "d2b48c", "00ced1", "8b8878", "0000cd", "e6e6fa", "f0e68c", "6495ed", "f0fff0", "ffe4c4", "ff7f50", "d3d3d3", "00bfff", "b03060", "6a5acd", "ffa07a", "8b7765", "20b2aa", "ff0000", "f5deb3", "eee8dc", "ba55d3", "faf0e6", "6b8e23", "8a2be2", "7cfc00", "ffdab9", "8b4513", "87cefa", "cdc9c9", "9932cc", "deb887", "eedd82", "228b22", "838b83", "b0c4de", "f08080", "696969", "ffdead", "da70d6", "bdb76b", "fdf5e6"};
			int flags[]=new int[]{65,73,81,83,97,99,113,121,129,137,145,147,161,163,177,185,321,323,329,337,339,353,355,369,371,377,385,387,393,401,403,417,419,433,435,1089,1097,1105,1107,1121,1123,1137,1145,1153,1161,1169,1171,1185,1187,1201,1209};
			for(int flag=0;flag< flags.length;++flag)
				{
				printGradientDef(w, "f"+flags[flag],
						"stop-color:#"+pastels[flag%pastels.length]+
						";stop-opacity:1;", "stop-color:white;stop-opacity:1;");
				}
			
			w.writeEndElement();//defs
			
			w.writeStartElement("g");
			w.writeAttribute("transform", "translate("+insets.left+","+insets.top+")");
			w.writeAttribute("class","maing");
			
			int y=insets.top;
			/* write title */
			w.writeStartElement("path");
			w.writeAttribute("class","maintitle");
			w.writeAttribute("d", this.hershey.svgPath(
					this.interval.getContig()+":"+StringUtils.niceInt(this.interval.getStart())+"-"+StringUtils.niceInt(this.interval.getEnd()),
					scaleRect(new Rectangle2D.Double(
							0,
							y,
							this.drawinAreaWidth,
							HEIGHT_MAIN_TITLE
							)))
					);
			writeTitle(w, this.intervalStr);
			w.writeEndElement();
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
			
			w.writeStartElement(SVG.NS,"g");
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
							w.writeAttribute("style","fill:url(#f"+record.getFlags()+");");
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
				w.writeEndElement();
				}
			w.writeEndElement();//g
			
			}
		
		private void plotSVG(final ArchiveFactory archive,final Path path) throws IOException,XMLStreamException {
			XMLStreamWriter w=null;
			
			final SamReaderFactory sfrf= SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT).
					referenceSequence(this.referenceFile);

			final Context context = new Context();
			
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
				final XMLOutputFactory xof=XMLOutputFactory.newFactory();
				try(OutputStream fout= archive.openOuputStream(filename)) {
					w=xof.createXMLStreamWriter(fout, "UTF-8");
					w.writeStartDocument("UTF-8", "1.0");
					printDocument(w,context);
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
					for(final Path filename: inputFiles)
						{
						LOG.info("Reading from "+filename);
						plotSVG(archive,filename);
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
