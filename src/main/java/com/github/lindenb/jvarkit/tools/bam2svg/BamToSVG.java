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

import java.awt.Dimension;
import java.awt.Insets;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.hershey.Hershey;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.net.UrlSupplier;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.rdf.ns.XLINK;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.readers.TabixReader;
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
modificationDate="20210728",
jvarkit_amalgamion =  true,
menu="BAM Visualization"
)
public class BamToSVG extends Launcher {
	private static final int HEIGHT_MAIN_TITLE=30;
	private static final int HEIGHT_SAMPLE_NAME=20;
	
	
	private static final Logger LOG = Logger.of(BamToSVG.class);

	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private File outputFile = null;
	@Parameter(names={"-i","--interval","--region"},description=IntervalParser.OPT_DESC,required=true)
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
	@Parameter(names={"--bases"},description="print bases in read")
	private boolean print_bases = false;
    @Parameter(names= {"--gff","--gff3"},description="Optional Tabix indexed GFF3 file.")
    private Path gtfFile = null;

	@DynamicParameter(names = "-D", description = "custom parameters. '-Dkey=value'. Undocumented.")
	private Map<String, String> dynaParams = new HashMap<>();


	private final Hershey hershey=new Hershey();
	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private SAMSequenceDictionary referenceDict=null;
	private CharSequence genomicSequence=null;
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
			dim_height+= HEIGHT_RULER;
			dim_height+= this.featureHeight;//ref seq
			dim_height+= this.featureHeight;//consensus
			if(!pos2variant.isEmpty())dim_height+=  this.featureHeight;//variants
			dim_height+= this.lines.size()* this.featureHeight;//consensus
			if(gtfFile!=null) dim_height+=this.featureHeight;
			return dim_height;
			}
		
		private void readVariantFile(final Path vcf)  throws IOException
			{
			try(VCFReader r= VCFReaderFactory.makeDefault().open(vcf, true)) {
				r.query(BamToSVG.this.interval).stream().
					filter(V->V.getGenotype(this.sampleName)!=null).
					forEach(ctx->this.pos2variant.add(ctx));
				}
			LOG.info("num variants "+this.pos2variant.size());
			}
		
		}
		
	private class SvgOutput {
		String sampleName;
		String svgFilename;
		Dimension dimension;
		final Set<String> rsIds = new HashSet<>();
		final Set<String> geneNames = new HashSet<>();
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
		
		if(rec.getReadPairedFlag()) {
			if(rec.getFirstOfPairFlag()) sb.append(" 1");
		if(rec.getSecondOfPairFlag()) sb.append(" 2");
		sb.append("/2");
		sb.append("\n");
		sb.append("pos: ").append(rec.getContig()+":"+StringUtils.niceInt(rec.getStart())+"-"+StringUtils.niceInt(rec.getEnd())).append("\n");
		if(rec.getStart()!=rec.getUnclippedStart() || rec.getEnd()!=rec.getUnclippedEnd()) {
			sb.append("unclipped: ").append(rec.getContig()+":"+StringUtils.niceInt(rec.getUnclippedStart())+"-"+StringUtils.niceInt(rec.getUnclippedEnd())).append("\n");
		}

		sb.append("strand:").append(rec.getReadNegativeStrandFlag()?" -":" +").append("\n"); 
		sb.append("flag:").append(rec.getFlags()).append("\n");
		sb.append("MAPQ:"+rec.getMappingQuality()).append("\n");
		if(rec.getMateUnmappedFlag())  {
			sb.append(" mate:unmapped\n");
			}
		else
			{
			sb.append("mate:"+rec.getMateReferenceName()+":"+rec.getMateAlignmentStart()+" "+(rec.getMateNegativeStrandFlag()?" -":" +")).append("\n");
			sb.append("tlen:"+StringUtils.niceInt(rec.getInferredInsertSize())).append("\n");
			}
		}
		if(rec.getReadFailsVendorQualityCheckFlag())  sb.append("fails-quality\n");
		if(rec.getDuplicateReadFlag()) sb.append("duplicate\n");
		if(rec.isSecondaryAlignment())  sb.append("secondary\n");
		if(rec.getSupplementaryAlignmentFlag())  sb.append("supplementary\n");
		for(SAMTagAndValue stav : rec.getAttributes()) {
			sb.append(stav.tag).append(":").append(stav.value).append("\n");
			}
		sb.append("\n");
		return sb.toString().trim();
		}
		
	
	private SvgOutput printDocument(
		final XMLStreamWriter w,
		final SAMSequenceDictionary dict,	
		final Context context
		) 
		throws XMLStreamException
		{
		final UrlSupplier urlSupplier=new UrlSupplier(dict);
		final String IDT_XNS="https://umr1087.univ-nantes.fr/";
		final Insets insets=new Insets(20, 20, 20, 20);
		final Dimension dim=new Dimension(insets.left+insets.right+this.drawinAreaWidth,0);
		dim.height += (int)context.getHeight();
		//dim.height+=100;
		
		
		
		w.writeStartElement("svg");
		w.writeAttribute("width", format(dim.width));
		w.writeAttribute("height", format(dim.height+insets.top+insets.bottom));
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
				".frame {stroke:darkgray;fill:none;}\n"+
				".maintitle {stroke:blue;fill:none;}\n"+
				".nocall {stroke:orange;fill:none;stroke-dasharray:2;}\n"+
				".homref {stroke:blue;fill:none;}\n"+
				".homvar {stroke:blue;fill:blue;}\n"+
				".het {stroke:blue;fill:blue;}\n"+
				"line.insert {stroke:red;stroke-width:4px;}\n"+
				".rulerline {stroke:lightgray;stroke-width:1px;}\n"+
				".rulerlabel {fill:gray;font-size:7px;}\n"+
				"text.maintitle {fill:darkgray;stroke:black;text-anchor:middle;font-size:25px;}\n"+
				"path.samplename {stroke:black;stroke-width:3px;}\n"+
				"circle.mutation {stroke:red;fill:none;stroke-width:"+(context.featureWidth<5?0.5:3.0)+"px;}"+
				"line.gene {stroke:blue;stroke-width:4px;}"+
				"rect.exon {stroke:blue;fill:orange;stroke-width:4px;}"+
				""
				);
		for(char base : new char[] {'A','C','G','T','N'}) {
			w.writeCharacters(".b"+base+" {stroke:"+AcidNucleics.cssColor(base)+";stroke-width:"+(context.featureWidth<5?0.5:1.0)+"px}\n");
			}
		w.writeEndElement();//style
		
		w.writeStartElement("defs");
		
		//mutations
		if(!context.pos2variant.isEmpty())
			{
			
			//nocall
			w.writeEmptyElement("rect");
			w.writeAttribute("id","nocall");
			w.writeAttribute("class","nocall");
			w.writeAttribute("x", "0");
			w.writeAttribute("y", "0");
			w.writeAttribute("width", format(context.featureWidth));
			w.writeAttribute("height",  format(context.featureHeight));
			
			//hom ref
			w.writeEmptyElement("rect");
			w.writeAttribute("id","homref");
			w.writeAttribute("class","homref");
			w.writeAttribute("x", "0");
			w.writeAttribute("y", "0");
			w.writeAttribute("width", format(context.featureWidth));
			w.writeAttribute("height",  format(context.featureHeight));
			
			
			
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
		
		w.writeEndElement();//defs
		
		w.writeStartElement("g");
		w.writeAttribute("transform", "translate("+insets.left+","+insets.top+")");
		w.writeAttribute("class","maing");
		
		double y=0;//not insets top because translate above
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
		w.writeCharacters(context.sampleName+" "+ SequenceDictionaryUtils.getBuildName(dict).orElse("") +" "+this.interval.toNiceString());
		w.writeEndElement();//text
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
			
			w.writeStartElement("rect");
			w.writeAttribute("class","rulerline");
			w.writeAttribute("x",format(x));
			w.writeAttribute("width",format(0.5));
			w.writeAttribute("y",format(y));
			w.writeAttribute("height",format(dim.height-y));
			writeTitle(w,StringUtils.niceInt(pos));
			w.writeEndElement();//line
			
			if(prev_ruler_printed==-1 ||  this.baseToPixel(prev_ruler_printed)+context.featureHeight < x)
				{
				x= (this.baseToPixel(pos)+this.baseToPixel(pos+1))/2.0;
				w.writeStartElement("text");
				w.writeAttribute("class", "rulerlabel");
				w.writeAttribute("x","0");
				w.writeAttribute("y","0");
				w.writeAttribute("transform", "translate("+format(x)+","+y+") rotate(90) ");
				w.writeCharacters(StringUtils.niceInt(pos));
				w.writeEndElement();//text
				prev_ruler_printed=pos;
				}
			
			}
		w.writeEndElement();//g for ruler
		y+=context.HEIGHT_RULER;
		
		
		
		/* write REFERENCE */
		w.writeStartElement("g");
		w.writeComment("REFERENCE");
		for(int pos=this.interval.getStart();
				context.featureWidth>5 && //ignore if too small
				pos<=this.interval.getEnd();++pos)
			{
			final char c= this.genomicSequence.charAt(pos-1);
			w.writeStartElement("use");
			w.writeAttribute("x",format( baseToPixel(pos)));
			w.writeAttribute("y",format(y));
			w.writeAttribute("xlink", XLINK.NS, "href", "#b"+c);
			writeTitle(w, String.valueOf(c)+" "+StringUtils.niceInt(pos));
			w.writeEndElement();//use
			}
		w.writeEndElement();//g
		y+= context.featureHeight;
		
		/* skip line for later consensus */
		final double consensus_y=y;
		y+=context.featureHeight;
		final Map<Integer,Counter<Character>> pos2consensus=new HashMap<Integer,Counter<Character>>();
		
		/* print variants */
		if(!context.pos2variant.isEmpty())
			{
			for(VariantContext ctx:context.pos2variant)
				{
				final Genotype g=ctx.getGenotype(context.sampleName);
				if(g==null) continue;
				
				
				if(ctx.getEnd()< this.interval.getStart()) continue;
				if(ctx.getStart()> this.interval.getEnd()) continue;
				
				String url = null;
				if(ctx.hasID()) {
					url = urlSupplier.of(ctx.getID()).stream().map(G->G.getUrl()).findFirst().orElse(null);
					}
				if(url==null) {
					url = hyperlink.apply(ctx).orElse(null);
					}
				if(url!=null) {
					w.writeStartElement("a");
					w.writeAttribute("href", url);
				}
				
				final double x0  = baseToPixel(ctx.getStart());
				w.writeStartElement("use");
				w.writeAttribute("x",format(x0));
				w.writeAttribute("y",format(y));
				final String svType;
				switch(g.getType()) {
					case HET:
					case HOM_REF:
					case HOM_VAR:
					case NO_CALL: svType = g.getType().name().toLowerCase().replace("_", "");break;
					default: svType="";
					}
				
				w.writeAttribute("xlink", XLINK.NS, "href", "#"+svType);
				final StringBuilder sb = new StringBuilder();
				sb.append(ctx.getContig()+":"+ctx.getStart()+" ");
				sb.append(ctx.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/")));
				sb.append("\n");
				sb.append(g.getType().name()).append("\n");
				if(g.isFiltered()) sb.append("FILTERED:"+g.getFilters()+"\n");
				if(g.hasDP()) sb.append("DP:"+g.getDP()).append("\n");
				if(g.hasGQ()) sb.append("GQ:"+g.getGQ()).append("\n");
				if(g.hasAD()) sb.append("AD:"+IntStream.of(g.getAD()).mapToObj(V->String.valueOf(V)).collect(Collectors.joining(","))).append("\n");
				if(g.hasPL()) sb.append("PL:"+IntStream.of(g.getPL()).mapToObj(V->String.valueOf(V)).collect(Collectors.joining(","))).append("\n");
				writeTitle(w,sb.toString());
				w.writeEndElement();//use
				if(url!=null) {
					w.writeEndElement();//a
					}
				}
			y+=context.featureHeight;
			}
		
		/* print all lines */
		w.writeComment("Alignments");
		for(int nLine=0;nLine< context.lines.size();++nLine)
			{
			w.writeStartElement("g");
			w.writeAttribute("transform", "translate(0,"+format(y)+")");
			List<SAMRecord> line= context.lines.get(nLine);
			//loop over records on that line
			for(SAMRecord record: line)
				{
				printSamRecord(w, context,record,pos2consensus);
				}
			w.writeEndElement();//g
			y+=context.featureHeight;
			}
		
		
		/* write consensus */
		w.writeComment("Consensus");
		final int max_depth = (int) pos2consensus.entrySet().stream().mapToLong(KV->KV.getValue().getTotal()).max().orElse(0L);
		for(int pos=this.interval.getStart();
				max_depth > 0L &&
				context.featureWidth>5 && //ignore if too small
				pos<=this.interval.getEnd();++pos)
			{
			final char rc= this.genomicSequence.charAt(pos-1);
			final Counter<Character> cons = pos2consensus.get(pos);
			if(cons==null) continue;
			final char ac=cons.getMostFrequent();
			final double cov_height = (cons.getTotal()/(double)max_depth)*(context.featureHeight-1);
			if(cons.getCountCategories()==1 && ac==rc) {
				final double x0  = baseToPixel(pos);
				w.writeStartElement("rect");
				w.writeAttribute("x", format(x0));
				w.writeAttribute("y", format(consensus_y+context.featureHeight-cov_height));
				w.writeAttribute("width", format(context.featureWidth-1));
				w.writeAttribute("height", format(cov_height));
				w.writeAttribute("style", "stroke:none;fill:gray;"); 
				writeTitle(w,String.valueOf(ac)+" "+cons.count(ac));
				w.writeEndElement();
				}
			else
				{
				double x0  = baseToPixel(pos);
				for(Character a : cons.keySetDecreasing()) {
					final double aw = (context.featureWidth-1)*((double)cons.count(a)/cons.getTotal());
					w.writeStartElement("rect");
					w.writeAttribute("x", format(x0));
					w.writeAttribute("y", format(consensus_y+context.featureHeight-cov_height));
					w.writeAttribute("width", format(aw));
					w.writeAttribute("height", format(cov_height));
					w.writeAttribute("style", "stroke:none;fill:"+AcidNucleics.cssColor(a)); 
					writeTitle(w, ""+a+" "+cons.count(a)+"/"+cons.getTotal()+" "+(int)(100.0*(double)cons.count(a)/cons.getTotal())+"%");
					w.writeEndElement();
					x0+=aw;
					}
				}
			}

		
		final Set<String> gene_names  = new HashSet<>();
		/* draw genes */
		 if(this.gtfFile!=null) {
			 w.writeStartElement("g");
			 try(TabixReader tbr= new TabixReader(this.gtfFile.toString())) {
				 final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
				 final String ctg = cvt.apply(this.interval.getContig());
				 if(!StringUtils.isBlank(ctg)) {
					 final List<GTFLine> genes  = new ArrayList<>();
					 final List<GTFLine> exons  = new ArrayList<>();
					 final GTFCodec codec = new GTFCodec();
					 final TabixReader.Iterator iter=tbr.query(ctg,this.interval.getStart(),this.interval.getEnd());
					 for(;;) {
						 final String line = iter.next();
						 if(line==null) break;
						 if(StringUtils.isBlank(line) ||  line.startsWith("#")) continue;
						 final GTFLine gtfline = codec.decode(line);
						 if(gtfline==null) continue;
						 if(gtfline.getType().equals("gene")) {
							 genes.add(gtfline);
							 
							 gtfline.getAttributes().entrySet().stream().
							 	filter(KV->KV.getKey().equals("gene_name")).
							 	map(KV->KV.getValue()).
							 	forEach(G->gene_names.add(G));
						 }
						 else if(gtfline.getType().equals("exon")) {
							 exons.add(gtfline);
						 }
					 	}/* end for */
				 for(GTFLine gtf:genes) {
					 w.writeStartElement("line");
	                 w.writeAttribute("class", "gene");
	                 w.writeAttribute("x1", format(baseToPixel(trim(gtf.getStart()))));
	                 w.writeAttribute("x2", format(baseToPixel(trim(gtf.getEnd()+1))));
	                 w.writeAttribute("y1", format(y+context.featureHeight*0.5));
	                 w.writeAttribute("y2", format(y+context.featureHeight*0.5));
	                 writeTitle(w,gtf.getAttributes().entrySet().stream().filter(KV->KV.getKey().equals("gene_name")).map(KV->KV.getValue()).findFirst().orElse(null));
	                 w.writeEndElement();
				 	 }
				 for(GTFLine gtf:exons) {
					 final double x1 = baseToPixel(trim(gtf.getStart()));
					 final double x2 = baseToPixel(trim(gtf.getEnd()+1));
					 w.writeStartElement("rect");
	                 w.writeAttribute("class", "exon");
	                 w.writeAttribute("x", format(x1));
	                 w.writeAttribute("width", format(x2-x1));
	                 w.writeAttribute("y1", format(y));
	                 w.writeAttribute("y2", format(y+context.featureHeight-1));
	                 writeTitle(w,gtf.getAttributes().entrySet().stream().filter(KV->KV.getKey().equals("gene_name")).map(KV->KV.getValue()).findFirst().orElse(null));
	                 w.writeEndElement();
				 	 }
				 }/* end if contig!null */
			 } /* end tabix open */
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		w.writeEndElement();
		y+=context.featureHeight;
		 } /* end gtf */
		
		w.writeEndElement();//g for sample
		w.writeComment("END SAMPLE: "+context.sampleName);

		y+= insets.bottom;
		/* surrounding frame for that sample */
		w.writeEmptyElement("rect");
		w.writeAttribute("class","frame");
		w.writeAttribute("x",format(0));
		w.writeAttribute("y",format(0));
		w.writeAttribute("width",format(drawinAreaWidth-1+insets.left+insets.right));
		w.writeAttribute("height",format(y-1));

		
		

		//w.writeEndElement();//g
		w.writeEndElement();//svg
		
		final SvgOutput svgOutput = new SvgOutput();
		svgOutput.sampleName = context.sampleName;
		dim.height+=insets.top+insets.bottom;
		svgOutput.dimension = dim;
		svgOutput.rsIds.addAll(context.pos2variant.stream().filter(V->V.hasID()).map(V->V.getID()).collect(Collectors.toSet()));
		svgOutput.geneNames.addAll(gene_names);
		return svgOutput;
		}
		@SuppressWarnings("fallthrough")
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
			if(bases==null || bases.equals(SAMRecord.NULL_SEQUENCE)) {
				return;
				}
			byte qualities[]=record.getBaseQualities();
			if(qualities==null||qualities.equals(SAMRecord.NULL_QUALS)) return;
			
			w.writeStartElement("g");
			writeTitle(w,getReadNameForPopup(record));

			
			int readPos=0;
			final Map<Integer,String> pos2insertions=new HashMap<Integer,String>();
			final List<CigarElement> cigarElements= cigar.getCigarElements();
			
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
							final int[] rgb = new int[3];
							if(record.getReadFailsVendorQualityCheckFlag()) {
								rgb[0]=255; rgb[1]=102; rgb[2]=255;
								}
							else if(record.getReadPairedFlag() && !record.getMateUnmappedFlag() && !record.getMateReferenceName().equals(this.interval.getContig())) {
								rgb[0]=102; rgb[1]=178; rgb[2]=255;
								}
							else if(record.getReadPairedFlag() && record.getMateUnmappedFlag()) {
								rgb[0]=51; rgb[1]=255; rgb[2]=51;
								}
							else if(record.getReadPairedFlag() && !record.getProperPairFlag()) {
								rgb[0]=153; rgb[1]=76; rgb[2]=0;
								}
							else
								{
								rgb[0]=192; rgb[1]=192; rgb[2]=192;
								}
							final double f= Math.max(record.getMappingQuality(),60)/60.0;
							for(int i=0;i< rgb.length;i++) rgb[i] = (int)(rgb[i]*f);
							w.writeAttribute("style","fill:rgb("+rgb[0]+","+rgb[1]+","+rgb[2]+")");
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
								char cb= genomicSequence.charAt(refPos+i-1);
								
								/* mutation ? */
								if(this.print_bases && cb!='N' && ca!='N' && cb!=ca && !in_clip && this.interval.contains(refPos+i))
									{
									w.writeEmptyElement("use");
									w.writeAttribute("x",format( baseToPixel(refPos+i)));
									//w.writeAttribute("y",format(y_top));never needed
									w.writeAttribute("xlink",XLINK.NS,"href", "#mut");
									}
								
								if( this.interval.contains(refPos+i) && 
									((this.print_bases&&context.featureWidth>5) || (cb!=ca ))
									 )
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
		
		private SvgOutput plotSVG(
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
				}
				
			if(this.vcf!=null)
				{
				context.readVariantFile(vcf);
				}
			
			
			
			context.featureWidth=  this.drawinAreaWidth/(double)(this.interval.getLengthOnReference()); 
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

			final SvgOutput svgOutput;
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			try(OutputStream fout= archive.openOuputStream(filename)) {
				w=xof.createXMLStreamWriter(fout, "UTF-8");
				w.writeStartDocument("UTF-8", "1.0");
				svgOutput = printDocument(w,this.referenceDict,context);
				w.writeEndDocument();
				w.flush();
				w.close();
				fout.flush();
				}
			svgOutput.svgFilename = filename;
			return svgOutput;
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
				this.interval  =  new IntervalParser(referenceDict).
					apply(this.intervalStr).
					orElseThrow(IntervalParser.exception(this.intervalStr))
					;
				
				final GenomicSequence gSequence=new GenomicSequence(this.indexedFastaSequenceFile, this.interval.getContig());
				this.genomicSequence = new DelegateCharSequence(gSequence) {
					@Override
					public char charAt(int index) {
						if(index<0 || index >= gSequence.length()) return 'N';
						return Character.toUpperCase(gSequence.charAt(index));
						}
					};
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
						final List<SvgOutput> svgOutputs = new ArrayList<>(inputFiles.size());
						for(final Path filename: inputFiles)
							{
							LOG.info("Reading from "+filename);
							final SvgOutput svgOutput = plotSVG(archive,mf,filename);
							svgOutputs.add(svgOutput);
							}
						mf.flush();
						
						Collections.sort(svgOutputs,(A,B)->A.sampleName.compareToIgnoreCase(B.sampleName));
						
						/* write HTML page */
						final XMLOutputFactory xof=XMLOutputFactory.newFactory();
						try(OutputStream fout= archive.openOuputStream(this.prefix+"index.html")) {
							final XMLStreamWriter w =xof.createXMLStreamWriter(fout, "UTF-8");
							w.writeStartElement("html");
							w.writeStartElement("head");
							
							w.writeEmptyElement("meta");
								w.writeAttribute("http-equiv", "Content-Type");
								w.writeAttribute("content", "text/html; charset=utf-8");
							w.writeEmptyElement("meta");
								w.writeAttribute("http-equiv", "author");
								w.writeAttribute("content", "Pierre Lindenbaum Phd");
	
							
							w.writeStartElement("title");
							w.writeCharacters(this.interval.toNiceString());
							w.writeEndElement();//title
							
							
							w.writeEndElement();//head
							
							w.writeStartElement("body");
							w.writeStartElement("h2");
							w.writeCharacters(SequenceDictionaryUtils.getBuildName(this.referenceDict).orElse("") + " "+ this.interval.toNiceString());
							w.writeEndElement();//h2
							w.writeEmptyElement("hr");
							final UrlSupplier urlSupplier = new UrlSupplier(this.referenceDict);
							w.writeCharacters("Hyperlinks: ");
							final Set<UrlSupplier.LabelledUrl> urls = new HashSet<>();
							urls.addAll( urlSupplier.of(this.interval));
							svgOutputs.stream().flatMap(SO->SO.rsIds.stream()).forEach(ID->urls.addAll(urlSupplier.of(ID)));
							svgOutputs.stream().flatMap(SO->SO.geneNames.stream()).forEach(ID->urls.addAll(urlSupplier.of("gene_name",ID)));

							
							
							for(UrlSupplier.LabelledUrl url : urls) {
								w.writeCharacters(" [ ");
								w.writeStartElement("a");
								w.writeAttribute("title", url.getUrl());
								w.writeAttribute("href",  url.getUrl());
								w.writeCharacters(url.getLabel());
								w.writeEndElement();
								w.writeCharacters(" ] ");
								}
							
							w.writeEmptyElement("hr");
							w.writeCharacters("Samples: ");
							for(SvgOutput out: svgOutputs) {
								w.writeCharacters(" [ ");
								w.writeStartElement("a");
								w.writeAttribute("title", out.sampleName);
								w.writeAttribute("href", "#"+out.sampleName);
								w.writeCharacters(out.sampleName);
								w.writeEndElement();
								w.writeCharacters(" ] ");
								
								}
							
							
							for(SvgOutput out: svgOutputs) {
								w.writeEmptyElement("hr");
								w.writeStartElement("div");
								w.writeEmptyElement("a");
								w.writeAttribute("name", out.sampleName);
								w.writeStartElement("h3");
								w.writeCharacters(out.sampleName);
								w.writeEndElement();//h3
								w.writeStartElement("img");
								w.writeAttribute("width", String.valueOf(out.dimension.width));
								w.writeAttribute("height", String.valueOf(out.dimension.height));
								w.writeAttribute("src", out.svgFilename);
								w.writeEndElement();//img
								w.writeEndElement();//div
								}
							w.writeEmptyElement("hr");

	                        w.writeStartElement("div");
	                        w.writeCharacters("Made with "+getClass().getSimpleName()+" version:"+getVersion()+". Pierre Lindenbaum PhD. 2022.");
	                        w.writeEndElement();//div
	                        
	                        
	                        w.writeEndElement();//body
	                        w.writeEndElement();//html
	                        w.flush();
	                        w.close();
							fout.flush();
							}

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
