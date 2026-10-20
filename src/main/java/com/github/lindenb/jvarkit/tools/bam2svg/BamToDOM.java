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

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.w3c.dom.Document;
import org.w3c.dom.DocumentFragment;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.hershey.Hershey;
import com.github.lindenb.jvarkit.illumina.IlluminaConstants;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.lang.JvarkitException;
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
import com.github.lindenb.jvarkit.svg.SVGDocument;
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
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;

public class BamToDOM extends SVGDocument {
	private static final int HEIGHT_MAIN_TITLE=30;
	private static final int HEIGHT_SAMPLE_NAME=20;
	
	
	private static final Logger LOG = Logger.of(BamToDOM.class);

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
	private Path fastaReference = null;
	private ReferenceSequenceFile indexedFastaSequenceFile=null;
	private SAMSequenceDictionary referenceDict=null;
	private CharSequence genomicSequence=null;
	private Locatable interval=null;

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
				r.query(BamToDOM.this.interval).stream().
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
		return  ((pos - this.interval.getStart())/(double)(this.interval.getLengthOnReference()))*(this.drawinAreaWidth);
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
		
		
		
		super.svgElement.setAttribute("width", format(dim.width));
		w.setAttribute("height", format(dim.height+insets.top+insets.bottom));
		w.writeDefaultNamespace(SVG.NS);
		w.writeNamespace("u",IDT_XNS);
		w.writeNamespace("rdf",RDF.NS);
		w.writeNamespace("xlink", XLINK.NS);

		
		super.setTitle(context.sampleName+" : "+this.interval.toNiceString());
		
		
		
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
			Element rect = super.rect(0, 0, context.featureWidth, context.featureWidth);
			rect.setAttribute("id", "nocall");
			rect.setAttribute("class", "nocall");
			
			super.defsElement.appendChild(rect);
			
			rect = super.rect(0, 0, context.featureWidth, context.featureWidth);
			rect.setAttribute("id", "homref");
			rect.setAttribute("class", "homref");
			super.defsElement.appendChild(rect);
			
			rect = super.rect(0, 0, context.featureWidth, context.featureWidth);
			rect.setAttribute("id", "homvar");
			rect.setAttribute("class", "homvar");
			super.defsElement.appendChild(rect);
			
			
			Element group = super.group();
			group.setAttribute("id","het");
			group.setAttribute("class","het");
			rect = super.rect(0, 0, context.featureWidth/2.0, context.featureWidth/2.0);
			group.appendChild(rect);
			rect = super.rect(context.featureWidth/2.0, context.featureWidth/2.0, context.featureWidth/2.0, context.featureWidth/2.0);
			group.appendChild(rect);
			rect = super.rect(0, 0, context.featureWidth, context.featureWidth);
			rect.setAttribute("style", "fill:none");
			group.appendChild(rect);
			super.defsElement.appendChild(group);

			
			}
		
		

		//base
		for(String base:new String[]{"a","A","t","T","g","G","c","C","n","N"})
			{
			Element path = super.path();
			super.defsElement.appendChild(path);
			path.setAttribute("id","b"+base);
			path.setAttribute("class","b"+base.toUpperCase());
			path.setAttribute("d",this.hershey.svgPath(
					base,
					0,
					0,
					context.featureWidth*0.95,
					context.featureHeight*0.95
					));
			setTitle(path,base);
			w.writeEndElement();
			}
		//mutation
		
		Element circle = super.circle((context.featureWidth/2.0), (context.featureWidth/2.0), (context.featureWidth/2.0));
		circle.setAttribute("id","mut");
		circle.setAttribute("class","mutation");
		setTitle(circle,"mutation");
		//TODO add 
		
		w.writeEndElement();//defs
		
		Element g = super.group();
		g.setAttribute("transform", "translate("+insets.left+","+insets.top+")");
		g.setAttribute("class","maing");
		
		double y=0;//not insets top because translate above
		/* write title */
		
		final Hyperlink hyperlink = Hyperlink.compile(dict);
		if(hyperlink.apply(this.interval).isPresent()) {
			w.writeStartElement("a");
			w.setAttribute("href",hyperlink.apply(this.interval).orElse(""));
			}
		
		Element text = super.text(
				(this.drawinAreaWidth/2.0), 
				y, 
				context.sampleName+" "+ SequenceDictionaryUtils.getBuildName(dict).orElse("") +" "+this.interval.toNiceString()
				);
		text.setAttribute("class","maintitle");
		//TODO append text
		
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
			
			Element rect = super.rect(x, y, 0.5, (dim.height-y));
			rect.setAttribute("class", "rulerline");
			setTitle(rect,StringUtils.niceInt(pos));
			//TODO append
			
			if(prev_ruler_printed==-1 ||  this.baseToPixel(prev_ruler_printed)+context.featureHeight < x)
				{
				x= (this.baseToPixel(pos)+this.baseToPixel(pos+1))/2.0;
				Element txt = super.text(0, 0,StringUtils.niceInt(pos));
				txt.setAttribute("class", "rulerlabel");
				txt.setAttribute("transform", "translate("+format(x)+","+y+") rotate(90) ");
				//TODO append
				
				
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
			final Element use = super.element("use");
			use.setAttribute("x",format( baseToPixel(pos)));
			use.setAttribute("y",format(y));
			use.setAttribute("href", "#b"+c);
			setTitle(use, String.valueOf(c)+" "+StringUtils.niceInt(pos));
			//TODO add 
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
				final Genotype gt=ctx.getGenotype(context.sampleName);
				if(gt==null) continue;
				
				
				if(ctx.getEnd()< this.interval.getStart()) continue;
				if(ctx.getStart()> this.interval.getEnd()) continue;
				
				String url = null;
				if(ctx.hasID()) {
					url = urlSupplier.of(ctx.getID()).stream().map(G->G.getUrl()).findFirst().orElse(null);
					}
				if(url==null) {
					url = hyperlink.apply(ctx).orElse(null);
					}
				
				
				final double x0  = baseToPixel(ctx.getStart());
				Element use = super.element("use");
				use.setAttribute("x",format(x0));
				use.setAttribute("y",format(y));
				final String svType;
				switch(gt.getType()) {
					case HET:
					case HOM_REF:
					case HOM_VAR:
					case NO_CALL: svType = gt.getType().name().toLowerCase().replace("_", "");break;
					default: svType="";
					}
				
				use.setAttribute("href", "#"+svType);
				final StringBuilder sb = new StringBuilder();
				sb.append(ctx.getContig()+":"+ctx.getStart()+" ");
				sb.append(ctx.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/")));
				sb.append("\n");
				sb.append(gt.getType().name()).append("\n");
				if(gt.isFiltered()) sb.append("FILTERED:"+gt.getFilters()+"\n");
				if(gt.hasDP()) sb.append("DP:"+gt.getDP()).append("\n");
				if(gt.hasGQ()) sb.append("GQ:"+gt.getGQ()).append("\n");
				if(gt.hasAD()) sb.append("AD:"+IntStream.of(gt.getAD()).mapToObj(V->String.valueOf(V)).collect(Collectors.joining(","))).append("\n");
				if(gt.hasPL()) sb.append("PL:"+IntStream.of(gt.getPL()).mapToObj(V->String.valueOf(V)).collect(Collectors.joining(","))).append("\n");
				setTitle(use,sb.toString());
				w.writeEndElement();//use
				use = super.anchor(use, url);
				
				//TODO append use
				}
			y+=context.featureHeight;
			}
		
		/* print all lines */
		super.comment("Alignments");
		for(int nLine=0;nLine< context.lines.size();++nLine)
			{
			Element g_align = super.group();
			g_align.setAttribute("transform", "translate(0,"+format(y)+")");
			List<SAMRecord> line= context.lines.get(nLine);
			//loop over records on that line
			for(SAMRecord record: line)
				{
				printSamRecord(wcontext,record,pos2consensus);
				}
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
				final Element rect = super.rect(
						x0,
						(consensus_y+context.featureHeight-cov_height),
						(context.featureWidth-1),
						(cov_height)
						);
				rect.setAttribute("style", "stroke:none;fill:gray;");
				setTitle(rect,String.valueOf(ac)+" "+cons.count(ac));
				}
			else
				{
				double x0  = baseToPixel(pos);
				for(Character a : cons.keySetDecreasing()) {
					final double aw = (context.featureWidth-1)*((double)cons.count(a)/cons.getTotal());
					final Element rect = super.rect(
						x0,
						(consensus_y+context.featureHeight-cov_height),
						aw,
						cov_height
						);
					rect.setAttribute("style",  "stroke:none;fill:"+AcidNucleics.cssColor(a));
					setTitle(rect, ""+a+" "+cons.count(a)+"/"+cons.getTotal()+" "+(int)(100.0*(double)cons.count(a)/cons.getTotal())+"%");
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
	                 w.setAttribute("class", "gene");
	                 w.setAttribute("x1", format(baseToPixel(trim(gtf.getStart()))));
	                 w.setAttribute("x2", format(baseToPixel(trim(gtf.getEnd()+1))));
	                 w.setAttribute("y1", format(y+context.featureHeight*0.5));
	                 w.setAttribute("y2", format(y+context.featureHeight*0.5));
	                 writeTitle(w,gtf.getAttributes().entrySet().stream().filter(KV->KV.getKey().equals("gene_name")).map(KV->KV.getValue()).findFirst().orElse(null));
	                 w.writeEndElement();
				 	 }
				 for(GTFLine gtf:exons) {
					 final double x1 = baseToPixel(trim(gtf.getStart()));
					 final double x2 = baseToPixel(trim(gtf.getEnd()+1));
					 w.writeStartElement("rect");
	                 w.setAttribute("class", "exon");
	                 w.setAttribute("x", format(x1));
	                 w.setAttribute("width", format(x2-x1));
	                 w.setAttribute("y1", format(y));
	                 w.setAttribute("y2", format(y+context.featureHeight-1));
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
		w.setAttribute("class","frame");
		w.setAttribute("x",format(0));
		w.setAttribute("y",format(0));
		w.setAttribute("width",format(drawinAreaWidth-1+insets.left+insets.right));
		w.setAttribute("height",format(y-1));

		
		

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
		private Node printSamRecord(
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
			final Cigar cigar = record.getCigar();
			if(cigar==null) return fragment();
			final byte bases[]=record.getReadBases();
			if(bases==null || bases.equals(SAMRecord.NULL_SEQUENCE)) {
				return fragment();
				}
			final byte qualities[]=record.getBaseQualities();
			if(qualities==null||qualities.equals(SAMRecord.NULL_QUALS)) return fragment();
			
			Element g = super.group();
			setTitle(g,getReadNameForPopup(record));
			

			
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
							final Element line = super.line(
									baseToPixel(c_start),
									mid_y,
									baseToPixel(c_end),
									mid_y
									);
							setTitle(line, op.name()+ce.getLength());
							line.setAttribute("class", "indel");
							g.appendChild(line);
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
						final StringBuilder sb=new StringBuilder();
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
						final Element path = super.path();
						path.setAttribute("d",sb.toString().trim());
						if( ce.getOperator()==CigarOperator.H || ce.getOperator()==CigarOperator.S)
							{
							path.setAttribute("style","fill:yellow;");
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
							path.setAttribute("style","fill:rgb("+rgb[0]+","+rgb[1]+","+rgb[2]+")");
							}
						
						if(op.consumesReadBases())
							{
							for(int i=0;i< ce.getLength();++i)
								{
								final int pos_i = refPos+i;
								boolean pos_in_interval =  (this.interval.getStart() <= pos_i && pos_i<=this.interval.getEnd() );
								char ca=(char)bases[readPos];
								if(!in_clip)
									{
									Counter<Character> count= consensus.get(pos_i) ;
									if(count==null)
										{
										count=new Counter<>();
										consensus.put(refPos+i,count) ;
										}
									count.incr(ca);
									}
								char cb= genomicSequence.charAt(refPos+i-1);
								
								/* mutation ? */
								if(this.print_bases && cb!='N' && ca!='N' && cb!=ca && !in_clip && pos_in_interval)
									{
									final Element use = super.element("use");
									use.setAttribute("x",format( baseToPixel(refPos+i)));
									use.setAttribute("href", "#mut");
									
									}
								
								if( pos_in_interval && 
									((this.print_bases&&context.featureWidth>5) || (cb!=ca ))
									 )
									{
									final Element use = super.element("use");
									use.setAttribute("x",format( baseToPixel(refPos+i)));
									//w.setAttribute("y",format(y_top));never needed
									use.setAttribute("href", "#b"+ca);
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
				
				
				final double x= baseToPixel(pos);
				final Element line = super.line(x, x, y_top5, y_bot5);
				setTitle(line, "Insertion "+insertion);
				line.setAttribute("class", "insert");
				g.appendChild(line);
				}

			
			return g;
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
				final boolean is_phased = IlluminaConstants.isDragenProximityEnabled(header);
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
						if(is_phased) {
							
							}
						
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
			final String filename= this.prefix+buildName+(buildName.isEmpty()?"":".")+ IOUtils.escapePath( this.interval.getContig())+"_"+this.interval.getStart()+"_"+this.interval.getEnd()+"."+context.sampleName+".svg";
			

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
				
				return 0;
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				return -1;
				}
			
			}
		
	public void setDrawinAreaWidth(int drawinAreaWidth) {
		this.drawinAreaWidth = drawinAreaWidth;
		}
		
	public void build(final SamReader samReader,final Path fastaReference, final Locatable loc) {
		this.fastaReference = Objects.requireNonNull(fastaReference, "null reference");
		try 
			{
			this.indexedFastaSequenceFile = indexedFastaSequenceFile;
			this.referenceDict = SequenceDictionaryUtils.extractRequired(indexedFastaSequenceFile);
			Objects.requireNonNull(this.indexedFastaSequenceFile, "null reference");
			Objects.requireNonNull(this.referenceDict, "null dict");
			super.removeAllChildren(super.rootElement);
			this.interval =Objects.requireNonNull(loc, "null interval");
			if(this.referenceDict.getSequence(loc.getContig())==null) {
				throw new JvarkitException.ContigNotFoundInDictionary(loc.getContig(), this.referenceDict);
				}
			
			}
		finally
			{
			
			}
	}
