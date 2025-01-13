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
package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.IntToDoubleFunction;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.transform.stream.StreamResult;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.ArrayResizer;
import com.github.lindenb.jvarkit.math.Average;
import com.github.lindenb.jvarkit.math.Median;
import com.github.lindenb.jvarkit.math.RunMedian;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.net.UrlSupplier;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalExtender;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.iterator.LineIterators;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
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
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC
## input

input a set of bam/cram files or one file with the suffix '.list' containing the path to the bams

## output

output is a HTML+SVG file

## example:

```
find dir -type f -name "*bam" > in.list 
java -jar dist/jvarkit.jar coverageplotter -R src/test/resources/rotavirus_rf.fa --region "RF01:100-200" in.list
```

## Screenshot

!(https://pbs.twimg.com/media/Fac3XR3aAAEJoXu?format=jpg&name=medium)[https://twitter.com/yokofakun/status/1560276675614887937]


END_DOC 
 */
@Program(
	name="coverageplotter",
	description="Display an image of depth to display any anomaly an intervals+bams",
	keywords={"cnv","bam","depth","coverage","svg"},
	creationDate="20200605",
	modificationDate="20241009",
	biostars = 9536274,
	jvarkit_amalgamion =  true,
	menu="CNV/SV"
	)
public class CoveragePlotter extends Launcher {
	private static final Logger LOG = Logger.build( CoveragePlotter.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refPath = null;
	@Parameter(names={"--regions","--region","--interval"},description = "Interval region",required=true)
	private String intervalStr=null;
	@Parameter(names={"--mapq"},description = "min mapping quality")
	private int min_mapq=1;
	@Parameter(names={"--max-depth"},description = "ignore position if depth > 'x'")
	private int max_depth=500;
	@Parameter(names={"--dimension","--dim"},description = "Image Dimension. " + DimensionConverter.OPT_DESC, converter=DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,300);
	@Parameter(names={"--extend","-x"},description = IntervalExtender.OPT_DESC)
	private String extendStr="3.0";
	@Parameter(names= {"--gff","--gff3"},description="Optional Tabix indexed GFF3 file. Will be used to retrieve an interval by gene name, or to display gene names in a region.")
	private Path gtfFile = null;
	@Parameter(names= {"--known"},description="Optional Tabix indexed Bed or VCF file containing known CNV. Both types must be indexed.")
	private Path knownCnvFile = null;
	@DynamicParameter(names = "-D", description = "style. Undocumented.",hidden=true)
	private Map<String, String> dynaParams = new HashMap<>();
	@Parameter(names= {"--include-center"},description="When calculating the median depth, also consider the original user's region, not only the extended interval.")
	private boolean include_original_interval_for_median = false;
	@Parameter(names= {"--ignore-known-containing"},description="Ignore known CNV containing the whole region (prevent large known CNV to be displayed) ")
	private boolean ignore_cnv_overlapping = false;
	@Parameter(names= {"--svg-only"},description="Force SVG-only output (default is HTML+SVG).")
	private boolean force_svg_output = false;
	@Parameter(names= {"--max-y"},description="Max normalized Y")
	private double max_normalized_y = 3.0;
	@Parameter(names= {"--css"},description="Custom CSS file. format <sample> <css>. One per line. eg. \"sample1 stroke:red;")
	private Path customCss = null;
	@Parameter(names= {"--smooth"},description="Run median smooth on this number of pixels. (ignore if <=1)")
	private int smooth_pixel_size = 10;
	@Parameter(names= {"--loess"},description="Run Loess smoothing on GC%. Experimental. For now, I find the smooting is too strong.")
	private boolean run_loess = false;
	@Parameter(names= {"--use-average"},description="Calculating the median depth can be memory consumming for large regions. If the region is larger than 'x', use 'average' instead of 'median'. "+ DistanceParser.OPT_DESCRIPTION,splitter=NoSplitter.class,converter = DistanceParser.StringConverter.class )
	private int use_average_instead_of_median_length_treshold = 2_000_000;


	
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");

	private static class SampleInfo {
		String sample;
		Path path;
		long nReads = 0L;
		Point2D maxPosition=null;
		String rgb=null;
		double meanDP=0.0;
		
		double medianAll=0.0;
		double medianInner=0.0;
		double medianOuter=0.0;
		}
	
	private List<Interval> getKnownCNVs(final Locatable region) {
		if(this.knownCnvFile==null) return Collections.emptyList();
		final String fname=this.knownCnvFile.getFileName().toString();

		final Predicate<Interval> rejectCnv = cnv->(this.ignore_cnv_overlapping && cnv.getStart() < region.getStart() && cnv.getEnd() > region.getEnd());
		final List<Interval> knowns = new ArrayList<>();
		if(fname.endsWith(".bed.gz")) {
			try(TabixReader tbr = new TabixReader(this.knownCnvFile.toString())) {
				final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
				final String ctg = cvt.apply(region.getContig());
					if(!StringUtils.isBlank(ctg)) {
						final BedLineCodec codec = new BedLineCodec();
						final TabixReader.Iterator iter = tbr.query(ctg,region.getStart(), region.getEnd());
						for(;;) {
							final String line = iter.next();
							if(line==null) break;
							final BedLine bed = codec.decode(line);
							if(bed==null) continue;
							final Interval rgn = new Interval(region.getContig(),bed.getStart(),bed.getEnd(),false,bed.getOrDefault(3, ""));
							if(rejectCnv.test(rgn)) continue;
							knowns.add(rgn);
							}
					}
				}
			catch(final Throwable err) {
				LOG.error(err);
				}
			}
		else if(FileExtensions.VCF_LIST.stream().anyMatch(X->fname.endsWith(X))) {
			try(VCFReader vcfFileReader= VCFReaderFactory.makeDefault().open(this.knownCnvFile,true)) {
				final ContigNameConverter cvt = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(vcfFileReader.getHeader()));
				final String ctg = cvt.apply(region.getContig());
				if(!StringUtils.isBlank(ctg)) {
					vcfFileReader.query(ctg, region.getStart(), region.getEnd()).
							stream().
							filter(VC->!VC.isSNP()).
							forEach(VC->{
								final List<String> list = new ArrayList<>();
								if(VC.hasID()) list.add(VC.getID());
								if(VC.hasAttribute(VCFConstants.SVTYPE))  list.add(VC.getAttributeAsString(VCFConstants.SVTYPE,"."));
								final Interval rgn= new Interval(region.getContig(),VC.getStart(),VC.getEnd(),false,String.join(";",list));
								if(rejectCnv.test(rgn)) return;
								knowns.add(rgn);
								});
					}
				}
			catch(final Throwable err) {
				LOG.error(err);
				}
			}
		else
			{
			LOG.warn("not a vcf of bed.gz file "+this.knownCnvFile);
			}
		return knowns;
		}
	
	private void drawKnownCnv(final XMLStreamWriter w,final Rectangle rectangle,final Locatable region) {
		final List<Interval> knows = getKnownCNVs(region);

		final Pileup<Interval> pileup = new Pileup<>();
		pileup.addAll(knows);

		if(!pileup.isEmpty() ) {
			final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*rectangle.getWidth();
			final double featureHeight = 10.0/pileup.getRowCount();
			for(int row=0;row< pileup.getRowCount();++row) {
				for(final Interval cnv:pileup.getRow(row)) {
					final double y= rectangle.getHeight()-12.0 + row*featureHeight;
					final double x1 = Math.max(rectangle.getX(),position2pixel.applyAsDouble(cnv.getStart()));
					final double x2 = Math.min(rectangle.getMaxX(),position2pixel.applyAsDouble(cnv.getEnd()));
					try {
						w.writeStartElement("rect");
						w.writeAttribute("class", "cnv");
						w.writeAttribute("x", format(x1));
						w.writeAttribute("y", format(y-1));
						w.writeAttribute("width", format(Math.max(0.5, x2-x1)));
						w.writeAttribute("height", format(featureHeight*0.9));
						title(w, "Known: "+new SimpleInterval(cnv).toNiceString()+" "+cnv.getName());
						w.writeEndElement();
						}
					catch(final Throwable err) {
						LOG.warn(err);
						}
					}
				}
			}
		}

	

	private Stream<Gff3Feature> getGenes(final Locatable region,final Predicate<Gff3Feature> filter) {
		if(this.gtfFile==null) return Stream.empty();
		TabixReader tbr = null;
		try {
			tbr= new TabixReader(this.gtfFile.toString());
			
			final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
			final String ctg = cvt.apply(region.getContig());
			if(StringUtils.isBlank(ctg)) {
				tbr.close();
				return Stream.empty();
				}
			
			final Gff3Codec codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
			final TabixReader.Iterator iter=tbr.query(ctg,region.getStart(),region.getEnd());
			final TabixReader tbrfinal = tbr;
			final AbstractIterator<Gff3Feature> iter2= new AbstractIterator<Gff3Feature>() {
				@Override
				protected Gff3Feature advance() {
					try {
						for(;;) {
							final String line = iter.next();
							if(line==null) return null;
							if(StringUtils.isBlank(line) ||  line.startsWith("#")) continue;
							final String tokens[]= CharSplitter.TAB.split(line);
							if(tokens.length<9 ) continue;
							tokens[0]=region.getContig();
							final LineIterator lr = LineIterators.of(new String[] {line});
							final Gff3Feature gffline = codec.decode(lr);
							if(gffline==null) continue;
							if(filter!=null && !filter.test(gffline)) continue;
							return gffline;
							}
						} catch (final IOException e) {
						LOG.error(e);
						return null;
						}
					}
				};
			return StreamSupport.stream(new IterableAdapter<Gff3Feature>(iter2).spliterator(),false).
					onClose(()->{ tbrfinal.close(); });
			}
		catch(Throwable err) {
			if(tbr!=null) tbr.close();
			return Stream.empty();
			}
	}
	
	
	private void drawGenes(final XMLStreamWriter w,final Rectangle rect,final Locatable region) {
		final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*rect.getWidth();
		final double y= rect.getMaxY()-4.0;
			final int geneSize = 10;
	
			getGenes(region,G->true).forEach(gtfline->{
					final double x1 = Math.max(rect.getX(),position2pixel.applyAsDouble(gtfline.getStart()));
					final double x2 = Math.min(rect.getMaxX(),position2pixel.applyAsDouble(gtfline.getEnd()));
					try {
					final String geneName= gtfline.getAttribute("gene_name").stream().findFirst().orElse("");
					if(gtfline.getType().equals("gene") ) {
						w.writeEmptyElement("text");
						w.writeAttribute("class","gene");
						w.writeAttribute("x", format(Math.max(x1,1)));
						w.writeAttribute("y", format(y-(geneSize+3)));
						w.writeCharacters(geneName);
						}
					else if(gtfline.getType().equals("exon") ) {
						w.writeStartElement("rect");
						w.writeAttribute("class","gene");
						w.writeAttribute("x", format(x1));
						w.writeAttribute("y", format(y-1));
						w.writeAttribute("width", format(x2-x1));
						w.writeAttribute("height", format(3));
						title(w,"exon "+geneName);
						w.writeEndElement();
						}
					else if(gtfline.getType().equals("transcript") ) {
						w.writeStartElement("line");
						w.writeAttribute("class","gene");
						w.writeAttribute("x1", format(x1));
						w.writeAttribute("y1", format(y));
						w.writeAttribute("x2", format(x2));
						w.writeAttribute("y2", format(y));
						title(w,"transcript "+geneName);
						w.writeEndElement();
						}
					} catch(final Throwable err) {
						LOG.warn(err);
					}
				});
		}
	
	
private void title(XMLStreamWriter w,final Object o) throws XMLStreamException{
	if(o==null) return ;
	final String s= String.valueOf(o);
	if(StringUtils.isBlank(s)) return ;
	w.writeStartElement("title");
	w.writeCharacters(s);
	w.writeEndElement();
	}
	
private String format(double v) {
	return this.decimalFormater.format(v);
	}
@Override
public int doWork(final List<String> args) {
	try
		{
		if(this.max_normalized_y < 2.0) {
			System.err.println("max Y should be >= 2.0");
			return -1;
			}
		if(this.smooth_pixel_size>1 && this.run_loess) {
			LOG.warn("loess: ignoring smoothing using running median size:"+this.smooth_pixel_size);
			}
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refPath);
		final Hyperlink hyperlink = Hyperlink.compile(dict);
		final SamReaderFactory samReaderFactory = SamReaderFactory.
					makeDefault().
					referenceSequence(CoveragePlotter.this.refPath).
					validationStringency(ValidationStringency.LENIENT)
					;
		final Map<String,String> sample2css;
		if(this.customCss!=null) {
			sample2css  = Files.lines(this.customCss).
					filter(S->!StringUtils.isBlank(S)).
					filter(S->!S.startsWith("#")).
					map(S->S.trim().split("\\s+", 2)).
					collect(Collectors.toMap(T->T[0],T->T[1]));
			}
		else
			{
			sample2css = Collections.emptyMap();
			}
		
		final List<Path> inputBams =  IOUtils.unrollPaths(args);
		
		if(inputBams.isEmpty()) {
			LOG.error("input bam file missing.");
			return -1;
			}
		final List<SampleInfo> samples = new ArrayList<>(inputBams.size());
		final Locatable the_locatable = new IntervalParser(dict).
			apply(this.intervalStr).orElseThrow(IntervalParser.exception("Cannot parse interval:\""+this.intervalStr+"\""));
		
		final SimpleInterval rawRegion = new SimpleInterval(the_locatable);
		final SimpleInterval extendedRegion = IntervalExtender.of(dict, this.extendStr).apply(rawRegion);
		if(!extendedRegion.contains(rawRegion)) {
			LOG.error("extended interval "+extendedRegion+" should contain the user region "+rawRegion);
			return -1;
			}
			
		if(!include_original_interval_for_median && extendedRegion.equals(rawRegion)) {
				LOG.error("extended region "+extendedRegion+" is same as raw region "+ rawRegion +" but median is only using extend regions. See options.");
			return -1;
			}
		
		
		final double[] gc_percent;
		if(run_loess) {
			/** extract Fastq sequence, and calculate GC% using sliding window */
			try {
				final Average average = new Average();
				final byte[] atgc;
				gc_percent = new double[extendedRegion.getLengthOnReference()];
				try(ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.refPath)) {
					atgc = ref.getSubsequenceAt(extendedRegion.getContig(), extendedRegion.getStart(), extendedRegion.getEnd()).getBases();
					if(atgc.length!=extendedRegion.getLengthOnReference()) {
						LOG.error("cannot retrieve ATGC sequence for "+extendedRegion);
						return -1;
						}
					}
				final int window_gc_size = 30;
				average.reset();
				for(int i=0;i< atgc.length;i++) {
					for(int j=Math.max(0, i-window_gc_size);j < Math.min(i+window_gc_size, atgc.length);++j) {
						switch(atgc[j]) {
							case 'g':case 'G': case 'c':case 'C': case 's':case 'S': average.accept(1.0); break;
							default:average.accept(0.0); break;
							}
						}
					gc_percent[i] = average.getAsDouble();
					}
				}
			catch(final IOException err) {
				LOG.error(err);
				return -1;
				}
			}
		else
			{
			gc_percent = null;
			}
		
		final int sampleFontSize=7;
		final double max_y = this.max_normalized_y;
		final ToDoubleFunction<Double> normToPixelY = NORM->  this.dimension.getHeight() * (1.0 -  (NORM/max_y));
		final double y_mid = normToPixelY.applyAsDouble(1.0);
		
		final ToDoubleFunction<Integer> pos2pixel = POS-> (POS - extendedRegion.getStart())/(double)extendedRegion.getLengthOnReference() * this.dimension.getWidth();
		final XMLOutputFactory xof = XMLOutputFactory.newFactory();
		final XMLStreamWriter w;
		boolean output_html = true;
		
		if(this.outputFile==null) {
			w = xof.createXMLStreamWriter(stdout(), "UTF-8");
			}
		else
			{
			w = xof.createXMLStreamWriter(new StreamResult(this.outputFile));
			output_html = !this.outputFile.getName().toLowerCase().endsWith(".svg");
			}
		
		if(this.force_svg_output || (this.outputFile!=null && this.outputFile.getName().toLowerCase().endsWith(".svg"))) {
			output_html = false;
			}
		
		if(output_html) {
			w.writeStartElement("html");
			w.writeStartElement("head");
			
			w.writeEmptyElement("meta");
				w.writeAttribute("http-equiv", "Content-Type");
				w.writeAttribute("content", "text/html; charset=utf-8");
			w.writeEmptyElement("meta");
				w.writeAttribute("http-equiv", "author");
				w.writeAttribute("content", "Pierre Lindenbaum Phd ");

			
			w.writeStartElement("title");
			w.writeCharacters(rawRegion.toNiceString()+" extended to "+extendedRegion.toNiceString());
			w.writeEndElement();//title
			w.writeStartElement("style");
			w.writeCharacters(""
				+ "table{font-family:Arial,Verdana,sans-serif;color:darkgray;font-size:14px;border-collapse:collapse;border:1px solid black;padding:5px;margin:5px;}"
				);
			w.writeEndElement();//title
			
			w.writeStartElement("script");
			// https://stackoverflow.com/questions/40108280/how-to-hide-svg-elements-with-javascript-on-a-button-click
			w.writeCharacters(""
					+"function toggleSample(s) {var E=document.getElementById(s);if(E!=null) {var D=E.style.display ;if(D == \"block\") {E.style.display=\"none\";} else {E.style.display=\"block\";} } return true;}");
			w.writeEndElement();//script
			
			w.writeEndElement();//head
			
			w.writeStartElement("body");
			w.writeStartElement("h2");
			w.writeCharacters(SequenceDictionaryUtils.getBuildName(dict).orElse("") + " "+ rawRegion.toNiceString()+" extended to "+extendedRegion.toNiceString());
			w.writeEndElement();//h2
			w.writeStartElement("div");
			//placeholder for xslt
			w.writeStartElement("div");
			w.writeAttribute("id", "form");
			w.writeEndElement();
			}
		else {
			w.writeStartDocument("UTF-8","1.0");
			}
		
		final Insets margin = new Insets(100,100, 100, 50);
		
		w.writeStartElement("svg");
		w.writeDefaultNamespace(SVG.NS);
		w.writeAttribute("width", format(this.dimension.width + margin.left + margin.right));
		w.writeAttribute("height", format(this.dimension.height + margin.top + margin.bottom));
		
		/* BEGIN metadata */
		w.writeStartElement("metadata");
		w.writeAttribute("id", "metadata");

		w.writeEndElement();//metadata
		/* EBND metadata */
		
		w.writeStartElement("style");
		w.writeCharacters(""
				+ "svg {stroke-linecap:round;}\n"
				+ ".geneh {fill:rgb(240,240,240);stroke:none;}\n"
				+ ".ymedian0 {fill:none;stroke:blue;stroke-dasharray:5,5;}\n"
				+ ".ymedian1 {fill:none;stroke:cyan;stroke-dasharray:5,5;}\n"
				+ ".vline {fill:none;stroke:gray;stroke-dasharray:5,5;}\n"
				+ ".frame {fill:none;stroke:darkgray;}\n"
				+ ".rgnbound {fill:none;stroke:green;}\n"
				+ ".hist {fill:none;stroke:darkgray;display:block;}\n"
				+ ".title1 {fill:darkgray;stroke:none;dominant-baseline:middle;text-anchor:middle;}\n"
				+ ".title2 {fill:darkgray;stroke:none;}\n"
				+ ".labely {fill:darkgray;stroke:none;dominant-baseline:middle;text-anchor:end;}\n"
				+ ".cnv {fill:blue;stroke:orange;opacity:0.5;}\n"
				+ ".gene {fill:green;stroke:orange;opacity:0.5;}\n" +
				sample2css.entrySet().stream().map(KV->"."+KV.getKey()+" {"+KV.getValue()+"}\n").collect(Collectors.joining("\n"))
				);
		
		w.writeEndElement();//style
		
		w.writeStartElement("g");
		w.writeStartElement("g");
		w.writeAttribute("transform", "translate("+format(margin.left)+","+format(margin.top)+")");

		//draw exons
		getGenes(extendedRegion,G->G.getType().equals("exon")).forEach(EX->{
			try {
				final double x1 = Math.max(0, pos2pixel.applyAsDouble(EX.getStart()));
				final double x2 = Math.min(this.dimension.getWidth(), pos2pixel.applyAsDouble(EX.getEnd()));

				w.writeStartElement("rect");
				w.writeAttribute("class", "geneh");
				w.writeAttribute("x", format(x1));
				w.writeAttribute("y", format(0));
				w.writeAttribute("width", format(Math.max(0.5, x2-x1)));
				w.writeAttribute("height", format(this.dimension.height));
				title(w,EX.getAttribute("gene_name").stream().map(GN->"exon "+GN).findFirst().orElse(""));
				w.writeEndElement();
				}
			catch(final Throwable err)  {
				LOG.warn(err);
				}
			});
		
		// draw mid -y
		double y = normToPixelY.applyAsDouble(1.0);
		w.writeStartElement("line");
		w.writeAttribute("class", "ymedian0");
		w.writeAttribute("x1", format(0));
		w.writeAttribute("x2", format(this.dimension.width));
		w.writeAttribute("y1", format(y));
		w.writeAttribute("y2", format(y));
		title(w,"median");
		w.writeEndElement();

		// draw 0.75 y line
		y = normToPixelY.applyAsDouble(1.5);
		w.writeStartElement("line");
		w.writeAttribute("class", "ymedian1");
		w.writeAttribute("x1", format(0));
		w.writeAttribute("x2", format(this.dimension.width));
		w.writeAttribute("y1", format(y));
		w.writeAttribute("y2", format(y));
		title(w,"median * 1.5");
		w.writeEndElement();
		
		// draw 0.25 y line
		y = normToPixelY.applyAsDouble(0.5);
		w.writeStartElement("line");
		w.writeAttribute("class", "ymedian1");
		w.writeAttribute("x1", format(0));
		w.writeAttribute("x2", format(this.dimension.width));
		w.writeAttribute("y1", format(y));
		w.writeAttribute("y2", format(y));
		title(w,"median * 0.5");
		w.writeEndElement();

		// draw frame
		w.writeStartElement("rect");
		w.writeAttribute("class", "frame");
		w.writeAttribute("x", format(0));
		w.writeAttribute("y", format(0));
		w.writeAttribute("width", format(this.dimension.width-1));
		w.writeAttribute("height", format(this.dimension.height-1));
		w.writeEndElement();

		// draw y axis
		w.writeStartElement("g");
		for(double v=0; v <= max_normalized_y;v+=0.25) {
			y= normToPixelY.applyAsDouble(v);
			w.writeStartElement("g");
			w.writeStartElement("line");
			w.writeAttribute("class", "ymedian1");
			w.writeAttribute("x1", format(-5));
			w.writeAttribute("x2", format(0));
			w.writeAttribute("y1", format(y));
			w.writeAttribute("y2", format(y));
			title(w,format(v));
			w.writeEndElement();
			
			w.writeStartElement("text");
			w.writeAttribute("class","labely");
			w.writeAttribute("x", format(-6));
			w.writeAttribute("y", format(y));
			w.writeCharacters(format(v));
			w.writeEndElement();//end
			
			w.writeEndElement();//g
			}
		w.writeStartElement("text");
		w.writeAttribute("style", "text-anchor:middle;");
		w.writeAttribute("x","0");
		w.writeAttribute("y","0");
		w.writeAttribute("transform", "translate("+format(-margin.left/2)+","+format(dimension.height/2.0)+") rotate(90)");
		w.writeCharacters("Normalized Depth.");
		w.writeEndElement();//end
		
		w.writeEndElement(); //g
		
		// draw x axis
		w.writeStartElement("g");
		double prev_x=0;
		for(int i=extendedRegion.getStart();i<=extendedRegion.getEnd();i++) {
			final double x0 = ((i-extendedRegion.getStart())/(double)extendedRegion.getLengthOnReference())*dimension.width;
			if(i>extendedRegion.getStart() && x0-prev_x > 100) {
				
				w.writeStartElement("line");
				w.writeAttribute("class", "vline");
				w.writeAttribute("x1", format(x0));
				w.writeAttribute("x2", format(x0));
				w.writeAttribute("y1", format(0));
				w.writeAttribute("y2", format(dimension.height));
				title(w,StringUtils.niceInt(i));
				w.writeEndElement();
				
				
				w.writeStartElement("text");
				w.writeAttribute("class","title2");
				w.writeAttribute("x","0");
				w.writeAttribute("y","0");
				w.writeAttribute("transform", "translate("+format(x0)+","+format(dimension.height+5)+") rotate(45)");
				w.writeCharacters(StringUtils.niceInt(i));
				w.writeEndElement();//end

				
				prev_x = x0;
			}
		}
		
		w.writeStartElement("text");
		w.writeAttribute("style", "text-anchor:middle;");
		w.writeAttribute("x","0");
		w.writeAttribute("y","0");
		w.writeAttribute("transform", "translate("+format(dimension.width/2)+","+format(dimension.height + margin.bottom/2.0)+")");
		w.writeCharacters("Genomic location");
		w.writeEndElement();//end
		w.writeEndElement();//g
		
		//draw genes
		drawGenes(w, new Rectangle(0,0,dimension.width,dimension.height), extendedRegion);
		// draw knowns
		drawKnownCnv(w, new Rectangle(0,0,dimension.width,dimension.height), extendedRegion);
		
		// draw original region vertical bounds 
		if(!extendedRegion.equals(rawRegion)) {
			
			for(int i=0;i<2;i++) {
				final int p = (i==0?rawRegion.getStart():rawRegion.getEnd());
				final double x = pos2pixel.applyAsDouble(p);
				w.writeStartElement("line");
				w.writeAttribute("class", "rgnbound");
				w.writeAttribute("x1", String.valueOf(x));
				w.writeAttribute("x2", String.valueOf(x));
				w.writeAttribute("y1", String.valueOf(0));
				w.writeAttribute("y2", String.valueOf(dimension.height));
				title(w,p);
				w.writeEndElement();
				}
			}
		
		final double depth[]= new double[extendedRegion.getLengthOnReference()];		
		for(final Path path: inputBams) {
			
			try(SamReader sr = samReaderFactory.open(path)) {
				if(!sr.hasIndex()) {
					LOG.error("index missng for "+path);
					return -1;
					}
				final SAMFileHeader header= sr.getFileHeader();
				final SampleInfo sampleInfo = new SampleInfo();
				sampleInfo.sample = header.getReadGroups().stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtil.isBlank(S)).
						findFirst().
						orElse(IOUtils.getFilenameWithoutCommonSuffixes(path));
				sampleInfo.path = path;
				
				LOG.info(sampleInfo.sample+" "+(1+samples.size())+"/"+inputBams.size()+ " "+extendedRegion);
				
				
				SequenceUtil.assertSequenceDictionariesEqual(dict,header.getSequenceDictionary());
				Arrays.fill(depth, 0);
				try(CloseableIterator<SAMRecord> siter = sr.queryOverlapping(extendedRegion.getContig(), extendedRegion.getStart(), extendedRegion.getEnd())) {
					while(siter.hasNext()) {
						final SAMRecord rec= siter.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(!SAMRecordDefaultFilter.accept(rec, this.min_mapq)) continue;
						int ref=rec.getStart();
						final Cigar cigar = rec.getCigar();
						if(cigar==null) continue;
						
						sampleInfo.nReads++;
						
						int max_end = extendedRegion.getEnd();
						// overlaping paired reads
						if(rec.getReadPairedFlag() &&
							!rec.getMateUnmappedFlag() &&
							rec.getReferenceIndex().equals(rec.getMateReferenceIndex()) &&
							rec.getStart() < rec.getMateAlignmentStart() &&
							rec.getEnd() > rec.getMateAlignmentStart())
							{
							max_end = Math.min(max_end, rec.getMateAlignmentStart());
							}
						
						for(CigarElement ce:cigar) {
							final CigarOperator op = ce.getOperator();
							final int len = ce.getLength();
							if(op.consumesReferenceBases()) {
								if(op.consumesReadBases()) {
									for(int i=0;i< len;i++) {
										final int pos = ref+i;
										if(pos < extendedRegion.getStart()) continue;
										if(pos > max_end) break;
										depth[pos-extendedRegion.getStart()]++;
										}
									}
								ref+=len;
								}
							}// loop cigar
						} // end samIter
					}// 
			
			
			if(this.run_loess) {
				final Map<Double, Average> gc2depth = new HashMap<>();
				for(int x=0;x < gc_percent.length;++x) {
					Average avg = gc2depth.get(gc_percent[x]);
					if(avg==null) {
						avg = new Average();
						gc2depth.put(gc_percent[x], avg);
						}
					avg.accept(depth[x]);
					}
				final List<Map.Entry<Double, Double>> L = gc2depth.entrySet().
						stream().
						map(KV->new AbstractMap.SimpleEntry<Double,Double>(KV.getKey(), KV.getValue().getAsDouble())).
						sorted((A,B)->A.getKey().compareTo(B.getKey())).
						collect(Collectors.toList());
				final double[] xval = L.stream().mapToDouble(KV->KV.getKey()).toArray();
				final double[] yval = L.stream().mapToDouble(KV->KV.getValue()).toArray();
				
				final LoessInterpolator loessInterpolator = new LoessInterpolator();
				final PolynomialSplineFunction pf= loessInterpolator.interpolate(xval, yval);
				for(int i=0;i< depth.length;i++) {
					final double v= pf.value(gc_percent[i]);
					if(v==0.0) continue;
					depth[i] = (depth[i]/v);
					}
				}
			else if(smooth_pixel_size>1) {
				int smooth_n_bases= (int)((smooth_pixel_size/dimension.getWidth())*extendedRegion.getLengthOnReference());
				if(smooth_n_bases%2==0) smooth_n_bases++;
				if(smooth_n_bases>1) {
					final double[] copy = new RunMedian(smooth_n_bases).apply(depth);
					for(int i=0;i< copy.length;i++) depth[i]=copy[i];
					}
				}
			
		/** get median */
			{
			final Median myMedianAll = new Median(extendedRegion.getLengthOnReference());
			final Median myMedianOuter = new Median(Math.max(100, extendedRegion.getLengthOnReference()-rawRegion.getLengthOnReference()));
			final Median myMedianInner = new Median(rawRegion.getLengthOnReference());
			final Average averageDepth = new Average();
			
			for(int i=0;i< depth.length;i++) {
				if(depth[i]>this.max_depth) continue;
				final double dp = depth[i];
				myMedianAll.accept(dp);
				final int pos = extendedRegion.getStart() + i;
				boolean use_for_dp=false;
				if(rawRegion.getStart()<=pos && pos<=rawRegion.getEnd()) {
					myMedianInner.accept(dp);
					if(this.include_original_interval_for_median) use_for_dp=true;
					}
				else
					{
					myMedianOuter.accept(dp);
					use_for_dp = true;
					}
				if(use_for_dp) {
					averageDepth.accept(depth[i]);
					}
				}
			if(averageDepth.getCount()==0L) {
				final String msg = "Skipping "+sampleInfo.sample +" "+extendedRegion+" because I cannot find any depth";
				LOG.warning(msg);
				continue;
				}
			sampleInfo.meanDP = averageDepth.getAsDouble();
			/** prevent long calculations */
			if(depth.length < this.use_average_instead_of_median_length_treshold ) {
				sampleInfo.medianInner = myMedianInner.get().orElse(0.0);
				sampleInfo.medianOuter = myMedianOuter.get().orElse(0.0);
				sampleInfo.medianAll = myMedianAll.get().orElse(0.0);
				}
			else
				{
				LOG.warn("big array ("+depth.length+") using average instead of median.");
				sampleInfo.medianInner = myMedianInner.getAverage().orElse(0.0);
				sampleInfo.medianOuter = myMedianOuter.getAverage().orElse(0.0);
				sampleInfo.medianAll = myMedianAll.getAverage().orElse(0.0);
				}
			}
		/** end median */

			
			final double medianDP = (this.include_original_interval_for_median?sampleInfo.medianAll:sampleInfo.medianOuter);
			if(medianDP <= 0.0) {
				final String msg = "Skipping "+sampleInfo.sample +" "+extendedRegion+" because median is 0";
				LOG.warning(msg);
				continue;
				}
			
			
			// normalize on median
			double[] normArray = new double[depth.length];
			for(int i=0;i< depth.length;i++) {
				normArray[i] = depth[i]/medianDP;
				}
			normArray= new ArrayResizer().resizeToDouble(normArray, dimension.width);

			w.writeCharacters("\n");
			w.writeComment(sampleInfo.sample);
			w.writeStartElement("g");
			
			Point2D max_position=null;
			double max_distance_to_1=0.0;
			final List<Point2D> points = new ArrayList<>(normArray.length);
			for(int i=0;i< normArray.length;i++) {
				final double normDepth = normArray[i];
				final double y2 = normToPixelY.applyAsDouble(normDepth);
				double distance_to_1 = Math.abs(normDepth-1.0);
				if(distance_to_1 > 0.3 && (max_position==null || distance_to_1 > max_distance_to_1)) {
					max_distance_to_1 = distance_to_1;
					max_position = new Point2D.Double(i,y2);
					}
				points.add(new Point2D.Double( i,y2));
				}

			//simplify array
			int x=1;
			while(x+1<points.size()) {
				if(points.get(x-1).getY()==points.get(x).getY() && points.get(x).getY()==points.get(x+1).getY()) {
					points.remove(x);
					}
				else
					{
					x++;
					}
				}

			
			w.writeStartElement("polyline");
			w.writeAttribute("id", sampleInfo.sample);
			
			final Color sampleRGB =  Color.getHSBColor((float) (samples.size()) / (float)inputBams.size(), 0.85f, 1.0f);
			sampleInfo.rgb = "rgb("+sampleRGB.getRed()+","+sampleRGB.getGreen()+","+sampleRGB.getBlue()+")";

			
			if(sample2css.containsKey(sampleInfo.sample)) {
				w.writeAttribute("class",sampleInfo.sample);
				w.writeAttribute("style", "display:block;");
				
				final String rgb2 = Arrays.stream(sample2css.get(sampleInfo.sample).split("[{;}]")).
						map(S->S.trim()).
						filter(S->S.startsWith("stroke:")).
						map(S->S.substring(7).trim()).
						findFirst().orElse(null);
				if(!StringUtil.isBlank(rgb2)) {
					sampleInfo.rgb = rgb2;
					}
				}
			else
				{
				w.writeAttribute("class","hist");
				w.writeAttribute("style", "display:block;stroke:" + sampleInfo.rgb +";");
				}
			w.writeAttribute("points",
					points.stream().map(P->format(P.getX())+","+format(P.getY())).collect(Collectors.joining(" "))
					);
			title(w,sampleInfo.sample);
			w.writeEndElement(); //polyline

			samples.add(sampleInfo);
			
			sampleInfo.maxPosition = max_position;
			w.writeEndElement(); //g
			}//end samReader
		} // end loop over paths/samples
				
				
				
		
		final String label_samples = samples.size()>10?"N="+StringUtils.niceInt(samples.size()):samples.stream().map(S->S.sample).collect(Collectors.joining(";"));
		final Set<String> all_genes = getGenes(extendedRegion,G->G.getType().equals("gene")).
				flatMap(G->G.getAttribute("gene_name").stream()).
				filter(F->!StringUtils.isBlank(F)).
				collect(Collectors.toCollection(TreeSet::new))
				;


		if(Boolean.parseBoolean(this.dynaParams.getOrDefault("min-max-samples", "true"))) {
			/** draw sample names */
			for(final SampleInfo si:samples) {
				final Point2D pt = si.maxPosition;
				if(pt==null) continue;
				double sny = pt.getY();
				if(sny>y_mid) sny+=sampleFontSize;
				
				w.writeStartElement("text");
				w.writeAttribute("class","title2");
				w.writeAttribute("x", format(pt.getX()));
				w.writeAttribute("y", format(Math.min(this.dimension.height-sampleFontSize,Math.max(sampleFontSize,sny))));
				w.writeCharacters(si.sample);
				w.writeEndElement();//end
				}
			}
		
		
		
		
		w.writeEndElement();//g
		
		w.writeStartElement("text");
		w.writeAttribute("class","title1");
		w.writeAttribute("x", format(margin.left+dimension.width/2.0));
		w.writeAttribute("y", format(margin.top/2));
		w.writeCharacters(
				extendedRegion.toNiceString()+" Length:"+StringUtils.niceInt(extendedRegion.getLengthOnReference())+
				" Sample(s):"+label_samples+" "+" Gene(s):"+
				(all_genes.size()>20?"N="+StringUtils.niceInt(all_genes.size()):String.join(";", all_genes))
				);
		w.writeEndElement();//end
		
		
		w.writeEndElement();//g
		w.writeEndElement();//svg
		if(output_html) {
			w.writeEndElement();//div
			
			w.writeStartElement("div");
			w.writeAttribute("style","display:none;");
			w.writeAttribute("id", "__PLACEHOLDER__");
			w.writeComment("Use this div to insert things later. ");
			w.writeCharacters("__PLACEHOLDER__");
			w.writeEndElement();
			
			
			w.writeStartElement("div");
			w.writeAttribute("id", "samples");

			
			w.writeStartElement("table");
			
			w.writeStartElement("caption");
			w.writeCharacters(extendedRegion.toNiceString());
			w.writeEndElement();
			
			w.writeStartElement("thead");
			
			w.writeStartElement("tr");
			
			w.writeStartElement("th");
			w.writeCharacters("Sample");
			w.writeEndElement();//th
			
			w.writeStartElement("th");
			w.writeCharacters("Mean Cov");
			w.writeEndElement();//th

			w.writeStartElement("th");
			w.writeCharacters("Median Cov in User Region");
			w.writeEndElement();//th
			
			w.writeStartElement("th");
			w.writeCharacters("Median Cov in Outer Region");
			w.writeEndElement();//th

			w.writeStartElement("th");
			w.writeCharacters("Median Cov in whole Region");
			w.writeEndElement();//th

			w.writeStartElement("th");
			w.writeCharacters("Ratio Inner/outer");
			w.writeEndElement();//th

			
			w.writeStartElement("th");
			w.writeCharacters("N. reads");
			w.writeEndElement();//th
			
			w.writeStartElement("th");
			w.writeCharacters("Path");
			w.writeEndElement();//th
			
			w.writeEndElement();//tr
			w.writeEndElement();//thead
			
			w.writeStartElement("tbody");
			for(SampleInfo si: samples) {
				w.writeStartElement("tr");

				w.writeStartElement("td");
				w.writeAttribute("style", "background-color:"+si.rgb);
				w.writeStartElement("a");
				w.writeAttribute("title",si.sample);
				w.writeAttribute("onclick","toggleSample(\""+si.sample+"\");");
				w.writeAttribute("href","#");
				w.writeCharacters(si.sample);
				w.writeEndElement();//a
				w.writeEndElement();//td
				
				w.writeStartElement("td");
				w.writeCharacters(format(si.meanDP));
				w.writeEndElement();//th

				w.writeStartElement("td");
				w.writeCharacters(format(si.medianInner));
				w.writeEndElement();//th
				
				w.writeStartElement("td");
				w.writeCharacters(format(si.medianOuter));
				w.writeEndElement();//th

				w.writeStartElement("td");
				w.writeCharacters(format(si.medianAll));
				w.writeEndElement();//th

				
				w.writeStartElement("td");
				if(si.medianOuter>0) {
					w.writeCharacters(String.valueOf((int)(100.0*(si.medianInner/si.medianOuter)))+"%");
					}
				
				w.writeEndElement();//th

				
				w.writeStartElement("td");
				w.writeCharacters(format(si.nReads));
				w.writeEndElement();//th

				w.writeStartElement("td");
				w.writeStartElement("span");
				w.writeAttribute("class", "path");
				w.writeAttribute("title",si.path.toString());
				w.writeCharacters(si.path.getFileName().toString());
				w.writeEndElement();//span
				w.writeEndElement();//th
				
				w.writeEndElement();//tr
				}
			w.writeEndElement();//tbody
			w.writeEndElement();//table
			w.writeEndElement();//div

			/** URLs **************/
			final Set<UrlSupplier.LabelledUrl> urls = new TreeSet<>();
			final UrlSupplier urlSupplier = new UrlSupplier(dict);
			urls.addAll(urlSupplier.of(rawRegion));
			getGenes(rawRegion,G->G.getType().equals("gene")).forEach(GFF->urls.addAll(urlSupplier.of(GFF)));
				if(!urls.isEmpty()) {
				w.writeStartElement("div");
				w.writeAttribute("id", "hyperlinks");
				
				w.writeStartElement("table");
				
				w.writeStartElement("caption");
				w.writeCharacters("Hyperlinks");
				w.writeEndElement();
				
				w.writeStartElement("thead");
				
				w.writeStartElement("tr");
				
				w.writeStartElement("th");
				w.writeCharacters("Source");
				w.writeEndElement();//th
				
				w.writeStartElement("th");
				w.writeCharacters("Target");
				w.writeEndElement();//th

				
				w.writeStartElement("th");
				w.writeCharacters("URL");
				w.writeEndElement();//th
	
				
				w.writeEndElement();//tr
				w.writeEndElement();//thead
				
				w.writeStartElement("tbody");
				for(UrlSupplier.LabelledUrl url: urls) {
					w.writeStartElement("tr");
	
					w.writeStartElement("td");
					w.writeCharacters(url.getLabel());
					w.writeEndElement();//th
	
					w.writeStartElement("td");
					w.writeCharacters(url.getTarget());
					w.writeEndElement();//th

					
					w.writeStartElement("td");
					w.writeStartElement("a");
					w.writeAttribute("title",url.getUrl());
					w.writeAttribute("target","_blank");
					w.writeAttribute("href",url.getUrl());
					w.writeCharacters(url.getDomain());
					w.writeEndElement();//a
					w.writeEndElement();//td
					
	
					
					w.writeEndElement();//tr
					}
				w.writeEndElement();//tbody
				w.writeEndElement();//table
				w.writeEndElement();//div
				} // end URLs
				
			/** KNOWN CNV */
			if(knownCnvFile!=null) {
				final List<Interval> knows = getKnownCNVs(extendedRegion).stream().
						filter(R->R.overlaps(rawRegion)).
						collect(Collectors.toList());
				if(!knows.isEmpty()) {

				w.writeStartElement("div");
				w.writeAttribute("id", "knowncnvs");

				w.writeStartElement("table");
				
				w.writeStartElement("caption");
				w.writeCharacters("Known CNV");
				w.writeEndElement();
				
				w.writeStartElement("thead");
				
				w.writeStartElement("tr");
				
				w.writeStartElement("th");
				w.writeCharacters("Position");
				w.writeEndElement();//th
				
				w.writeStartElement("th");
				w.writeCharacters("Length");
				w.writeEndElement();//th

				
				w.writeStartElement("th");
				w.writeCharacters("Name");
				w.writeEndElement();//th
	
				w.writeStartElement("th");
				w.writeCharacters("Fraction Region");
				w.writeEndElement();//th

				w.writeStartElement("th");
				w.writeCharacters("Fraction Known");
				w.writeEndElement();//th

				
				w.writeEndElement();//tr
				w.writeEndElement();//thead
				
				w.writeStartElement("tbody");
				for(Interval k:knows) {
					w.writeStartElement("tr");
	
					w.writeStartElement("td");
					
					w.writeStartElement("a");
					w.writeAttribute("title",new SimpleInterval(k).toNiceString());
					w.writeAttribute("href",hyperlink.apply(k).orElse(""));
					w.writeAttribute("target","_blank");
					w.writeCharacters(new SimpleInterval(k).toNiceString());
					w.writeEndElement();//a
					
					w.writeEndElement();//td
					
					int len = k.getLengthOnReference();
					w.writeStartElement("td");
					w.writeCharacters(StringUtils.niceInt(len));
					w.writeEndElement();//td
					
					w.writeStartElement("td");
					w.writeCharacters(k.getName());
					w.writeEndElement();//td
					
					
					if(k.overlaps(rawRegion)) {
						final double n2= rawRegion.getIntersectionLength(k);
						
						w.writeStartElement("td");
						w.writeCharacters((int)(100.0*(n2/rawRegion.getLengthOnReference()))+"%");
						w.writeEndElement();//td
						
						w.writeStartElement("td");
						w.writeCharacters((int)(100.0*(n2/len))+"%");
						w.writeEndElement();//td
						}
					else
						{
						w.writeEmptyElement("td");
						w.writeEmptyElement("td");
						}

					w.writeEndElement();//tr
					}
				w.writeEndElement();//tbody
				w.writeEndElement();//table
				w.writeEndElement();//div
				}
				
			}
				

			w.writeStartElement("div");
			w.writeAttribute("id", "about");

			w.writeCharacters("User interval:" + rawRegion.toNiceString()+" (length:"+StringUtils.niceInt(rawRegion.getLengthOnReference())+"). ");
			w.writeCharacters("Extended interval: " + extendedRegion.toNiceString()+" (length:"+StringUtils.niceInt(extendedRegion.getLengthOnReference())+"). ");
			w.writeCharacters("Mapping quality: " + this.min_mapq +". ");
			w.writeCharacters("Loess: " + this.run_loess +". ");
			w.writeCharacters("Percentil: " + (depth.length< this.use_average_instead_of_median_length_treshold?"median":"average") +". ");
			w.writeCharacters("Include original interval for mean/median depth calculation: " + this.include_original_interval_for_median +". ");
			w.writeCharacters("Date: " + StringUtils.now() +". ");
			w.writeEndElement();//div

			
			w.writeEmptyElement("hr");

			w.writeStartElement("div");
			w.writeCharacters("Made with "+getClass().getSimpleName()+" version:"+getVersion()+". Pierre Lindenbaum PhD. 2022.");
			w.writeEndElement();//div
			
			
			w.writeEndElement();//body
			w.writeEndElement();//html
			}
		else {
			w.writeEndDocument();
			}
		w.flush();
		w.close();
		return 0;
		}
	catch(final Throwable err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		}
	}

public static void main(final String[] args) {
	new CoveragePlotter().instanceMainWithExit(args);
	}

}
