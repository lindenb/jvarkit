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
package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
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

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.ArrayResizer;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalExtender;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
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
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC
input a set of bam/cram files or one file with the suffix '.list' containing the path to the bams

output is a HTML+SVG file

```
java -jar dist/coverageplotter.jar -R src/test/resources/rotavirus_rf.fa --region "RF01:100-200" src/test/resources/*.bam 
```


END_DOC 
 */
@Program(
	name="coverageplotter",
	description="Display an image of depth to display any anomaly an intervals+bams",
	keywords={"cnv","bam","depth","coverage","svg"},
	creationDate="20200605",
	modificationDate="20220808"
	)
public class CoveragePlotter extends Launcher {
	private static final Logger LOG = Logger.build( CoveragePlotter.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refPath = null;
	@Parameter(names={"--region","--interval"},description = "Interval region",required=true)
	private String intervalStr=null;
	@Parameter(names={"--mapq"},description = "min mapping quality")
	private int min_mapq=1;
	@Parameter(names={"--max-depth"},description = "ignore position if depth > 'x'")
	private int max_depth=500;
	@Parameter(names={"--dimension","--dim"},description = "Image Dimension. " + DimensionConverter.OPT_DESC, converter=DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,300);
	@Parameter(names={"--extend","-x"},description = IntervalExtender.OPT_DESC)
	private String extendStr="3.0";
	@Parameter(names= {"--gtf"},description="Optional Tabix indexed GTF file. Will be used to retrieve an interval by gene name, or to display gene names in a region.")
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
	
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");

	private static class SampleInfo {
		String sample;
		Path path;
		long nReads = 0L;
		Point2D maxPosition=null;
		String rgb=null;
		double medianDP=0.0;
		double meanDP=0.0;
		}
	
	private void drawKnownCnv(final XMLStreamWriter w,final Rectangle rectangle,final Locatable region) {
		if(this.knownCnvFile==null) return;
		final String fname=this.knownCnvFile.getFileName().toString();

		final Pileup<Interval> pileup = new Pileup<>();
		final Predicate<Interval> rejectCnv = cnv->(this.ignore_cnv_overlapping && cnv.getStart() < region.getStart() && cnv.getEnd() > region.getEnd());

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
							if(rejectCnv.test(rgn)) return;
							pileup.add(rgn);
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
								pileup.add(rgn);
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
		if(!pileup.isEmpty() ) {
			final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*rectangle.getWidth();
			final double featureHeight = 4.0/pileup.getRowCount();
			for(int row=0;row< pileup.getRowCount();++row) {
				for(final Interval cnv:pileup.getRow(row)) {
					final double y= rectangle.getHeight()-8.0 + row*featureHeight;
					final double x1 = position2pixel.applyAsDouble(cnv.getStart());
					final double x2 = position2pixel.applyAsDouble(cnv.getEnd());
					try {
						w.writeEmptyElement("rect");
						w.writeAttribute("class", "cnv");
						w.writeAttribute("x", format(x1));
						w.writeAttribute("y", format(y-1));
						w.writeAttribute("width", format(Math.max(0.5, x2-x1)));
						w.writeAttribute("height", format(featureHeight*0.9));
						}
					catch(final Throwable err) {
						LOG.warn(err);
						}
					}
				}
			}
		}

	

	private Stream<GTFLine> getGenes(final Locatable region) {
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
			
			final GTFCodec codec = new GTFCodec();
			final TabixReader.Iterator iter=tbr.query(ctg,region.getStart(),region.getEnd());
			final TabixReader tbrfinal = tbr;
			final AbstractIterator<GTFLine> iter2= new AbstractIterator<GTFLine>() {
				@Override
				protected GTFLine advance() {
					try {
						for(;;) {
							final String line = iter.next();
							if(line==null) return null;
							if(StringUtils.isBlank(line) ||  line.startsWith("#")) continue;
							final String tokens[]= CharSplitter.TAB.split(line);
							if(tokens.length<9 ) continue;
							tokens[0]=region.getContig();
							final GTFLine gtfline = codec.decode(line);
							if(gtfline==null) continue;
							return gtfline;
							}
						} catch (final IOException e) {
						LOG.error(e);
						return null;
						}
					}
				};
			return StreamSupport.stream(new IterableAdapter<GTFLine>(iter2).spliterator(),false).
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
	
			getGenes(region).forEach(gtfline->{
					final double x1 = position2pixel.applyAsDouble(gtfline.getStart());
					final double x2 = position2pixel.applyAsDouble(gtfline.getEnd());
					try {
					if(gtfline.getType().equals("gene") ) {
						w.writeEmptyElement("text");
						w.writeAttribute("class","gene");
						w.writeAttribute("x", format(Math.max(x1,1)));
						w.writeAttribute("y", format(y-(geneSize+3)));
						w.writeCharacters(gtfline.getAttribute("gene_name"));
						}
					else if(gtfline.getType().equals("exon") ) {
						w.writeStartElement("rect");
						w.writeAttribute("class","gene");
						w.writeAttribute("x", format(x1));
						w.writeAttribute("y", format(y-1));
						w.writeAttribute("width", format(x2-x1));
						w.writeAttribute("height", format(3));
						title(w,"exon");
						w.writeEndElement();
						}
					else if(gtfline.getType().equals("transcript") ) {
						w.writeStartElement("line");
						w.writeAttribute("class","gene");
						w.writeAttribute("x1", format(x1));
						w.writeAttribute("y1", format(y));
						w.writeAttribute("x2", format(x2));
						w.writeAttribute("y2", format(y));
						title(w,"transcript");
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
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refPath);
		final SamReaderFactory samReaderFactory = SamReaderFactory.
					makeDefault().
					referenceSequence(CoveragePlotter.this.refPath).
					validationStringency(ValidationStringency.LENIENT)
					;
		
		final List<Path> inputBams =  IOUtils.unrollPaths(args);
		
		if(inputBams.isEmpty()) {
			LOG.error("input bam file missing.");
			return -1;
			}
		final List<SampleInfo> samples = new ArrayList<>(inputBams.size());
		final Locatable the_locatable = IntervalParserFactory.newInstance(dict).make().
			apply(this.intervalStr).orElseThrow(()->new IllegalArgumentException("Cannot parse interval:\""+this.intervalStr+"\""));
		
		
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
			
		final int sampleFontSize=7;
		final BufferedImage image = new BufferedImage(this.dimension.width,this.dimension.height,BufferedImage.TYPE_INT_ARGB);
		final double max_y = this.max_normalized_y;
		final ToDoubleFunction<Double> normToPixelY = NORM->  this.dimension.getHeight() * (1.0 -  (NORM/max_y));
		final double y_mid = normToPixelY.applyAsDouble(1.0);
		final DiscreteMedian<Integer> discreteMedian = new DiscreteMedian<>();
		
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
		
		if(this.force_svg_output) {
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
			w.writeCharacters(extendedRegion.toNiceString());
			w.writeEndElement();//h2
			w.writeStartElement("div");
			}
		else {
			w.writeStartDocument("UTF-8","1.0");
			}
		w.writeStartElement("svg");
		w.writeDefaultNamespace(SVG.NS);
		w.writeAttribute("width", format(this.dimension.width));
		w.writeAttribute("height", format(this.dimension.height));
		w.writeStartElement("style");
		w.writeCharacters(""
				+ ".geneh {fill:rgb(240,240,240);stroke:none;}\n"
				+ ".ymedian0 {fill:none;stroke:blue;}\n"
				+ ".ymedian1 {fill:none;stroke:cyan;}\n"
				+ ".frame {fill:none;stroke:darkgray;}\n"
				+ ".rgnbound {fill:none;stroke:green;}\n"
				+ ".hist {fill:none;stroke:darkgray;display:block;}\n"
				+ ".title1 {fill:darkgray;stroke:none;}\n"
				+ ".title2 {fill:darkgray;stroke:none;}\n"
				+ ".cnv {fill:blue;stroke:orange;opacity:0.5;}\n"
				+ ".gene {fill:green;stroke:orange;opacity:0.5;}\n"
				);
		w.writeEndElement();
		
		w.writeStartElement("g");
		//draw exons
		getGenes(extendedRegion).filter(G->G.getType().equals("exon")).forEach(EX->{
			try {
				final double x1 = pos2pixel.applyAsDouble(EX.getStart());
				final double x2 = pos2pixel.applyAsDouble(EX.getEnd());

				w.writeEmptyElement("rect");
				w.writeAttribute("class", "geneh");
				w.writeAttribute("x", format(x1));
				w.writeAttribute("y", format(0));
				w.writeAttribute("width", format(Math.max(0.5, x2-x1)));
				w.writeAttribute("height", format(this.dimension.height));
				}
			catch(Throwable err)  {
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

		//draw genes
		drawGenes(w, new Rectangle(0,0,image.getWidth(),image.getHeight()), extendedRegion);
		// draw knowns
		drawKnownCnv(w, new Rectangle(0,0,image.getWidth(),image.getHeight()), extendedRegion);
		
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
		
		final int depth[]= new int[extendedRegion.getLengthOnReference()];		
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

				
			discreteMedian.clear();
			double sumDepth=0;
			int countDepth=0;
			for(int i=0;i< depth.length;i++) {
				if(depth[i]>this.max_depth) continue;
				if(!include_original_interval_for_median) {
						final int pos = extendedRegion.getStart() + i;
						if(rawRegion.getStart()<=pos && pos<=rawRegion.getEnd()) continue;
						}
				discreteMedian.add(depth[i]);
				sumDepth+=depth[i];
				countDepth++;
				}
			if(countDepth<=0) {
				final String msg = "Skipping "+sampleInfo.sample +" "+extendedRegion+" because I cannot find any depth";
				LOG.warning(msg);
				continue;
				}
			sampleInfo.meanDP = sumDepth/countDepth;
			
			sampleInfo.medianDP = discreteMedian.getMedian().orElse(0);
			if(sampleInfo.medianDP <= 0.0) {
				final String msg = "Skipping "+sampleInfo.sample +" "+extendedRegion+" because median is 0";
				LOG.warning(msg);
				continue;
				}
			
			
			// normalize on median
			double[] normArray = new double[depth.length];
			for(int i=0;i< depth.length;i++) {
				normArray[i] = depth[i]/sampleInfo.medianDP;
				}
			normArray= new ArrayResizer().resizeToDouble(normArray, image.getWidth());

			
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

			final Color sampleRGB =  Color.getHSBColor((float) (samples.size()) / (float)inputBams.size(), 0.85f, 1.0f);
			sampleInfo.rgb="rgb("+sampleRGB.getRed()+","+sampleRGB.getGreen()+","+sampleRGB.getBlue()+")";
			
			w.writeStartElement("polyline");
			w.writeAttribute("id", sampleInfo.sample);
			w.writeAttribute("class", "hist");
			w.writeAttribute("style", "display:block;stroke:" + sampleInfo.rgb +";");
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
		final Set<String> all_genes = getGenes(extendedRegion).
				filter(G->G.getType().equals("gene")).
				map(G->G.getAttribute("gene_name")).
				filter(F->!StringUtils.isBlank(F)).
				collect(Collectors.toCollection(TreeSet::new))
				;

		w.writeStartElement("text");
		w.writeAttribute("class","title1");
		w.writeAttribute("x", format(10));
		w.writeAttribute("y", format(10));
		w.writeCharacters(
				extendedRegion.toNiceString()+" Length:"+StringUtils.niceInt(extendedRegion.getLengthOnReference())+
				" Sample(s):"+label_samples+" "+" Gene(s):"+
				(all_genes.size()>20?"N="+StringUtils.niceInt(all_genes.size()):String.join(";", all_genes))
				);
		w.writeEndElement();//end

		
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
		
		
		w.writeEndElement();//g
		w.writeEndElement();//svg
		if(output_html) {
			w.writeEndElement();//div
			w.writeStartElement("div");
			
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
			w.writeCharacters("Median Cov");
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
				w.writeCharacters(format(si.medianDP));
				w.writeEndElement();//th
				
				w.writeStartElement("td");
				w.writeCharacters(format(si.nReads));
				w.writeEndElement();//th

				w.writeStartElement("td");
				w.writeCharacters(si.path.toString());
				w.writeEndElement();//th
				
				w.writeEndElement();//tr
				}
			w.writeEndElement();//tbody
			w.writeEndElement();//table

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
