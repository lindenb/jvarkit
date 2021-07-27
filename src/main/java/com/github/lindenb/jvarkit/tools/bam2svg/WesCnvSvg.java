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
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.CoverageFactory;
import com.github.lindenb.jvarkit.samtools.util.IntervalExtender;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.stream.HtsCollectors;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.XLINK;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;


/**
BEGIN_DOC

## Input

input is a set of bam file or a file with the '*.list' suffix containing the path to the bam files.

## Example

```
$ find dir -name "*.bam"  > bam.list
$ java -jar dist/wescnvsvg.jar -R ref.fasta -B cnv.bed bam.list > cnv.svg 
```

## Screenshots

https://twitter.com/yokofakun/status/1022503372669300738 : 

![ScreenShot](https://pbs.twimg.com/media/DjCpKcYXgAAq4fw.jpg:large)

https://twitter.com/yokofakun/status/1022805656905150464

![ScreenShot](https://pbs.twimg.com/media/DjG8Do0XsAA4U46.jpg:large)

https://twitter.com/yokofakun/status/1023953180927963137

![ScreenShot](https://pbs.twimg.com/media/DjXP8_4X0AAtSQZ.jpg)

https://twitter.com/yokofakun/status/1024315948847849472

![ScreenShot](https://pbs.twimg.com/media/DjcZ6vNXcAEYjEt.jpg)

https://twitter.com/yokofakun/status/1024315948847849472

![ScreenShot](https://pbs.twimg.com/media/DjcZ6vNXcAEYjEt.jpg)

https://twitter.com/yokofakun/status/1025330308193779712

![ScreenShot](https://pbs.twimg.com/media/Djq0Se-W0AEAbyR.jpg)

https://twitter.com/yokofakun/status/1040592885786263554

![ScreenShot](https://pbs.twimg.com/media/DnDttNgX4AAtxax.jpg)

https://twitter.com/yokofakun/status/1040577235856580608

![ScreenShot](https://pbs.twimg.com/media/DnDfaGLXcAArg0P.jpg)

https://twitter.com/yokofakun/status/1057625407913111557

![ScreenShot](https://pbs.twimg.com/media/Dq1whOTX0AAzkZc.jpg)

https://twitter.com/yokofakun/status/1180046139502059521

![ScreenShot](https://pbs.twimg.com/media/EGBdtO7WoAEHjEE?format=jpg&name=small)


END_DOC
 */
@Program(name="wescnvsvg",
description="SVG visualization of bam DEPTH for multiple regions",
keywords={"bam","alignment","graphics","visualization","svg","wes","bed","capture","exome"},
modificationDate="20210726",
creationDate="20180726"
)
public class WesCnvSvg  extends Launcher {
	private static final Logger LOG = Logger.build(WesCnvSvg.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-B","--bed","-b","--capture","-r","-rgn","--region","--interval"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class)
	private IntervalListProvider intervalListProvider = IntervalListProvider.unspecified();
	@Parameter(names={"-R","--ref","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"-w","--width"},description="Page width")
	private int drawinAreaWidth = 1000 ;
	@Parameter(names={"-height","--height"},description="Sample Track height")
	private int sampleTrackHeight = 100 ;
	@Parameter(names={"-smooth","--smooth"},description= "how to smooth data")
	private CoverageFactory.ScaleType scaleType = CoverageFactory.ScaleType.AVERAGE;
	@Parameter(names={"--title"},description="document title")
	private String domSvgTitle=WesCnvSvg.class.getSimpleName();
	@Parameter(names={"-u","--url","--hyperlink"},description= "creates a hyperlink an area is 'clicked'. " + Hyperlink.OPT_DESC,converter=Hyperlink.StringConverter.class,splitter=NoSplitter.class)
	private Hyperlink hyperlinkType = Hyperlink.empty();
	@Parameter(names={"-css","--css"},description="custom svg css stylesheet")
	private File cssFile = null;
	@Parameter(names={"-x","--extend"},description= IntervalExtender.OPT_DESC)
	private String extendWhat= "0";
	@Parameter(names={"-Q","--mapq"},description="Min mapping quality")
	private int minMappingQuality = 1;
	@Parameter(names={"--normalize"},description="normalize on median",hidden=true)
	private boolean normalize_on_median_flag =false;
	@Parameter(names={"--vcf"},description="plot VCF data")
	private Path vcfFile = null;
	@DynamicParameter(names = "-D", description = "other parameters. '-Dkey=value'. Undocumented.")
	private Map<String, String> dynaParams = new HashMap<>();
	
	
	private class BamInput
		{
		Path bamPath;
		String sample;
		final List<double[]> coverages = new ArrayList<>();
		double getPixelHeight() {
			return WesCnvSvg.this.sampleTrackHeight;
			}
		}
	
	
	private class CaptureInterval extends SimpleInterval
		{
		double pixelx=0.0;
		CaptureInterval(final Locatable locatable) {
			super(locatable);
			}
		
		public double getPixelX1() {
			return pixelx;
			}
		
		public double getPixelWidth(){
			return (this.getLengthOnReference()/(double)WesCnvSvg.this.countBasesToBeDisplayed)*WesCnvSvg.this.drawinAreaWidth;
		}
		
		String getName() {
			return this.getContig()+":"+StringUtils.niceInt(this.getStart())+"-"+StringUtils.niceInt(this.getEnd());
			}
		String getId() {
			return getContig()+"_"+this.getStart()+"_"+this.getEnd();
			}
		}
	
	//private final List<CaptureInterval> intervals = new ArrayList<>();
	private final List<BamInput> bamInputs = new ArrayList<>();
	private ReferenceSequenceFile indexedFastaSequenceFile;
	private DecimalFormat decimalFormater = new DecimalFormat("##.##");
	//private double globalMaxDepth = 0.0;
	private int countBasesToBeDisplayed = 0;
	private final int gc_win=100;

		 
	
	/** convert double to string */
	private String format(double v)
		{
		return this.decimalFormater.format(v);
		}
	
	private double getGcPercent(final GenomicSequence seq,int chromStart,int chromEnd)
		{
		final SAMSequenceRecord ssr = seq.getSAMSequenceRecord();
		while(chromEnd-chromStart+1< this.gc_win)
			{
			chromEnd++;
			chromStart--;
			}
		chromStart = Math.max(chromStart,1);
		chromEnd = Math.min(ssr.getSequenceLength(),chromEnd);
		if(chromStart>chromEnd) return 0.0;
		
		final GenomicSequence.GCPercent gcpercent=seq.getGCPercent(chromStart-1, chromEnd);
		if(gcpercent.isEmpty()) return 0.0;
		return gcpercent.getGCPercent();
		}	
	
	@Override
	public int doWork(final List<String> args) {
		XMLStreamWriter w = null;
		BufferedReader r = null;
		OutputStream fout=null;
		VCFReader vcfReader = null;
		try
			{
			
			this.indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);			
			final IntervalExtender extender = IntervalExtender.of(refDict,this.extendWhat);
			if(extender.isShriking()) {
				LOG.error("extend factor <1.0");
				return -1;
				}
			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(refDict);
			final ContigDictComparator contigDictCompare = new ContigDictComparator(refDict);
			
			final List<CaptureInterval> userIntervals = this.intervalListProvider.
					stream().
					/* https://stackoverflow.com/questions/25523375 */
					map(loc->contigNameConverter.convertToSimpleInterval(loc).<RuntimeException>orElseThrow(()->new RuntimeException(new JvarkitException.ContigNotFoundInDictionary(loc.getContig(),refDict)))).
					map(extender).
					collect(HtsCollectors.mergeIntervals()).
					map(L->new CaptureInterval(L)).
					sorted(contigDictCompare.createLocatableComparator()).
					collect(Collectors.toCollection(ArrayList::new))
					;
		
			if(userIntervals.isEmpty())
				{
				LOG.error("no interval or bed defined");
				return -1;
				}
			
			this.countBasesToBeDisplayed = userIntervals.stream().
					mapToInt(R->R.getLengthOnReference()).
					sum();
			if(this.countBasesToBeDisplayed<1) {
				LOG.error("Nothing to display. BED count==0");
				return -1;
				}
			else
				{
				double x1=0;
				for(int i=0;i< userIntervals.size();++i) {
					final CaptureInterval ci = userIntervals.get(i);
					ci.pixelx = x1;
					x1+=ci.getPixelWidth();
					}
				}
			
			/* distinct ordered contigs */
			final List<String> distinctContigs = userIntervals.stream().
					map(R->R.getContig()).
					collect(Collectors.toSet()).
					stream().
					sorted(contigDictCompare).
					collect(Collectors.toList());
			
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT).
					referenceSequence(this.faidx);
			for(final Path bamFile:IOUtils.unrollPaths(args)) {
				final BamInput bi = new BamInput();
				bi.bamPath = bamFile;
				try(SamReader samReader = srf.open(bamFile)) {
					final SAMSequenceDictionary samDict = SequenceDictionaryUtils.extractRequired(samReader.getFileHeader());
					if(!SequenceUtil.areSequenceDictionariesEqual(refDict, samDict)) {
						throw new JvarkitException.DictionariesAreNotTheSame(refDict, samDict);
						}
					bi.sample = samReader.getFileHeader().getReadGroups().stream().
							map(V->V.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamFile));
					final CoverageFactory covFactory = new CoverageFactory().setMappingQuality(this.minMappingQuality);
					for(final String contig: distinctContigs) {
						final List<CaptureInterval> contig_intervals = userIntervals.stream().
								filter(R->R.getContig().equals(contig)).
								collect(Collectors.toList());
						
						final CoverageFactory.SimpleCoverage coverage = covFactory.getSimpleCoverage(samReader,
							contig_intervals,
							null);
						for(CaptureInterval rgn:contig_intervals)  {
							final double[] array = coverage.getSubCoverage(rgn).scale(this.scaleType, (int)rgn.getPixelWidth());
							bi.coverages.add(array);
							}
						}
					}
				this.bamInputs.add(bi);
				}
			if(this.bamInputs.isEmpty()) {
				LOG.error("no bam input");
				return -1;
				}
			
			if(this.vcfFile!=null) {
				vcfReader = VCFReaderFactory.makeDefault().open(this.vcfFile,true);
				final SAMSequenceDictionary vcfDict = vcfReader.getHeader().getSequenceDictionary();
				if(vcfDict!=null) SequenceUtil.assertSequenceDictionariesEqual(refDict, vcfDict);
				}
			
			
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			if(this.outputFile==null)
				{
				w=xof.createXMLStreamWriter(stdout(), "UTF-8");
				}
			else
				{
				fout =Files.newOutputStream(this.outputFile);
				w=xof.createXMLStreamWriter(fout, "UTF-8");
				}
			
			final Function<List<Point2D.Double>,String> points2str = (L)->
				L.stream().map(S->format(S.getX())+","+format(S.getY())).
				collect(Collectors.joining(" "));
				
			final Consumer<List<Point2D.Double>> simplifyPoints = (L)->{
				for(int z=0;z+1< L.size();++z) {
					if(L.get(z).getY()==L.get(z+1).getY())
						{
						L.get(z).x = L.get(z+1).x;
						L.remove(z+1);
						}
					}
			};
			
			w.writeStartDocument("UTF-8", "1.0");
			
			
			final Dimension dim=new Dimension(this.drawinAreaWidth,0);
			final int bed_header_height=20;
			dim.height += bed_header_height;
			dim.height+=(int)this.bamInputs.stream().
					mapToDouble(B->B.getPixelHeight()).
					sum();
			
			LOG.debug("drawing area: "+dim);
			
			w.writeStartElement("svg");
			w.writeAttribute("width", String.valueOf(dim.width));
			w.writeAttribute("height", String.valueOf(dim.height));
			w.writeDefaultNamespace(SVG.NS);
			w.writeNamespace("xlink", XLINK.NS);
			
			// https://stackoverflow.com/questions/15717970
			w.writeStartElement("style");
			if(this.cssFile!=null) {
				w.writeCharacters(IOUtil.slurp(this.cssFile));
				}
			else {
				w.writeCharacters(
						"g.maing {stroke:black;stroke-width:0.5px;fill:whitesmoke;font-size:10pt;}\n"+
						"text.sampleLabel {stroke:none;stroke-width:0.5px;fill:blue;}" +
						"text.captureLabel {stroke:none;stroke-width:0.5px;fill:slategrey;text-anchor:middle;}" +
						"polygon.area {stroke:darkgray;stroke-width:0.5px;fill:url('#grad01');}" +
						"line.linedp {stroke:darkcyan;stroke-width:0.3px;opacity:0.4;}" +
						"text.linedp {fill-opacity:0.6;font-size:7px;stroke:none;stroke-width:0.5px;fill:darkcyan;}" +
						"rect.sampleFrame { fill:none;stroke:slategray;stroke-width:0.3px;}" +
						"rect.clickRgn {fill:none;stroke:none;pointer-events:all;}" +
						"polyline.gc {stroke:lightcoral;stroke-width:0.3px;fill:none;}"+
						"polyline.clipping {stroke:orange;stroke-width:0.8px;fill:none;}"+
						"circle.ar {fill:orange;stroke:none;}"+
						"circle.aa {fill:red;stroke:none;}"+
						"circle.rr {fill:green;stroke:none;}"
						);
				}
			w.writeEndElement();//style
			
			w.writeStartElement("title");
			w.writeCharacters(this.domSvgTitle);
			w.writeEndElement();

			w.writeStartElement("defs");
			// alleles
			final double genotype_radius = Double.parseDouble(this.dynaParams.getOrDefault("gt.radius", "1.5"));
			w.writeEmptyElement("circle");
			w.writeAttribute("r", format(genotype_radius));
			w.writeAttribute("id","rr");
			w.writeAttribute("class","rr");
			
			w.writeEmptyElement("circle");
			w.writeAttribute("r", format(genotype_radius));
			w.writeAttribute("id","ar");
			w.writeAttribute("class","ar");
			
			w.writeEmptyElement("circle");
			w.writeAttribute("r", format(genotype_radius));
			w.writeAttribute("id","aa");
			w.writeAttribute("class","aa");
			
			//gradient
			w.writeStartElement("linearGradient");
			w.writeAttribute("id","grad01");
			w.writeAttribute("x1","50%");
			w.writeAttribute("x2","50%");
			w.writeAttribute("y1","0%");
			w.writeAttribute("y2","100%");
			w.writeEmptyElement("stop");
				w.writeAttribute("offset","0%");
				w.writeAttribute("style","stop-color:lightgray;stop-opacity:1;");
			w.writeEmptyElement("stop");
				w.writeAttribute("offset","100%");
				w.writeAttribute("style","stop-color:gray;stop-opacity:1;");
			w.writeEndElement();

			
			// gc percent
			for(final CaptureInterval ci: userIntervals)
				{
				final GenomicSequence genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile,ci.getContig());
				final int gc_percent_width= (int)ci.getPixelWidth();
				final List<Point2D.Double> points= new ArrayList<>(gc_percent_width);
				for(int x=0;x< gc_percent_width;++x)
					{
					int pos1= ci.getStart()+(int)(((x+0)/ci.getPixelWidth())*ci.getLengthOnReference());
					int pos2= ci.getStart()+(int)(((x+1)/ci.getPixelWidth())*ci.getLengthOnReference());
					double gc_percent = getGcPercent(genomicSequence,pos1,pos2);
					double y = this.sampleTrackHeight - this.sampleTrackHeight*gc_percent;
					
					points.add(new Point2D.Double(x, y));
					}
				simplifyPoints.accept(points);
				
				w.writeStartElement("polyline");
				w.writeAttribute("class","gc");
				w.writeAttribute("id","z"+ci.getId());
				w.writeAttribute("points",points2str.apply(points));
				w.writeStartElement("title");
				w.writeCharacters("GC %");
				w.writeEndElement();
				w.writeEndElement();
				}
			
			w.writeEndElement();//defs
			
			
			w.writeStartElement("script");
			
			final StringBuilder openBrowserFunction = new StringBuilder(
					"function openGenomeBrowser(contig,chromStart,chromEnd) {\n"
					);
			if(!this.hyperlinkType.isEmpty())
				{
				openBrowserFunction.append("var url=\""+this.hyperlinkType.getPattern()+"\".replace(/__CHROM__/g,contig).replace(/__START__/g,chromStart).replace(/__END__/g,chromEnd);\n");
				openBrowserFunction.append("window.open(url,'_blank');\n");

				}
			else
				{
				//nothing
				}
			openBrowserFunction.append("}\n");
			
			w.writeCData(
				openBrowserFunction.toString() +
				"function clicked(evt,contig,chromStart,chromEnd){\n" +
			    "    var e = evt.target;\n" +
			    "    var dim = e.getBoundingClientRect();\n" +
			    "    var x = 1.0 * evt.clientX - dim.left;\n" + 
			    "    var cLen = 1.0* (chromEnd - chromStart); if(cLen<1) cLen=1.0;\n" + 
			    "    var pos1 = chromStart + parseInt(((x+0)/dim.width)*cLen);\n" +
			    "    var pos2 = chromStart + parseInt(((x+1)/dim.width)*cLen);\n" +
			    "   openGenomeBrowser(contig,pos1,pos2);\n" +
			    "}\n");                
			w.writeEndElement();//script
			
			
			w.writeStartElement("g");
			w.writeAttribute("class", "maing");
			
			int y=0;
			w.writeStartElement("g");
			w.writeComment("interval background");
			for(final CaptureInterval ci:userIntervals)
				{
				w.writeStartElement("text");
					w.writeAttribute("class", "captureLabel");
					w.writeAttribute("textLength",String.valueOf(ci.getPixelWidth()*0.8));
					w.writeAttribute("lengthAdjust", "spacing");
					w.writeAttribute("x",String.valueOf(ci.getPixelX1()+ci.getPixelWidth()/2.0));
					w.writeAttribute("y",String.valueOf(bed_header_height-2));
					w.writeCharacters(ci.getName());
					w.writeStartElement("title");
					w.writeCharacters(ci.toNiceString());
					w.writeEndElement();//title
				w.writeEndElement();//text

				
				w.writeStartElement("rect");
					w.writeAttribute("style","stroke:black;fill:none;");
					w.writeAttribute("x",String.valueOf(ci.getPixelX1()));
					w.writeAttribute("y", "0");
					w.writeAttribute("width",String.valueOf(ci.getPixelWidth()));
					w.writeAttribute("height",String.valueOf(dim.height));
					w.writeStartElement("title");
						w.writeCharacters(ci.getName());
					w.writeEndElement();
				w.writeEndElement();
				}
			w.writeEndElement();//interval background
			
			y+= bed_header_height;
			
			for(final BamInput bi:this.bamInputs)
				{	
				w.writeComment(bi.bamPath.toString());
				w.writeStartElement("g");
				w.writeAttribute("transform","translate(0,"+y+")");
				
				if(normalize_on_median_flag) {
					final double medianDepth = Math.max(1.0,Percentile.median().evaluate(bi.coverages.stream().flatMapToDouble(B->Arrays.stream(B)).toArray()).orElse(1.0));
					LOG.info("median"+medianDepth);
					for(final double[] coverage_array: bi.coverages) {
						for(int px=0;px< coverage_array.length;px++) {
							coverage_array[px]/=medianDepth;
							}
						}
					}
				
				final double maxDepth = bi.coverages.stream().flatMapToDouble(B->Arrays.stream(B)).max().orElse(1.0);
				
				for(int ridx=0;ridx<userIntervals.size();ridx++) {
					final CaptureInterval ci = userIntervals.get(ridx);
					final String clickedAttribute = "clicked(evt,\""+ci.getContig()+"\","+ci.getStart()+","+ci.getEnd()+")";
					final double[] coverage_array = bi.coverages.get(ridx);
					final double leftx =ci.getPixelX1();
					w.writeStartElement("g");
					w.writeAttribute("transform","translate("+leftx+",0)");
					
					final int segment_width = (int)ci.getPixelWidth();

					//coverage
					{
					final List<Point2D.Double> points = new ArrayList<>(segment_width);
					points.add(new Point2D.Double(0,bi.getPixelHeight()));

					for(int px=0;px< coverage_array.length;px++) {
						final double y_avg_cov=  coverage_array[px];
						final double new_y = bi.getPixelHeight()-(y_avg_cov/maxDepth)*bi.getPixelHeight();
						points.add(new Point2D.Double(px,new_y));
						}
					points.add(new Point2D.Double(ci.getPixelWidth(),bi.getPixelHeight()));
					simplifyPoints.accept(points);				
					points.add(new Point2D.Double(leftx,bi.getPixelHeight()));//close
					w.writeEmptyElement("polygon");
					w.writeAttribute("class","area");
					w.writeAttribute("title",ci.toNiceString());
					//w.writeAttribute("onclick", clickedAttribute);
					w.writeAttribute("points",points2str.apply(points));
					}
					//w.writeEndElement();//g
					
					
					
					
					int depthshift=1;
					for(;;) {
						final int numdiv =(int) (maxDepth/depthshift);
						if(numdiv<=10) break;
						depthshift*=10;
						}
					
					
					int depth=depthshift;
					while(depth<  maxDepth)
						{
						double new_y = bi.getPixelHeight()-(depth/maxDepth)*bi.getPixelHeight();
						
						w.writeStartElement("text");
							w.writeAttribute("class", "linedp");
							w.writeAttribute("x","1");
							w.writeAttribute("y",String.valueOf(new_y+1));
							w.writeCharacters(String.valueOf(depth));
						w.writeEndElement();//text
						
						w.writeStartElement("line");
							w.writeAttribute("class","linedp");
							w.writeAttribute("stroke-dasharray","4");
							w.writeAttribute("x1", "0");
							w.writeAttribute("x2", String.valueOf(ci.getPixelWidth()));
							w.writeAttribute("y1", String.valueOf(new_y));
							w.writeAttribute("y2", String.valueOf(new_y));
							w.writeStartElement("title");
								w.writeCharacters(String.valueOf(depth));
							w.writeEndElement();
						w.writeEndElement();//line
						depth+=depthshift;
						
						

						}
					// polyline
					w.writeEmptyElement("use");
					w.writeAttribute("xlink",XLINK.NS,"href","#z"+ci.getId());
					w.writeAttribute("x","0");
					w.writeAttribute("y","0");

					
					//click
					w.writeStartElement("rect");
					w.writeAttribute("class","clickRgn");
					w.writeAttribute("onclick", clickedAttribute);
					w.writeAttribute("x", "0");
					w.writeAttribute("y", "0");
					w.writeAttribute("width", String.valueOf(ci.getPixelWidth()));
					w.writeAttribute("height", String.valueOf(bi.getPixelHeight()));
					w.writeEndElement();
					
					
					//genotype
					if(vcfReader!=null) {
						try(CloseableIterator<VariantContext> iter = vcfReader.query(ci)) {
							while(iter.hasNext()) {
								final VariantContext ctx = iter.next();
								final Genotype gt = ctx.getGenotype(bi.sample);
								if(gt==null) break;
								String allele_id = null;
								switch(gt.getType()) {
									case HET: allele_id = "ar"; break;
									case HOM_REF: allele_id = "rr"; break;
									case HOM_VAR: allele_id = "aa"; break;
									default: allele_id=null;break;
									}
								if(allele_id!=null) {
									w.writeEmptyElement("use");
									w.writeAttribute("xlink",XLINK.NS,"href","#"+allele_id);
									w.writeAttribute("x",format(((ctx.getStart()-ci.getStart())/(double)ci.getLengthOnReference())*ci.getPixelWidth()));
									w.writeAttribute("y",format(bi.getPixelHeight()-2*genotype_radius));
									}
								}
							}
						}
					
					
					w.writeEndElement();//g
					}
				
				//frame for this sample
				w.writeStartElement("rect");
					w.writeAttribute("class","sampleFrame");
					w.writeAttribute("x","0");
					w.writeAttribute("y", "0");
					w.writeAttribute("width",String.valueOf(dim.width));
					w.writeAttribute("height",String.valueOf(bi.getPixelHeight()));
				w.writeEndElement();//rect

				
				w.writeStartElement("text");
					w.writeAttribute("class", "sampleLabel");
					w.writeAttribute("x","5");
					w.writeAttribute("y","12");
					w.writeStartElement("title");
					w.writeCharacters(bi.bamPath.toString());
					w.writeEndElement();
					w.writeCharacters(bi.sample);
				w.writeEndElement();//text

				
				w.writeEndElement();//g
				y+=bi.getPixelHeight();
				}
			
			w.writeStartElement("g");
			w.writeComment("interval lines");
			for(int n=0;n<= userIntervals.size();n++)
				{
				w.writeEmptyElement("line");
				String x1=
						n< userIntervals.size()?
						String.valueOf(userIntervals.get(n).getPixelX1()):
						String.valueOf(dim.width)
						;
				w.writeAttribute("x1", x1);
				w.writeAttribute("y1", "0");
				w.writeAttribute("x2", x1);
				w.writeAttribute("y2", String.valueOf(dim.height));
				}
			w.writeEndElement();//interval lines

			
			w.writeEndElement();//g
			w.writeEndElement();//svg
			w.writeEndDocument();
			w.flush();
			w.close();
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfReader);
			CloserUtil.close(w);
			CloserUtil.close(fout);
			CloserUtil.close(r);
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(this.bamInputs);
			}
		}

	
public static void main(final String[] args)
	{
	new WesCnvSvg().instanceMainWithExit(args);
	}
}
