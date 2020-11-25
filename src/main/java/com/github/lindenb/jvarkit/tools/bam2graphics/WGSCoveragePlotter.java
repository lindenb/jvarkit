/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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

import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Result;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.SequenceUtil;
/**
BEGIN_DOC
input is an interval of a file source of interval (bed, vcf, gtf, interval_list , ,etc...)


```
java -jar dist/coverageplotter.jar -R src/test/resources/rotavirus_rf.fa -B src/test/resources/S1.bam -B src/test/resources/S2.bam "RF01:1-4000" -w 50 | less -r
```


END_DOC 
 */
@Program(
	name="wgscoverageplotter",
	description="Whole genome coverage plotter",
	keywords={"cnv","bam","depth","coverage"},
	creationDate="20201125",
	modificationDate="20201125"
	)
public class WGSCoveragePlotter extends Launcher {
	private static final Logger LOG = Logger.build( WGSCoveragePlotter.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refPath = null;
	@Parameter(names={"--min-contig-length"},description="Skip chromosome with length < 'x'. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_contig_length = 0;
	@Parameter(names={"--skip-contig-regex"},description="Skip chromosome matching this regular expression")
	private String skipContigExpr = "(NC_007605|hs37d5)";

	@Parameter(names={"--percentile"},description="How to we bin values")
	private DiscreteMedian.Tendency percentile = DiscreteMedian.Tendency.median;

	@Parameter(names={"--mapq"},description = "min mapping quality")
	private int min_mapq=1;
	@Parameter(names={"-C","--max-depth"},description = "Max depth to display")
	private int max_depth=100;
	@Parameter(names={"--clip","--cap"},description = "Don't allow coverage to be greater than 'max-depth' in the SVG file.")
	private boolean cap_depth= false;

	@Parameter(names={"--dimension"},description = "Image Dimension. " + DimensionConverter.OPT_DESC, converter=DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,500);
	@DynamicParameter(names = "-D", description = "style",hidden=true)
	private Map<String, String> dynaParams = new HashMap<>();
	@Parameter(names={"--disable-paired-overlap"},description="Disable: Count overlapping bases with mate for paired-end")
	private boolean disable_paired_overlap_flag=false;

	
	private Document document = null;
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");


	
	
private class ChromInfo {
	final SAMSequenceRecord ssr;
	double x = 0;
	double width=0;
	ChromInfo(final SAMSequenceRecord ssr) {
		this.ssr=ssr;
		}
	}
private String format(double v) {
return this.decimalFormater.format(v);
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
	try
		{
		final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		dbf.setNamespaceAware(true);
		final DocumentBuilder db = dbf.newDocumentBuilder();
		this.document=db.newDocument();
		
		
		final Element svgRoot = element("svg");
		this.document.appendChild(svgRoot);
		svgRoot.setAttribute("width",format(dimension.width+1));
		svgRoot.setAttribute("height",format(dimension.height+1));
			
		final Element maintitle = element("title");
		svgRoot.appendChild(maintitle);
		
		final Element descr = element("desc");
		svgRoot.appendChild(descr);
		descr.appendChild(text("Author: Pierre Lindenbaum"));
		
		final Element style = element("style");
		svgRoot.appendChild(style);
		style.appendChild(text(
				"g.maing {stroke:black;stroke-width:0.5px;fill:none;}\n"+
				"text.title {stroke:none;fill:black;stroke-width:1px;text-anchor:middle;}\n"+
				"rect.frame {stroke:darkgray;fill:none;stroke-width:0.5px;}\n" + 
				"text.ruler {stroke:none;fill:black;stroke-width:1px;text-anchor:end;}\n"+
				"text.chromName {stroke:none;fill:black;stroke-width:1px;text-anchor:middle;}\n"+
				"polygon.cov0 {stroke:lightgray;fill:antiquewhite;stroke-width:0.5px;}\n"+
				"polygon.cov1 {stroke:lightgray;fill:beige;stroke-width:0.5px;}\n"+
				"rect.average {stroke:green;fill:green;stroke-width:0.5px;}\n"+
				"rect.median {stroke:red;fill:red;stroke-width:0.5px;}\n"+
				"line.ruler {stroke:darkgray;stroke-dasharray:4;fill:none;stroke-width:1;}\n"+
				""
				));

		final Element mainG = element("g");
		mainG.setAttribute("class","maing");
		svgRoot.appendChild(mainG);
		
		

		
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refPath);
		final List<ChromInfo> chromInfos = dict.getSequences().stream().
				filter(SR->SR.getSequenceLength()>=this.min_contig_length).
				filter(SR->StringUtils.isBlank(this.skipContigExpr) || !SR.getSequenceName().matches(this.skipContigExpr)).
				map(SR->new ChromInfo(SR)).
				collect(Collectors.toList());
		if(chromInfos.isEmpty()) {
			LOG.info("no valid chromosome was found in "+this.refPath);
			return -1;
			}
		final double pixelsBetweenChromosomes = Double.parseDouble(dynaParams.getOrDefault("distanceBetweenChromosomes", "1"));
		final double marginLeft = Double.parseDouble(dynaParams.getOrDefault("margin-left", "100"));
		final double marginRight = Double.parseDouble(dynaParams.getOrDefault("margin-right", "10"));
		final double marginTop = Double.parseDouble(dynaParams.getOrDefault("margin-top", "80"));
		final double marginBottom = Double.parseDouble(dynaParams.getOrDefault("margin-bottom", "10"));
		final double drawingWidth= this.dimension.width - (marginLeft+marginRight);
		final double drawingHeight = this.dimension.height - (marginTop+marginBottom);
		
		if(drawingWidth<=0) {
			LOG.error("Bad image dimension");
			return -1;
		}
		
		final Element g_chroms = element("g");
		g_chroms.setAttribute("transform", "translate("+marginLeft+","+marginTop+")");
		g_chroms.setAttribute("id", "chromosomes");
		mainG.appendChild(g_chroms);
		
		final long sumLength = chromInfos.stream().mapToLong(CI->CI.ssr.getSequenceLength()).sum();
		final double pixelsPerBase = (drawingWidth - (pixelsBetweenChromosomes*(chromInfos.size()-1)))/(double)sumLength;
		double x = marginRight;
		for(int i=0;i< chromInfos.size();i++) {
			final ChromInfo ci = chromInfos.get(i);
			ci.x = x;
			ci.width = ci.ssr.getSequenceLength() * pixelsPerBase;
			x=(ci.x+ci.width)+pixelsBetweenChromosomes;
			}
		
		
		final SamReaderFactory samReaderFactory = SamReaderFactory.
					makeDefault().
					referenceSequence(this.refPath).
					validationStringency(ValidationStringency.LENIENT)
					;
		final Path input = Paths.get(oneAndOnlyOneFile(args));
		maintitle.appendChild(text(input.getFileName().toString()));
		
		
		final Element g_label = element("text",input.getFileName().toString());
		g_label.setAttribute("class", "title");
		g_label.setAttribute("x", format(marginLeft+drawingWidth/2.0));
		g_label.setAttribute("y", format(marginTop/2.0));
		g_label.appendChild(element("title",input.toString()));
		mainG.appendChild(g_label);
		
		final ProgressFactory.Watcher<SAMSequenceRecord> progress= ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
		try(SamReader sr = samReaderFactory.open(input)) {
			if(!sr.hasIndex()) {
				LOG.error("input bam "+input+" is not indexed.");
				return -1;
				}
			final SAMFileHeader header = sr.getFileHeader();
						
			SequenceUtil.assertSequenceDictionariesEqual(dict, SequenceDictionaryUtils.extractRequired(header));
			for(final ChromInfo ci: chromInfos) {
				progress.apply(ci.ssr);
				int coverage[]  = new int[ci.ssr.getSequenceLength()];
				Arrays.fill(coverage, 0);
				try(SAMRecordIterator iter= sr.queryOverlapping(ci.ssr.getSequenceName(), 1, ci.ssr.getSequenceLength())) {
						while(iter.hasNext()) {
						final SAMRecord rec= iter.next();
						
						if(!SAMRecordDefaultFilter.accept(rec, this.min_mapq)) continue;
						
						int max_end1 = coverage.length;
						
						if(!this.disable_paired_overlap_flag && 
							rec.getReadPairedFlag() && 
							!rec.getMateUnmappedFlag() &&
							rec.getReferenceIndex().equals(rec.getMateReferenceIndex()) &&
							rec.getAlignmentStart() < rec.getMateAlignmentStart() &&
							rec.getAlignmentEnd() > rec.getMateAlignmentStart()
							) {
							max_end1 = rec.getMateAlignmentStart() - 1;
							}
						for(final AlignmentBlock block:rec.getAlignmentBlocks()) {
							final int pos1=block.getReferenceStart();
							final int len = block.getLength();
							for(int i=0;i< len;i++) {
								if(pos1+i-1>=0 && pos1 +i <= max_end1) {
									coverage[pos1 + i -1]++;
									}
								}
							}
						}
					}
				final Element g = element("g");
				g.setAttribute("transform", "translate("+ci.x+",0)");
				g.setAttribute("id", "chrom-"+ci.ssr.getSequenceName());
				g_chroms.appendChild(g);
				
				final DiscreteMedian<Integer> ctg_median = new DiscreteMedian<>();

				final Element polyline = element("polygon");
				g.appendChild(polyline);
				final List<Point2D> points  = new ArrayList<>((int)ci.width);
				polyline.setAttribute("class", "cov"+(ci.ssr.getSequenceIndex()%2));
				points.add(new Point2D.Double(0,drawingHeight));
				
				for(int i=0;i< coverage.length;++i) {
					ctg_median.add(coverage[i]);
					}
				
				for(x=0;x< ci.width;x++) {
					int x1 = (int)((x+0)*(1.0/pixelsPerBase));
					final int x2 = (int)Math.min(ci.ssr.getSequenceLength(),(x+1)*(1.0/pixelsPerBase));
					final DiscreteMedian<Integer> median = new DiscreteMedian<>();
					
					while(x1<x2 && x1 < coverage.length) {
						median.add(coverage[x1]);
						
						++x1;
					    }
					final OptionalDouble dbl = median.getTendency(this.percentile);
					if(dbl.isPresent()) {
						double dbl2 = dbl.getAsDouble();
						if(dbl2 > this.max_depth && this.cap_depth) dbl2=this.max_depth;
						points.add(new Point2D.Double(x,drawingHeight - (dbl2/this.max_depth)*drawingHeight));
						}
					
					}
				points.add(new Point2D.Double(ci.width,drawingHeight));
				
				polyline.setAttribute("points",points.stream().map(P->format(P.getX())+","+format(P.getY())).collect(Collectors.joining(" ")));
				polyline.appendChild(element("title", ci.ssr.getSequenceName()));
				
				for(int side=0;side<2;++side) {
					final OptionalDouble g_dbl = (side==0?ctg_median.getAverage():ctg_median.getMedian());
					if(g_dbl.isPresent() && g_dbl.orElse(0) < this.max_depth) {
						double h = drawingHeight - ((g_dbl.getAsDouble()/this.max_depth)*drawingHeight);
						final Element hline = element("rect");
						hline.appendChild(element("title",(side==0?"average:":"median:")+format(g_dbl.orElse(0))));
						g.appendChild(hline);
						hline.setAttribute("class",(side==0?"average":"median"));
						hline.setAttribute("x", "0");
						hline.setAttribute("y", format(h));
						hline.setAttribute("width", format(ci.width));
						hline.setAttribute("height", "1");
						}
					}
					

				
				final Element title = element("text", ci.ssr.getSequenceName());
				title.setAttribute("x", format(ci.width/2.0));
				title.setAttribute("y", format(-12.0));
				title.setAttribute("class", "chromName");
				g.appendChild(title);
				
				final Element rect= element("rect");
				rect.setAttribute("x", "0");
				rect.setAttribute("y", "0");
				rect.setAttribute("width", format(ci.width));
				rect.setAttribute("height", format(drawingHeight));
				rect.setAttribute("class", "frame");
				g.appendChild(rect);
				
				coverage=null;
				System.gc();
				}
			progress.close();
			
			final Element rulers= element("g");
			g_chroms.appendChild(rulers);
			int prev_y=-1;
			for(int i=1;i <= 10.0;i++) {
				final int y = (int)((i/10.0)*this.max_depth);
				if(y==prev_y ||y==0) continue;
				double h = drawingHeight - ((i/10.0)*drawingHeight);
				final Element hline = element("line");
				rulers.appendChild(hline);
				hline.setAttribute("class", "ruler");
				hline.setAttribute("x1", "-2");
				hline.setAttribute("y1", format(h));
				hline.setAttribute("x2", format(drawingWidth));
				hline.setAttribute("y2", format(h));

				final Element label = element("text",""+y);
				label.setAttribute("class", "ruler");
				label.setAttribute("x", "-5");
				label.setAttribute("y", format(h));
				rulers.appendChild(label);
				
				prev_y=y;
				}
			
			}
		final Transformer tr = TransformerFactory.newInstance().newTransformer();
		final Result result;
		
		if(this.outputFile!=null)
			{
			result = new StreamResult(this.outputFile.toFile());
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
		}
	}

public static void main(final String[] args) {
	new WGSCoveragePlotter().instanceMainWithExit(args);
	}

}
