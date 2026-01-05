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
package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import javax.xml.XMLConstants;
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
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.log.ProgressFactory;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.rdf.ns.RDF;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.svg.SVG;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.SequenceUtil;
/**
 BEGIN_DOC

## Files

Input is an indexed BAM or CRAM file

Output is a SVG file

## Example
```
java -jar dist/jvarkit.jar wgscoverageplotter --dimension 1500x500 -C -1 --clip -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam --include-contig-regex "RF.*" --percentile median  > ~/jeter.svg
```

## Cited in

 * Vinicius A.C & al. Comparative analyses of Theobroma cacao and T. grandiflorum mitogenomes reveal conserved gene content embedded within complex and plastic structures". Gene Volume 849 2023 .https://doi.org/10.1016/j.gene.2022.146904.
 * Varani AM, Silva SR, Lopes S, Barbosa JBF, Oliveira D, Correa MA, Moraes AP, Miranda VFO, Prosdocimi F. 2022. The complete organellar genomes of the entheogenic plant Psychotria viridis (Rubiaceae), a main component of the ayahuasca brew. PeerJ 10:e14114 https://doi.org/10.7717/peerj.14114
 * Cho, C.H., Park, S.I., Huang, TY. et al. Genome-wide signatures of adaptation to extreme environments in red algae. Nat Commun 14, 10 (2023). https://doi.org/10.1038/s41467-022-35566-x
 * Rhys T. White, Martina Jelocnik, Natalie Klukowski, Md. Hakimul Haque, Subir Sarker. The first genomic insight into Chlamydia psittaci sequence type (ST)24 from a healthy captive psittacine host in Australia demonstrates evolutionary proximity with strains from psittacine, human, and equine hosts. Veterinary Microbiology, Volume 280, 2023, 109704, ISSN 0378-1135. https://doi.org/10.1016/j.vetmic.2023.109704.
 * Rhys T. White. A discovery down under: decoding the draft genome sequence of Pantoea stewartii from Australia's Critically Endangered western ground parrot/kyloring (Pezoporus flaviventris). https://doi.org/10.1099/mgen.0.001101 . Microbial Genomics Vol9 Issue 9
 * genomic and transcriptomic insight into the evolution and breeding of turfgrass species, poa annua and bouteloua dactyloides. Thesis. Christopher W. Benson. https://etda.libraries.psu.edu/files/final_submissions/29057 **WGSCoveragePlotter**
 * Jo, J., Park, JS., Won, H. et al. The first Chromosomal-level genome assembly of Sageretia thea using Nanopore long reads and Pore-C technology. Sci Data 11, 959 (2024). https://doi.org/10.1038/s41597-024-03798-9
 * Tim P. Bean, Hannah Farley, Jennifer Nascimento-Schulze, Tim Regan,Scottish oyster mortality event and association with Vibrio aestuarianus, Aquaculture Reports, Volume 39, 2024, 102480, ISSN 2352-5134, https://doi.org/10.1016/j.aqrep.2024.102480. **WGSCoveragePlotter**
 * del Olmo V, Redondo-Rio A, Garcia AB, Limtong S, Saus E, Gabaldon T (2025) Insights into the origin, hybridisation and adaptation of Candida metapsilosis hybrid pathogens. PLoS Pathog 21(1): e1012864. https://doi.org/10.1371/journal.ppat.1012864
 * Survival of the "unfittest": clinical emergence of hyper-multidrug-resistant Nakaseomyces glabratus with rare nonfunctional Erg3 and Erg11 and severely impaired fitness. Hans Carolus & al. bioRxiv 2025.02.05.636719; doi: https://doi.org/10.1101/2025.02.05.636719 
 * Kinnby, A., Robert, C., Havenhand, J.N., Brostrom, G., Bussiere, L. and De Wit, P. (2025), The Roles of Plasticity and Selection in Rapid Phenotypic Changes at the Pacific Oyster Invasion Front in Europe. Mol Ecol e17684. https://doi.org/10.1111/mec.17684
 * Dijkstra, J., A. C. van Westerhoven, L. Gomez-Gil, et al. 2025. "Accessory Chromosome Contributes to Virulence of Banana Infecting Fusarium oxysporum Tropical Race 4." Molecular Plant Pathology 26, no. 9: e70146. https://doi.org/10.1111/mpp.70146.
 * Wan J, Weldon E, Ganser G, Morris ERA, Hughes EV, Bordin AI, Heine PA, Hust M, Cohen ND, Gill JJ, Liu M. 0. Immunogenic Streptococcus equi cell surface proteins identified by ORFeome phage display. mSphere 0:e00626-25. https://doi.org/10.1128/msphere.00626-25 

## Screenshot

https://twitter.com/yokofakun/status/1331898068002861056

![twitter](https://pbs.twimg.com/media/EnvaOnNW4AAkGTz?format=jpg&name=medium "Screenshot")


END_DOC
 */
@Program(
	name="wgscoverageplotter",
	description="Whole genome coverage plotter",
	keywords={"svg","bam","depth","coverage"},
	creationDate="20201125",
	modificationDate="20230505",
	biostars={104063,475162,9536274},
	jvarkit_amalgamion =  true,
	nfcore = "https://nf-co.re/modules/jvarkit_wgscoverageplotter.html",
	menu="BAM Visualization"
	)
public class WGSCoveragePlotter extends Launcher {
	private static final Logger LOG = Logger.of( WGSCoveragePlotter.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refPath = null;
	@Parameter(names={"--min-contig-length"},description="Skip chromosome with length < 'x'. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_contig_length = 0;
	@Parameter(names={"-X","--skip-contig-regex"},description="Skip chromosomes matching this regular expression. Ignore if blank.")
	private String skipContigExpr = "";
	@Parameter(names={"-I","--include-contig-regex"},description="Only keep chromosomes matching this regular expression. Ignore if blank.")
	private String includeContigExpr = "";
	@Parameter(names={"--percentile"},description="How to we bin the coverage under one pixel.")
	private DiscreteMedian.Tendency percentile = DiscreteMedian.Tendency.median;
	@Parameter(names={"--mapq"},description = "min mapping quality")
	private int min_mapq=1;
	@Parameter(names={"-C","--max-depth"},description = "Max depth to display. The special value '-1' will first compute the average depth and the set the max depth to 2*average")
	private int max_depth=100;
	@Parameter(names={"--clip","--cap"},description = "Don't show coverage to be greater than 'max-depth' in the SVG file.")
	private boolean cap_depth= false;

	@Parameter(names={"--dimension"},description = "Image Dimension. " + DimensionConverter.OPT_DESC, converter=DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,500);
	@DynamicParameter(names = "-D", description = "set some css style elements. '-Dkey=value'. Undocumented.")
	private Map<String, String> dynaParams = new HashMap<>();
	@Parameter(names={"--disable-paired-overlap"},description="Count overlapping bases with mate for paired-end")
	private boolean disable_paired_overlap_flag=false;
	@Parameter(names={"--points"},description="Plot the coverage using points instead of areas.")
	private boolean plot_using_points =false;
	@Parameter(names={"--partition"},description="When using the option --samples, use this partition "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;
	@Parameter(names={"--samples"},description="Limit to those groups. See also --partition. Multiple separated with commas.")
	private String limitSamples = "";

	
	
	private Document document = null;
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");


	
	
private class ChromInfo {
	final SAMSequenceRecord ssr;
	int idx = 0;
	double x = 0;
	double width = 0;
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
private Element rdf(final String tag) {
	return this.document.createElementNS(RDF.NS, "rdf:"+tag);
	}
private Text text(final Object o) {
	return this.document.createTextNode(o==null?"":String.valueOf(o));
	}
private Element element(final String tag,final Object content) {
	final Element E = element(tag);
	E.appendChild(text(content));
	return E;
	}
/** get coverage in given chromosome for all bases of the chromosome */
private int[] getContigCoverage(final SamReader sr, final SAMSequenceRecord ssr,final SamRecordFilter samReadFilter) throws IOException {
	final int[] coverage;
	try {
		coverage = new int[ssr.getSequenceLength()];
	} catch(final Throwable err) {
		throw new OutOfMemoryError("Error: cannot alloc memory for sizeof(int)*"+ssr.getSequenceLength()+" (chrom "+ssr.getSequenceName()+")");
		}
	Arrays.fill(coverage, 0);
	try(SAMRecordIterator iter= sr.queryOverlapping(ssr.getSequenceName(), 1, ssr.getSequenceLength())) {
			while(iter.hasNext()) {
			final SAMRecord rec= iter.next();
			if(samReadFilter.filterOut(rec)) continue;
			
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
	return coverage;
	}

@Override
public int doWork(final List<String> args) {
	try
		{
		if(!StringUtils.isBlank(this.skipContigExpr) && !StringUtils.isBlank(this.includeContigExpr)) {
			LOG.error("Both include/exclude patterns are not blank.");
			return -1;
			}
		
		final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		dbf.setNamespaceAware(true);
		final DocumentBuilder db = dbf.newDocumentBuilder();
		this.document=db.newDocument();
		
		
		final Element svgRoot = element("svg");
		this.document.appendChild(svgRoot);
		svgRoot.setAttribute("width",format(dimension.width));
		svgRoot.setAttribute("height",format(dimension.height));
		svgRoot.setAttributeNS(XMLConstants.XML_NS_URI, XMLConstants.XML_NS_PREFIX + ":base","https://umr1087.univ-nantes.fr/");
			
		final Element maintitle = element("title");
		svgRoot.appendChild(maintitle);
		
		final Element metaE = element("metadata");
		metaE.setAttribute("id", "metadata");
		svgRoot.appendChild(metaE);
		final Element rdfRDF = rdf("RDF");
		metaE.appendChild(rdfRDF);
		final Element rdfDesc = rdf("Description");
		rdfDesc.setAttributeNS(RDF.NS, "rdf:ID", "#");
		rdfRDF.appendChild(rdfDesc);
		
		
		final Element descr = element("desc");
		svgRoot.appendChild(descr);
		descr.appendChild(text("Made with "+ getProgramName()+" "+getVersion()+" Author: Pierre Lindenbaum"));
		
		final Element style = element("style");
		svgRoot.appendChild(style);
		style.appendChild(text(
				"g.maing {stroke:black;stroke-width:0.5px;fill:none;}\n"+
				"text.title {stroke:none;fill:black;stroke-width:1px;text-anchor:middle;}\n"+
				"rect.frame {stroke:darkgray;fill:none;stroke-width:0.5px;}\n" + 
				"text.ruler {stroke:none;fill:darkgray;stroke-width:1px;text-anchor:end;}\n"+
				"text.chromName {stroke:none;fill:darkgray;stroke-width:1px;text-anchor:middle;}\n"+
				"text.xLabel {stroke:none;fill:black;stroke-width:1px;text-anchor:middle;}\n"+
				"text.yLabel {stroke:none;fill:black;stroke-width:1px;text-anchor:middle;}\n"+
				".cov0 {stroke:"+dynaParams.getOrDefault("stroke-cov0",plot_using_points?"slategray":"gray")+";fill:"+dynaParams.getOrDefault("fill-cov0","antiquewhite")+";stroke-width:0.5px;}\n"+
				".cov1 {stroke:"+dynaParams.getOrDefault("stroke-cov1",plot_using_points?"#61666B":"gray")+";fill:"+dynaParams.getOrDefault("fill-cov0","beige")+";stroke-width:0.5px;}\n"+
				"rect.average {stroke:green;fill:green;stroke-width:0.5px;opacity:"+dynaParams.getOrDefault("average-opacity","0.3")+";}\n"+
				"rect.median {stroke:red;fill:red;stroke-width:0.5px;opacity:"+dynaParams.getOrDefault("median-opacity","0.3")+";}\n"+
				"line.ruler {stroke:darkgray;stroke-dasharray:4;fill:none;stroke-width:1;}\n"+
				""
				));

		if(this.plot_using_points) {
			final int dz = Integer.parseInt(dynaParams.getOrDefault("point-size","5"));
			final Element defsE = element("defs");
			svgRoot.appendChild(defsE);
			final Element g = element("g");
			defsE.appendChild(g);

			g.setAttribute("id", "cross");
			Element L= element("line");
			L.setAttribute("x1",format(-dz));L.setAttribute("x2",format(dz));
			L.setAttribute("y1",format(0));L.setAttribute("y2",format(0));
			g.appendChild(L);
			L= element("line");
			L.setAttribute("x1",format(0));L.setAttribute("x2",format(0));
			L.setAttribute("y1",format(-dz));L.setAttribute("y2",format(dz));
			g.appendChild(L);
			}
				
		final Element mainG = element("g");
		mainG.setAttribute("class","maing");
		svgRoot.appendChild(mainG);
		
		

		
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refPath);
		final List<ChromInfo> chromInfos = dict.getSequences().stream().
				filter(SR->SR.getSequenceLength()>=this.min_contig_length).
				filter(SR->StringUtils.isBlank(this.skipContigExpr) || !SR.getSequenceName().matches(this.skipContigExpr)).
				filter(SR->StringUtils.isBlank(this.includeContigExpr) || SR.getSequenceName().matches(this.includeContigExpr)).
				map(SR->new ChromInfo(SR)).
				collect(Collectors.toList());
		if(chromInfos.isEmpty()) {
			LOG.error("no valid chromosome was found in "+this.refPath);
			return -1;
			}
		else
			{
			LOG.info("Will plot the following contigs:" + chromInfos.stream().
					map(CI->CI.ssr.getContig()).
					collect(Collectors.joining(", ")));
			}
		final double pixelsBetweenChromosomes = Double.parseDouble(dynaParams.getOrDefault("distanceBetweenChromosomes", "1"));
		final double marginLeft = Double.parseDouble(dynaParams.getOrDefault("margin-left", "100"));
		final double marginRight = Double.parseDouble(dynaParams.getOrDefault("margin-right", "10"));
		final double marginTop = Double.parseDouble(dynaParams.getOrDefault("margin-top", "80"));
		final double marginBottom = Double.parseDouble(dynaParams.getOrDefault("margin-bottom", "50"));
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
		
		final long genomeLength = chromInfos.stream().mapToLong(CI->CI.ssr.getSequenceLength()).sum();
		if(genomeLength<=0L) {
			throw new IllegalStateException("genome-length<=0 :"+genomeLength);
			}
		
		final double pixelsPerBase = (drawingWidth - (pixelsBetweenChromosomes*(chromInfos.size()-1)))/(double)genomeLength;
		if ( pixelsPerBase <=0) {
			throw new IllegalStateException("pixel per base <=0 with drawing-width:"+drawingWidth+" distance-between-contigs:"+pixelsBetweenChromosomes+" n-chromosomes:"+chromInfos.size()+" genome-length:"+genomeLength);
			}
		double x = marginRight;
		for(int i=0;i< chromInfos.size();i++) {
			final ChromInfo ci = chromInfos.get(i);
			ci.idx = i;
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
		
		
		
		try(SamReader sr = samReaderFactory.open(input)) {
			if(!sr.hasIndex()) {
				LOG.error("input bam "+input+" is not indexed.");
				return -1;
				}
			final SAMFileHeader header = sr.getFileHeader();
			
			final Set<String> limitRGIds;
			if(StringUtils.isBlank(this.limitSamples)) {
				limitRGIds = null;
				}
			else
				{
				final Set<String> sampleSet = Arrays.stream(this.limitSamples.split("[ ,]+")).
						filter(S->!StringUtils.isBlank(S)).
						collect(Collectors.toSet());
				limitRGIds = header.getReadGroups().stream().
					filter(RG->sampleSet.contains(this.groupBy.apply(RG, ""))).
					map(RG->RG.getId()).
					collect(Collectors.toSet());
				if(limitRGIds.isEmpty()) {
					LOG.error("samples are specified "+groupBy.name()+":("+this.limitSamples+") but I cannot find any in the read groups of "+input);
					return -1;
					}
				}
			
			final SamRecordFilter samReadFilter = new SamRecordFilter() {
				@Override
				public boolean filterOut(final SAMRecord rec) {
					if(!SAMRecordDefaultFilter.accept(rec, min_mapq)) return true;
					if(limitRGIds!=null) {
						final SAMReadGroupRecord rg = rec.getReadGroup();
						if(rg==null) return true;
						if(!limitRGIds.contains(rg.getId())) return true;
						}
					return false;
					}
				@Override
				public boolean filterOut(SAMRecord first, SAMRecord second) {
					return filterOut(first) && filterOut(second);
					}
				};
			
			//title is path + sample(s)
			final Element g_label = element("text",input.getFileName().toString() + " " + header.getReadGroups().stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toCollection(TreeSet::new)).stream().
					collect(Collectors.joining(" "))
					);
			mainG.appendChild(g_label);
			g_label.setAttribute("class", "title");
			g_label.setAttribute("x", format(marginLeft+drawingWidth/2.0));
			g_label.setAttribute("y", format(marginTop/2.0));
			g_label.appendChild(element("title",input.toString()));
			
			SequenceUtil.assertSequenceDictionariesEqual(dict, SequenceDictionaryUtils.extractRequired(header));
			
			
			if(this.max_depth==-1) {
				LOG.info("Computing "+ this.percentile.name() +" depth...");
				final ProgressFactory.Watcher<SAMSequenceRecord> progress= ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
				final DiscreteMedian<Integer> median = new DiscreteMedian<>();
				for(final ChromInfo ci: chromInfos) {
					int coverage[]  = getContigCoverage(sr, progress.apply(ci.ssr), samReadFilter);
					for(int i=0;i<coverage.length;i++) {
						median.add(coverage[i]);
						}
					coverage=null;
					System.gc();
					}
				progress.close();
				final DiscreteMedian.Tendency t;
				
				if(this.percentile.equals(DiscreteMedian.Tendency.median)) {
					t = this.percentile;
				} else {
					// if percentile is min, everything will be at 0...
					t = DiscreteMedian.Tendency.average;
					}
				max_depth=median.isEmpty()?0:(int)(Double.parseDouble(dynaParams.getOrDefault("factor-y", "2.0"))*median.getTendency(t).orElse(0.0));
				if(max_depth<=0) max_depth = 1;
				LOG.info("Now using max depth for y axis="+this.max_depth);
				}
			
			
			if(max_depth<1) {
				LOG.error("Bad user 'max-depth':"+max_depth);
				return -1;
				}
			
			final ProgressFactory.Watcher<SAMSequenceRecord> progress= ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
			for(final ChromInfo ci: chromInfos) {
				progress.apply(ci.ssr);
				int coverage[]  = getContigCoverage(sr, ci.ssr, samReadFilter);

				final Element g = element("g");
				g.setAttribute("transform", "translate("+ci.x+",0)");
				g.setAttribute("id", "chrom-"+ci.ssr.getSequenceName());
				g_chroms.appendChild(g);
				
				final DiscreteMedian<Integer> ctg_median = new DiscreteMedian<>();

				final Element polyline = element(this.plot_using_points?"g":"polygon");
				g.appendChild(polyline);
				final List<Point2D> points  = new ArrayList<>((int)ci.width);
				polyline.setAttribute("class", "cov"+(ci.idx%2));
				if(!plot_using_points) points.add(new Point2D.Double(0,drawingHeight));
				
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
						if(dbl2 > this.max_depth && this.cap_depth) {
							if(this.plot_using_points) continue;
							dbl2=this.max_depth;
						}
						points.add(new Point2D.Double(x,drawingHeight - (dbl2/this.max_depth)*drawingHeight));
						}
					
					}
				if(!plot_using_points)  points.add(new Point2D.Double(ci.width,drawingHeight));
				
				if(this.plot_using_points) {
					for(final Point2D pt: points) {
						final Element use = element("use");
						use.setAttribute("href", "#cross");
						use.setAttribute("x",format(pt.getX()));
						use.setAttribute("y",format(pt.getY()));
						polyline.appendChild(use);
						}
					}
				else {
					//simplify, remove points on the same line y
					int i=1;
					while(i+3< points.size())/* ignore first and last */ {
						final Point2D p1 = points.get(i+0);
						final Point2D p2 = points.get(i+1);
						final Point2D p3 = points.get(i+2);
						if(p1.getY()==p2.getY() && p2.getY()==p3.getY()) {
							points.remove(i+1);
							}
						else
							{
							i++;
							}
						}
					
					polyline.setAttribute("points",points.stream().map(P->format(P.getX())+","+format(P.getY())).collect(Collectors.joining(" ")));
					polyline.appendChild(element("title", ci.ssr.getSequenceName()));
					}
				
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
				title.appendChild(element("title",ci.ssr.getSequenceName()+" "+StringUtils.niceInt(ci.ssr.getSequenceLength())+" bp"));
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
				final double h = drawingHeight - ((i/10.0)*drawingHeight);
				final Element hline = element("line");
				rulers.appendChild(hline);
				hline.setAttribute("class", "ruler");
				hline.setAttribute("x1", "-2");
				hline.setAttribute("y1", format(h));
				hline.setAttribute("x2", format(drawingWidth));
				hline.setAttribute("y2", format(h));
				hline.appendChild(element("title",StringUtils.niceInt(y)));

				final Element label = element("text",StringUtils.niceInt(y));
				label.setAttribute("class", "ruler");
				label.setAttribute("x", "-5");
				label.setAttribute("y", format(h));
				rulers.appendChild(label);
				
				prev_y=y;
				}
			final Element xLabel = element("text",this.refPath.getFileName().toString()+" "+
					SequenceDictionaryUtils.getBuildName(dict).orElse("")+
					" ("+ StringUtils.niceInt(genomeLength)+" bp)");
			xLabel.setAttribute("class", "xLabel");
			xLabel.setAttribute("x", format(marginLeft+drawingWidth/2.0));
			xLabel.setAttribute("y", format(this.dimension.height - marginBottom/2.0));
			xLabel.appendChild(element("title",this.refPath.toString()));
			mainG.appendChild(xLabel);
			
			
			final Element yLabel = element("text","Coverage bp ("+this.percentile.name()+")");
			yLabel.setAttribute("class", "yLabel");
			yLabel.setAttribute("x","0");
			yLabel.setAttribute("y","0");
			yLabel.setAttribute("transform", "translate("+format(marginLeft/4.0)+","+format(marginTop+drawingHeight/2.0)+") rotate(90)");
			mainG.appendChild(yLabel);
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
		document=null;
		}
	}

public static void main(final String[] args) {
	new WGSCoveragePlotter().instanceMainWithExit(args);
	}

}
