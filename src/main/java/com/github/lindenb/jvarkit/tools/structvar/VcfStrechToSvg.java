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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.Text;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.CoverageFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Exon;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.util.vcf.AFExtractorFactory;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

## Motivation

another VCF to SVG

plot  allele ratio for each sample for each diallelic variant in a given region

colors define the genotype type : HOM_REF, HET, HOM_VAR

the size the dots represent the inverse internal allele frequency.

## Example

```
java -jar dist/vcfstrech2svg.jar --bed intervals.bed  -o TMP indexed.vcf.gz
```

## Screenshot

![twitter](https://pbs.twimg.com/media/EvtRyyXWEAEqpz7?format=jpg&name=small "Screenshot")

[https://twitter.com/yokofakun/status/1367778079813341185](https://twitter.com/yokofakun/status/1367778079813341185)


END_DOC
*/
@Program(name="vcfstrech2svg",
description="another VCF to SVG",
keywords={"vcf","deletion","cnv","svg"},
creationDate="20210304",
modificationDate="20210309"
)
public class VcfStrechToSvg extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfStrechToSvg.class).make();
	/** which format to use */
	private enum FORMAT {PL ,AD}
	
	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-r","--region","--bed"},description="BED File",required=true)
	private Path bedFile = null;
	@Parameter(names= {"--gzip"},description="Generate gzipped compressed svg files.")
	private boolean compressed_svg=false;
	@Parameter(names={"-w","--width"},description="image width.")
	private int image_width_pixel  = 1_000;
	@Parameter(names={"--gq"},description="minimum FORMAT/GQ")
	private int minGQ  = 1;
	@Parameter(names={"--dp"},description="minimum sum(FORMAT/AD[0]+FORMAT/AD[1]). Doesn't work with FORMAT/PL.")
	private int minDP  = 1;
	@Parameter(names={"--keep-filtered"},description="keep FILTERed variants")
	private boolean accept_filtered=false;
	@Parameter(names={"--pack-distance"},description="pack variant in the same area if they're close to 'x' bp. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int withinDistance  = 10_000;
	@Parameter(names={"--extend"},description="Extend each area with 'x' bp. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int extendSet  = 1_000;
	@Parameter(names={"--max-af"},description="Discard variant with an internal AF > 'x' "+FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private double maxAf  = 1.0;
	@Parameter(names={"--min-af"},description="Discard variant with an internal AF < 'x' "+FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private double minAF  = 0.0;
	@Parameter(names={"-u","--url","--hyperlink"},description= "creates a hyperlink an area is clicked. " + Hyperlink.OPT_DESC,converter=Hyperlink.StringConverter.class,splitter=NoSplitter.class)
	private Hyperlink hyperlinkType = Hyperlink.empty();
	@Parameter(names={"--samples"},description= "Limit to those samples. (comma or space separated). If the list starst with '^' then the samples are excluded.")
	private String sampleStr= null;
	@Parameter(names= {"--no-tooltip"},description="remove contextual tooltips (reduce the size of the svg)")
	private boolean remove_tooltip=false;
	@Parameter(names= {"--gtf"},description="Plot gene structure using this GTF file.")
	private Path gff3Path = null;
	@Parameter(names= {"--manifest"},description="Output BED manifest")
	private Path manifestPath = null;
	@Parameter(names= {"--color-tag"},description="specify the optional INFO/tag defining a named svg color. If defined for a variant, a vertical line with the color will be painted. ")
	private String colorTag = null;
	@Parameter(names= {"--reference","-R"},description=CRAM_INDEXED_REFENCE)
	private Path refFaidx=null;
	@Parameter(names= {"--bam-list","--bam","--bams"},description="Add BAM/CRAM for plotting the depth. A file with the suffix '.list' is a list of path to the bams/crams.")
	private List<String> bamListPath=new ArrayList<>();
	@Parameter(names= {"--format"},description="wich format to use to calculate the allele depth ratio.")
	private FORMAT useFormat = FORMAT.AD;
	@Parameter(names= {"-hr","--hom-ref"},description="Hide HOM_REF genotypes (0/0)")
	private boolean hide_hom_ref=false;
	@Parameter(names= {"--af"},description=AFExtractorFactory.OPT_DESC)
	private String afExtractorFactoryStr= AFExtractorFactory.FORMAT_GT;

	
	@SuppressWarnings("serial")
	@DynamicParameter(names = "--param", description = "Other parameters. Undocumented")
	public Map<String, String> dynamicParams = new HashMap<String,String>() {{{
		put("gt.r1","1");
		put("gt.r2","7");
		put("sample.height","50");
		}}};

	
	private final DecimalFormat decimalFormater = new DecimalFormat("##.#");
	private Document document = null;
	private final Map<String, Path> sample2bam = new HashMap<>();
	private AFExtractorFactory.AFExtractor afExtractor = null;
	private IntervalTreeMap<Gene> all_genes = new IntervalTreeMap<>();
	
	private class VariantSet implements Locatable {
		private final List<VariantContext> variants = new ArrayList<>();
		double x = 0.0;
		double width = 0.0;
		CoverageFactory.SimpleCoverage coverage = null;
		
		VariantSet(final VariantContext vc) {
			this.variants.add(vc);
			}
		@Override
		public String getContig() {
			return variants.get(0).getContig();
			}
		@Override
		public int getStart() {
			return Math.max(0,variants.get(0).getStart() - extendSet);
			}
		@Override
		public int getEnd() {
			return variants.stream().mapToInt(V->V.getEnd()).max().getAsInt() + extendSet ;
			}
		@Override
		public String toString() {
			return getContig()+":"+StringUtils.niceInt(getStart())+"-"+StringUtils.niceInt(getEnd()) +" N="+variants.size();
			}
		public ToDoubleFunction<Integer> createBaseToPixelFunction() {
			return POS-> ((POS-this.getStart())/(double)this.getLengthOnReference()) * this.width;
			}
		}
	/** wrape node into a genomic hyperlink if needed */
	private Node wrapLoc(final Element node,final Locatable loc) {
		if(loc==null) return node;
		final String url = this.hyperlinkType.apply(loc);
		if(StringUtils.isBlank(url)) return node;
		final Element a = element("a");
		a.setAttribute("target","_blank");
		a.setAttribute("href",url);
		a.appendChild(node);
		return a;
		}

	/** creates a SVG element */
	private Element element(final String tag) {
		return this.document.createElementNS(SVG.NS, tag);
		}
	/** creates a DOM Text element */
	private Text text(final Object o) {
		return this.document.createTextNode(o==null?"":String.valueOf(o));
		}
	/** creates a SVG element with a simple text */
	private Element element(final String tag,final Object content) {
		final Element E = element(tag);
		E.appendChild(text(content));
		return E;
		}
	/** format a double number */
	private String format(double v) {
		return this.decimalFormater.format(v);
	}
	
	/** get allele frequency for this variant */
	private OptionalDouble getAF(final VariantContext ctx) {
		return this.afExtractor.getMinAlleleFrequency(ctx);
		}
	
	/** shall we accept a variant */
	private boolean acceptVariant(final VariantContext V) {
		if(!accept_filtered && V.isFiltered()) return false;
		if(!V.isBiallelic()) return false;
		final OptionalDouble af = getAF(V);
		if(!af.isPresent()) return false;
		if(af.getAsDouble() > this.maxAf) return false;
		if(af.getAsDouble() < this.minAF) return false;
		return true;
		}
	
	private void run(
			final ArchiveFactory archive,
			final BedLine bed,
			final VCFHeader header,
			final VCFReader in,
			final PrintWriter manifest
			) {
		LOG.info("processing "+bed);
		final Set<String> limitSamples;
		if(StringUtils.isBlank(this.sampleStr)) {
			limitSamples = new TreeSet<>(header.getSampleNamesInOrder());
			}
		else if(this.sampleStr.startsWith("^")) {
			limitSamples = new TreeSet<>(header.getSampleNamesInOrder());
			limitSamples.removeAll(Arrays.stream(this.sampleStr.substring(1).split("[, ]+")).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toSet()));
			}
		else
			{
			limitSamples = Arrays.stream(this.sampleStr.split("[, ]+")).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toCollection(TreeSet::new));
			}
		final SAMSequenceDictionary dict=header.getSequenceDictionary();
		try(CloseableIterator<VariantContext> iter=in.query(bed)){
			final List<VariantSet> L = iter.stream().
					filter(V->acceptVariant(V)).
					map(V->new VariantSet(V)).
					collect(Collectors.toCollection(ArrayList::new));
			if(L.isEmpty()) {
				LOG.warn("No valid variant found for \""+bed+"\"");
				return;
				}
			int i=0;
			while(i +1 < L.size()) {
				if(L.get(i).withinDistanceOf(L.get(i+1), this.withinDistance)) {
					L.get(i).variants.addAll(L.get(i+1).variants);
					L.remove(i+1);
					}
				else
					{
					i++;
					}
			}
		
		
		final int margin_left= 50;
		final int margin_right= 10;
		final double drawingAreaWidth = image_width_pixel-(margin_left+margin_right);
		final int intervalLength = L.stream().mapToInt(V->V.getLengthOnReference()).sum();
		double x=0;
		for(i=0;i< L.size();i++) {
			L.get(i).x = x;
			x+= (L.get(i).getLengthOnReference()/(double)intervalLength)*drawingAreaWidth;
			}
		for(i=0;i< L.size();i++) {
			L.get(i).width = (i+1<L.size()?L.get(i+1).x:drawingAreaWidth)-L.get(i).x;
			}
		
		try {
			final DocumentBuilderFactory db = DocumentBuilderFactory.newInstance();
			final DocumentBuilder dom = db.newDocumentBuilder();
			this.document = dom.newDocument();
			final Element svgRoot = element("svg");
			this.document.appendChild(svgRoot);
			final String mainTitleStr = SequenceDictionaryUtils.getBuildName(dict).orElse("")+" "+
					new SimpleInterval(bed).toNiceString()+" length:"+StringUtils.niceInt(bed.getLengthOnReference());
			/* SVG title */
			{
			final Element title = element("title");
			svgRoot.appendChild(title);
			title.appendChild(text(mainTitleStr));
			}

			/* SVG style */
			{
			final String gtopacity = this.dynamicParams.getOrDefault("gt.opacity","0.7");
			final Element style = element("style");
			svgRoot.appendChild(style);
			style.appendChild(text(
					".maintitle {text-anchor:middle;fill:blue} "+
					".vc {stroke-width:0.5px;} "+
					".transcript {stroke:black;stroke-width:1px;opacity:1;}"+
					".exon {stroke:black;stroke-width:0.5px;fill:blue;opacity:1;}"+
					".sample {fill:blue;font-size:7px;} "+
					".samplelabel {stroke:gray;stroke-width:0.5px;font-size:"+this.dynamicParams.getOrDefault("sample.fontsize","7")+"px;}\n" +
					".coverage { fill:gray; stroke:yellow;opacity:0.2;} " +
					".frame { fill:none; stroke: darkgray;} " +
					".area0 {fill:white;}\n" +
					".area1 {fill:floralwhite;}\n" +
					"circle.HOM_REF {fill:green;opacity:"+gtopacity+";stroke-width:0.5px;}\n" +
					"circle.HET {fill:blue;opacity:"+gtopacity+";stroke-width:0.5px;}\n" +
					"circle.HOM_VAR {fill:red;opacity:"+gtopacity+";stroke-width:0.5px;}\n" +
					"a {cursor: pointer;}\n"
					));
			}
			/* desc */
			{
			final Element descr = element("desc");
			svgRoot.appendChild(descr);
			descr.appendChild(text("Author: Pierre Lindenbaum\n" +
					JVarkitVersion.getInstance().getCompilationDate()+"\n"+
					JVarkitVersion.getInstance().getGitHash()
					));
			
			}
			
			
			// main title
			{
			final Element gtitle= element("text",mainTitleStr);
			gtitle.setAttribute("class", "maintitle");
			gtitle.setAttribute("x", format(this.image_width_pixel/2.0));
			gtitle.setAttribute("y", "15");
			svgRoot.appendChild(wrapLoc(gtitle,bed));
			}
			
			
			
			int margin_top= 50;
			double y = margin_top;

			final double min_circle_radius = Double.parseDouble(this.dynamicParams.getOrDefault("gt.r1","1"));
			final double max_circle_radius = Double.parseDouble(this.dynamicParams.getOrDefault("gt.r2","7"));
			final Element main_g = element("g");
			svgRoot.appendChild(main_g);
			
			
			/** plot genes */
			if(!this.all_genes.isEmpty())
				{
				final double transcript_height = 5;
				final double exon_height = (transcript_height-1);
				final double save_y = y;
				final Element g_genes= element("g");
				g_genes.setAttribute("transform", "translate("+margin_left+",0)");
				main_g.appendChild(g_genes);

				/* loop over each vset */
				for(i=0;i< L.size();i++) {
					final VariantSet vset = L.get(i);
					// find transcript in this vset
					final List<Transcript> transcripts = this.all_genes.getOverlapping(vset).
							stream().
							flatMap(G->G.getTranscripts().stream()).
							filter(T->T.overlaps(vset)).
							collect(Collectors.toList());
					if(transcripts.isEmpty()) continue;
					final Element g_vset = element("g");
					g_vset.setAttribute("transform", "translate("+format(vset.x)+",0)");
					g_genes.appendChild(g_vset);
					// y in this vset
					double y2 = save_y;
					
					/* convert base to pixel */
					final ToDoubleFunction<Integer> base2pixel = vset.createBaseToPixelFunction();

					/* loop over transcripts */
					for(final Transcript tr: transcripts) {
						final Element g_tr = element("g");
						g_tr.setAttribute("transform", "translate(0,"+format(y2)+")");
						g_vset.appendChild(g_tr);
						final Element line = element("line");
						line.setAttribute("class", "transcript");
						line.setAttribute("x1", format(Math.max(0, base2pixel.applyAsDouble(tr.getStart()))));
						line.setAttribute("y1", format(transcript_height/2.0));
						line.setAttribute("x2", format(Math.min(vset.width, base2pixel.applyAsDouble(tr.getEnd()))));
						line.setAttribute("y2", format(transcript_height/2.0));
						line.appendChild(element("title",tr.getId()));
						g_tr.appendChild(wrapLoc(line,tr));
						
						/* loop over exons */
						for(final Exon exon:tr.getExons()) {
							if(!exon.overlaps(vset)) continue;
							final Element exRect = element("rect");
							exRect.setAttribute("class", "exon");
							final double x_start = Math.max(0, base2pixel.applyAsDouble(exon.getStart()));
							exRect.setAttribute("x", format(x_start));
							exRect.setAttribute("y", format(transcript_height/2.0 - exon_height/2.0));
							final double x_end = Math.min(vset.width,base2pixel.applyAsDouble(exon.getEnd()));
							exRect.setAttribute("width", format(x_end  - x_start));
							exRect.setAttribute("height", format(exon_height));
							exRect.appendChild(element("title",exon.getName()));
							g_tr.appendChild(wrapLoc(exRect,exon));
							}
						y2+= transcript_height+0.5;
						}
					y= Math.max(y, y2);
					}
				y++;
				}

			
			final double sample_height= Double.parseDouble(this.dynamicParams.getOrDefault("sample.height","25"));
			final double sample_height2 = sample_height - (max_circle_radius*2.0);
			
			int space_between_samples =2;
			int got_n_samples = 0;
			for(final String sn: header.getSampleNamesInOrder()) {
				if(!limitSamples.contains(sn)) continue;
				
				boolean got_this_sample = false;
				final Element g_sample = element("g");
				g_sample.setAttribute("transform", "translate("+margin_left+","+format(y)+")");
				
				/* get coverage */
				final int maxCoverage;
				if(this.sample2bam.containsKey(sn)) {
					final CoverageFactory coverageFactory = new CoverageFactory();
					
					
					try(SamReader sr= this.openSamReader(this.sample2bam.get(sn))) {
						/* loop over each variant set */
						for(final VariantSet vset: L) {
							vset.coverage = coverageFactory.getSimpleCoverage(sr, vset, sn);
							}
						}
					maxCoverage = L.stream().
							flatMapToInt(V->V.coverage.stream()).
							max().
							orElse(0);
					}
				else
					{
					maxCoverage = 0;
					for(final VariantSet vset: L) {
						vset.coverage=null;
						}
					}

				/* loop over each variant set */
				for(i=0;i< L.size();i++) {
					final VariantSet vset = L.get(i);
					final Element g_vset = element("g");
					g_vset.setAttribute("transform", "translate("+format(vset.x)+",0)");
					g_sample.appendChild(g_vset);
					
					/* convert base to pixel */
					final ToDoubleFunction<Integer> base2pixel = vset.createBaseToPixelFunction();
					
					// plot set length
					final Element rect  = element("rect");
					rect.setAttribute("class", "area"+(i%2));
					rect.setAttribute("x", "0");
					rect.setAttribute("y", "0");
					rect.setAttribute("width",format(vset.width));
					rect.setAttribute("height", format(sample_height));
					if(!remove_tooltip) rect.appendChild(element("title",vset.toString()));
					g_vset.appendChild(rect);
										
					// plot coverage
					if(maxCoverage>0 && this.sample2bam.containsKey(sn)) {
						final double[] scaled = vset.coverage.scaleAverage((int)vset.width);
						final StringBuilder sb=new StringBuilder();
						sb.append("0,"+sample_height);
						for(int t=0;t< scaled.length;t++) {
							if(t>1 && t+1 < scaled.length &&
									format(scaled[t-1]).equals(format(scaled[t+1])) &&
									format(scaled[t-1]).equals(format(scaled[t]))) continue;
							sb.append(" ").append(t).append(",");
							sb.append(format(sample_height*(1.0 - scaled[t]/maxCoverage)));
							}
						sb.append(" "+format(vset.width)+","+sample_height);
						final Element polyline  = element("polyline");
						polyline.setAttribute("class", "coverage");
						polyline.setAttribute("points",sb.toString());
						g_vset.appendChild(polyline);
						vset.coverage=null;
						}
					
					
					//plot vertical line if colorTag defined
					if(!StringUtils.isBlank(this.colorTag)) {
						for(final VariantContext vc: vset.variants) {
							if(!vc.hasAttribute(this.colorTag)) continue;
							final String cssColor = vc.getAttributeAsString(this.colorTag, "");
							if(StringUtils.isBlank(cssColor)) continue;
							final double x0 = base2pixel.applyAsDouble(vc.getStart());
							final Element line  = element("line");
							line.setAttribute("class", "vc");
							line.setAttribute("style", "stroke:"+cssColor);
							line.setAttribute("x1", format(x0));
							line.setAttribute("y1", "0");
							line.setAttribute("x2", format(x0));
							line.setAttribute("y2", format(sample_height));
							g_vset.appendChild(line);
							}
						}
					
					// print all variants in this vcfset for this sample
					for(final VariantContext vc: vset.variants) {
						final Genotype gt= vc.getGenotype(sn);
						if(gt.isNoCall()) continue;
						if(hide_hom_ref && gt.isHomRef()) continue;
						if(gt.hasGQ() && gt.getGQ() < this.minGQ) continue;
						final OptionalDouble alt_ratio  = getAltRatio(gt);
						if(!alt_ratio.isPresent()) continue;

						
						final OptionalDouble af = getAF(vc);
						final double circle_radius = min_circle_radius + (max_circle_radius - min_circle_radius)*(1.0-af.orElse(1.0));
						
						//  HOMREF=0;  HET =0.5; HOMVAR = 1;

						final double gtx = base2pixel.applyAsDouble(vc.getStart());
						final double gty= sample_height- ( sample_height2*alt_ratio.getAsDouble() + (sample_height-sample_height2)/2.0);
						final Element circle  = element("circle");
						circle.setAttribute("class",gt.getType().name());
						circle.setAttribute("cx", format(gtx));
						circle.setAttribute("cy", format(gty));
						circle.setAttribute("r", format(circle_radius));
						if(!remove_tooltip) circle.appendChild(element("title",vc.getStart()+" "+(vc.hasID()?vc.getID():"")+" "+
								vc.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/"))+" "+
								gt.getType().name()+" AF="+format(af.orElse(-1))));
						g_vset.appendChild(wrapLoc(circle,vc));
						got_this_sample = true;
						}
					}
				
				final Element frame_sample = element("rect");
				frame_sample.setAttribute("class", "frame");
				frame_sample.setAttribute("x", "0");
				frame_sample.setAttribute("y", "0");
				frame_sample.setAttribute("width",format(drawingAreaWidth));
				frame_sample.setAttribute("height", format(sample_height));
				g_sample.appendChild(frame_sample);
				
				final Element label = element("text",sn +( maxCoverage==0?"":" Max Cov. "+maxCoverage));
				label.setAttribute("class","samplelabel");
				label.setAttribute("x","0");
				label.setAttribute("y","0");
				//label.setAttribute("transform", "translate("+format(-10)+","+0+") rotate(90) ");
				label.setAttribute("transform", "translate(12,12)");
				if(!remove_tooltip) label.appendChild(element("title",sn));
				g_sample.appendChild(label);

				if(got_this_sample) {
					got_n_samples++;
					main_g.appendChild(g_sample);
					y+= sample_height + space_between_samples;
					}
				else
					{
					LOG.warn("no valid data for sample "+ sn+" in "+bed);
					}
				}
			// remove extra sample space
			y-=space_between_samples;
			
			svgRoot.setAttribute("width",format(this.image_width_pixel+1));
			svgRoot.setAttribute("height",format(y+1));

			
			if(got_n_samples==0) {
				LOG.info("no sample/genotype found for "+bed);
				return;
				}
			//save
			final Transformer tr = TransformerFactory.newInstance().newTransformer();
			final String filename =    bed.getContig()+"_"+bed.getStart()+"_"+bed.getEnd()+ ".svg"+(this.compressed_svg?".gz":"");
			LOG.info("writing "+filename);
			if(this.compressed_svg) {
				try(final OutputStream pw=archive.openOuputStream(filename)) {
					try(GZIPOutputStream gzout = new GZIPOutputStream(pw)) {
						tr.transform(new DOMSource(this.document),new StreamResult(gzout));
						gzout.finish();
						gzout.flush();
						}
					pw.flush();
					}
				}
			else
				{
				try(final PrintWriter pw=archive.openWriter(filename)) {
					tr.transform(new DOMSource(this.document),new StreamResult(pw));
					pw.flush();
					}
				}
			manifest.print(bed.getContig());
			manifest.print("\t");
			manifest.print(bed.getStart()-1);
			manifest.print("\t");
			manifest.print(bed.getEnd());
			manifest.print("\t");
			manifest.print(filename);
			manifest.println();
			}
		catch(final Throwable err) {
			throw new RuntimeException(err);
			}
		finally {
			this.document = null;
			}
		}
	}
	
	private OptionalDouble getAltRatio(final Genotype gt) {
	switch(this.useFormat) {
		case PL: {
			if(!gt.hasPL()) return OptionalDouble.empty();
			final int pl[]= gt.getPL();
			if( pl==null ||  pl.length!=3)  return OptionalDouble.empty();
			final double countR = pl[2 /* yes, inversion */]*2 + pl[1];
			final double countA = pl[0 /* yes, inversion */]*2 + pl[1];
			final double DP = countR + countA;
			final double alt_ratio  = countA/DP;
			return OptionalDouble.of(alt_ratio);
			}
		case AD:
		default:{
			if(!gt.hasAD()) return OptionalDouble.empty();
			final int ad[]= gt.getAD();
			if(ad==null || ad.length!=2)  return OptionalDouble.empty();
			
			final double countR = ad[0];
			final double countA = ad[1];
			final double DP = countR + countA;
			if(DP<=0 || DP <= this.minDP) return OptionalDouble.empty();;
			
			//  HOMREF=0;  HET =0.5; HOMVAR = 1;
			final double alt_ratio  = countA/DP;
			return OptionalDouble.of(alt_ratio);
			}
		}
	}
	

	/** open a sam file. Check it is indexed */
	private SamReader openSamReader(final Path path) throws IOException {
		final SamReader sr= super.createSamReaderFactory().
				referenceSequence(this.refFaidx).
				open(path);
		if(!sr.hasIndex()) {
			sr.close();
			throw new SAMException("File "+path+" is not indexed.");
			}
		return sr;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			this.afExtractor = new AFExtractorFactory().parseFieldExtractor(this.afExtractorFactoryStr);
			
			final String input = super.oneAndOnlyOneFile(args);
			
			if(!this.bamListPath.isEmpty()) {
				LOG.info("reading bam list");
				for(Path bamPath:IOUtils.unrollPaths(this.bamListPath)) {
					try(SamReader sr = openSamReader(bamPath)) {
						final SAMFileHeader hdr = sr.getFileHeader();
						for(final String sn:hdr.getReadGroups().
								stream().
								map(RG->RG.getSample()).
								filter(S->!StringUtils.isBlank(S)).collect(Collectors.toSet())) {
							if(this.sample2bam.containsKey(sn)) {
								LOG.info("duplicate sample in bam "+bamPath+" and "+ this.sample2bam.get(sn));
								return -1;
								}
							this.sample2bam.put(sn, bamPath);
							}
						}
					}
				}
			
			try (VCFReader r= VCFReaderFactory.makeDefault().open(input, true)) {
				final VCFHeader header=r.getHeader();
				this.afExtractor.validateHeader(header);
				
				final String searchFormat;
				switch(this.useFormat) {
					case PL: searchFormat = VCFConstants.GENOTYPE_PL_KEY;break;
					case AD://through
					default: searchFormat = VCFConstants.GENOTYPE_ALLELE_DEPTHS; break;
					}
				if(header.getFormatHeaderLine(searchFormat)==null) {
					LOG.error("FORMAT/"+searchFormat+" undefined in "+input);
					return -1;
					}
				
				if(!header.hasGenotypingData()) {
					LOG.error("No genotype in input");
					return -1;
					}
				final SAMSequenceDictionary dict = header.getSequenceDictionary();
				if(dict!=null && this.refFaidx!=null) {
					SequenceUtil.assertSequenceDictionariesEqual(dict, SequenceDictionaryUtils.extractRequired(this.refFaidx));
					}
				
				if(this.gff3Path!=null) {
					LOG.info("reading gtf" + this.gff3Path);
					try(GtfReader gtfReader = new GtfReader(this.gff3Path)) {
						if(dict!=null)gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
						gtfReader.getAllGenes().forEach(G->{
							this.all_genes.put(new Interval(G), G);
						});
					}
				}
				
				try(BufferedReader br = IOUtils.openPathForBufferedReading(this.bedFile)) {
					try(PrintWriter manifest = (this.manifestPath==null?
							new PrintWriter(new NullOuputStream()):
							IOUtils.openPathForPrintWriter(this.manifestPath))) {
						try(ArchiveFactory out=ArchiveFactory.open(this.outputFile)) {
							final BedLineCodec codec = new BedLineCodec();
							br.lines().
								map(L->codec.decode(L)).
								filter(B->B!=null).
								forEach(B->{
								run(out,B,header,r,manifest);
								});
							}
						manifest.flush();
						}
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			}
		}

	public static void main(final String[] args)
		{
		new VcfStrechToSvg().instanceMainWithExit(args);
		}
	}
