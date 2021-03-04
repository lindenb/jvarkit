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
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
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
import org.w3c.dom.Text;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

## Motivation

another VCF to SVG

## Example



## See also

```
bcftools roh
```

END_DOC
*/
@Program(name="vcfstrech2svg",
description="another VCF to SVG",
keywords={"vcf","deletion","cnv","svg"},
creationDate="20210304",
modificationDate="20210304",
generate_doc=false
)
public class VcfStrechToSvg extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfStrechToSvg.class).make();

	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-r","--region","--bed"},description="BED File",required=true)
	private Path bedFile = null;
	@Parameter(names= {"--gzip"},description="Generate gzipped compressed svg files.")
	private boolean compressed_svg=false;
	@Parameter(names={"-w","--width"},description="image width.")
	private int image_width_pixel  = 1_000;
	@Parameter(names={"-GQ"},description="minimum FORMAT/GQ")
	private int minGQ  = 1;
	@Parameter(names={"-DP"},description="minimum sum(FORMAT/AD[0]+FORMAT/AD[1])")
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


	@SuppressWarnings("serial")
	@DynamicParameter(names = "--param", description = "Other parameters. Undocumented")
	public Map<String, String> dynamicParams = new HashMap<String,String>() {{{
		put("gt.r","2");
		put("sample.height","50");
		}}};

	
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private Document document = null;
	
	private class VariantSet implements Locatable {
		private final List<VariantContext> variants = new ArrayList<>();
		double x = 0.0;
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
	
	private String format(double v) {
		return this.decimalFormater.format(v);
	}
	
	private void run(final ArchiveFactory archive,final BedLine bed,final VCFHeader header,final VCFReader in) {
		LOG.info("processing "+bed);
		final Predicate<VariantContext> acceptVariant= V->{
			if(!accept_filtered && V.isFiltered()) return false;
			if(!V.isBiallelic()) return false;
			if(this.maxAf<1.0 || this.minAF>0.0) {
				int ac=0;
				int an=0;
				for(final Genotype gt: V.getGenotypes())  {
					if(gt.isNoCall()) continue;
					for(final Allele a:gt.getAlleles()) {
						an++;
						if(!a.isReference()) ac++;
						}
					}
				if(an==0) return false;
				final double af = ac/(double)an;
				if(af > this.maxAf) return false;
				if(af < this.minAF) return false;
				}
			return true;
			};
		
		
		final SAMSequenceDictionary dict=header.getSequenceDictionary();
		try(CloseableIterator<VariantContext> iter=in.query(bed)){
			final List<VariantSet> L = iter.stream().
					filter(acceptVariant).
					map(V->new VariantSet(V)).
					collect(Collectors.toCollection(ArrayList::new));
			if(L.isEmpty()) return;
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
		int margin_left= 50;
		int margin_right= 10;
		final double drawingAreaWidth = image_width_pixel-(margin_left+margin_right);
		final int intervalLength = L.stream().mapToInt(V->V.getLengthOnReference()).sum();
		double x=0;
		for(i=0;i< L.size();i++) {
			L.get(i).x = x;
			x+= (L.get(i).getLengthOnReference()/(double)intervalLength)*drawingAreaWidth;
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
			final Element style = element("style");
			svgRoot.appendChild(style);
			style.appendChild(text(
					".maintitle {text-anchor:middle;fill:blue} "+
					".sample {fill:blue;font-size:7px;} "+
					".samplelabel {stroke:gray;stroke-width:0.5px;font-size:"+this.dynamicParams.getOrDefault("sample.fontsize","7")+"px;}\n" +
					".frame { fill:none; stroke: darkgray;} " +
					".area0 {fill:white;}\n" +
					".area1 {fill:floralwhite;}\n" +
					"circle.HOM_REF {fill:green;stroke-width:0.5px;}\n" +
					"circle.HET {fill:blue;stroke-width:0.5px;}\n" +
					"circle.HOM_VAR {fill:red;stroke-width:0.5px;}\n" +
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
			Element gtitle= element("text",mainTitleStr);
			gtitle.setAttribute("class", "maintitle");
			gtitle.setAttribute("x", format(this.image_width_pixel/2.0));
			gtitle.setAttribute("y", "15");
			svgRoot.appendChild(gtitle);
			}
			
			int margin_top= 50;
			double y = margin_top;
			final double circle_radius = Double.parseDouble(this.dynamicParams.getOrDefault("gt.r","2"));
			final Element main_g = element("g");
			svgRoot.appendChild(main_g);
			
			
			final double sample_height= Double.parseDouble(this.dynamicParams.getOrDefault("sample.height","25"));
			final double sample_height2 = sample_height - (circle_radius*2.0);
			
			int space_between_samples =2;
			boolean got_one_sample = false;
			for(final String sn:header.getSampleNamesInOrder()) {
				if(L.stream().flatMap(L2->L2.variants.stream()).map(V->V.getGenotype(sn)).noneMatch(G->G.hasAD() && (G.isHomRef() ||G.isHet() || G.isHomVar()))) {
					LOG.info("no data for "+sn);
					continue;
				}
				
				final Element g_sample = element("g");
				g_sample.setAttribute("transform", "translate("+margin_left+","+format(y)+")");
				main_g.appendChild(g_sample);
				
				final Element label = element("text",sn);
				label.setAttribute("class","samplelabel");
				label.setAttribute("x","0");
				label.setAttribute("y","0");
				label.setAttribute("transform", "translate("+format(-10)+","+0+") rotate(90) ");
				g_sample.appendChild(label);

				
				
				final Text title = text(sn);
				g_sample.appendChild(title);

				for(i=0;i< L.size();i++) {
					final VariantSet vset = L.get(i);
					final double vsetwidth=(i+1<L.size()?L.get(i+1).x:drawingAreaWidth)-vset.x;
					
					final Element rect  = element("rect");
					rect.setAttribute("class", "area"+(i%2));
					rect.setAttribute("x", format(vset.x));
					rect.setAttribute("y", "0");
					rect.setAttribute("width",format(vsetwidth));
					rect.setAttribute("height", format(sample_height));
					rect.appendChild(element("title",vset.toString()));
					g_sample.appendChild(rect);
					
					for(VariantContext vc: vset.variants) {
						final Genotype gt= vc.getGenotype(sn);
						if(gt.isNoCall() || !gt.hasAD()) continue;
						final int ad[]= gt.getAD();
						if(ad==null || ad.length!=2) continue;
						if(gt.hasGQ() && gt.getGQ() < this.minGQ) continue;
						
						final double countR = ad[0];
						final double countA = ad[1];
						final double DP = countR + countA;
						if(DP<=0 || DP <= this.minDP) continue;
						//  HOMREF=0;  HET =0.5; HOMVAR = 1;
						double alt_ratio  = countA/DP;

						double gtx = vset.x + ((vc.getStart()-vset.getStart())/(double)vset.getLengthOnReference())*vsetwidth;
						double gty= sample_height- ( sample_height2*alt_ratio + (sample_height-sample_height2)/2.0);
						final Element circle  = element("circle");
						circle.setAttribute("class",gt.getType().name());
						circle.setAttribute("cx", format(gtx));
						circle.setAttribute("cy", format(gty));
						circle.setAttribute("r", format(circle_radius));
						circle.appendChild(element("title",vc.getStart()+" "+(vc.hasID()?vc.getID():"")+" "+vc.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/"))+" "+gt.getType().name()));
						g_sample.appendChild(circle);
						got_one_sample =  true;
						}
					
					}
				
				
				final Element frame_sample = element("rect");
				frame_sample.setAttribute("class", "frame");
				frame_sample.setAttribute("x", "0");
				frame_sample.setAttribute("y", "0");
				frame_sample.setAttribute("width",format(drawingAreaWidth));
				frame_sample.setAttribute("height", format(sample_height));

				
				
				g_sample.appendChild(frame_sample);

				y+= sample_height + space_between_samples;
				}
			// remove extra sample space
			y-=space_between_samples;
			
			svgRoot.setAttribute("width",format(this.image_width_pixel+1));
			svgRoot.setAttribute("height",format(y+1));

			
			if(!got_one_sample) {
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

			}
		catch(final Exception err) {
			throw new RuntimeException(err);
			}
		finally {
			this.document = null;
			}
		}
	}
	
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = super.oneAndOnlyOneFile(args);
			
			try (VCFReader r= VCFReaderFactory.makeDefault().open(input, true)) {
				final VCFHeader header=r.getHeader();
				if(header.getFormatHeaderLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS)==null) {
					LOG.error("FORMAT/"+VCFConstants.GENOTYPE_ALLELE_DEPTHS+" undefined in "+input);
					return -1;
					}
				if(!header.hasGenotypingData()) {
					LOG.error("No genotype in input");
					return -1;
					}
				try(BufferedReader br = IOUtils.openPathForBufferedReading(this.bedFile)) {
					try(ArchiveFactory out=ArchiveFactory.open(this.outputFile)) {
						final BedLineCodec codec = new BedLineCodec();
						br.lines().
							map(L->codec.decode(L)).
							filter(B->B!=null).
							forEach(B->{
							run(out,B,header,r);
							});
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
