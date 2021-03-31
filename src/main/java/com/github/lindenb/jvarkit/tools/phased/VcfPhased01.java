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
package com.github.lindenb.jvarkit.tools.phased;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.IntFunction;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Text;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

# Motivation

For @154sns , 10 July 2019

Displays phased genotypes from X10 genomics into SVG.

# Input

input is a list of indexed vcf files or a file with the suffix '.list' containing the full path to the vcfs

# Example

```
find . -type f -name "*.vcf.gz" > in.list
java -jar dist/vcfphased01.jar -r "chr1:1000-2000" -xp '1001,1010' in.list 
```

# Screenshots

https://twitter.com/yokofakun/status/1148964221482414080

![https://pbs.twimg.com/media/D_Hwd2dXoAAzx8g.jpg](https://pbs.twimg.com/media/D_Hwd2dXoAAzx8g.jpg)

END_DOC
*/

@Program(name="vcfphased01",
description="X10 Phased SVG to Scalar Vector Graphics (SVG)",
keywords={"x10","phased","genotypes","svg"},
creationDate="20190710",
modificationDate="20190711",
biostars=9462569
)
public class VcfPhased01 extends Launcher {
	private static final Logger LOG=Logger.build(VcfPhased01.class).make();
	private static final String X10_FORMAT_PS="PS";

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;

	@Parameter(names={"-r","--interval","--region"},description="interval CHROM:START-END",required=true)
	private String intervalStr = null;
	@Parameter(names={"-p","--pos","--highligth"},description="Highligth positions in the VCFs. (comma separated)")
	private String hihlight = "";
	@Parameter(names={"-xp","--xpos","--extra-highligth"},description="Extra Highligth positions that are not always in the vcfs. (comma separated)")
	private String extrahihlight = "";
	@Parameter(names={"-k","--knownGenes"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private String knownGeneUri=KnownGene.getDefaultUri();

	private SAMSequenceDictionary theDict= null;
	private SimpleInterval theInterval = null;
	private final List<PhasedVcf> phasedVcfs = new ArrayList<>();
	private Document dom = null;
	private int sample_height=99;
	private int sample_width=1000;
	
	
	
	private class PhasedVcf
		{
		@SuppressWarnings("unused")
		final Path path;
		String sample;
		final List<VariantContext> variants = new ArrayList<>(10_000);
		
		PhasedVcf(final Path path)
			{
			this.path= path;
			}
		
		  Element createTrack()
		 	{
			final Set<Integer> highlightpos= Arrays.stream(hihlight.split("[ ,;/]")).filter(S->!StringUtils.isBlank(S)).map(S->Integer.parseInt(S)).collect(Collectors.toCollection(HashSet::new));
			final Set<Integer> extrahighlightpos= Arrays.stream(extrahihlight.split("[ ,;/]")).filter(S->!StringUtils.isBlank(S)).map(S->Integer.parseInt(S)).collect(Collectors.toCollection(HashSet::new));
			  
			highlightpos.addAll(extrahighlightpos);
			
			
			variants.stream().forEach(V->extrahighlightpos.remove(V.getStart()));
			  
			final Element track = element("g");
			//label
			final Element label = element("text");
			label.appendChild(text(this.sample));
			label.appendChild(title(this.sample));
			label.setAttribute("x", "-2");
			label.setAttribute("y", fmt(sample_height/2));
			label.setAttribute("class", "sn");
			track.appendChild(label);
			
			final Element g =  element("g");
			track.appendChild(g);
			
			this.variants.stream().forEach(V->extrahighlightpos.remove(V.getStart()));
			
			// unphased section see https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf
			final List<Integer> switch_haploblocks = new ArrayList<>();
			if(this.variants.stream().map(V->V.getGenotype(0)).anyMatch(G->G.hasAnyAttribute(X10_FORMAT_PS))) {
				int i=0;
				while(i< this.variants.size())
					{
					final VariantContext ctxi = this.variants.get(i);
					@SuppressWarnings("deprecation")
					int psi = ctxi.getGenotype(0).getAttributeAsInt(X10_FORMAT_PS, -1);
					int j=i+1;
					while(j<this.variants.size())
						{
						VariantContext ctxj = this.variants.get(j);
						@SuppressWarnings("deprecation")
						int psj = ctxj.getGenotype(0).getAttributeAsInt(X10_FORMAT_PS, -1);
						if(psi!=psj) {
							switch_haploblocks.add(ctxj.getStart());
							break;
							}
						j++;
						}
					i=j;
					}
				LOG.info("number transitions for haplo-blocks "+switch_haploblocks.size());
				}
			
			int i=0;
			boolean side_flip = true;
			int block_start =  theInterval.getStart();
			for(;;)
				{
				double block_x1 = positionToPixel(block_start);
				double block_x2 = positionToPixel(i>=switch_haploblocks.size()?theInterval.getEnd():switch_haploblocks.get(i));
				double block_width = 1+(block_x2-block_x1);
				side_flip=!side_flip;
				
				// phase left
				Element rect = element("rect");
				rect.setAttribute("x",fmt(block_x1));
				rect.setAttribute("y", "0");
				rect.setAttribute("width",fmt(block_width));
				rect.setAttribute("height",String.valueOf(sample_height/3));
				rect.setAttribute("class",(side_flip?"phaseL":"phaseR"));
				g.appendChild(rect);
				
				// phase mid
				rect = element("rect");
				rect.setAttribute("x",fmt(block_x1));
				rect.setAttribute("y",String.valueOf(sample_height/3));
				rect.setAttribute("width",fmt(block_width));
				rect.setAttribute("height",String.valueOf(sample_height/3));
				rect.setAttribute("class",(side_flip?"phaseM":"phaseN"));
				g.appendChild(rect);
				
				// phase right
				rect = element("rect");
				rect.setAttribute("x",fmt(block_x1));
				rect.setAttribute("y",String.valueOf(2*sample_height/3));
				rect.setAttribute("width",fmt(block_width));
				rect.setAttribute("height",String.valueOf(sample_height/3));
				rect.setAttribute("class",(side_flip?"phaseR":"phaseL"));
				g.appendChild(rect);

				
				
				if(i>=switch_haploblocks.size()) break;
				block_start = switch_haploblocks.get(i);
				i++;
				}
			
			// hightlight
			for(final VariantContext ctx:this.variants) {
				if(!highlightpos.contains(ctx.getStart())) continue;
				Element rect = element("rect");
				rect.setAttribute("x",fmt(positionToPixel(ctx.getStart())));
				rect.setAttribute("y",String.valueOf(-5));
				rect.setAttribute("width",fmt(1));
				rect.setAttribute("height",String.valueOf(sample_height+10));
				rect.setAttribute("class","highlight");
				rect.appendChild(title(StringUtils.niceInt(ctx.getStart())));
				g.appendChild(rect);

			}



			
			for(final VariantContext ctx: this.variants)
				{
				final double gt_radius=1;
				double x = positionToPixel(ctx.getStart());
				
				String title= ctx.getContig()+":"+ctx.getStart()+" "+ctx.getType().name()+" "+ctx.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/"));
				final Genotype gt = ctx.getGenotype(0);
				if(gt.isNoCall()) continue;
				

				if(ctx.isFiltered())
					{
					title+=" FILTERED";
					}
				final List<Allele> alleles = gt.getAlleles();
				if(alleles.size()!=2) LOG.warn("No a diploid ??");
				for(int idx=0;idx< alleles.size();++idx) {
					final Allele a  = alleles.get(idx);
					if(a.isNoCall()) continue;
					
					final Element circle = element("circle");
					circle.setAttribute("cx", fmt(x));
					circle.setAttribute("r", fmt(gt_radius *  (highlightpos.contains(ctx.getStart())?2:1)));

					if(a.isReference()) {
						circle.setAttribute("class", "gtR");
						}
					else if(a.isSymbolic()) {
						circle.setAttribute("class", "gtS");
						}
					else if(ctx.isIndel()) {
						circle.setAttribute("class", "gtD");
						}
					else
						{
						circle.setAttribute("class", "gtA");
						}
					
					double y;
					if(!gt.isPhased() || ctx.isFiltered())
						{
						y = sample_height/2 +  gt_radius*(idx==0?-1:1);
						}
					else if(idx==0) {
						y = sample_height/6;
						
						}
					else {
						y = 2*sample_height/3+ sample_height/6 ;
						}
					
					if(highlightpos.contains(ctx.getStart())) y+=10;
					
					circle.setAttribute("cy",fmt(y));
					circle.appendChild(title(title+"  ("+a.getDisplayString()+")"));
					if(!a.equals(Allele.SPAN_DEL)) {
						g.appendChild(circle);
						}
					}
				
				
				}
			// extra hightlight
			for(final Integer pos: extrahighlightpos) {
				if(pos< theInterval.getStart() || pos> theInterval.getEnd()) continue;
				Element rect = element("rect");
				rect.setAttribute("x",fmt(positionToPixel(pos)));
				rect.setAttribute("y",String.valueOf(-5));
				rect.setAttribute("width",fmt(1));
				rect.setAttribute("height",String.valueOf(sample_height+10));
				rect.setAttribute("class","xhighlight");
				rect.appendChild(title(StringUtils.niceInt(pos)));
				g.appendChild(rect);
				}
			
			// frame
			Element rect = element("rect");
			rect.setAttribute("x", "0");
			rect.setAttribute("y","0");
			rect.setAttribute("width",String.valueOf(sample_width));
			rect.setAttribute("height",String.valueOf(sample_height));
			rect.setAttribute("class","frame");
			g.appendChild(rect);

			return track;
		 	}
		
		}
	private double positionToPixel(int pos) {
		return sample_width*((pos-theInterval.getStart())/(double)theInterval.getLengthOnReference());
	}
	private void openPhasedVcf(final Path path) throws IOException {
		final PhasedVcf phased = new PhasedVcf(path);
		try(final VCFReader vcfFileReader  = VCFReaderFactory.makeDefault().open(path,true)) {
			final VCFHeader header = vcfFileReader.getHeader();
			final List<String> samples = header.getSampleNamesInOrder();
			if(samples.isEmpty()) throw new IOException("no sample in "+path);
			if(samples.size()>1) throw new IOException("more than one sample in "+path);

			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			if(theDict==null) {
				theDict = dict;
				}
			else if(!SequenceUtil.areSequenceDictionariesEqual(dict, theDict))
				{
				throw new JvarkitException.DictionariesAreNotTheSame(dict, theDict);
				}
			if(this.theInterval==null) 
				{
				this.theInterval = IntervalParserFactory.
						newInstance().
						dictionary(dict).
						make().
						apply(this.intervalStr).
						orElseThrow(IntervalParserFactory.exception(this.intervalStr));
				}
			phased.sample = samples.get(0);
			phased.variants.addAll(vcfFileReader.query(this.theInterval).
					stream().
					collect(Collectors.toList()));
			}
		for(int i=0;i< phased.variants.size();i++) {
			
			}
		this.phasedVcfs.add(phased);
		}
	
	private String fmt(double f) {
		return String.format("%.2f",f);
	}
	
	private Element element(final String tagName) {
		return this.dom.createElementNS(SVG.NS,tagName);
		}
	private Text text(final String text) {
		return this.dom.createTextNode(text);
		}
	private Element title(final String s)
		{
		final Element e = element("title");
		e.appendChild(text(s));
		return e;
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final int margin_top=30;
			final int margin_left=100;
			for(final Path path:IOUtils.unrollPaths(args)) {
				this.openPhasedVcf(path);
				}
			if(this.phasedVcfs.isEmpty()) {
				LOG.error("no vcf defined");
				return -1;
				}
			final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			final DocumentBuilder db = dbf.newDocumentBuilder();
			this.dom = db.newDocument();
			
			Element svgRoot = this.dom.createElementNS(SVG.NS, "svg");
			svgRoot.setAttribute("style", "stroke-width:0.5");
			this.dom.appendChild(svgRoot);
			svgRoot.appendChild(title("vcfphased01"));
			
			svgRoot.appendChild(title(theInterval.getContig()+":"+StringUtils.niceInt(theInterval.getStart())+"-"+StringUtils.niceInt(theInterval.getEnd())));
			
			Element style= element("style");
			style.appendChild(text(
					  ".frame {fill:none;stroke:darkgray;} "
					+ ".phaseL {fill:#d5f4e6;stroke:lightray;} "
					+ ".phaseM {fill:gainsboro;stroke:lightray;} "
					+ ".phaseN {fill:#c6c6c6;stroke:lightray;} "
					+ ".phaseR {fill:#f1e3dd;stroke:lightray;} "
					+ ".gtM {fill:gray ;stroke:none;} "
					+ ".gtR {fill:green;stroke:lightseagreen ;} "
					+ ".gtA {fill:red;stroke:maroon ;} "
					+ ".gtS {fill:yellow;stroke:orange ;} "
					+ ".gtD {fill:blue;stroke:orange ;} "
					+ ".sn {text-anchor:end;} "
    				+ ".kgexon {fill:lightray;stroke:black;}"
    				+ ".kgcds {fill:yellow;stroke:black;opacity:0.7;}"
    				+ ".kgname {fill:black;stroke:none;font-size:6px;text-anchor:end;}"
    				+ ".kgtr {fill:none;stroke:black;}"
					+ ".unphased {fill:yellow; stroke:red;fill-opacity:0.2;}"
					+ ".highlight {fill:orange;stroke:red;fill-opacity:0.3;}"
					+ ".xhighlight {fill:darkorchid;stroke:none;fill-opacity:0.3;}"
					
					));
			svgRoot.appendChild(style);
			Element main=element("g");
			svgRoot.appendChild(main);

			double y = margin_top;
			
			for(int i=0;i< this.phasedVcfs.size();i++)
				{
				final Element g = this.phasedVcfs.get(i).createTrack();
				g.setAttribute("transform", "translate("+margin_left+","+y+")");
				y+=sample_height+10;
				main.appendChild(g);
				}
			
			final Element mainTitle=element("text");
			mainTitle.appendChild(text(theInterval.getContig()+":"+StringUtils.niceInt(theInterval.getStart())+"-"+StringUtils.niceInt(theInterval.getEnd())));
			mainTitle.setAttribute("x", "5");
			mainTitle.setAttribute("y", fmt(margin_top/2));
			main.appendChild(mainTitle);

			
			
			final List<KnownGene> genes = new ArrayList<>();
			if(!StringUtils.isBlank(this.knownGeneUri)) {
				final ContigNameConverter convert = ContigNameConverter.fromOneDictionary(theDict);
				try(BufferedReader br = IOUtils.openURIForBufferedReading(this.knownGeneUri)) {
					br.lines().
						filter(S->!StringUtils.isBlank(S)).
						filter(S->!S.startsWith("#")).
						map(S->CharSplitter.TAB.split(S)).
						map(T->new KnownGene(T)).
						filter(K->theInterval.getContig().equals(convert.apply(K.getContig()))).
						filter(K->CoordMath.overlaps(theInterval.getStart(), theInterval.getEnd(), K.getTxStart()+1,  K.getTxEnd())).
						forEach(K->genes.add(K));

					}
				}
			final double featureHeight=12;
			final double TRANSCRIPT_HEIGHT=10;
			for(final KnownGene kg:genes)
				{
				final IntFunction<Integer> trim = P-> Math.max(theInterval.getStart(),Math.min(theInterval.getEnd(), P));
				//int cdsHeigh= 5;
				double exonHeight=TRANSCRIPT_HEIGHT-5;
				double midY=TRANSCRIPT_HEIGHT/2;
		
				Element g = element("g");
				main.appendChild(g);
				
				g.setAttribute("transform", "translate("+margin_left+","+y+")");
				
				g.appendChild(title(kg.getName()));
				
				Element label = element("text");
				label.setAttribute("x", fmt(-10));
				label.setAttribute("y", fmt(featureHeight-featureHeight/2.0));
				label.setAttribute("class", "kgname");
				label.appendChild(text(kg.getName()));
				g.appendChild(label);
				
				/* transcript line */
				Element line = element("line");
				
				line.setAttribute("class","kgtr");
				line.setAttribute("x1",fmt(positionToPixel(trim.apply(kg.getTxStart()+1))));
				line.setAttribute("y1",fmt(midY));
				line.setAttribute("x2",fmt(positionToPixel(trim.apply(kg.getTxEnd()))));
				line.setAttribute("y2",fmt(midY));
				g.appendChild(line);
				
			
				/* exons */
				for(final KnownGene.Exon exon:kg.getExons())
					{
					if(!CoordMath.overlaps(theInterval.getStart(), theInterval.getEnd(), exon.getStart()+1, exon.getEnd())) continue;

					Element r = element("rect");
					r.setAttribute("class","kgexon");
					
					r.setAttribute("x",fmt(positionToPixel(trim.apply(exon.getStart()+1))));
					r.setAttribute("y",String.valueOf(midY-exonHeight/2));
					r.setAttribute("width",fmt(positionToPixel(trim.apply(exon.getEnd()))-positionToPixel((trim.apply(exon.getStart())))));
					r.setAttribute("height",fmt(exonHeight));
					r.appendChild(title(exon.getName()));
					g.appendChild(r);
					}
			
				
				//String label=String.format("%15s", g.getName());
				//w.writeEmptyElement("path");
				//double fontHeight=Math.min(10,0.8*TRANSCRIPT_HEIGHT);
				//w.writeAttribute("d",this.hershey.svgPath(label,-insets.left,midY-fontHeight/2,insets.left*0.9,fontHeight));
	
				
				y+=featureHeight;
				}

			
			
			svgRoot.setAttribute("width", fmt(sample_width+margin_left+1));
			svgRoot.setAttribute("height",fmt(y+10));
			
			// dump XML
			final TransformerFactory transformerFactory = TransformerFactory.newInstance();
            final Transformer transformer = transformerFactory.newTransformer();
            final Writer w = super.openPathOrStdoutAsPrintWriter(this.outputFile);
            transformer.transform(new DOMSource(this.dom), new StreamResult(w));
            w.flush();
            w.close();
			return 0;
			} 
		catch (final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	public static void main(final String[] args) {
		new VcfPhased01().instanceMainWithExit(args);
	}
}
