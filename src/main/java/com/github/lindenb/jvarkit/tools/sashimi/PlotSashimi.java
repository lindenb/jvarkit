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
package com.github.lindenb.jvarkit.tools.sashimi;

import java.io.File;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
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
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Exon;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

/**
BEGIN_DOC

## Input

input is a set of indexed BAM files or a file with the '.list' suffix containing the path to the bams.

## Example

```
 java -jar dist/plotsashimi.jar -r "chr3:38597400-38599300" -m jeter.mf --gtf jeter.gtf ENCFF331CGL.rnaseq.bam -o TMP

$ find TMP/ -name "*.svg"
TMP/86/f9362065ddce1af4b31b47be80fff6/chr3_38599300_38599300.svg

$ cat jeter.mf 
#chrom	start	end	bam	Sample	Genes	svg
chr3	38597399	38599300	ENCFF331CGL.rnaseq.bam	SCN5A	.	86/f9362065ddce1af4b31b47be80fff6/chr3_38599300_38599300.svg

```


```
$ cat jeter.list
ENCFF331CGL.rnaseq.bam


$ cat jeter.bed 
chr3	38595150	38599347
chr3	38595350	38599500



$ java -jar dist/plotsashimi.jar -r jeter.bed -m jeter.mf --gtf jeter.gtf jeter.list -o jeter.zip

$ unzip -l jeter.zip 
Archive:  jeter.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
    36566  2019-12-05 15:52   f6/d0bf04098aa24eb3111666ed171a11/chr3_38599347_38599347.svg
    36364  2019-12-05 15:52   82/fd0c7ad61ecc0fa22eda93e75a6943/chr3_38599500_38599500.svg
---------                     -------
    72930                     2 files


$ cat jeter.mf | column -t
#chrom  start     end       bam                     Genes   Samples  svg
chr3    38595150  38599347  ENCFF331CGL.rnaseq.bam  SCN5A   .        f6/d0bf04098aa24eb3111666ed171a11/chr3_38599347_38599347.svg
chr3    38595350  38599500  ENCFF331CGL.rnaseq.bam  SCN5A   .        82/fd0c7ad61ecc0fa22eda93e75a6943/chr3_38599500_38599500.svg

```


## Screenshot

* https://twitter.com/yokofakun/status/1202587424725127168

![https://twitter.com/yokofakun/status/1202587424725127168](https://pbs.twimg.com/media/ELBy4vAX0AABSzF?format=jpg&name=small)

* https://twitter.com/yokofakun/status/1202905778140712960

![https://twitter.com/yokofakun/status/1202905778140712960](https://pbs.twimg.com/media/ELGUbZ_W4AAMKxj?format=png&name=small)

* https://twitter.com/yokofakun/status/1207337424935936001

![https://twitter.com/yokofakun/status/1207337424935936001](https://pbs.twimg.com/media/EMFS-xGXsAAVR20?format=jpg&name=small)

END_DOC

**/
@Program(name="plotsashimi",
description="Print Sashimi plots from Bam",
keywords={"bam","visualization","svg","rna","exon","rnaseq"},
modificationDate="20191104",
creationDate="20191117",
biostars=497894
)
public class PlotSashimi extends Launcher {
private static final Logger LOG = Logger.build(PlotSashimi.class).make();

@Parameter(names={"-g","--gtf"},description=GtfReader.OPT_DESC)
private Path gtfPath = null;

@Parameter(names={"-r","--region","--interval"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,required=true)
private IntervalListProvider intervalListProvider= IntervalListProvider.empty();

@Parameter(names={"-R","--reference"},description="For Reading CRAM. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
private Path faidx = null;

@Parameter(names={"-o","--out"},description=ArchiveFactory.OPT_DESC,required=true)
private Path outputFile=null;

@Parameter(names={"-m","--manifest"},description="Manifest Bed file output containing chrom/start/end of each gene")
private Path manifestFile = null;
@Parameter(names={"-D","--use-deletion"},description="also use the D operator in the cigar string (default is use only 'N').")
private boolean use_D_operator= false;
@Parameter(names={"-w","--width"},description="image width.")
private int image_width_pixel  = 1_000;
@Parameter(names={"--skip-empty"},description="Do not generate a SVG file if there is no read in the interval")
private boolean skip_region_without_read = false;
@Parameter(names={"--mapq"},description="Min mapping quality")
private int min_mapq = 0;
@Parameter(names={"--force-max-coverage"},description="Force the maximum coverage to this value. ignored if <=0 ")
private int force_max_coverage = 0;
@Parameter(names={"--css"},description="Custom CSS stylesheet")
private Path cssPath = null;
@Parameter(names={"-u","--url","--hyperlink"},description= "creates a hyperlink an area is clicked. " + Hyperlink.OPT_DESC,converter=Hyperlink.StringConverter.class,splitter=NoSplitter.class)
private Hyperlink hyperlinkType = Hyperlink.empty();

@Parameter(names={"--partition"},description=SAMRecordPartition.OPT_DESC)
private SAMRecordPartition partition= SAMRecordPartition.sample;
@Parameter(names= {"--gzip"},description="Generate gzipped compressed svg files.")
private boolean compressed_svg=false;

@SuppressWarnings("serial")
@DynamicParameter(names = "--param", description = "Other parameters.",hidden=true)
public Map<String, String> dynamicParams = new HashMap<String,String>() {{{
	put("coverage.height","300");
	}}};

private final IntervalTreeMap<Gene> geneMap = new IntervalTreeMap<>();
private Document document;
private final DecimalFormat decimalFormater = new DecimalFormat("##.##");

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
/** convert double to string */
private String format(double v)
	{
	return this.decimalFormater.format(v);
	}
/** best ticks separation */
private int bestTicks(final int max) {
	if(max<=10) return 1;
	final int ndigit=(int)Math.ceil(Math.log10(max-1));
	return Math.max(1,(int)Math.pow(10, ndigit-2));
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

/** create the SVG itself */
private void plotSashimi(
	final ArchiveFactory archive,
	final SamReader samReader,
	final Locatable interval,
	final Path bamPath,
	final PrintWriter manifest
	) {
	final int drawing_width = Math.max(100,this.image_width_pixel);
	final int coverageHeight =  Math.max(100,Integer.parseInt(this.dynamicParams.getOrDefault("coverage.height","300")));
	final double pixelperbase =  drawing_width/(double)interval.getLengthOnReference();
	final SAMFileHeader header= samReader.getFileHeader();
	final Collection<Gene> genes= this.geneMap.getOverlapping(interval);
	final Set<String> geneNames = genes.stream().map(G->G.getGeneName()).filter(S->!StringUtils.isBlank(S)).collect(Collectors.toCollection(TreeSet::new));
	
	/** extract the sample name or just use the filename */
	final String sampleName = StringUtils.ifBlank(
		header.getReadGroups().
			stream().
			map(G->this.partition.apply(G)).
			filter(S->!StringUtils.isBlank(S)).
			sorted().
			collect(Collectors.joining(";"))
			,
			bamPath.getFileName().toString()
			);
	
	final Function<Integer, Double> pos2pixel = POS-> (POS - interval.getStart())/(double)interval.getLengthOnReference() * drawing_width;
	
	final Counter<SimpleInterval> gaps = new  Counter<>();
	final int coverage[] = new int[interval.getLengthOnReference()];
	try(SAMRecordIterator iter=samReader.queryOverlapping(interval.getContig(), interval.getStart(), interval.getEnd())) {
		/** no read here, skip */
		boolean got_one = false;
		while(iter.hasNext()) {
			final SAMRecord rec = iter.next();
			if(rec.getReadUnmappedFlag()) continue;
			if(rec.getReadFailsVendorQualityCheckFlag()) continue;
			if(rec.isSecondaryOrSupplementary()) continue;
			if(rec.getDuplicateReadFlag()) continue;
			if(rec.getMappingQuality()< this.min_mapq) continue;
			final Cigar cigar =rec.getCigar();
			if(cigar==null || cigar.isEmpty()) continue;
			got_one = true;
			int ref= rec.getAlignmentStart();
			for(final CigarElement ce:cigar) {
				if(ref> interval.getEnd()) break;
				final CigarOperator op = ce.getOperator();
				
				if(op.equals(CigarOperator.N) || (use_D_operator && op.equals(CigarOperator.D) ))  {
					gaps.incr(new SimpleInterval(rec.getContig(),ref,ref+ce.getLength()-1));
					}

				
				if(op.consumesReferenceBases()) {
					if(op.consumesReadBases()) {
						for(int x=0;x<ce.getLength();++x) {
							final int pos1 = ref+x;
							if(pos1< interval.getStart()) continue;
							if(pos1> interval.getEnd()) break;
							coverage[pos1-interval.getStart()]++;
							}
						}
					ref+=ce.getLength();
					}
				}
			}

		if(!got_one && this.skip_region_without_read) return;
		}
	
		final int max_coverage;
		
		if(this.force_max_coverage>0) {
			max_coverage = this.force_max_coverage;
		} else
			{
			max_coverage = Math.max(1,Arrays.stream(coverage).max().orElse(0));
			}

	
		while(this.document.hasChildNodes()) {
			this.document.removeChild(this.document.getFirstChild());
			}
	
		final Element svgRoot = element("svg");
		this.document.appendChild(svgRoot);
		
		/* SVG title */
		{
		final Element title = element("title");
		svgRoot.appendChild(title);
		title.appendChild(text(interval.toString() +(!geneNames.isEmpty() && geneNames.size()<3?" "+String.join(" ", geneNames):"")));
		}
		
		/* SVG style */
		{
		final Element style = element("style");
		svgRoot.appendChild(style);
		style.appendChild(text(
				this.cssPath==null?
				".coverage { fill:red;fill:url('#grad01')} " +
				".maintitle {text-anchor:middle;fill:blue} "+
				".sample {fill:blue;font-size:7px;} "+
				".frame { fill:none; stroke: darkgray;} " +
				".arcK { fill:none; stroke: blue; stroke-linecap:round;opacity:0.8;} " +
				".arcU { fill:none; stroke: red; stroke-linecap:round;opacity:0.8;} " +
				".transcript { fill:darkgray; stroke: darkgray;} " +
				".exon { fill:green; stroke: darkgray;} " +
				".frame { fill:none; stroke: darkgray;} " +
				".rulerline {stroke:lightgray;stroke-width:0.5px;}\n"+
				".exonline {stroke:green;stroke-width:0.5px;opacity:0.5;}\n"+
				".rulerlabel {stroke:gray;stroke-width:0.5px;font-size:7px;}\n" +
				"a {cursor: pointer;}\n"
				: IOUtils.slurpPath(this.cssPath)
				));
		}
		
		
		// SVG def
		{
		final Element defs = element("defs");
		svgRoot.appendChild(defs);
		//linear gradient	
			{
			Element grad= element("linearGradient");
			defs.appendChild(grad);
			grad.setAttribute("id","grad01");
			grad.setAttribute("gradientTransform","rotate(90)");
			
			Element stop = element("stop");
			grad.appendChild(stop);
			stop.setAttribute("offset", "0%");
			stop.setAttribute("stop-color",(max_coverage>50?"red":max_coverage>20?"green":"blue"));
			
			stop = element("stop");
			grad.appendChild(stop);
			stop.setAttribute("offset", "100%");
			stop.setAttribute("stop-color", "darkblue");
			
			}
		}		

		
		final Element descr = element("desc");
		svgRoot.appendChild(descr);
		descr.appendChild(text("Author: Pierre Lindenbaum\n" +
				JVarkitVersion.getInstance().getCompilationDate()+"\n"+
				JVarkitVersion.getInstance().getGitHash()
				));
	
		final Element maing = element("g");
		svgRoot.appendChild(maing);
		int y=0;
		
		// main title
		Element gtitle= element("text",new SimpleInterval(interval).toNiceString()+(StringUtils.isBlank(sampleName)?"":" "+sampleName)+ (geneNames.isEmpty()?"": " "+String.join(" ", geneNames)));
		gtitle.setAttribute("class", "maintitle");
		gtitle.setAttribute("x", format(drawing_width/2));
		gtitle.setAttribute("y", "15");
		svgRoot.appendChild(gtitle);
		y+=20;
		// sample name
		if(!StringUtils.isBlank(sampleName))
			{	
			gtitle= element("text",sampleName);
			gtitle.setAttribute("class", "sample");
			gtitle.setAttribute("x", "5");
			gtitle.setAttribute("y", "20");
			svgRoot.appendChild(gtitle);
			}
		y+=50;
		final int prev_y=y;
		
		/** horizontal ruler */
		{
		final Element ruler_gh = element("g");
		maing.appendChild(ruler_gh);
		
		
		final int sep= bestTicks(interval.getLengthOnReference());
		
		for(int pos =interval.getStart(); pos<= interval.getEnd();++pos)
			{
			if(pos%sep!=0) continue;
			double x= pos2pixel.apply(pos);
			final Element line = element("line");
			ruler_gh.appendChild(line);
			line.setAttribute("class","rulerline");
			line.appendChild(element("title",StringUtils.niceInt(pos)));
			line.setAttribute("x1",format(x));
			line.setAttribute("x2",format(x));
			line.setAttribute("y1",format(y));
			line.setAttribute("y2",format(y+coverageHeight));
			
			final Element label = element("text",StringUtils.niceInt(pos));
			label.setAttribute("class","rulerlabel");
			label.setAttribute("x","0");
			label.setAttribute("y","0");
			label.setAttribute("transform", "translate("+format(x)+","+y+") rotate(90) ");
			ruler_gh.appendChild(label);
			}
		}

		/** vertical ruler */
		{
		final Element ruler_gv = element("g");
		maing.appendChild(ruler_gv);
		
		
		final int sep= bestTicks(max_coverage);
		
		for(int pos =0; pos<= max_coverage;++pos)
			{
			if(pos%sep!=0) continue;
			
			double ry= (int)(y + coverageHeight - (pos/(double) max_coverage)*coverageHeight);
			final Element line = element("line");
			ruler_gv.appendChild(line);
			line.setAttribute("class","rulerline");
			line.appendChild(element("title",StringUtils.niceInt(pos)));
			line.setAttribute("x1",format(0));
			line.setAttribute("x2",format(drawing_width));
			line.setAttribute("y1",format(ry));
			line.setAttribute("y2",format(ry));
			
			final Element label = element("text",StringUtils.niceInt(pos));
			label.setAttribute("class","rulerlabel");
			label.setAttribute("x","1");
			label.setAttribute("y", format(ry));
			ruler_gv.appendChild(label);
			}
		}

		/** vertical lines of exons */
		final Element exon_v = element("g");
		
		
		final Element covPath = element("path");
		covPath.setAttribute("class", "coverage");
		maing.appendChild(covPath);
		
		
		final StringBuilder sb = new StringBuilder();
		sb.append( "M 0 "+format(y+coverageHeight));
		for(int k=0;k< coverage.length;k++)
			{
			//if(k+1< coverage.length && coverage[k]==coverage[k+1]) continue; 
			final double dpy= y+ coverageHeight - coverageHeight*(coverage[k]/(double)max_coverage);
			sb.append(" L "+format(pixelperbase*k)+" "+format(dpy));
			}
		sb.append(" L "+format(drawing_width)+" "+format(y+ coverageHeight));
		sb.append(" Z");
		covPath.setAttribute("d", sb.toString());
		covPath.appendChild(element("title","Coverage. Max:"+StringUtils.niceInt(max_coverage)));
		
		int next_y=y+coverageHeight;
		/* plot arc */
		if(!gaps.isEmpty()) {
			 boolean drawAbove  = true;
			int max_occurence = (int)gaps.count(gaps.getMostFrequent());
			for(final SimpleInterval intron:gaps.keySet()) {
				final int occurence = (int)gaps.count(intron);
				
				boolean is_known_intron = genes.stream().
					flatMap(G->G.getTranscripts().stream()).
					flatMap(T->T.getIntrons().stream()).
					filter(E->E.overlaps(intron)).anyMatch(
						I->I.getStart()==intron.getStart() &&
						   I.getEnd() == intron.getEnd()
							);
				
				
				final int junctionStart = intron.getStart()-1;
                final int junctionEnd = intron.getEnd()+1;
                
                
				if(!CoordMath.encloses(interval.getStart(),interval.getEnd(), junctionStart, junctionEnd)) continue;
                
				final double xstart = pos2pixel.apply(junctionStart);
				final double xend = pos2pixel.apply(junctionEnd);
				double ystart=  y+ coverageHeight - coverageHeight*(coverage[junctionStart - interval.getStart()]/(double)max_coverage);
				double yend=  y+ coverageHeight - coverageHeight*(coverage[junctionEnd- interval.getStart()]/(double)max_coverage);
				
				final Element arc= element("path");
				sb.setLength(0);
				
				double x_mid = (xend - xstart)/2.0;
				double x2 = xstart + x_mid;
				final double y2;
				
				// small gap: always print it under xaxis
				if(xend - xstart < 30) drawAbove=true;
				
				if(drawAbove)
					{
					ystart = y+ coverageHeight;
					yend = y+ coverageHeight;
					y2= y+ coverageHeight + x_mid;
					next_y = (int)Math.max(next_y, y+ coverageHeight + x_mid/2 + 10);
					}
				else
					{
					y2= Math.max(0,ystart+(yend-ystart)/2.0 - x_mid);
					}
				sb.append("M "+format(xstart)+" "+format(ystart));
				sb.append(" Q "+format(x2)+" "+format(y2)+" "+format(xend)+" "+format(yend));
				
				arc.setAttribute("d", sb.toString());
				arc.setAttribute("class","arc"+(is_known_intron?"K":"U"));
				
				final double stroke_width = 1 /* show very small one */ + (occurence/(double)max_occurence)*10;
				
				arc.setAttribute("style", "stroke-width:"+format(stroke_width)+"px;");
				arc.appendChild(element("title",new SimpleInterval(interval.getContig(),junctionStart,junctionEnd).toNiceString()+" ("+StringUtils.niceInt(occurence)+") "+(is_known_intron?"known":"unknown")));
				maing.appendChild(wrapLoc(arc,intron));
				
			
				
                drawAbove = !drawAbove;
			}
		}
		
		y = next_y;
		y+=2;


		/** pileup transcripts */	
		final List<List<Transcript>> transcriptRows = new ArrayList<>();
		for(final Transcript transcript: genes.
				stream().
				flatMap(L->L.getTranscripts().stream()).
				filter(T->T.overlaps(interval)).
				sorted((A,B)->Integer.compare(A.getStart(), B.getStart())).
				collect(Collectors.toList())) {
			int rowidx=0;
			while(rowidx< transcriptRows.size()) {
				final List<Transcript> row =  transcriptRows.get(rowidx);
				final Transcript last = row.get(row.size()-1);
				if(!last.overlaps(transcript)) {
					row.add(transcript);
					break;
					}
				rowidx++;
				}
			if(rowidx==transcriptRows.size())
				{
				final List<Transcript> row = new ArrayList<>();
				row.add(transcript);
				transcriptRows.add(row);
				}
			}
		/** plot transcripts */
		final Element transcripts_g = element("g");
		maing.appendChild(transcripts_g);
		
		final int transcript_height = Math.max(10, Integer.parseInt(this.dynamicParams.getOrDefault("transcript.height", "12")));
		for(final List<Transcript> row: transcriptRows) {
			final Element grow = element("g");
			transcripts_g.appendChild(grow);
			
			for(final Transcript transcript: row) {
				final Element transcript_g = element("g");
				grow.appendChild(transcript_g);
				final Element tr = element("line");
				final double midy= y+transcript_height/2.0;
				tr.setAttribute("class", "transcript");
				double tr_x1 = Math.max(1.0,pos2pixel.apply(transcript.getStart()));
				tr.setAttribute("x1", format(tr_x1));
				tr.setAttribute("y1", format(midy));
				tr.setAttribute("x2", format(Math.min(drawing_width,pos2pixel.apply(transcript.getEnd()))));
				tr.setAttribute("y2", format(midy));
				tr.appendChild(element("title",transcript.getId()+" "+transcript.getGene().getGeneName()));

				transcript_g.appendChild(wrapLoc(tr,transcript));
				
				final Element label=element("text",transcript.getId()+" "+transcript.getGene().getGeneName());
				label.setAttribute("class", "rulerlabel");
				label.setAttribute("x", format(tr_x1));
				label.setAttribute("y", format(midy));
				transcript_g.appendChild(label);

				
				
				for(final Exon exon:transcript.getExons()) {
					if(!exon.overlaps(interval)) continue;
					final Element exon_rect = element("rect");
					
					
					exon_rect.setAttribute("class", "exon");
					final double exonx1= Math.max(1.0,pos2pixel.apply(exon.getStart()));
					final double exonx2= Math.min(drawing_width,pos2pixel.apply(exon.getEnd()));
					exon_rect.setAttribute("x", format(exonx1));
					exon_rect.setAttribute("y", format(y));
					exon_rect.setAttribute("height", format(transcript_height));
					exon_rect.setAttribute("width", format(exonx2-exonx1));
					exon_rect.appendChild(element("title",exon.getName()+" "+transcript.getId()+" "+transcript.getGene().getGeneName()));
					transcript_g.appendChild(wrapLoc(exon_rect,exon));
					
					
					
					for(int side=0;side<2;++side) {
						final double x = pos2pixel.apply(side==0?exon.getStart():exon.getEnd());
						final Element line = element("line");
						exon_v.appendChild(line);
						
						line.setAttribute("class", "exonline");
						line.setAttribute("x1", format(x));
						line.setAttribute("x2", format(x));
						line.setAttribute("y2", format(y+transcript_height));
						line.setAttribute("y1", format(prev_y));
							
						}
					
					}
				}
			y+=transcript_height+2;
			}
		maing.appendChild(exon_v);
		
		
		/* final frame */
		final Element frame_rect = element("rect");
		frame_rect.setAttribute("class", "frame");
		frame_rect.setAttribute("x", "0");
		frame_rect.setAttribute("y", "0");
		frame_rect.setAttribute("width", format(drawing_width));
		frame_rect.setAttribute("height",format(y));
		svgRoot.appendChild(frame_rect);
		
		svgRoot.setAttribute("width",format(drawing_width+1));
		svgRoot.setAttribute("height",format(y+1));
		
		
		try {
		
		final Transformer tr = TransformerFactory.newInstance().newTransformer();
		
		final String md5 = StringUtils.md5(interval.getContig()+":"+interval.getStart()+":"+interval.getEnd()+":"+bamPath.toString());
		final String filename =  md5.substring(0,2) + File.separatorChar + md5.substring(2) + 
					File.separator+ interval.getContig()+"_"+interval.getStart()+"_"+interval.getEnd()+
					(StringUtils.isBlank(sampleName)?"":"."+sampleName.replaceAll("[/\\:]", "_")) +
					".svg"+(this.compressed_svg?".gz":"");
		

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
		manifest.print(interval.getContig());
		manifest.print('\t');
		manifest.print(interval.getStart()-1);
		manifest.print('\t');
		manifest.print(interval.getEnd());
		manifest.print('\t');
		manifest.print(bamPath.toString());
		manifest.print('\t');
		manifest.print(geneNames.isEmpty()?".":String.join(",",geneNames));
		manifest.print('\t');
		manifest.print(StringUtils.isBlank(sampleName)?".":sampleName);
		manifest.print('\t');
		manifest.print((archive.isTarOrZipArchive()?"":this.outputFile.toString()+File.separator)+filename);
		manifest.println();
		} catch(final Exception err) {
			throw new RuntimeException(err);
		}
	}

@Override
public int doWork(final List<String> args) {
	ArchiveFactory archive=null;
	PrintWriter manifest = null;
	try
		{
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		dbf.setNamespaceAware(true);
		DocumentBuilder db = dbf.newDocumentBuilder();
		this.document = db.newDocument();

		
		final SamReaderFactory srf = super.createSamReaderFactory();
		if(faidx!=null) {
			srf.referenceSequence(this.faidx);
			}
		if(this.gtfPath!=null) {
			try(GtfReader gtfReader=new GtfReader(this.gtfPath)) {
				if(this.faidx!=null) gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(this.faidx)));
				gtfReader.
					getAllGenes().
					stream().
					forEach(G-> this.geneMap.put(new Interval(G), G));
				}
			}
		
		archive = ArchiveFactory.open(this.outputFile);
		
		manifest = new PrintWriter(this.manifestFile==null?new NullOuputStream():IOUtils.openPathForWriting(manifestFile));
		manifest.println("#chrom\tstart\tend\tbam\tGenes\tSamples\tsvg");

		
		for(final Path bam: IOUtils.unrollPaths(args)) {
			try(SamReader sr = srf.open(bam)) {
				if(!sr.hasIndex()) {
					LOG.error("Bam is not indexed "+bam);
					return -1;
					}
				final SAMFileHeader header= sr.getFileHeader();
				final ArchiveFactory final_archive = archive;
				final PrintWriter final_manifest = manifest;
				this.intervalListProvider.
					dictionary(header.getSequenceDictionary()).
					stream().
					forEach(R->{
						plotSashimi(final_archive,sr,R,bam,final_manifest);
					});
					
				}
			}
		return 0;
		} 
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	finally {
		CloserUtil.close(manifest);
		CloserUtil.close(archive);
		}
	}

public static void main(final String[] args) {
	new PlotSashimi().instanceMainWithExit(args);
	}
}
