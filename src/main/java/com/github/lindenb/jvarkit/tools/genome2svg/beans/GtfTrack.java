package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;

public class GtfTrack extends Track {
	private String tabix;
	
	private static class Transcript implements Locatable {
		GTFLine transcript = null;
		final List<GTFLine> exons = new ArrayList<>();
		Transcript() {
			}
		@Override
		public String getContig() {
			return transcript.getContig();
			}
		@Override
		public int getStart() {
			return transcript.getStart();
			}
		@Override
		public int getEnd() {
			return transcript.getEnd();
			}
		}
	
	public GtfTrack() {
		
		}
	public void setTabix(String tabix) {
		this.tabix = tabix;
		}
	public String getTabix() {
		return tabix;
		}
	
	
	@Override
	public void paint(SVGContext ctx) {
		ctx.styleNode.appendChild(ctx.text("."+getId()+"strand {fill:none;stroke:gray;stroke-width:0.5px;}\n"));
		ctx.styleNode.appendChild(ctx.text("."+getId()+"kgtr {fill:none;stroke:gray;stroke-width:1px;}\n"));
		ctx.styleNode.appendChild(ctx.text("."+getId()+"kgexon {fill:lightgray;stroke:gray;stroke-width:1px;}\n"));

		

		Element polyline = ctx.element("polyline");
		ctx.defsNode.appendChild(polyline);
		polyline.setAttribute("id", "strandF");
		polyline.setAttribute("class", getId()+"strand");
		polyline.setAttribute("points", "-5,-5 0,0 -5,5");
		
		polyline = ctx.element("polyline");
		ctx.defsNode.appendChild(polyline);
		polyline.setAttribute("id", "strandR");
		polyline.setAttribute("class", getId()+"strand");
		polyline.setAttribute("points", "5,-5 0,0 5,5");
		
		
		final Element g0 = ctx.element("g");
		ctx.tracksNode.appendChild(g0);
		insertTitle(g0,ctx);
		
		final Pileup<Transcript> pileup = new Pileup<>(ctx.createCollisionPredicate());
		final GTFCodec codec = new GTFCodec();
		try(TabixReader tbx = new TabixReader(getTabix())) {
			final Map<String,Transcript> id2transcript= new HashMap<>();
			TabixReader.Iterator iter = tbx.query(ctx.loc.getContig(), ctx.loc.getStart(), ctx.loc.getEnd());
			for(;;) {
				final String line = iter.next();
				if(line==null) break;
				final GTFLine gtf = codec.decode(line);
				if(gtf==null) continue;
				if(!gtf.hasAttribute("transcript_id")) continue;
				final String transcript_id = gtf.getAttribute("transcript_id");
				if(StringUtils.isBlank(transcript_id)) continue;
				if(!(gtf.getType().equals("transcript") || gtf.getType().equals("exon"))) {
					continue;
					}
				Transcript tr = id2transcript.get(transcript_id);
				if(tr==null) {
					tr =new Transcript();
					id2transcript.put(transcript_id,tr);
					}
				if(gtf.getType().equals("transcript")) {
					tr.transcript = gtf;
					}
				else if(gtf.getType().equals("exon")) {
					tr.exons.add(gtf);
					}
				}
			id2transcript.values().
				stream().
				filter(T->T.transcript!=null && !T.exons.isEmpty()).
				sorted((A,B)->Integer.compare(A.getStart(),B.getStart())).
				forEach(T->pileup.add(T));
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}

		final double exonHeight = getFeatureHeight();
		for(List<Transcript> row:pileup.getRows()) {
			double midY= ctx.y + exonHeight/2.0;
			for(Transcript tr:row) {
				final Element g = ctx.element("g");
				g0.appendChild(g);
				
				/* transcript line */
				Element line = ctx.element("line");
				g.appendChild(line);
				line.setAttribute("class",getId()+"kgtr");
				line.setAttribute("x1",String.valueOf(ctx.pos2pixel(ctx.trimpos(tr.getStart()))));
				line.setAttribute("y1",String.valueOf(midY));
				line.setAttribute("x2",String.valueOf(ctx.pos2pixel(ctx.trimpos(tr.getEnd()))));
				line.setAttribute("y2",String.valueOf(midY));
				

				
				/* strand symbols */
				for(double pixX=0;
					pixX< ctx.image_width;
					pixX+=30)
					{
					double pos1= ctx.pixel2genomic(pixX);
					if(pos1< tr.getStart()) continue;
					if(pos1> tr.getEnd()) break;
					Element use = ctx.element("use");
					g.appendChild(use);
					use.setAttribute("href", "#strand"+(tr.transcript.isPostiveStrand()?"F":"R"));
					use.setAttribute("x",String.valueOf(pixX));
					use.setAttribute("y",String.valueOf(midY));
					}
				
				/* exons */
				for(final GTFLine exon: tr.exons)
					{
					if(!exon.overlaps(ctx.loc)) continue;
					Element rect = ctx.element("rect");
					g.appendChild(rect);
					rect.setAttribute("class",getId()+"kgexon");
					
					rect.setAttribute("x",String.valueOf(ctx.pos2pixel(ctx.trimpos(exon.getStart()))));
					rect.setAttribute("y",String.valueOf(midY-exonHeight/2));
					rect.setAttribute("width",String.valueOf(
							ctx.pos2pixel(ctx.trimpos(exon.getEnd()))-ctx.pos2pixel(ctx.trimpos(exon.getStart()))));
					rect.setAttribute("height",String.valueOf(exonHeight));
					}
				
				}
			
			// legend
			{
			final Set<String> geneNames = row.stream().map(KG->KG.transcript.getAttribute("gene_name")).collect(Collectors.toCollection(LinkedHashSet::new));
			if(geneNames.size()<10) {
				final Element legend = ctx.element("text",String.join(" ", geneNames));
				legend.setAttribute("style", "text-anchor:end;font-size:10px;");
				legend.setAttribute("x", format(-5));
				legend.setAttribute("y",String.valueOf(ctx.y+10));
				g0.appendChild(legend);
				}
			}

			
			ctx.y += exonHeight;
			ctx.y += 2;
			}
		
	}
}